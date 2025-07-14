import bpy
import bpy_extras
from bpy_extras import bmesh_utils
import bmesh
import math
import mathutils
import bl_math
import sys

class SloppyQuadUVUnfold(bpy.types.Operator):
    bl_idname = "operator.sloppy_quad_unfold"
    bl_label = "Quad Unfold UVs"
    bl_description = "Attempt to unwrap quad mesh UVs as if from flat piece of cloth, with the assumption that all face corners should be right-angled originally"
    bl_options = {'REGISTER', 'UNDO'}

    iP = bpy.props.IntProperty
    fP = bpy.props.FloatProperty
    fvP = bpy.props.FloatVectorProperty
    bP = bpy.props.BoolProperty
    eP = bpy.props.EnumProperty
    bvP = bpy.props.BoolVectorProperty
    sP = bpy.props.StringProperty

    def update_unfold_mode(self, context):
        context.scene.sloppy_props.procuv_unfold_mode = self.unfold_mode
        return None
    
    def update_initial_quad(self, context):
        context.scene.sloppy_props.procuv_initial_quad = self.initial_quad
        return None
        
    # region ProcUV Properties
    unfold_mode : eP(
        name = "Mode",
        description = "Unfolding Mode",
        items = [
            ("A", "All Directions", "Work outwards from initial quad"),
            ("B", "X First", "Work outwards in (local) X first until no more faces are found, calculating edge lengths to use in next stage, then work outwards in Y. repeating first step when more faces are found in X"),
            ("C", "Y First", "Work outwards in (local) Y first until no more faces are found, calculating edge lengths to use in next stage, then work outwards in X. repeating first step when more faces are found in Y")
            ],
        update=update_unfold_mode
        ) # type: ignore

    initial_quad : eP(
        name = "Initial Quad(s)",
        description = "Which quad (for each UV island, if there are more than one,) to start from",
        items = [
            ("A", "Most Regular", "Start with quad (or virtual quad, formed by combining most regular face with adjacent triangle as calculated from face corner angles) that is flattest and/or that has normal most similar to island average normals"),
            ("B", "Selected", "First selected quad (or virtual quad, formed from two first triangles selected, or, if only 1 triangle is selected in island, formed by combining with adjacent triangle as calculated from face corner angles) in each island (if an island has no faces selected, initial quad falls back to Most Regular)")
            ],
        update=update_initial_quad
        ) # type: ignore

    reg_flat_fac : fP(
        name = "Flatness Factor",
        description = "How much flatness influences initial quad selection",
        default = 1.0,
        min = 0.0,
        max = 2.0
        ) # type: ignore

    reg_norm_fac : fP(
        name = "Normal Factor",
        description = "How much island average normal influences initial quad selection",
        default = 1.0,
        min = 0.0,
        max = 2.0
        ) # type: ignore

    pre_calc_edge_lengths : bP(
        name = "Pre-Calculate Edge Lengths",
        description = "Calculate average edge lengths of columns and rows before setting UV coordinates",
        default = False
        ) # type: ignore

    offset_per_island : fvP(
        name = "Offset/Island",
        description = "Offset UVs per island",
        size = 2
        ) # type: ignore

    quant_avg_norm : bP(
        name = "Quantize Average Normal",
        description = "Quantize island's averaged normal to up, down, right, left",
        default = False
        ) # type: ignore

    only_move_loops_in_face : bP(
        name = "Only Edit Current Loops",
        description = "Avoid moving loops of other faces than current face when editing UV loops",
        default = False
        ) # type: ignore
    #endregion

    def execute(self, context):
        props = context.scene.sloppy_props
        bm = bmesh.from_edit_mesh(bpy.context.active_object.data)
        uv_layer = bm.loops.layers.uv.verify()
        islands = bmesh_utils.bmesh_linked_uv_islands(bm, uv_layer)
        
        self.unfold_mode = props.procuv_unfold_mode
        self.initial_quad = props.procuv_initial_quad

        # region init vars
        current_dir_array = []
        max_avg_angle = 0.0
        min_avg_angle = 360.0
        current_avg_normal = mathutils.Vector((0,0,0))
        current_normal = mathutils.Vector((0,0,0))
        current_face_center = mathutils.Vector((0,0,0))

        prev_co_dict = [
                {
                    "xy_str": "x0y0",
                    "co_ur": mathutils.Vector((0,0)),
                    "co_ul": mathutils.Vector((0,0)),
                    "co_ll": mathutils.Vector((0,0)),
                    "co_lr": mathutils.Vector((0,0)),
                },
            ]

        # region Attributes
        attr_dict = [
            {"name": "edge_u", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "edge_l", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "edge_d", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "edge_r", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "corner_ur", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "corner_ul", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "corner_ll", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "corner_lr", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "co_ur_x", "type": "FLOAT", "domain": "FACE", "layer": None},
            {"name": "co_ul_x", "type": "FLOAT", "domain": "FACE", "layer": None},
            {"name": "co_ll_x", "type": "FLOAT", "domain": "FACE", "layer": None},
            {"name": "co_lr_x", "type": "FLOAT", "domain": "FACE", "layer": None},
            {"name": "co_ur_y", "type": "FLOAT", "domain": "FACE", "layer": None},
            {"name": "co_ul_y", "type": "FLOAT", "domain": "FACE", "layer": None},
            {"name": "co_ll_y", "type": "FLOAT", "domain": "FACE", "layer": None},
            {"name": "co_lr_y", "type": "FLOAT", "domain": "FACE", "layer": None},
            {"name": "loop_corner_index", "type": "INT", "domain": "CORNER", "layer": None},
            {"name": "ind_face_uv", "type": "FLOAT_VECTOR", "domain": "CORNER", "layer": None},
            {"name": "ind_face_uv_orig", "type": "FLOAT_VECTOR", "domain": "CORNER", "layer": None},
            {"name": "virtual_quad", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "vq_other_tri", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "vq_diagonal", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "trailing_tri", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "process_sequence", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "process_sequence_edge", "type": "INT", "domain": "EDGE", "layer": None},
            {"name": "process_sequence_max", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "xi", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "yi", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "island_index", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "flatness", "type": "FLOAT", "domain": "FACE", "layer": None},
            {"name": "normal_regularity", "type": "FLOAT", "domain": "FACE", "layer": None},
            {"name": "nom", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "nom_dir", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "nom_dir_edge", "type": "INT", "domain": "EDGE", "layer": None},
            {"name": "avg_length_x", "type": "FLOAT", "domain": "FACE", "layer": None},
            {"name": "avg_length_y", "type": "FLOAT", "domain": "FACE", "layer": None},
            ]
        
        for nam in attr_dict:
            attr = props.find_or_add_attribute(nam["name"], nam["type"], nam["domain"])
            nam["layer"] = props.get_attribute_layer(nam["name"], nam["type"], nam["domain"], bm)
        
        flatness = props.get_dict_layer("flatness", attr_dict)
        normal_regularity = props.get_dict_layer("normal_regularity", attr_dict)
        island_index = props.get_dict_layer("island_index", attr_dict)
        edge_u = props.get_dict_layer("edge_u", attr_dict)
        edge_l = props.get_dict_layer("edge_l", attr_dict)
        edge_d = props.get_dict_layer("edge_d", attr_dict)
        edge_r = props.get_dict_layer("edge_r", attr_dict)
        corner_ur = props.get_dict_layer("corner_ur", attr_dict)
        corner_ul = props.get_dict_layer("corner_ul", attr_dict)
        corner_ll = props.get_dict_layer("corner_ll", attr_dict)
        corner_lr = props.get_dict_layer("corner_lr", attr_dict)
        co_ur_x = props.get_dict_layer("co_ur_x", attr_dict)
        co_ul_x = props.get_dict_layer("co_ul_x", attr_dict)
        co_ll_x = props.get_dict_layer("co_ll_x", attr_dict)
        co_lr_x = props.get_dict_layer("co_lr_x", attr_dict)
        co_ur_y = props.get_dict_layer("co_ur_y", attr_dict)
        co_ul_y = props.get_dict_layer("co_ul_y", attr_dict)
        co_ll_y = props.get_dict_layer("co_ll_y", attr_dict)
        co_lr_y = props.get_dict_layer("co_lr_y", attr_dict)
        virtual_quad_attr = props.get_dict_layer("virtual_quad", attr_dict)
        diagonal_index = props.get_dict_layer("vq_diagonal", attr_dict)
        vq_other_tri = props.get_dict_layer("vq_other_tri", attr_dict)
        trailing_tri = props.get_dict_layer("trailing_tri", attr_dict)
        process_sequence_attr = props.get_dict_layer("process_sequence", attr_dict)
        process_sequence_edge = props.get_dict_layer("process_sequence_edge", attr_dict)
        nom = props.get_dict_layer("nom", attr_dict)
        nom_dir = props.get_dict_layer("nom_dir", attr_dict)
        nom_dir_edge = props.get_dict_layer("nom_dir_edge", attr_dict)
        process_sequence_max = props.get_dict_layer("process_sequence_max", attr_dict)
        loop_corner_index = props.get_dict_layer("loop_corner_index", attr_dict)
        ind_face_uv = props.get_dict_layer("ind_face_uv", attr_dict)
        ind_face_uv_orig = props.get_dict_layer("ind_face_uv_orig", attr_dict)
        xi = props.get_dict_layer("xi", attr_dict)
        yi = props.get_dict_layer("yi", attr_dict)
        avg_length_x = props.get_dict_layer("avg_length_x", attr_dict)
        avg_length_y = props.get_dict_layer("avg_length_y", attr_dict)
        
        co_attr_arr = [[co_ur_x, co_ur_y], [co_ul_x, co_ul_y], [co_ll_x, co_ll_y], [co_lr_x, co_lr_y]]
        # endregion

        quant_arr = [
            [mathutils.Vector((1,0,0)), 
                [
                    mathutils.Vector((0,0,-1)),
                    mathutils.Vector((0,-1,0)),
                    mathutils.Vector((0,0,1)),
                    mathutils.Vector((0,1,0))
                ],
                "Left"
             ],
             [mathutils.Vector((-1,0,0)),
                [
                    mathutils.Vector((0,0,-1)),
                    mathutils.Vector((0,1,0)),
                    mathutils.Vector((0,0,1)),
                    mathutils.Vector((0,-1,0))
                ],
                "Right"
             ],
             [mathutils.Vector((0,1,0)),
                [
                    mathutils.Vector((0,0,-1)),
                    mathutils.Vector((1,0,0)),
                    mathutils.Vector((0,0,1)),
                    mathutils.Vector((-1,0,0))
                ],
                "Towards"
             ],
             [mathutils.Vector((0,-1,0)),
                [
                    mathutils.Vector((0,0,-1)),
                    mathutils.Vector((-1,0,0)),
                    mathutils.Vector((0,0,1)),
                    mathutils.Vector((1,0,0))
                ],
                "Away"
             ],
             [mathutils.Vector((0,0,-1)),
                [
                    mathutils.Vector((0,-1,0)),
                    mathutils.Vector((1,0,0)),
                    mathutils.Vector((0,1,0)),
                    mathutils.Vector((-1,0,0))
                ],
                "Up"
             ],
             [mathutils.Vector((0,0,1)),
                [
                    mathutils.Vector((0,1,0)),
                    mathutils.Vector((1,0,0)),
                    mathutils.Vector((0,-1,0)),
                    mathutils.Vector((-1,0,0))
                ],
                "Down"
             ],
        ]
        # endregion

        # region functions
        def check_list_for_duplicate(the_list):
            seen = set()
            for the_item in the_list:
                if the_item in seen: return True
                seen.add(the_item)
            return False
        
        def quant_sort(e):
            quant_dot = e[0].dot(current_normal)
            return quant_dot

        def map_sort_by_coord(e):
            return e[0]

        def map_sort_by_edge_length(e):
            return e[2]

        def regularity_sort(e):
            flat = (1.0 - e[flatness]) * -self.reg_flat_fac
            norm_reg = e[normal_regularity] * -self.reg_norm_fac
            return flat + norm_reg
        
        def adj_vert_angle_sort(e):
            return e[1]
        
        def calc_edge_center(e):
            return e.verts[0].co.lerp(e.verts[1].co, 0.5)

        def dot_sort_u(e):
            evec = calc_edge_center(e) - current_face_center
            edir = evec.normalized()
            return edir.dot(current_dir_array[0])
        
        def dot_sort_l(e):
            evec = calc_edge_center(e) - current_face_center
            edir = evec.normalized()
            return edir.dot(current_dir_array[1])
        
        def dot_sort_d(e):
            evec = calc_edge_center(e) - current_face_center
            edir = evec.normalized()
            return edir.dot(current_dir_array[2])
        
        def dot_sort_r(e):
            evec = calc_edge_center(e) - current_face_center
            edir = evec.normalized()
            return edir.dot(current_dir_array[3])
        
        def get_other_dir_edges(edge_arr, corner_verts_arr, source_face, other_face, dir):
            out_array = [None, None, None, None]
            virtual_edges = [None, None, None, None]
            corner_verts = corner_verts_arr[0]
            adj_edge = None
            adjacency_verts = []
            if dir == 0:
                other_face[edge_d] = edge_arr[dir].index
                out_array[2] = edge_arr[dir]
                adj_edge = edge_arr[dir]
                adjacency_verts = [corner_verts[0], corner_verts[1]]
                virtual_edges[2] = adj_edge
            if dir == 1:
                other_face[edge_r] = edge_arr[dir].index
                out_array[3] = edge_arr[dir]
                adj_edge = edge_arr[dir]
                adjacency_verts = [corner_verts[1], corner_verts[2]]
                virtual_edges[3] = adj_edge
            if dir == 2:
                other_face[edge_u] = edge_arr[dir].index
                out_array[0] = edge_arr[dir]
                adj_edge = edge_arr[dir]
                adjacency_verts = [corner_verts[2], corner_verts[3]]
                virtual_edges[0] = adj_edge
            if dir == 3:
                other_face[edge_l] = edge_arr[dir].index
                out_array[1] = edge_arr[dir]
                adj_edge = edge_arr[dir]
                adjacency_verts = [corner_verts[3], corner_verts[0]]
                virtual_edges[1] = adj_edge
            sides = len(other_face.edges)

            if sides == 3:
                main_tri = other_face
                tri_trailing = False
                diagonal_edge = None
                far_tri = None
                adj_verts_arr = []
                for loop in main_tri.loops:
                    if loop.vert in adjacency_verts:
                        adj_vert_arr_app = [loop.vert]
                        adj_vert_angle = math.degrees(loop.calc_angle())
                        adj_vert_arr_app.append(adj_vert_angle)
                        adj_verts_arr.append(adj_vert_arr_app)
                adj_verts_arr.sort(key=adj_vert_angle_sort)
                for avf in adj_verts_arr[0][0].link_faces:
                    if avf != main_tri and avf != source_face:
                        main_tri_adjacent = False
                        for avfe in avf.edges:
                            if avfe not in edge_arr:
                                if main_tri in avfe.link_faces:
                                    if avfe.seam == False:
                                        main_tri_adjacent = True
                                        diagonal_edge = avfe
                                    if avfe.seam == True:
                                        tri_trailing = True
                        if main_tri_adjacent == True:
                            far_tri = avf
                            if len(far_tri.edges) > 3:
                                tri_trailing = True
                if far_tri == None:
                    tri_trailing = True
                if tri_trailing == False:
                    if far_tri[island_index] == main_tri[island_index]:
                        main_tri[virtual_quad_attr] = 1
                        far_tri[virtual_quad_attr] = 1
                        main_tri[trailing_tri] = 0
                        far_tri[trailing_tri] = 0
                        main_tri[vq_other_tri] = far_tri.index
                        far_tri[vq_other_tri] = main_tri.index
                        main_tri[diagonal_index] = diagonal_edge.index
                        avf[diagonal_index] = diagonal_edge.index
                    if far_tri[island_index] != main_tri[island_index]:
                        main_tri[trailing_tri] = 1
                        main_tri[virtual_quad_attr] = 0
                        tri_trailing = True
                if tri_trailing == True:
                    main_tri[trailing_tri] = 1
                    main_tri.select = True
                    main_tri[virtual_quad_attr] = 0
                    out_array = []

                if tri_trailing == False:
                    for mte in main_tri.edges:
                        if mte != adj_edge and mte != diagonal_edge:
                            if dir == 0:
                                if mte in corner_verts[0].link_edges:
                                    virtual_edges[3] = mte
                                    out_array[3] = mte
                                if mte in corner_verts[1].link_edges:
                                    virtual_edges[1] = mte
                                    out_array[1] = mte
                            if dir == 1:
                                if mte in corner_verts[1].link_edges:
                                    virtual_edges[0] = mte
                                    out_array[0] = mte
                                if mte in corner_verts[2].link_edges:
                                    virtual_edges[2] = mte
                                    out_array[2] = mte
                            if dir == 2:
                                if mte in corner_verts[2].link_edges:
                                    virtual_edges[1] = mte
                                    out_array[1] = mte
                                if mte in corner_verts[3].link_edges:
                                    virtual_edges[3] = mte
                                    out_array[3] = mte
                            if dir == 3:
                                if mte in corner_verts[3].link_edges:
                                    virtual_edges[2] = mte
                                    out_array[2] = mte
                                if mte in corner_verts[0].link_edges:
                                    virtual_edges[0] = mte
                                    out_array[0] = mte
                    for fte in far_tri.edges:
                        if fte != diagonal_edge:
                            if dir == 0:
                                if fte in corner_verts[0].link_edges:
                                    virtual_edges[3] = fte
                                    out_array[3] = fte
                                elif fte in corner_verts[1].link_edges:
                                    virtual_edges[1] = fte
                                    out_array[1] = fte
                                else:
                                    virtual_edges[0] = fte
                                    out_array[0] = fte
                            if dir == 1:
                                if fte in corner_verts[1].link_edges:
                                    virtual_edges[0] = fte
                                    out_array[0] = fte
                                elif fte in corner_verts[2].link_edges:
                                    virtual_edges[2] = fte
                                    out_array[2] = fte
                                else:
                                    virtual_edges[1] = fte
                                    out_array[1] = mte
                            if dir == 2:
                                if fte in corner_verts[2].link_edges:
                                    virtual_edges[1] = fte
                                    out_array[1] = fte
                                elif fte in corner_verts[3].link_edges:
                                    virtual_edges[3] = fte
                                    out_array[3] = fte
                                else:
                                    virtual_edges[2] = fte
                                    out_array[2] = fte
                            if dir == 3:
                                if fte in corner_verts[3].link_edges:
                                    virtual_edges[2] = fte
                                    out_array[2] = fte
                                elif fte in corner_verts[0].link_edges:
                                    virtual_edges[0] = fte
                                    out_array[0] = fte
                                else:
                                    virtual_edges[3] = fte
                                    out_array[3] = fte
                    main_tri[edge_u] = virtual_edges[0].index
                    far_tri[edge_u] = virtual_edges[0].index
                    main_tri[edge_l] = virtual_edges[1].index
                    far_tri[edge_l] = virtual_edges[1].index
                    main_tri[edge_d] = virtual_edges[2].index
                    far_tri[edge_d] = virtual_edges[2].index
                    main_tri[edge_r] = virtual_edges[3].index
                    far_tri[edge_r] = virtual_edges[3].index

            if sides == 4:
                other_face[virtual_quad_attr] = 0
                for ofe in other_face.edges:
                    if ofe not in source_face.edges:
                        newdir = None
                        if dir == 0:
                            if ofe in corner_verts[0].link_edges:
                                out_array[3] = ofe
                                newdir = 3
                            if ofe in corner_verts[1].link_edges:
                                out_array[1] = ofe
                                newdir = 1
                            if ofe not in corner_verts[0].link_edges and ofe not in corner_verts[1].link_edges:
                                out_array[0] = ofe
                                newdir = 0
                        if dir == 1:
                            if ofe in corner_verts[1].link_edges:
                                out_array[0] = ofe
                                newdir = 0
                            if ofe in corner_verts[2].link_edges:
                                out_array[2] = ofe
                                newdir = 2
                            if ofe not in corner_verts[1].link_edges and ofe not in corner_verts[2].link_edges:
                                out_array[1] = ofe
                                newdir = 1
                        if dir == 2:
                            if ofe in corner_verts[2].link_edges:
                                out_array[1] = ofe
                                newdir = 1
                            if ofe in corner_verts[3].link_edges:
                                out_array[3] = ofe
                                newdir = 3
                            if ofe not in corner_verts[2].link_edges and ofe not in corner_verts[3].link_edges:
                                out_array[2] = ofe
                                newdir = 2
                        if dir == 3:
                            if ofe in corner_verts[3].link_edges:
                                out_array[2] = ofe
                                newdir = 2
                            if ofe in corner_verts[0].link_edges:
                                out_array[0] = ofe
                                newdir = 0
                            if ofe not in corner_verts[3].link_edges and ofe not in corner_verts[0].link_edges:
                                out_array[3] = ofe
                                newdir = 3
                        # print(f"New direction {newdir}: edge {ofe.index}")
                    if ofe in source_face.edges:
                        if dir == 0:
                            out_array[2] = ofe
                            newdir = 2
                        elif dir == 1:
                            out_array[3] = ofe
                            newdir = 3
                        elif dir == 2:
                            out_array[0] = ofe
                            newdir = 0
                        elif dir == 3:
                            out_array[1] = ofe
                            newdir = 1
                        # print(f"New direction {newdir}: edge {ofe.index}. Adjacent edge.")
                # for edgi, edge in enumerate(out_array):
                    # print(f"{edgi} - Direction:{dir}: {edge}")
                other_face[edge_u] = out_array[0].index
                other_face[edge_l] = out_array[1].index
                other_face[edge_d] = out_array[2].index
                other_face[edge_r] = out_array[3].index

            return out_array
        
        def check_winding(face, corner_verts_arr):
            poss = 0
            negs = 0
            nones = 0
            out_str = ""
            last_i = corner_verts_arr[1][-1]
            for vi in corner_verts_arr[1]:
                if vi == None:
                    nones += 1
                else:
                    result = vi - last_i
                    if result > 0:
                        poss += 1
                    if result < 0:
                        negs += 1
                    last_i = vi
            if poss > negs:
                out_str = "CCW"
            if poss < negs:
                out_str = "CW"
            if nones > 0:
                out_str = "Winding check failed"
            return out_str
        
        def tri_angle_edge_sort(e):
            angles = 0.0
            for vert in e[0].verts:
                for loop in vert.link_loops:
                    if loop in e[1].loops:
                        angles += math.degrees(loop.calc_angle())
            return angles
        
        def get_adjacent_triangle(tri):
            tri_edges = []
            out_tri = [None, None]
            for te in tri.edges:
                for tef in te.link_faces:
                    if tef != tri and len(tef.edges) == 3:
                        tri_edge = []
                        tri_edge.append(te)
                        tri_edge.append(tri)
                        tri_edge.append(tef)
                        tri_edges.append(tri_edge)
            if len(tri_edges) > 0:
                tri_edges.sort(key=tri_angle_edge_sort)
                out_tri[0] = tri_edges[0][2]
                out_tri[1] = tri_edges[0][0]
            else:
                pass
            return out_tri
        
        def get_any_adjacent_quad(source_face):
            out_quad = None
            for qe in source_face.edges:
                for qef in qe.link_faces:
                    if qef != source_face and len(qef.edges) == 4:
                        out_quad = qef
                        return out_quad
            else:
                pass
            return out_quad
        
        def get_init_dir_edges(arr):
            out_array = []
            u_arr = arr.copy()
            l_arr = arr.copy()
            d_arr = arr.copy()
            r_arr = arr.copy()
            u_arr.sort(key=dot_sort_u)
            l_arr.sort(key=dot_sort_l)
            d_arr.sort(key=dot_sort_d)
            r_arr.sort(key=dot_sort_r)
            out_array.append(u_arr[0])
            out_array.append(l_arr[0])
            out_array.append(d_arr[0])
            out_array.append(r_arr[0])
            return out_array
        
        def get_corner_verts(edge_arr, v_arr):
            out_arrays = []
            out_array = [None, None, None, None]
            out_array_is = [None, None, None, None]
            v_iter = 0
            for v in v_arr:
                ul = 0
                ll = 0
                lr = 0
                ur = 0
                if v in edge_arr[0].verts:
                    ur += 1
                    ul += 1
                if v in edge_arr[1].verts:
                    ul += 1
                    ll += 1
                if v in edge_arr[2].verts:
                    ll += 1
                    lr += 1
                if v in edge_arr[3].verts:
                    lr += 1
                    ur += 1
                if ur == 2:
                    out_array[0] = v
                    out_array_is[0] = v_iter
                    # print(f"Upper right corner: {v.index} (seq. index: {v_iter})")
                if ul == 2:
                    out_array[1] = v
                    out_array_is[1] = v_iter
                    # print(f"Upper left corner: {v.index} (seq. index: {v_iter})")
                if ll == 2:
                    out_array[2] = v
                    out_array_is[2] = v_iter
                    # print(f"Lower left corner: {v.index} (seq. index: {v_iter})")
                if lr == 2:
                    out_array[3] = v
                    out_array_is[3] = v_iter
                    # print(f"Lower right corner: {v.index} (seq. index: {v_iter})")
                v_iter += 1
            out_arrays.append(out_array)
            out_arrays.append(out_array_is)
            return out_arrays

        def viz_quad(qe, qc, qf, vqf, virt, wnd):
            fnamlen = 2 + len(str(qf.index))
            if virt == True:
                fnamlen = 3 + len(str(vqf[0].index)) + len(str(vqf[1].index))
            max_l = max(len(str(qc[1].index)), len(str(qc[2].index)), len(str(qe[1].index)))
            max_c = 2 + max(fnamlen, len(str(qe[0].index)), len(str(qe[2].index)))
            max_r = max(len(str(qc[0].index)), len(str(qe[3].index)), len(str(qc[3].index)))
            max_vline = max(max_l, max_r)
            ulc = str(qc[1].index).rjust(max_vline)
            uline = str(qe[0].index).center(max_c,"-")
            urc = str(qc[0].index).ljust(max_vline)
            vline = "|".center(max_vline)
            vspace = " ".center(max_c)
            vu = "^".center(max_vline)
            vd = "v".center(max_vline)
            larrow = "<"
            rarrow = ">"
            lline = str(qe[1].index).rjust(max_vline)
            fnam = str(qf.index).center(max_c)
            if virt == True:
                fnama = str(vqf[0].index) + "/" + str(vqf[1].index)
                fnam = fnama.center(max_c)
            rline = str(qe[3].index).ljust(max_vline)
            llc = str(qc[2].index).rjust(max_vline)
            dline = str(qe[2].index).center(max_c, "-")
            lrc = str(qc[3].index).ljust(max_vline)
            quad_msg = "\n{0}{1}{2}\n{3}{4}{3}\n{3}{4}{3}\n{3}{4}{3}\n{5}{6}{7}\n{3}{4}{3}\n{3}{4}{3}\n{3}{4}{3}\n{8}{9}{10}\n".format(ulc,uline,urc,vline,vspace,lline,fnam,rline, llc, dline, lrc)
            if wnd == "CCW":
                quad_msg = "\n{0}{1}{13}{2}\n{12}{4}{3}\n{3}{4}{3}\n{3}{4}{3}\n{5}{6}{7}\n{3}{4}{3}\n{3}{4}{3}\n{3}{4}{11}\n{8}{14}{9}{10}\n".format(ulc,uline,urc,vline,vspace,lline,fnam,rline, llc, dline, lrc, vu, vd, larrow, rarrow)
            if wnd == "CW":
                quad_msg = "\n{0}{14}{1}{2}\n{3}{4}{11}\n{3}{4}{3}\n{3}{4}{3}\n{5}{6}{7}\n{3}{4}{3}\n{3}{4}{3}\n{12}{4}{3}\n{8}{9}{13}{10}\n".format(ulc,uline,urc,vline,vspace,lline,fnam,rline, llc, dline, lrc, vu, vd, larrow, rarrow)
            return quad_msg

        def viz_add_quad(qe, oqe, qc, oqc, qf, oqf, vqf, ovqf, virt, ovirt, add_dir):
            fnamlen = 2 + len(str(qf.index))
            fnamolen = 2 + len(str(oqf.index))
            if virt == True:
                fnamlen = 3 + len(str(vqf[0].index)) + len(str(vqf[1].index))
            if ovirt == True:
                fnamolen = 3 + len(str(ovqf[0].index)) + len(str(ovqf[1].index))
            max_l = max(len(str(qc[1].index)), len(str(oqc[1].index)), len(str(qc[2].index)), len(str(oqc[2].index)), len(str(qe[1].index)), len(str(oqe[1].index)))
            max_c = 6 + max(fnamlen, fnamolen, len(str(qe[0].index)), len(str(oqe[0].index)), len(str(qe[2].index)), len(str(oqe[2].index)))
            max_r = max(len(str(qc[0].index)), len(str(oqc[0].index)), len(str(qe[3].index)), len(str(oqe[3].index)), len(str(qc[3].index)), len(str(oqc[3].index)))
            max_vline = max(max_l, max_r)
            vspace_sum = max_c + (2 * max_vline)
            ulc = str(qc[1].index).rjust(max_vline)
            oulc = str(oqc[1].index).rjust(max_vline)
            uline = str(qe[0].index).center(max_c,"=")
            ouline = str(oqe[0].index).center(max_c,"-")
            urc = str(qc[0].index).ljust(max_vline)
            ourc = str(oqc[0].index).ljust(max_vline)
            vline = "||".center(max_vline)
            ovline = "|".center(max_vline)
            vspace = " ".center(max_c)
            ovspace = " ".center(vspace_sum)
            ohspace = " ".center(max_c + 2)
            ovu = "^".center(vspace_sum)
            ovd = "v".center(vspace_sum)
            olarrow = "<".center(max_c + 2)
            orarrow = ">".center(max_c + 2)
            lline = str(qe[1].index).rjust(max_vline)
            olline = str(oqe[1].index).rjust(max_vline)
            fnam = str(qf.index).center(max_c)
            if virt == True:
                fnama = str(vqf[0].index) + "/" + str(vqf[1].index)
                fnam = fnama.center(max_c)
            ofnam = str(oqf.index).center(max_c)
            if ovirt == True:
                ofnama = str(ovqf[0].index) + "/" + str(ovqf[1].index)
                ofnam = ofnama.center(max_c)
            rline = str(qe[3].index).ljust(max_vline)
            orline = str(oqe[3].index).ljust(max_vline)
            llc = str(qc[2].index).rjust(max_vline)
            ollc = str(oqc[2].index).rjust(max_vline)
            dline = str(qe[2].index).center(max_c, "=")
            odline = str(oqe[2].index).center(max_c, "-")
            lrc = str(qc[3].index).ljust(max_vline)
            olrc = str(oqc[3].index).ljust(max_vline)

            quad_msg_aa = "{0}{1}{2}".format(ulc,uline,urc)
            quad_msg_avsp = "{0}{1}{0}".format(vline,vspace)
            quad_msg_ab = "{0}{1}{2}".format(lline,fnam,rline)
            quad_msg_ac = "{0}{1}{2}".format(llc,dline,lrc)
            quad_msg_ba = "{0}{1}{2}".format(oulc,ouline,ourc)
            quad_msg_bvsp = "{0}{1}{0}".format(ovline,vspace)
            quad_msg_bb = "{0}{1}{2}".format(olline,ofnam,orline)
            quad_msg_bc = "{0}{1}{2}".format(ollc,odline,olrc)

            quad_msg = ""
            if add_dir == 0:
                quad_msg = "\n{0}\n{3}\n{3}\n{3}\n{1}\n{3}\n{3}\n{3}\n{2}\n{8}\n{9}\n{8}\n{4}\n{7}\n{7}\n{7}\n{5}\n{7}\n{7}\n{7}\n{6}\n".format(quad_msg_ba, quad_msg_bb, quad_msg_bc, quad_msg_bvsp, quad_msg_aa, quad_msg_ab, quad_msg_ac, quad_msg_avsp, ovspace, ovu)
            if add_dir == 1:
                quad_msg = "\n{0}{8}{4}\n{3}{8}{7}\n{3}{8}{7}\n{3}{8}{7}\n{1}{9}{5}\n{3}{8}{7}\n{3}{8}{7}\n{3}{8}{7}\n{2}{8}{6}\n".format(quad_msg_ba, quad_msg_bb, quad_msg_bc, quad_msg_bvsp, quad_msg_aa, quad_msg_ab, quad_msg_ac, quad_msg_avsp, ohspace, olarrow)
            if add_dir == 2:
                quad_msg = "\n{0}\n{3}\n{3}\n{3}\n{1}\n{3}\n{3}\n{3}\n{2}\n{8}\n{9}\n{8}\n{4}\n{7}\n{7}\n{7}\n{5}\n{7}\n{7}\n{7}\n{6}\n".format(quad_msg_aa, quad_msg_ab, quad_msg_ac, quad_msg_avsp, quad_msg_ba, quad_msg_bb, quad_msg_bc, quad_msg_bvsp, ovspace, ovd)
            if add_dir == 3:
                quad_msg = "\n{0}{8}{4}\n{3}{8}{7}\n{3}{8}{7}\n{3}{8}{7}\n{1}{9}{5}\n{3}{8}{7}\n{3}{8}{7}\n{3}{8}{7}\n{2}{8}{6}\n".format(quad_msg_aa, quad_msg_ab, quad_msg_ac, quad_msg_avsp, quad_msg_ba, quad_msg_bb, quad_msg_bc, quad_msg_bvsp, ohspace, orarrow)
            return quad_msg
        # endregion

        previous_island_end_x = 0.0
        previous_island_end_y = 0.0

        for ii, island in enumerate(islands):
            for iif in island:
                iif[island_index] = ii

        for ii, island in enumerate(islands):
            bm.verts.ensure_lookup_table()
            bm.edges.ensure_lookup_table()
            bm.faces.ensure_lookup_table()
            edge_length_map_x = [[],[]]
            edge_length_map_y = [[],[]]
            island_max_x = 0.0
            island_max_y = 0.0
            island_min_x = 999.0
            island_min_y = 999.0
            process_sequence = 0
            faces_remain = []
            faces_done = []
            verts_done = []
            faces_selected = []
            island_loops = []
            init_face = None
            init_virtual_face = []
            avg_norm = mathutils.Vector((0,0,0))
            if props.verbose: print(f"Island {ii}:")

            # region prep init search
            for face in island:
                avg_norm += face.normal
                avg_angle = 0.0
                angle_count = 0
                for fe in face.edges:
                    angle_count += 1
                    if fe.is_boundary == True or len(face.edges) < 2 or len(face.edges) > 2:
                        avg_angle += 90
                    else:
                        avg_angle += math.degrees(fe.calc_face_angle())
                for fl in face.loops:
                    island_loops.append(fl.index)
                avg_angle /= angle_count
                new_max_angle = max(avg_angle, max_avg_angle)
                new_min_angle = min(avg_angle, min_avg_angle)
                max_avg_angle = new_max_angle
                min_avg_angle = new_min_angle
                face[flatness] = avg_angle
                faces_remain.append(face)
                if face.select == True:
                    faces_selected.append(face)
            
            current_avg_normal = avg_norm.normalized()

            if self.quant_avg_norm == True:
                current_normal = current_avg_normal
                quant_arr.sort(key=quant_sort)
                current_avg_normal = quant_arr[0][0] * -1.0
                if props.verbose == True: print(f"Average island normals quantized to {quant_arr[0][2]}")

            bpy.context.active_object.data.update()

            for face in island:
                new_flatness = props.remap_val(face[flatness], min_avg_angle, max_avg_angle, 0.0, 1.0)
                face[flatness] = new_flatness
                n_dot_avg = face.normal.dot(current_avg_normal)
                face[normal_regularity] = props.remap_val(n_dot_avg, -1.0, 1.0, 0.0, 1.0)
            
            bpy.context.active_object.data.update()
            # endregion

            # region find init face(s)
            initial_face_found = False
            
            while initial_face_found == False and len(faces_remain) > 0:
                is_selected_face = False
                init_face_is_virtual = False
                init_virtual_face_diagonal = None
                if self.initial_quad == "B":
                    if len(faces_selected) > 0:
                        init_face = faces_selected[0]
                        init_virtual_face = [init_face]
                        if len(init_face.edges) < 4:
                            if len(faces_selected) > 1:
                                if len(faces_selected[1].edges) == 3:
                                    init_face_is_virtual = True
                                    init_virtual_face = [init_face, faces_selected[1]]
                                    is_selected_face = True
                                if len(faces_selected[1].edges) == 4:
                                    init_face_is_virtual = False
                                    init_face = faces_selected[1]
                                    init_virtual_face = [faces_selected[1]]
                                    is_selected_face = True
                                else:
                                    if props.verbose: print("Too inconsistent mesh! Aborting!")
                                    return  {"CANCELLED"}
                            else:
                                other_tri = get_adjacent_triangle(init_face)
                                if other_tri[0] != None:
                                    init_face_is_virtual = True
                                    init_virtual_face = [init_face, other_tri[0]]
                                    is_selected_face = other_tri[0].select
                                    init_virtual_face_diagonal = other_tri[1]
                                else:
                                    other_quad = get_any_adjacent_quad(init_face)
                                    if other_quad != None:
                                        init_face = other_quad
                                        init_virtual_face = [init_face]
                                        is_selected_face = init_face.select
                                        init_face_is_virtual = False
                                    else:
                                        if props.verbose: print("Too inconsistent mesh! Aborting!")
                                        return  {"CANCELLED"}
                    if len(faces_selected) == 0:
                        faces_remain.sort(key=regularity_sort)
                        init_face = faces_remain[0]
                        init_virtual_face = [init_face]
                        if len(init_face.edges) == 3:
                            other_tri = get_adjacent_triangle(init_face)
                            if other_tri[0] != None:
                                init_face_is_virtual = True
                                init_virtual_face = [init_face, other_tri[0]]
                                init_virtual_face_diagonal = other_tri[1]
                            else:
                                other_quad = get_any_adjacent_quad(init_face)
                                if other_quad != None:
                                    init_face = other_quad
                                    init_virtual_face = [init_face]
                                else:
                                    if props.verbose: print("Too inconsistent mesh! Aborting!")
                                    return  {"CANCELLED"}
                else:
                    faces_remain.sort(key=regularity_sort)
                    init_face = faces_remain[0]
                    init_virtual_face = [init_face]
                    if len(init_face.edges) == 3:
                        other_tri = get_adjacent_triangle(init_face)
                        if other_tri[0] != None:
                            init_face_is_virtual = True
                            init_virtual_face = [init_face, other_tri[0]]
                            init_virtual_face_diagonal = other_tri[1]
                        else:
                            other_quad = get_any_adjacent_quad(init_face)
                            if other_quad != None:
                                init_face = other_quad
                                init_virtual_face = [init_face]
                            else:
                                if props.verbose: print("Too inconsistent mesh! Aborting!")
                                return  {"CANCELLED"}

                current_normal = init_face.normal

                if init_face_is_virtual == True:
                    virtual_normal = mathutils.Vector((0,0,0))
                    for vf in init_virtual_face:
                        virtual_normal += vf.normal
                    current_normal = virtual_normal.normalized()

                quant_arr.sort(key=quant_sort)
                current_quant_index = 0
                current_dir_array = quant_arr[current_quant_index][1]
                init_face_dir_arr = [None, None, None, None]

                if props.verbose == True:
                    msg_start = "Most regular face for island "
                    if is_selected_face == True:
                        msg_start = "Selected face for island "
                    msg = "{:}{:}: {:} with a flatness of {:.3f} and a normal regularity of {:.3f}.\nFacing quantized to {:}\n".format(msg_start, ii, init_face.index, init_face[flatness], init_face[normal_regularity], quant_arr[0][2])
                    if init_face_is_virtual == True:
                        msg_start = "Most regular virtual face for island "
                        if is_selected_face == True:
                            msg_start = "Selected virtual face for island "
                        msg = "{:}{:}: {:}/{:} with a flatness of {:.3f} and a normal regularity of {:.3f}.\nFacing quantized to {:}\n".format(msg_start, ii, init_virtual_face[0].index, init_virtual_face[1].index, init_face[flatness], init_face[normal_regularity], quant_arr[0][2])
                    print(msg)
                
                if init_face_is_virtual == True:
                    vif_edges = []
                    for vif in init_virtual_face:
                        vif.select = True
                    for avife in init_virtual_face[0].edges:
                        if avife != init_virtual_face_diagonal:
                            vif_edges.append(avife)
                    for bvife in init_virtual_face[1].edges:
                        if bvife != init_virtual_face_diagonal:
                            vif_edges.append(bvife)
                    current_face_center = calc_edge_center(init_virtual_face_diagonal)
                    init_face_dir_arr = get_init_dir_edges(vif_edges)
                    while check_list_for_duplicate(init_face_dir_arr) == True and current_quant_index < 4:
                        current_quant_index += 1
                        try:
                            current_dir_array = quant_arr[current_quant_index][1]
                        except:
                            if props.verbose: print(f"Failed to get directions for virtual face {init_virtual_face[0].index}/{init_virtual_face[1].index}! Marking faces as trailing triangles...")
                            else:
                                pass
                    
                    if current_quant_index > 3:
                        init_face[process_sequence_attr] = process_sequence
                        if init_face in faces_remain:
                            faces_remain.remove(init_face)
                        if init_face not in faces_done:
                            faces_done.append(init_face)
                        init_face[trailing_tri] = 1
                        init_face.select = True
                        init_face[virtual_quad_attr] = 0
                        for vif in init_virtual_face:
                            vif[process_sequence_attr] = process_sequence
                            if vif in faces_remain:
                                faces_remain.remove(vif)
                            if vif not in faces_done:
                                faces_done.append(vif)
                            vif[trailing_tri] = 1
                            vif.select = True
                            vif[virtual_quad_attr] = 0
                        init_virtual_face = []
                    
                    if current_quant_index <= 3:
                        vif_verts = []
                        for vif in init_virtual_face:
                            for vifv in vif.verts:
                                if vifv not in vif_verts:
                                    vif_verts.append(vifv)
                            vif[edge_u] = init_face_dir_arr[0].index
                            vif[edge_l] = init_face_dir_arr[1].index
                            vif[edge_d] = init_face_dir_arr[2].index
                            vif[edge_r] = init_face_dir_arr[3].index
                            vif[process_sequence_attr] = process_sequence
                        these_corner_verts = get_corner_verts(init_face_dir_arr, vif_verts)
                        initial_face_found = True
                        for tcv in these_corner_verts[0]:
                            if tcv == None:
                                initial_face_found = False
                                if props.verbose == True: print(f"Failed to calculate corners for virtual face {init_virtual_face[0].index}/{init_virtual_face[1].index}. Will pick another initial face.")
                        if initial_face_found == False:
                            for vif in init_virtual_face:
                                if vif in faces_remain:
                                    faces_remain.remove(vif)
                                if vif not in faces_done:
                                    faces_done.append(vif)
                                vif[trailing_tri] = 1
                                vif.select = True
                                vif[virtual_quad_attr] = 0
                            init_virtual_face = []
                else:
                    if_verts = [ifv for ifv in init_face.verts]
                    init_face.select = True
                    ife_edges = []
                    for ife in init_face.edges:
                        ife_edges.append(ife)
                    current_face_center = init_face.calc_center_median()
                    init_face_dir_arr = get_init_dir_edges(ife_edges)
                    init_face[edge_u] = init_face_dir_arr[0].index
                    init_face[edge_l] = init_face_dir_arr[1].index
                    init_face[edge_d] = init_face_dir_arr[2].index
                    init_face[edge_r] = init_face_dir_arr[3].index
                    init_face[process_sequence_attr] = process_sequence
                    these_corner_verts = get_corner_verts(init_face_dir_arr, if_verts)
                    initial_face_found = True
                    for tcv in these_corner_verts[0]:
                        if tcv == None:
                            initial_face_found = False
                    if initial_face_found == False:
                        if init_face in faces_remain:
                            faces_remain.remove(init_face)
                        if init_face not in faces_done:
                            faces_done.append(init_face)
                        init_face[trailing_tri] = 1
                        init_face.select = True
                        init_face[virtual_quad_attr] = 0
                        init_virtual_face = []
            
            bpy.context.active_object.data.update()

            next_round = []
            init_face_corner_verts = []
            init_virtual_loops = []
            # endregion

            # region proc init face(s)
            if initial_face_found == True:
                vif_verts = []
                for vif in init_virtual_face:
                    vif[xi] = 0
                    vif[yi] = 0
                    if vif in faces_remain:
                        faces_remain.remove(vif)
                    for vifv in vif.verts:
                        if vifv not in vif_verts:
                            vif_verts.append(vifv)
                    for vifl in vif.loops:
                        if vifl not in init_virtual_loops:
                            init_virtual_loops.append(vifl)
                init_face_corner_verts = get_corner_verts(init_face_dir_arr, vif_verts)
                if props.verbose == True:
                    print(f"There are {len(init_virtual_face)} faces in the initial virtual face. Their indices:")
                    for ivf in init_virtual_face:
                        print(ivf.index)
                for vif in init_virtual_face:
                    try:
                        vif[corner_ur] = init_face_corner_verts[0][0].index
                        vif[corner_ul] = init_face_corner_verts[0][1].index
                        vif[corner_ll] = init_face_corner_verts[0][2].index
                        vif[corner_lr] = init_face_corner_verts[0][3].index
                    except:
                        if props.verbose == True: print(f"Failed to calculate corner vertices for face {vif.index}. This should not happen...")
                        if vif in init_virtual_face:
                            init_virtual_face.remove(vif)
                        pass
                    vif[xi] = 0
                    vif[yi] = 0

                avg_edge_length_x = (init_face_dir_arr[0].calc_length() + init_face_dir_arr[2].calc_length()) / 2
                avg_edge_length_y = (init_face_dir_arr[1].calc_length() + init_face_dir_arr[3].calc_length()) / 2
                for vif in init_virtual_face:
                    vif[avg_length_x] = avg_edge_length_x
                    vif[avg_length_y] = avg_edge_length_y
                
                # region init uv edit
                if self.pre_calc_edge_lengths == False:
                    avg_edge_length_x = (init_face_dir_arr[0].calc_length() + init_face_dir_arr[2].calc_length()) / 4
                    avg_edge_length_y = (init_face_dir_arr[1].calc_length() + init_face_dir_arr[3].calc_length()) / 4
                    
                    offset_co = mathutils.Vector((previous_island_end_x, previous_island_end_y))

                    init_face_corner_vert_cos = [mathutils.Vector((-avg_edge_length_x, avg_edge_length_y)) + offset_co, mathutils.Vector((-avg_edge_length_x, -avg_edge_length_y)) + offset_co, mathutils.Vector((avg_edge_length_x, -avg_edge_length_y)) + offset_co, mathutils.Vector((avg_edge_length_x, avg_edge_length_y)) + offset_co]
                    for ifvco, ifva in zip(init_face_corner_vert_cos, co_attr_arr):
                        for vif in init_virtual_face:
                            vif[ifva[0]] = ifvco.x
                            vif[ifva[1]] = ifvco.y

                    for iv, ivco in zip(init_face_corner_verts[0], init_face_corner_vert_cos):
                        for ivl in iv.link_loops:
                            if ivl.index in island_loops:
                                should_edit = True
                                if self.only_move_loops_in_face == True:
                                    should_edit = ivl in init_virtual_loops
                                if should_edit == True:
                                    ivl[uv_layer].uv = ivco
                                    new_island_max_x = max(island_max_x, ivco.x)
                                    island_max_x = new_island_max_x
                                    new_island_min_x = min(island_min_x, ivco.x)
                                    island_min_x = new_island_min_x
                                    new_island_max_y = max(island_max_y, ivco.y)
                                    island_max_y = new_island_max_y
                                    new_island_min_y = min(island_min_y, ivco.y)
                                    island_min_y = new_island_min_y
                                    ivl[uv_layer].pin_uv = True
                        if iv not in verts_done:
                            verts_done.append(iv)
                # endregion
                
                # region init pre-calc
                if self.pre_calc_edge_lengths == True:
                    this_edge_length_map_x_item = [init_face, init_face_is_virtual, init_virtual_face, [init_face_corner_verts[0][1],init_face_corner_verts[0][2]], [init_face_corner_verts[0][0],init_face_corner_verts[0][3]]]
                    this_edge_length_map_y_item = [init_face, init_face_is_virtual, init_virtual_face, [init_face_corner_verts[0][2],init_face_corner_verts[0][3]], [init_face_corner_verts[0][0],init_face_corner_verts[0][1]]]
                    map_x_index = 0
                    map_y_index = 0
                    try:
                        map_x_index = edge_length_map_x[0].index(0)
                    except:
                        map_x_index = len(edge_length_map_x[0])
                        edge_length_map_x[0].append(0)
                        edge_length_map_x[1].append([0, [], 0.0, 0])
                    edge_length_map_x[1][map_x_index][1].append(this_edge_length_map_x_item)
                    edge_length_map_x[1][map_x_index][2] += avg_edge_length_x
                    edge_length_map_x[1][map_x_index][3] += 1

                    try:
                        map_y_index = edge_length_map_y[0].index(0)
                    except:
                        map_y_index = len(edge_length_map_y[0])
                        edge_length_map_y[0].append(0)
                        edge_length_map_y[1].append([0, [], 0.0, 0])
                    edge_length_map_y[1][map_y_index][1].append(this_edge_length_map_y_item)
                    edge_length_map_y[1][map_y_index][2] += avg_edge_length_y
                    edge_length_map_y[1][map_y_index][3] += 1
                # endregion
                
                # region init ind uv attr
                for ci, corner in enumerate(init_face_corner_verts[0]):
                    for corner_loop in corner.link_loops:
                        if corner_loop in init_virtual_loops:
                            corner_loop[loop_corner_index] = ci
                            if ci == 0:
                                corner_loop[ind_face_uv][0] = 1.0
                                corner_loop[ind_face_uv][1] = 1.0
                            if ci == 1:
                                corner_loop[ind_face_uv][0] = 0.0
                                corner_loop[ind_face_uv][1] = 1.0
                            if ci == 2:
                                corner_loop[ind_face_uv][0] = 0.0
                                corner_loop[ind_face_uv][1] = 0.0
                            if ci == 3:
                                corner_loop[ind_face_uv][0] = 1.0
                                corner_loop[ind_face_uv][1] = 0.0
                # endregion
                
                # region init ind uv orig
                for ci, corner in zip(init_face_corner_verts[0], init_face_corner_verts[0]):
                    for corner_loop in corner.link_loops:
                        if corner_loop in init_virtual_loops:
                            if ci == 0:
                                corner_loop[ind_face_uv_orig][0] = 1.0
                                corner_loop[ind_face_uv_orig][1] = 1.0
                            if ci == 1:
                                corner_loop[ind_face_uv_orig][0] = 0.0
                                corner_loop[ind_face_uv_orig][1] = 1.0
                            if ci == 2:
                                corner_loop[ind_face_uv_orig][0] = 0.0
                                corner_loop[ind_face_uv_orig][1] = 0.0
                            if ci == 3:
                                corner_loop[ind_face_uv_orig][0] = 1.0
                                corner_loop[ind_face_uv_orig][1] = 0.0
                # endregion
                
                # region init face info
                init_face_winding = check_winding(init_face, init_face_corner_verts)
                if props.verbose == True:
                    if init_face_is_virtual == False:
                        print(f"Processing initial face of island {ii}: {init_face.index}. Not virtual.")
                    if init_face_is_virtual == True:
                        print(f"Processing initial face of island {ii}: {init_virtual_face[0].index}/{init_virtual_face[1].index}. Virtual quad.")
                    try:
                        print(viz_quad(init_face_dir_arr, init_face_corner_verts[0], init_face, init_virtual_face, init_face_is_virtual, init_face_winding))
                    except:
                        print(f"Failed to visualize initial face(s)! Operation will probably fail...")
                # endregion

                if init_face not in faces_done:
                    faces_done.append(init_face)

                # region fill round 0
                for di, ife in enumerate(init_face_dir_arr):
                    dirs_to_do = [0,1,2,3]
                    if self.unfold_mode == "B":
                        dirs_to_do = [1,3]
                    elif self.unfold_mode == "C":
                        dirs_to_do = [0,2]
                    if di in dirs_to_do:
                        for ifef in ife.link_faces:
                            if ifef in faces_remain and ifef not in next_round and ifef != init_face and ifef not in init_virtual_face:
                                adj_face = init_face
                                other_virtual_face = [ifef]
                                other_verts = [ov for ov in ifef.verts]
                                if init_face_is_virtual == True:
                                    for vif in init_virtual_face:
                                        if ife in vif.edges:
                                            adj_face = vif
                                other_dir_edges = get_other_dir_edges(init_face_dir_arr, init_face_corner_verts, adj_face, ifef, di)
                                is_other_face_virtual = ifef[virtual_quad_attr]
                                if is_other_face_virtual == True:
                                    other_virtual_face.append(bm.faces[ifef[vq_other_tri]])
                                    for otv in other_virtual_face[1].verts:
                                        if otv not in other_verts:
                                            other_verts.append(otv)
                                ifef[nom_dir] = di
                                ifef[nom] = adj_face.index

                                if ifef[trailing_tri] == 1:
                                    if props.verbose: print(f"Initial face of island {ii}: skipped face index {ifef.index} in direction {di} because it is a trailing triangle.")
                                    ifef[process_sequence_attr] = process_sequence
                                    ife[process_sequence_edge] = process_sequence
                                    ife[nom_dir_edge] = di
                                    for ifefl in ifef.loops:
                                        ifefl[uv_layer].pin_uv = False
                                    process_sequence += 1
                                    faces_remain.remove(ifef)
                                    ifef.select = True
                                if ifef[trailing_tri] == 0:
                                    other_corner_verts = get_corner_verts(other_dir_edges, other_verts)
                                    skip_append = False

                                    if props.verbose:
                                        print(f"Initial face of island {ii}: Added face index {ifef.index} in direction {di} to next round. Virtual face: {init_face_is_virtual}.")
                                        try:
                                            print(viz_add_quad(init_face_dir_arr, other_dir_edges, init_face_corner_verts[0], other_corner_verts[0], init_face, ifef, init_virtual_face, other_virtual_face, init_face_is_virtual, is_other_face_virtual, di))
                                        except:
                                            print(f"Face index {ifef.index} failed to generate valid corners! Face skipped.")
                                            skip_append = True
                                    
                                    if skip_append == True:
                                        ifef[trailing_tri] = 1
                                        if ifef in faces_remain:
                                            faces_remain.remove(ifef)
                                        if ifef not in faces_done:
                                            faces_done.append(ifef)
                                    if skip_append == False:
                                        next_round.append(ifef)
                                    ifef[process_sequence_attr] = process_sequence
                                    if di == 1 or di == 3:
                                        if di == 1: ifef[xi] = init_face[xi] - 1
                                        if di == 3: ifef[xi] = init_face[xi] + 1
                                        ifef[yi] = init_face[yi]
                                    if di == 0 or di == 2:
                                        ifef[xi] = init_face[xi]
                                        if di == 0: ifef[yi] = init_face[yi] + 1
                                        if di == 2: ifef[yi] = init_face[yi] - 1
                                    ife[process_sequence_edge] = process_sequence
                                    ife[nom_dir_edge] = di
                                    process_sequence += 1
                # endregion

            if props.verbose:
                if initial_face_found == False:
                    if props.verbose: print(f"No usable initial face found for island {ii}.")
            # endregion

            bpy.context.active_object.data.update()
            round = 0
            do_a = True
            retries = 0

            # region traverse mesh
            tododo = len(faces_remain)
            while tododo > 0:
                if props.verbose: print(f"{len(faces_remain)} remaining faces.")
                round += 1
                this_round = next_round.copy()
                next_round = []

                for face in faces_remain:
                    if face[trailing_tri] == 1:
                        faces_remain.remove(face)
                        face.select = True
                        for loop in face.loops:
                            loop[uv_layer].pin_uv = False
                        if props.verbose: print(f"Removed face {face.index} from remaining faces because it is a trailing triangle.")
                for ttf in this_round:
                    if ttf[trailing_tri] == 1:
                        this_round.remove(ttf)
                        ttf.select = True
                        for ttfl in ttf.loops:
                            ttfl[uv_layer].pin_uv = False
                        if props.verbose: print(f"Removed face {ttf.index} from round because it is a trailing triangle.")

                for trf in this_round:
                    # region proc t face
                    if props.verbose: print(f"Island {ii} - Round {round}: Processing face {trf.index} with direction {trf[nom_dir]}.")
                    these_edges = [bm.edges[trf[edge_u]], bm.edges[trf[edge_l]], bm.edges[trf[edge_d]], bm.edges[trf[edge_r]]]

                    is_virtual_quad = False
                    virtual_quad = [trf, None]
                    virtual_verts = []
                    virtual_loops = []
                    corners = []
                    cos = [None, None, None, None]
                    nom_cos = [
                        mathutils.Vector((bm.faces[trf[nom]][co_ur_x], bm.faces[trf[nom]][co_ur_y])),
                        mathutils.Vector((bm.faces[trf[nom]][co_ul_x], bm.faces[trf[nom]][co_ul_y])),
                        mathutils.Vector((bm.faces[trf[nom]][co_ll_x], bm.faces[trf[nom]][co_ll_y])),
                        mathutils.Vector((bm.faces[trf[nom]][co_lr_x], bm.faces[trf[nom]][co_lr_y]))
                        ]

                    if trf in faces_remain:
                        faces_remain.remove(trf)
                    if trf not in faces_done:
                        faces_done.append(trf)

                    trf.select = True

                    if trf[virtual_quad_attr] == 1:
                        trof = bm.faces[trf[vq_other_tri]]
                        trof.select = True
                        virtual_quad[0] = trf
                        virtual_quad[1] = trof
                        is_virtual_quad = True
                        if trof in faces_remain:
                            faces_remain.remove(trof)
                    
                    avg_edge_length_x = (these_edges[0].calc_length() + these_edges[2].calc_length()) / 2
                    avg_edge_length_y = (these_edges[1].calc_length() + these_edges[3].calc_length()) / 2
                    
                    for vf in virtual_quad:
                        if vf != None:
                            vf[avg_length_x] = avg_edge_length_x
                            vf[avg_length_y] = avg_edge_length_y
                            for vfv in vf.verts:
                                if vfv not in virtual_verts:
                                    virtual_verts.append(vfv)
                            for vfl in vf.loops:
                                if vfl not in virtual_loops:
                                    virtual_loops.append(vfl)
                    corners = get_corner_verts(these_edges, virtual_verts)
                    corner_missing = False

                    for corner_i, corner in enumerate(corners[0]):
                        if corner == None:
                            for vf in virtual_quad:
                                if vf != None:
                                    vf[trailing_tri] = 1
                                    vq_is = [None, None]
                                    for vqfi, vqf in enumerate(virtual_quad):
                                        try:
                                            vq_is[vqfi] = vqf.index
                                        except:
                                            vq_is.pop()
                                    if props.verbose == True and corner_missing == False:
                                        print(f"Virtual quad corner {corner_i} failed to compute! Skipping virtual quad {vq_is}...")
                                    corner_missing = True
                    
                    if corner_missing == False:
                        for ci, corner in enumerate(corners[0]):
                            for corner_loop in corner.link_loops:
                                if corner_loop in virtual_loops:
                                    corner_loop[loop_corner_index] = ci
                                    if ci == 0:
                                        corner_loop[ind_face_uv][0] = 1.0
                                        corner_loop[ind_face_uv][1] = 1.0
                                    if ci == 1:
                                        corner_loop[ind_face_uv][0] = 0.0
                                        corner_loop[ind_face_uv][1] = 1.0
                                    if ci == 2:
                                        corner_loop[ind_face_uv][0] = 0.0
                                        corner_loop[ind_face_uv][1] = 0.0
                                    if ci == 3:
                                        corner_loop[ind_face_uv][0] = 1.0
                                        corner_loop[ind_face_uv][1] = 0.0
                        
                        for ci, corner in zip(corners[1], corners[0]):
                            for corner_loop in corner.link_loops:
                                if corner_loop in virtual_loops:
                                    if ci == 0:
                                        corner_loop[ind_face_uv_orig][0] = 1.0
                                        corner_loop[ind_face_uv_orig][1] = 1.0
                                    if ci == 1:
                                        corner_loop[ind_face_uv_orig][0] = 0.0
                                        corner_loop[ind_face_uv_orig][1] = 1.0
                                    if ci == 2:
                                        corner_loop[ind_face_uv_orig][0] = 0.0
                                        corner_loop[ind_face_uv_orig][1] = 0.0
                                    if ci == 3:
                                        corner_loop[ind_face_uv_orig][0] = 1.0
                                        corner_loop[ind_face_uv_orig][1] = 0.0
                        
                        if props.verbose:
                            if is_virtual_quad == False:
                                wnd = check_winding(trf, corners)
                                print(viz_quad(these_edges, corners[0], trf, virtual_quad, is_virtual_quad, wnd))
                            else:
                                print(viz_quad(these_edges, corners[0], trf, virtual_quad, is_virtual_quad, ""))

                        if trf[nom_dir] == 0:
                            cos[0] = nom_cos[0] + mathutils.Vector((0,avg_edge_length_y))
                            cos[1] = nom_cos[1] + mathutils.Vector((0,avg_edge_length_y))
                            cos[2] = nom_cos[1]
                            cos[3] = nom_cos[0]
                        elif trf[nom_dir] == 1:
                            cos[0] = nom_cos[1]
                            cos[1] = nom_cos[1] + mathutils.Vector((-avg_edge_length_x,0))
                            cos[2] = nom_cos[2] + mathutils.Vector((-avg_edge_length_x,0))
                            cos[3] = nom_cos[2]
                        elif trf[nom_dir] == 2:
                            cos[0] = nom_cos[3]
                            cos[1] = nom_cos[2]
                            cos[2] = nom_cos[2] + mathutils.Vector((0,-avg_edge_length_y))
                            cos[3] = nom_cos[3] + mathutils.Vector((0,-avg_edge_length_y))
                        elif trf[nom_dir] == 3:
                            cos[0] = nom_cos[1] + mathutils.Vector((avg_edge_length_x,0))
                            cos[1] = nom_cos[0]
                            cos[2] = nom_cos[3]
                            cos[3] = nom_cos[3] + mathutils.Vector((avg_edge_length_x,0))
                        # endregion

                        # region t uv edit
                        if self.pre_calc_edge_lengths == False:
                            for vi, cornerv in enumerate(corners[0]):
                                for vil in cornerv.link_loops:
                                    if vil.index in island_loops:
                                        should_edit = True
                                        if self.only_move_loops_in_face == True:
                                            should_edit = vil in virtual_loops
                                        if should_edit == True:
                                            vil[uv_layer].uv = cos[vi]
                                            new_island_max_x = max(island_max_x, cos[vi].x)
                                            island_max_x = new_island_max_x
                                            new_island_min_x = min(island_min_x, cos[vi].x)
                                            island_min_x = new_island_min_x
                                            new_island_max_y = max(island_max_y, cos[vi].y)
                                            island_max_y = new_island_max_y
                                            new_island_min_y = min(island_min_y, cos[vi].y)
                                            island_min_y = new_island_min_y
                                            vil[uv_layer].pin_uv = True
                                if cornerv not in verts_done:
                                    verts_done.append(cornerv)
                            for trfvco, trfva in zip(cos, co_attr_arr):
                                for vif in virtual_quad:
                                    if vif != None:
                                        vif[trfva[0]] = trfvco.x
                                        vif[trfva[1]] = trfvco.y
                        # endregion
                        
                        # region t pre-calc prep
                        if self.pre_calc_edge_lengths == True:
                            this_edge_length_map_x_item = [trf, is_virtual_quad, virtual_quad, [corners[0][1],corners[0][2]], [corners[0][0],corners[0][3]]]
                            this_edge_length_map_y_item = [trf, is_virtual_quad, virtual_quad, [corners[0][2],corners[0][3]], [corners[0][0],corners[0][1]]]
                            map_x_index = 0
                            map_y_index = 0
                            try:
                                map_x_index = edge_length_map_x[0].index(trf[xi])
                            except:
                                map_x_index = len(edge_length_map_x[0])
                                edge_length_map_x[0].append(trf[xi])
                                edge_length_map_x[1].append([trf[xi], [], 0.0, 0])
                            edge_length_map_x[1][map_x_index][1].append(this_edge_length_map_x_item)
                            edge_length_map_x[1][map_x_index][2] += avg_edge_length_x
                            edge_length_map_x[1][map_x_index][3] += 1

                            try:
                                map_y_index = edge_length_map_y[0].index(trf[yi])
                            except:
                                map_y_index = len(edge_length_map_y[0])
                                edge_length_map_y[0].append(trf[yi])
                                edge_length_map_y[1].append([trf[yi], [], 0.0, 0])
                            edge_length_map_y[1][map_y_index][1].append(this_edge_length_map_y_item)
                            edge_length_map_y[1][map_y_index][2] += avg_edge_length_y
                            edge_length_map_y[1][map_y_index][3] += 1
                        # endregion
                        
                        # region fill next round
                        if self.unfold_mode == "A":
                            dirs_to_do =  [0,1,2,3]
                            for di, ife in enumerate(these_edges):
                                if di in dirs_to_do:
                                    for ifef in ife.link_faces:
                                        if ifef in faces_remain and ifef not in next_round and ifef != trf:
                                            adj_face = trf
                                            other_verts = [ifefv for ifefv in ifef.verts]
                                            if is_virtual_quad == True:
                                                for vif in virtual_quad:
                                                    if ife in vif.edges:
                                                        adj_face = vif
                                            other_dir_arr = get_other_dir_edges(these_edges, corners, adj_face, ifef, di)
                                            ifef[nom_dir] = di
                                            ifef[nom] = adj_face.index
                                            process_sequence += 1
                                            ifef[process_sequence_attr] = process_sequence
                                            other_virtual_face = [ifef, bm.faces[ifef[vq_other_tri]]]
                                            other_face_is_virtual = ifef[virtual_quad_attr]
                                            if other_face_is_virtual == True:
                                                for ovf in other_virtual_face:
                                                    for ovfv in ovf.verts:
                                                        if ovfv not in other_verts:
                                                            other_verts.append(ovfv)
                                            if ifef[trailing_tri] == 1:
                                                if props.verbose: print(f"Island {ii} - Round {round}: Skipped face index {ifef.index} in direction {di} because it is a trailing triangle.")
                                                ifef[process_sequence_attr] = process_sequence
                                                ife[process_sequence_edge] = process_sequence
                                                ife[nom_dir_edge] = di
                                                process_sequence += 1
                                                faces_remain.remove(ifef)
                                                ifef.select = True
                                                for ifefl in ifef.loops:
                                                    ifefl[uv_layer].pin_uv = False
                                            if ifef[trailing_tri] == 0:
                                                skip_append = False
                                                try:
                                                    other_corner_verts = get_corner_verts(other_dir_arr, other_verts)
                                                except:
                                                    skip_append = True
                                                if props.verbose:
                                                    print(f"Island {ii} - Round {round}: Added face index {ifef.index} in direction {di} to next round of island {ii}. Virtual face: {other_face_is_virtual}. Current process sequence index: {process_sequence}")
                                                    try:
                                                        print(viz_add_quad(these_edges, other_dir_arr, corners[0], other_corner_verts[0], adj_face, ifef, virtual_quad, other_virtual_face, is_virtual_quad, other_face_is_virtual, di))
                                                    except:
                                                        print(f"Failed to calculate corners for face {ifef.index}! Operation will probably fail...")
                                                if skip_append == True:
                                                    ifef[trailing_tri] = 1
                                                    if ifef in faces_remain:
                                                        faces_remain.remove(ifef)
                                                    if ifef not in faces_done:
                                                        faces_done.append(ifef)
                                                if skip_append == False:
                                                    next_round.append(ifef)
                                                ifef[process_sequence_attr] = process_sequence
                                                ife[process_sequence_edge] = process_sequence
                                                ife[nom_dir_edge] = di
                                                if di == 1 or di == 3:
                                                    if di == 1: ifef[xi] = adj_face[xi] - 1
                                                    if di == 3: ifef[xi] = adj_face[xi] + 1
                                                    ifef[yi] = adj_face[yi]
                                                if di == 0 or di == 2:
                                                    ifef[xi] = adj_face[xi]
                                                    if di == 0: ifef[yi] = adj_face[yi] + 1
                                                    if di == 2: ifef[yi] = adj_face[yi] - 1
                                                process_sequence += 1

                        if self.unfold_mode == "B" or self.unfold_mode == "C":
                            dirs_to_do = [1, 3]
                            dirs_to_do_a = [1, 3]
                            dirs_to_do_b = [0, 2]
                            if self.unfold_mode == "C":
                                dirs_to_do_a = [0, 2]
                                dirs_to_do_b = [1, 3]

                            if do_a == True:
                                dirs_to_do = dirs_to_do_a
                            if do_a == False:
                                dirs_to_do = dirs_to_do_b

                            for di, ife in enumerate(these_edges):
                                if di in dirs_to_do:
                                    for ifef in ife.link_faces:
                                        if ifef in faces_remain and ifef not in next_round and ifef != trf:
                                            adj_face = trf
                                            other_verts = [ifefv for ifefv in ifef.verts]
                                            if is_virtual_quad == True:
                                                for vif in virtual_quad:
                                                    if ife in vif.edges:
                                                        adj_face = vif
                                            try:
                                                other_dir_arr = get_other_dir_edges(these_edges, corners, adj_face, ifef, di)
                                            except:
                                                ifef[trailing_tri] = 1
                                            ifef[nom_dir] = di
                                            ifef[nom] = adj_face.index
                                            ifef[process_sequence_attr] = process_sequence
                                            other_virtual_face = [ifef, bm.faces[ifef[vq_other_tri]]]
                                            other_face_is_virtual = ifef[virtual_quad_attr]
                                            if other_face_is_virtual == True:
                                                for ovf in other_virtual_face:
                                                    for ovfv in ovf.verts:
                                                        if ovfv not in other_verts:
                                                            other_verts.append(ovfv)
                                            if ifef[trailing_tri] == 1:
                                                if props.verbose: print(f"Island {ii} - Round {round}: Skipped face index {ifef.index} in direction {di} because it is a trailing triangle.")
                                                ifef[process_sequence_attr] = process_sequence
                                                ife[process_sequence_edge] = process_sequence
                                                ife[nom_dir_edge] = di
                                                process_sequence += 1
                                                faces_remain.remove(ifef)
                                                ifef.select = True
                                                for ifefl in ifef.loops:
                                                    ifefl[uv_layer].pin_uv = False
                                            if ifef[trailing_tri] == 0:
                                                skip_append = False
                                                try:
                                                    other_corner_verts = get_corner_verts(other_dir_arr, other_verts)
                                                except:
                                                    skip_append = True
                                                if props.verbose:
                                                    print(f"Island {ii} - Round {round}: Added face index {ifef.index} in direction {di} to next round of island {ii}. Virtual face: {other_face_is_virtual}. Current process sequence index: {process_sequence}")
                                                    try:
                                                        print(viz_add_quad(these_edges, other_dir_arr, corners[0], other_corner_verts[0], adj_face, ifef, virtual_quad, other_virtual_face, is_virtual_quad, other_face_is_virtual, di))
                                                    except:
                                                        print(f"Failed to calculate corners for face {ifef.index}! Operation will probably fail...")
                                                if skip_append == True:
                                                    ifef[trailing_tri] = 1
                                                    if ifef in faces_remain:
                                                        faces_remain.remove(ifef)
                                                    if ifef not in faces_done:
                                                        faces_done.append(ifef)
                                                if skip_append == False:
                                                    next_round.append(ifef)
                                                ifef[process_sequence_attr] = process_sequence
                                                ife[process_sequence_edge] = process_sequence
                                                ife[nom_dir_edge] = di
                                                if di == 1 or di == 3:
                                                    if di == 1: ifef[xi] = adj_face[xi] - 1
                                                    if di == 3: ifef[xi] = adj_face[xi] + 1
                                                    ifef[yi] = adj_face[yi]
                                                if di == 0 or di == 2:
                                                    ifef[xi] = adj_face[xi]
                                                    if di == 0: ifef[yi] = adj_face[yi] + 1
                                                    if di == 2: ifef[yi] = adj_face[yi] - 1
                                                process_sequence += 1
                    # endregion
                
                if len(next_round) == 0:
                    if self.unfold_mode == "B" or self.unfold_mode == "C":
                        if do_a == True:
                            do_a = False
                            if props.verbose: print("No more in main axis, trying secondary axis")
                        elif do_a == False:
                            do_a = True
                            if props.verbose: print("No more in secondary axis, returning to main axis")
                        next_round = faces_done.copy()
                    else:
                        if len(faces_remain) > 0:
                            next_round.append(faces_remain[0])
                        retries += 1

                tododo = len(faces_remain)

                bpy.context.active_object.data.update()
            
            # endregion

            if self.pre_calc_edge_lengths == False:
                offset_add_x = island_max_x - island_min_x
                if props.verbose == True: print(f"Island {ii} width: {offset_add_x}")
                offset_add_y = island_max_y - island_min_y
                if props.verbose == True: print(f"Island {ii} height: {offset_add_y}")
                previous_island_end_x += (self.offset_per_island[0] * offset_add_x)
                previous_island_end_y += (self.offset_per_island[1] * offset_add_y)

            # region apply pre-calc
            if self.pre_calc_edge_lengths == True and (len(edge_length_map_x[0]) + len(edge_length_map_y[0])) > 0:
                sorted_x = edge_length_map_x[1].copy()
                sorted_x.sort(key=map_sort_by_coord)
                sorted_y = edge_length_map_y[1].copy()
                sorted_y.sort(key=map_sort_by_coord)

                list_a = sorted_x
                list_b = sorted_y

                if self.unfold_mode == "C":
                    list_a = sorted_y
                    list_b = sorted_x

                last_a_len = previous_island_end_x
                last_b_len = previous_island_end_y
                if self.unfold_mode == "C":
                    last_a_len = previous_island_end_y
                    last_b_len = previous_island_end_x
                
                this_total_a_len = 0.0
                this_total_b_len = 0.0

                a_name = "column"
                b_name = "row"
                if self.unfold_mode == "C":
                    a_name = "row"
                    b_name = "column"
                
                if props.verbose: print(f"Island {ii} has {len(list_a)} {a_name}s and {len(list_b)} {b_name}s ")

                for list_a_item in list_a:
                    a_len = list_a_item[2]/list_a_item[3]
                    old_a_len = 0.0 + last_a_len
                    new_a_len = 0.0 + last_a_len + a_len
                    if props.verbose: print(f"Average edge length for {a_name} {list_a_item[0]}: {a_len}\nTotal previous length: {last_a_len}")
                    for face_data in list_a_item[1]:
                        for alc in face_data[3]:
                            for alcl in alc.link_loops:
                                if alcl.index in island_loops:
                                    should_edit = True
                                    if self.only_move_loops_in_face == True:
                                        should_edit = alcl in virtual_loops
                                    if should_edit == True:
                                        if self.unfold_mode == "C":
                                            alcl[uv_layer].uv.y = old_a_len
                                        else:
                                            alcl[uv_layer].uv.x = old_a_len
                        for amc in face_data[4]:
                            for amcl in amc.link_loops:
                                if amcl.index in island_loops:
                                    should_edit = True
                                    if self.only_move_loops_in_face == True:
                                        should_edit = amcl in virtual_loops
                                    if should_edit == True:
                                        if self.unfold_mode == "C":
                                            amcl[uv_layer].uv.y = new_a_len
                                        else:
                                            amcl[uv_layer].uv.x = new_a_len
                    last_a_len += a_len
                    this_total_a_len += a_len
                
                bpy.context.active_object.data.update()

                for list_b_item in list_b:
                    b_len = list_b_item[2]/list_b_item[3]
                    old_b_len = 0.0 + last_b_len
                    new_b_len = 0.0 + last_b_len + b_len
                    if props.verbose: print(f"Average edge length for {b_name} {list_b_item[0]}: {b_len}\nTotal previous length: {last_b_len}")
                    for face_data in list_b_item[1]:
                        for blc in face_data[3]:
                            for blcl in blc.link_loops:
                                if blcl.index in island_loops:
                                    should_edit = True
                                    if self.only_move_loops_in_face == True:
                                        should_edit = blcl in virtual_loops
                                    if should_edit == True:
                                        if self.unfold_mode == "C":
                                            blcl[uv_layer].uv.x = old_b_len
                                        else:
                                            blcl[uv_layer].uv.y = old_b_len
                                        blcl[uv_layer].pin_uv = True
                        for bmc in face_data[4]:
                            for bmcl in bmc.link_loops:
                                if bmcl.index in island_loops:
                                    should_edit = True
                                    if self.only_move_loops_in_face == True:
                                        should_edit = bmcl in virtual_loops
                                    if should_edit == True:
                                        if self.unfold_mode == "C":
                                            bmcl[uv_layer].uv.x = new_b_len
                                        else:
                                            bmcl[uv_layer].uv.y = new_b_len
                                        bmcl[uv_layer].pin_uv = True
                    last_b_len += b_len
                    this_total_b_len += b_len
                
                if self.unfold_mode == "C":
                    previous_island_end_x += (self.offset_per_island[0] * this_total_b_len)
                    previous_island_end_y += (self.offset_per_island[1] * this_total_a_len)
                else:
                    previous_island_end_x += (self.offset_per_island[0] * this_total_a_len)
                    previous_island_end_y += (self.offset_per_island[1] * this_total_b_len)
                bpy.context.active_object.data.update()
            # endregion

            for face in island:
                face[process_sequence_max] = process_sequence
        
        bpy.context.active_object.data.update()

        try:
            bpy.ops.uv.unwrap(method='MINIMUM_STRETCH', fill_holes=True, correct_aspect=True, use_subsurf_data=False, margin=0, no_flip=False, iterations=10, use_weights=False, weight_group="uv_importance", weight_factor=1)
        except:
            bpy.ops.uv.unwrap(method='ANGLE_BASED', fill_holes=True, correct_aspect=True, use_subsurf_data=False, margin=0, no_flip=False, iterations=10, use_weights=False, weight_group="uv_importance", weight_factor=1)

        bpy.context.active_object.data.update()

        for face in bm.faces:
            for loop in face.loops:
                loop[uv_layer].pin_uv = False
        return {"FINISHED"}

class SloppyDeTri(bpy.types.Operator):
    bl_idname = "operator.sloppy_detri"
    bl_label = "Detriangulate from Selected"
    bl_description = "Attempt to detriangulate mesh using either two selected triangles or selected diagonal edge"
    bl_options = {'REGISTER', 'UNDO'}

    iP = bpy.props.IntProperty
    fP = bpy.props.FloatProperty
    fvP = bpy.props.FloatVectorProperty
    bP = bpy.props.BoolProperty
    eP = bpy.props.EnumProperty
    bvP = bpy.props.BoolVectorProperty
    sP = bpy.props.StringProperty

    def update_unfold_mode(self, context):
        context.scene.sloppy_props.procuv_unfold_mode = self.unfold_mode
        return None
    
    def update_initial_quad(self, context):
        context.scene.sloppy_props.procuv_initial_quad = self.initial_quad
        return None
        
    # region ProcUV Properties
    unfold_mode : eP(
        name = "Mode",
        description = "Unfolding Mode",
        items = [
            ("A", "All Directions", "Work outwards from initial quad"),
            ("B", "X First", "Work outwards in (local) X first until no more faces are found, calculating edge lengths to use in next stage, then work outwards in Y. repeating first step when more faces are found in X"),
            ("C", "Y First", "Work outwards in (local) Y first until no more faces are found, calculating edge lengths to use in next stage, then work outwards in X. repeating first step when more faces are found in Y")
            ],
        update=update_unfold_mode
        ) # type: ignore

    respect_seams : bP(
        name = "Respect Seams",
        description = "Do not cross seams",
        default = False
        ) # type: ignore
    #endregion

    def execute(self, context):
        props = context.scene.sloppy_props
        bm = bmesh.from_edit_mesh(bpy.context.active_object.data)
        uv_layer = bm.loops.layers.uv.verify()
        islands = bmesh_utils.bmesh_linked_uv_islands(bm, uv_layer)
        
        faces_done = []
        edges_done = []

        suitable_selection = False
        init_tris = []
        init_diagonal = None
        for face in bm.faces:
            if face.select == True:
                if len(face.verts) == 3:
                    init_tris.append(face)
        if len(init_tris) == 2:
            for edge in init_tris[0].edges:
                if edge in init_tris[1].edges:
                    init_diagonal = edge
                    suitable_selection = True
        if suitable_selection == False:
            for edge in bm.edges:
                if edge.select == True:
                    init_diagonal = edge
        if init_diagonal != None:
            for ef in edge.link_faces:
                if len(ef.verts) != 3:
                    suitable_selection = False

        next_faces = bmesh.ops.dissolve_edges(bm, [init_diagonal], False)

        # region init vars
        current_dir_array = []
        max_avg_angle = 0.0
        min_avg_angle = 360.0
        current_avg_normal = mathutils.Vector((0,0,0))
        current_normal = mathutils.Vector((0,0,0))
        current_face_center = mathutils.Vector((0,0,0))

        prev_co_dict = [
                {
                    "xy_str": "x0y0",
                    "co_ur": mathutils.Vector((0,0)),
                    "co_ul": mathutils.Vector((0,0)),
                    "co_ll": mathutils.Vector((0,0)),
                    "co_lr": mathutils.Vector((0,0)),
                },
            ]

        # region Attributes
        attr_dict = [
            {"name": "edge_u", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "edge_l", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "edge_d", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "edge_r", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "corner_ur", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "corner_ul", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "corner_ll", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "corner_lr", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "co_ur_x", "type": "FLOAT", "domain": "FACE", "layer": None},
            {"name": "co_ul_x", "type": "FLOAT", "domain": "FACE", "layer": None},
            {"name": "co_ll_x", "type": "FLOAT", "domain": "FACE", "layer": None},
            {"name": "co_lr_x", "type": "FLOAT", "domain": "FACE", "layer": None},
            {"name": "co_ur_y", "type": "FLOAT", "domain": "FACE", "layer": None},
            {"name": "co_ul_y", "type": "FLOAT", "domain": "FACE", "layer": None},
            {"name": "co_ll_y", "type": "FLOAT", "domain": "FACE", "layer": None},
            {"name": "co_lr_y", "type": "FLOAT", "domain": "FACE", "layer": None},
            {"name": "loop_corner_index", "type": "INT", "domain": "CORNER", "layer": None},
            {"name": "ind_face_uv", "type": "FLOAT_VECTOR", "domain": "CORNER", "layer": None},
            {"name": "ind_face_uv_orig", "type": "FLOAT_VECTOR", "domain": "CORNER", "layer": None},
            {"name": "virtual_quad", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "vq_other_tri", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "vq_diagonal", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "trailing_tri", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "process_sequence", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "process_sequence_edge", "type": "INT", "domain": "EDGE", "layer": None},
            {"name": "process_sequence_max", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "xi", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "yi", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "island_index", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "flatness", "type": "FLOAT", "domain": "FACE", "layer": None},
            {"name": "normal_regularity", "type": "FLOAT", "domain": "FACE", "layer": None},
            {"name": "nom", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "nom_dir", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "nom_dir_edge", "type": "INT", "domain": "EDGE", "layer": None},
            {"name": "avg_length_x", "type": "FLOAT", "domain": "FACE", "layer": None},
            {"name": "avg_length_y", "type": "FLOAT", "domain": "FACE", "layer": None},
            ]
        
        for nam in attr_dict:
            attr = props.find_or_add_attribute(nam["name"], nam["type"], nam["domain"])
            nam["layer"] = props.get_attribute_layer(nam["name"], nam["type"], nam["domain"], bm)
        
        flatness = props.get_dict_layer("flatness", attr_dict)
        normal_regularity = props.get_dict_layer("normal_regularity", attr_dict)
        island_index = props.get_dict_layer("island_index", attr_dict)
        edge_u = props.get_dict_layer("edge_u", attr_dict)
        edge_l = props.get_dict_layer("edge_l", attr_dict)
        edge_d = props.get_dict_layer("edge_d", attr_dict)
        edge_r = props.get_dict_layer("edge_r", attr_dict)
        corner_ur = props.get_dict_layer("corner_ur", attr_dict)
        corner_ul = props.get_dict_layer("corner_ul", attr_dict)
        corner_ll = props.get_dict_layer("corner_ll", attr_dict)
        corner_lr = props.get_dict_layer("corner_lr", attr_dict)
        co_ur_x = props.get_dict_layer("co_ur_x", attr_dict)
        co_ul_x = props.get_dict_layer("co_ul_x", attr_dict)
        co_ll_x = props.get_dict_layer("co_ll_x", attr_dict)
        co_lr_x = props.get_dict_layer("co_lr_x", attr_dict)
        co_ur_y = props.get_dict_layer("co_ur_y", attr_dict)
        co_ul_y = props.get_dict_layer("co_ul_y", attr_dict)
        co_ll_y = props.get_dict_layer("co_ll_y", attr_dict)
        co_lr_y = props.get_dict_layer("co_lr_y", attr_dict)
        virtual_quad_attr = props.get_dict_layer("virtual_quad", attr_dict)
        diagonal_index = props.get_dict_layer("vq_diagonal", attr_dict)
        vq_other_tri = props.get_dict_layer("vq_other_tri", attr_dict)
        trailing_tri = props.get_dict_layer("trailing_tri", attr_dict)
        process_sequence_attr = props.get_dict_layer("process_sequence", attr_dict)
        process_sequence_edge = props.get_dict_layer("process_sequence_edge", attr_dict)
        nom = props.get_dict_layer("nom", attr_dict)
        nom_dir = props.get_dict_layer("nom_dir", attr_dict)
        nom_dir_edge = props.get_dict_layer("nom_dir_edge", attr_dict)
        process_sequence_max = props.get_dict_layer("process_sequence_max", attr_dict)
        loop_corner_index = props.get_dict_layer("loop_corner_index", attr_dict)
        ind_face_uv = props.get_dict_layer("ind_face_uv", attr_dict)
        ind_face_uv_orig = props.get_dict_layer("ind_face_uv_orig", attr_dict)
        xi = props.get_dict_layer("xi", attr_dict)
        yi = props.get_dict_layer("yi", attr_dict)
        avg_length_x = props.get_dict_layer("avg_length_x", attr_dict)
        avg_length_y = props.get_dict_layer("avg_length_y", attr_dict)
        
        co_attr_arr = [[co_ur_x, co_ur_y], [co_ul_x, co_ul_y], [co_ll_x, co_ll_y], [co_lr_x, co_lr_y]]
        # endregion

        quant_arr = [
            [mathutils.Vector((1,0,0)), 
                [
                    mathutils.Vector((0,0,-1)),
                    mathutils.Vector((0,-1,0)),
                    mathutils.Vector((0,0,1)),
                    mathutils.Vector((0,1,0))
                ],
                "Left"
             ],
             [mathutils.Vector((-1,0,0)),
                [
                    mathutils.Vector((0,0,-1)),
                    mathutils.Vector((0,1,0)),
                    mathutils.Vector((0,0,1)),
                    mathutils.Vector((0,-1,0))
                ],
                "Right"
             ],
             [mathutils.Vector((0,1,0)),
                [
                    mathutils.Vector((0,0,-1)),
                    mathutils.Vector((1,0,0)),
                    mathutils.Vector((0,0,1)),
                    mathutils.Vector((-1,0,0))
                ],
                "Towards"
             ],
             [mathutils.Vector((0,-1,0)),
                [
                    mathutils.Vector((0,0,-1)),
                    mathutils.Vector((-1,0,0)),
                    mathutils.Vector((0,0,1)),
                    mathutils.Vector((1,0,0))
                ],
                "Away"
             ],
             [mathutils.Vector((0,0,-1)),
                [
                    mathutils.Vector((0,-1,0)),
                    mathutils.Vector((1,0,0)),
                    mathutils.Vector((0,1,0)),
                    mathutils.Vector((-1,0,0))
                ],
                "Up"
             ],
             [mathutils.Vector((0,0,1)),
                [
                    mathutils.Vector((0,1,0)),
                    mathutils.Vector((1,0,0)),
                    mathutils.Vector((0,-1,0)),
                    mathutils.Vector((-1,0,0))
                ],
                "Down"
             ],
        ]
        # endregion

        # region functions
        def check_list_for_duplicate(the_list):
            seen = set()
            for the_item in the_list:
                if the_item in seen: return True
                seen.add(the_item)
            return False
        
        def quant_sort(e):
            quant_dot = e[0].dot(current_normal)
            return quant_dot

        def map_sort_by_coord(e):
            return e[0]

        def map_sort_by_edge_length(e):
            return e[2]

        def regularity_sort(e):
            flat = (1.0 - e[flatness]) * -self.reg_flat_fac
            norm_reg = e[normal_regularity] * -self.reg_norm_fac
            return flat + norm_reg
        
        def adj_vert_angle_sort(e):
            return e[1]
        
        def calc_edge_center(e):
            return e.verts[0].co.lerp(e.verts[1].co, 0.5)

        def dot_sort_u(e):
            evec = calc_edge_center(e) - current_face_center
            edir = evec.normalized()
            return edir.dot(current_dir_array[0])
        
        def dot_sort_l(e):
            evec = calc_edge_center(e) - current_face_center
            edir = evec.normalized()
            return edir.dot(current_dir_array[1])
        
        def dot_sort_d(e):
            evec = calc_edge_center(e) - current_face_center
            edir = evec.normalized()
            return edir.dot(current_dir_array[2])
        
        def dot_sort_r(e):
            evec = calc_edge_center(e) - current_face_center
            edir = evec.normalized()
            return edir.dot(current_dir_array[3])
        
        def get_other_dir_edges(edge_arr, corner_verts_arr, source_face, other_face, dir):
            out_array = [None, None, None, None]
            virtual_edges = [None, None, None, None]
            corner_verts = corner_verts_arr[0]
            adj_edge = None
            adjacency_verts = []
            if dir == 0:
                other_face[edge_d] = edge_arr[dir].index
                out_array[2] = edge_arr[dir]
                adj_edge = edge_arr[dir]
                adjacency_verts = [corner_verts[0], corner_verts[1]]
                virtual_edges[2] = adj_edge
            if dir == 1:
                other_face[edge_r] = edge_arr[dir].index
                out_array[3] = edge_arr[dir]
                adj_edge = edge_arr[dir]
                adjacency_verts = [corner_verts[1], corner_verts[2]]
                virtual_edges[3] = adj_edge
            if dir == 2:
                other_face[edge_u] = edge_arr[dir].index
                out_array[0] = edge_arr[dir]
                adj_edge = edge_arr[dir]
                adjacency_verts = [corner_verts[2], corner_verts[3]]
                virtual_edges[0] = adj_edge
            if dir == 3:
                other_face[edge_l] = edge_arr[dir].index
                out_array[1] = edge_arr[dir]
                adj_edge = edge_arr[dir]
                adjacency_verts = [corner_verts[3], corner_verts[0]]
                virtual_edges[1] = adj_edge
            sides = len(other_face.edges)

            if sides == 3:
                main_tri = other_face
                tri_trailing = False
                diagonal_edge = None
                far_tri = None
                adj_verts_arr = []
                for loop in main_tri.loops:
                    if loop.vert in adjacency_verts:
                        adj_vert_arr_app = [loop.vert]
                        adj_vert_angle = math.degrees(loop.calc_angle())
                        adj_vert_arr_app.append(adj_vert_angle)
                        adj_verts_arr.append(adj_vert_arr_app)
                adj_verts_arr.sort(key=adj_vert_angle_sort)
                for avf in adj_verts_arr[0][0].link_faces:
                    if avf != main_tri and avf != source_face:
                        main_tri_adjacent = False
                        for avfe in avf.edges:
                            if avfe not in edge_arr:
                                if main_tri in avfe.link_faces:
                                    if avfe.seam == False:
                                        main_tri_adjacent = True
                                        diagonal_edge = avfe
                                    if avfe.seam == True:
                                        tri_trailing = True
                        if main_tri_adjacent == True:
                            far_tri = avf
                            if len(far_tri.edges) > 3:
                                tri_trailing = True
                if far_tri == None:
                    tri_trailing = True
                if tri_trailing == False:
                    if far_tri[island_index] == main_tri[island_index]:
                        main_tri[virtual_quad_attr] = 1
                        far_tri[virtual_quad_attr] = 1
                        main_tri[trailing_tri] = 0
                        far_tri[trailing_tri] = 0
                        main_tri[vq_other_tri] = far_tri.index
                        far_tri[vq_other_tri] = main_tri.index
                        main_tri[diagonal_index] = diagonal_edge.index
                        avf[diagonal_index] = diagonal_edge.index
                    if far_tri[island_index] != main_tri[island_index]:
                        main_tri[trailing_tri] = 1
                        main_tri[virtual_quad_attr] = 0
                        tri_trailing = True
                if tri_trailing == True:
                    main_tri[trailing_tri] = 1
                    main_tri.select = True
                    main_tri[virtual_quad_attr] = 0
                    out_array = []

                if tri_trailing == False:
                    for mte in main_tri.edges:
                        if mte != adj_edge and mte != diagonal_edge:
                            if dir == 0:
                                if mte in corner_verts[0].link_edges:
                                    virtual_edges[3] = mte
                                    out_array[3] = mte
                                if mte in corner_verts[1].link_edges:
                                    virtual_edges[1] = mte
                                    out_array[1] = mte
                            if dir == 1:
                                if mte in corner_verts[1].link_edges:
                                    virtual_edges[0] = mte
                                    out_array[0] = mte
                                if mte in corner_verts[2].link_edges:
                                    virtual_edges[2] = mte
                                    out_array[2] = mte
                            if dir == 2:
                                if mte in corner_verts[2].link_edges:
                                    virtual_edges[1] = mte
                                    out_array[1] = mte
                                if mte in corner_verts[3].link_edges:
                                    virtual_edges[3] = mte
                                    out_array[3] = mte
                            if dir == 3:
                                if mte in corner_verts[3].link_edges:
                                    virtual_edges[2] = mte
                                    out_array[2] = mte
                                if mte in corner_verts[0].link_edges:
                                    virtual_edges[0] = mte
                                    out_array[0] = mte
                    for fte in far_tri.edges:
                        if fte != diagonal_edge:
                            if dir == 0:
                                if fte in corner_verts[0].link_edges:
                                    virtual_edges[3] = fte
                                    out_array[3] = fte
                                elif fte in corner_verts[1].link_edges:
                                    virtual_edges[1] = fte
                                    out_array[1] = fte
                                else:
                                    virtual_edges[0] = fte
                                    out_array[0] = fte
                            if dir == 1:
                                if fte in corner_verts[1].link_edges:
                                    virtual_edges[0] = fte
                                    out_array[0] = fte
                                elif fte in corner_verts[2].link_edges:
                                    virtual_edges[2] = fte
                                    out_array[2] = fte
                                else:
                                    virtual_edges[1] = fte
                                    out_array[1] = mte
                            if dir == 2:
                                if fte in corner_verts[2].link_edges:
                                    virtual_edges[1] = fte
                                    out_array[1] = fte
                                elif fte in corner_verts[3].link_edges:
                                    virtual_edges[3] = fte
                                    out_array[3] = fte
                                else:
                                    virtual_edges[2] = fte
                                    out_array[2] = fte
                            if dir == 3:
                                if fte in corner_verts[3].link_edges:
                                    virtual_edges[2] = fte
                                    out_array[2] = fte
                                elif fte in corner_verts[0].link_edges:
                                    virtual_edges[0] = fte
                                    out_array[0] = fte
                                else:
                                    virtual_edges[3] = fte
                                    out_array[3] = fte
                    main_tri[edge_u] = virtual_edges[0].index
                    far_tri[edge_u] = virtual_edges[0].index
                    main_tri[edge_l] = virtual_edges[1].index
                    far_tri[edge_l] = virtual_edges[1].index
                    main_tri[edge_d] = virtual_edges[2].index
                    far_tri[edge_d] = virtual_edges[2].index
                    main_tri[edge_r] = virtual_edges[3].index
                    far_tri[edge_r] = virtual_edges[3].index

            if sides == 4:
                other_face[virtual_quad_attr] = 0
                for ofe in other_face.edges:
                    if ofe not in source_face.edges:
                        newdir = None
                        if dir == 0:
                            if ofe in corner_verts[0].link_edges:
                                out_array[3] = ofe
                                newdir = 3
                            if ofe in corner_verts[1].link_edges:
                                out_array[1] = ofe
                                newdir = 1
                            if ofe not in corner_verts[0].link_edges and ofe not in corner_verts[1].link_edges:
                                out_array[0] = ofe
                                newdir = 0
                        if dir == 1:
                            if ofe in corner_verts[1].link_edges:
                                out_array[0] = ofe
                                newdir = 0
                            if ofe in corner_verts[2].link_edges:
                                out_array[2] = ofe
                                newdir = 2
                            if ofe not in corner_verts[1].link_edges and ofe not in corner_verts[2].link_edges:
                                out_array[1] = ofe
                                newdir = 1
                        if dir == 2:
                            if ofe in corner_verts[2].link_edges:
                                out_array[1] = ofe
                                newdir = 1
                            if ofe in corner_verts[3].link_edges:
                                out_array[3] = ofe
                                newdir = 3
                            if ofe not in corner_verts[2].link_edges and ofe not in corner_verts[3].link_edges:
                                out_array[2] = ofe
                                newdir = 2
                        if dir == 3:
                            if ofe in corner_verts[3].link_edges:
                                out_array[2] = ofe
                                newdir = 2
                            if ofe in corner_verts[0].link_edges:
                                out_array[0] = ofe
                                newdir = 0
                            if ofe not in corner_verts[3].link_edges and ofe not in corner_verts[0].link_edges:
                                out_array[3] = ofe
                                newdir = 3
                        # print(f"New direction {newdir}: edge {ofe.index}")
                    if ofe in source_face.edges:
                        if dir == 0:
                            out_array[2] = ofe
                            newdir = 2
                        elif dir == 1:
                            out_array[3] = ofe
                            newdir = 3
                        elif dir == 2:
                            out_array[0] = ofe
                            newdir = 0
                        elif dir == 3:
                            out_array[1] = ofe
                            newdir = 1
                        # print(f"New direction {newdir}: edge {ofe.index}. Adjacent edge.")
                # for edgi, edge in enumerate(out_array):
                    # print(f"{edgi} - Direction:{dir}: {edge}")
                other_face[edge_u] = out_array[0].index
                other_face[edge_l] = out_array[1].index
                other_face[edge_d] = out_array[2].index
                other_face[edge_r] = out_array[3].index

            return out_array
        
        def check_winding(face, corner_verts_arr):
            poss = 0
            negs = 0
            nones = 0
            out_str = ""
            last_i = corner_verts_arr[1][-1]
            for vi in corner_verts_arr[1]:
                if vi == None:
                    nones += 1
                else:
                    result = vi - last_i
                    if result > 0:
                        poss += 1
                    if result < 0:
                        negs += 1
                    last_i = vi
            if poss > negs:
                out_str = "CCW"
            if poss < negs:
                out_str = "CW"
            if nones > 0:
                out_str = "Winding check failed"
            return out_str
        
        def tri_angle_edge_sort(e):
            angles = 0.0
            for vert in e[0].verts:
                for loop in vert.link_loops:
                    if loop in e[1].loops:
                        angles += math.degrees(loop.calc_angle())
            return angles
        
        def get_adjacent_triangle(tri):
            tri_edges = []
            out_tri = [None, None]
            for te in tri.edges:
                for tef in te.link_faces:
                    if tef != tri and len(tef.edges) == 3:
                        tri_edge = []
                        tri_edge.append(te)
                        tri_edge.append(tri)
                        tri_edge.append(tef)
                        tri_edges.append(tri_edge)
            if len(tri_edges) > 0:
                tri_edges.sort(key=tri_angle_edge_sort)
                out_tri[0] = tri_edges[0][2]
                out_tri[1] = tri_edges[0][0]
            else:
                pass
            return out_tri
        
        def get_any_adjacent_quad(source_face):
            out_quad = None
            for qe in source_face.edges:
                for qef in qe.link_faces:
                    if qef != source_face and len(qef.edges) == 4:
                        out_quad = qef
                        return out_quad
            else:
                pass
            return out_quad
        
        def get_init_dir_edges(arr):
            out_array = []
            u_arr = arr.copy()
            l_arr = arr.copy()
            d_arr = arr.copy()
            r_arr = arr.copy()
            u_arr.sort(key=dot_sort_u)
            l_arr.sort(key=dot_sort_l)
            d_arr.sort(key=dot_sort_d)
            r_arr.sort(key=dot_sort_r)
            out_array.append(u_arr[0])
            out_array.append(l_arr[0])
            out_array.append(d_arr[0])
            out_array.append(r_arr[0])
            return out_array
        
        def get_corner_verts(edge_arr, v_arr):
            out_arrays = []
            out_array = [None, None, None, None]
            out_array_is = [None, None, None, None]
            v_iter = 0
            for v in v_arr:
                ul = 0
                ll = 0
                lr = 0
                ur = 0
                if v in edge_arr[0].verts:
                    ur += 1
                    ul += 1
                if v in edge_arr[1].verts:
                    ul += 1
                    ll += 1
                if v in edge_arr[2].verts:
                    ll += 1
                    lr += 1
                if v in edge_arr[3].verts:
                    lr += 1
                    ur += 1
                if ur == 2:
                    out_array[0] = v
                    out_array_is[0] = v_iter
                    # print(f"Upper right corner: {v.index} (seq. index: {v_iter})")
                if ul == 2:
                    out_array[1] = v
                    out_array_is[1] = v_iter
                    # print(f"Upper left corner: {v.index} (seq. index: {v_iter})")
                if ll == 2:
                    out_array[2] = v
                    out_array_is[2] = v_iter
                    # print(f"Lower left corner: {v.index} (seq. index: {v_iter})")
                if lr == 2:
                    out_array[3] = v
                    out_array_is[3] = v_iter
                    # print(f"Lower right corner: {v.index} (seq. index: {v_iter})")
                v_iter += 1
            out_arrays.append(out_array)
            out_arrays.append(out_array_is)
            return out_arrays

        def viz_quad(qe, qc, qf, vqf, virt, wnd):
            fnamlen = 2 + len(str(qf.index))
            if virt == True:
                fnamlen = 3 + len(str(vqf[0].index)) + len(str(vqf[1].index))
            max_l = max(len(str(qc[1].index)), len(str(qc[2].index)), len(str(qe[1].index)))
            max_c = 2 + max(fnamlen, len(str(qe[0].index)), len(str(qe[2].index)))
            max_r = max(len(str(qc[0].index)), len(str(qe[3].index)), len(str(qc[3].index)))
            max_vline = max(max_l, max_r)
            ulc = str(qc[1].index).rjust(max_vline)
            uline = str(qe[0].index).center(max_c,"-")
            urc = str(qc[0].index).ljust(max_vline)
            vline = "|".center(max_vline)
            vspace = " ".center(max_c)
            vu = "^".center(max_vline)
            vd = "v".center(max_vline)
            larrow = "<"
            rarrow = ">"
            lline = str(qe[1].index).rjust(max_vline)
            fnam = str(qf.index).center(max_c)
            if virt == True:
                fnama = str(vqf[0].index) + "/" + str(vqf[1].index)
                fnam = fnama.center(max_c)
            rline = str(qe[3].index).ljust(max_vline)
            llc = str(qc[2].index).rjust(max_vline)
            dline = str(qe[2].index).center(max_c, "-")
            lrc = str(qc[3].index).ljust(max_vline)
            quad_msg = "\n{0}{1}{2}\n{3}{4}{3}\n{3}{4}{3}\n{3}{4}{3}\n{5}{6}{7}\n{3}{4}{3}\n{3}{4}{3}\n{3}{4}{3}\n{8}{9}{10}\n".format(ulc,uline,urc,vline,vspace,lline,fnam,rline, llc, dline, lrc)
            if wnd == "CCW":
                quad_msg = "\n{0}{1}{13}{2}\n{12}{4}{3}\n{3}{4}{3}\n{3}{4}{3}\n{5}{6}{7}\n{3}{4}{3}\n{3}{4}{3}\n{3}{4}{11}\n{8}{14}{9}{10}\n".format(ulc,uline,urc,vline,vspace,lline,fnam,rline, llc, dline, lrc, vu, vd, larrow, rarrow)
            if wnd == "CW":
                quad_msg = "\n{0}{14}{1}{2}\n{3}{4}{11}\n{3}{4}{3}\n{3}{4}{3}\n{5}{6}{7}\n{3}{4}{3}\n{3}{4}{3}\n{12}{4}{3}\n{8}{9}{13}{10}\n".format(ulc,uline,urc,vline,vspace,lline,fnam,rline, llc, dline, lrc, vu, vd, larrow, rarrow)
            return quad_msg

        def viz_add_quad(qe, oqe, qc, oqc, qf, oqf, vqf, ovqf, virt, ovirt, add_dir):
            fnamlen = 2 + len(str(qf.index))
            fnamolen = 2 + len(str(oqf.index))
            if virt == True:
                fnamlen = 3 + len(str(vqf[0].index)) + len(str(vqf[1].index))
            if ovirt == True:
                fnamolen = 3 + len(str(ovqf[0].index)) + len(str(ovqf[1].index))
            max_l = max(len(str(qc[1].index)), len(str(oqc[1].index)), len(str(qc[2].index)), len(str(oqc[2].index)), len(str(qe[1].index)), len(str(oqe[1].index)))
            max_c = 6 + max(fnamlen, fnamolen, len(str(qe[0].index)), len(str(oqe[0].index)), len(str(qe[2].index)), len(str(oqe[2].index)))
            max_r = max(len(str(qc[0].index)), len(str(oqc[0].index)), len(str(qe[3].index)), len(str(oqe[3].index)), len(str(qc[3].index)), len(str(oqc[3].index)))
            max_vline = max(max_l, max_r)
            vspace_sum = max_c + (2 * max_vline)
            ulc = str(qc[1].index).rjust(max_vline)
            oulc = str(oqc[1].index).rjust(max_vline)
            uline = str(qe[0].index).center(max_c,"=")
            ouline = str(oqe[0].index).center(max_c,"-")
            urc = str(qc[0].index).ljust(max_vline)
            ourc = str(oqc[0].index).ljust(max_vline)
            vline = "||".center(max_vline)
            ovline = "|".center(max_vline)
            vspace = " ".center(max_c)
            ovspace = " ".center(vspace_sum)
            ohspace = " ".center(max_c + 2)
            ovu = "^".center(vspace_sum)
            ovd = "v".center(vspace_sum)
            olarrow = "<".center(max_c + 2)
            orarrow = ">".center(max_c + 2)
            lline = str(qe[1].index).rjust(max_vline)
            olline = str(oqe[1].index).rjust(max_vline)
            fnam = str(qf.index).center(max_c)
            if virt == True:
                fnama = str(vqf[0].index) + "/" + str(vqf[1].index)
                fnam = fnama.center(max_c)
            ofnam = str(oqf.index).center(max_c)
            if ovirt == True:
                ofnama = str(ovqf[0].index) + "/" + str(ovqf[1].index)
                ofnam = ofnama.center(max_c)
            rline = str(qe[3].index).ljust(max_vline)
            orline = str(oqe[3].index).ljust(max_vline)
            llc = str(qc[2].index).rjust(max_vline)
            ollc = str(oqc[2].index).rjust(max_vline)
            dline = str(qe[2].index).center(max_c, "=")
            odline = str(oqe[2].index).center(max_c, "-")
            lrc = str(qc[3].index).ljust(max_vline)
            olrc = str(oqc[3].index).ljust(max_vline)

            quad_msg_aa = "{0}{1}{2}".format(ulc,uline,urc)
            quad_msg_avsp = "{0}{1}{0}".format(vline,vspace)
            quad_msg_ab = "{0}{1}{2}".format(lline,fnam,rline)
            quad_msg_ac = "{0}{1}{2}".format(llc,dline,lrc)
            quad_msg_ba = "{0}{1}{2}".format(oulc,ouline,ourc)
            quad_msg_bvsp = "{0}{1}{0}".format(ovline,vspace)
            quad_msg_bb = "{0}{1}{2}".format(olline,ofnam,orline)
            quad_msg_bc = "{0}{1}{2}".format(ollc,odline,olrc)

            quad_msg = ""
            if add_dir == 0:
                quad_msg = "\n{0}\n{3}\n{3}\n{3}\n{1}\n{3}\n{3}\n{3}\n{2}\n{8}\n{9}\n{8}\n{4}\n{7}\n{7}\n{7}\n{5}\n{7}\n{7}\n{7}\n{6}\n".format(quad_msg_ba, quad_msg_bb, quad_msg_bc, quad_msg_bvsp, quad_msg_aa, quad_msg_ab, quad_msg_ac, quad_msg_avsp, ovspace, ovu)
            if add_dir == 1:
                quad_msg = "\n{0}{8}{4}\n{3}{8}{7}\n{3}{8}{7}\n{3}{8}{7}\n{1}{9}{5}\n{3}{8}{7}\n{3}{8}{7}\n{3}{8}{7}\n{2}{8}{6}\n".format(quad_msg_ba, quad_msg_bb, quad_msg_bc, quad_msg_bvsp, quad_msg_aa, quad_msg_ab, quad_msg_ac, quad_msg_avsp, ohspace, olarrow)
            if add_dir == 2:
                quad_msg = "\n{0}\n{3}\n{3}\n{3}\n{1}\n{3}\n{3}\n{3}\n{2}\n{8}\n{9}\n{8}\n{4}\n{7}\n{7}\n{7}\n{5}\n{7}\n{7}\n{7}\n{6}\n".format(quad_msg_aa, quad_msg_ab, quad_msg_ac, quad_msg_avsp, quad_msg_ba, quad_msg_bb, quad_msg_bc, quad_msg_bvsp, ovspace, ovd)
            if add_dir == 3:
                quad_msg = "\n{0}{8}{4}\n{3}{8}{7}\n{3}{8}{7}\n{3}{8}{7}\n{1}{9}{5}\n{3}{8}{7}\n{3}{8}{7}\n{3}{8}{7}\n{2}{8}{6}\n".format(quad_msg_aa, quad_msg_ab, quad_msg_ac, quad_msg_avsp, quad_msg_ba, quad_msg_bb, quad_msg_bc, quad_msg_bvsp, ohspace, orarrow)
            return quad_msg
        # endregion

        previous_island_end_x = 0.0
        previous_island_end_y = 0.0

        for ii, island in enumerate(islands):
            for iif in island:
                iif[island_index] = ii

        for ii, island in enumerate(islands):
            bm.verts.ensure_lookup_table()
            bm.edges.ensure_lookup_table()
            bm.faces.ensure_lookup_table()
            edge_length_map_x = [[],[]]
            edge_length_map_y = [[],[]]
            island_max_x = 0.0
            island_max_y = 0.0
            island_min_x = 999.0
            island_min_y = 999.0
            process_sequence = 0
            faces_remain = []
            faces_done = []
            verts_done = []
            faces_selected = []
            island_loops = []
            init_face = None
            init_virtual_face = []
            avg_norm = mathutils.Vector((0,0,0))
            if props.verbose: print(f"Island {ii}:")

            # region prep init search
            for face in island:
                avg_norm += face.normal
                avg_angle = 0.0
                angle_count = 0
                for fe in face.edges:
                    angle_count += 1
                    if fe.is_boundary == True or len(face.edges) < 2 or len(face.edges) > 2:
                        avg_angle += 90
                    else:
                        avg_angle += math.degrees(fe.calc_face_angle())
                for fl in face.loops:
                    island_loops.append(fl.index)
                avg_angle /= angle_count
                new_max_angle = max(avg_angle, max_avg_angle)
                new_min_angle = min(avg_angle, min_avg_angle)
                max_avg_angle = new_max_angle
                min_avg_angle = new_min_angle
                face[flatness] = avg_angle
                faces_remain.append(face)
                if face.select == True:
                    faces_selected.append(face)
            
            current_avg_normal = avg_norm.normalized()

            if self.quant_avg_norm == True:
                current_normal = current_avg_normal
                quant_arr.sort(key=quant_sort)
                current_avg_normal = quant_arr[0][0] * -1.0
                if props.verbose == True: print(f"Average island normals quantized to {quant_arr[0][2]}")

            bpy.context.active_object.data.update()

            for face in island:
                new_flatness = props.remap_val(face[flatness], min_avg_angle, max_avg_angle, 0.0, 1.0)
                face[flatness] = new_flatness
                n_dot_avg = face.normal.dot(current_avg_normal)
                face[normal_regularity] = props.remap_val(n_dot_avg, -1.0, 1.0, 0.0, 1.0)
            
            bpy.context.active_object.data.update()
            # endregion

            # region find init face(s)
            initial_face_found = False
            
            while initial_face_found == False and len(faces_remain) > 0:
                is_selected_face = False
                init_face_is_virtual = False
                init_virtual_face_diagonal = None
                if self.initial_quad == "B":
                    if len(faces_selected) > 0:
                        init_face = faces_selected[0]
                        init_virtual_face = [init_face]
                        if len(init_face.edges) < 4:
                            if len(faces_selected) > 1:
                                if len(faces_selected[1].edges) == 3:
                                    init_face_is_virtual = True
                                    init_virtual_face = [init_face, faces_selected[1]]
                                    is_selected_face = True
                                if len(faces_selected[1].edges) == 4:
                                    init_face_is_virtual = False
                                    init_face = faces_selected[1]
                                    init_virtual_face = [faces_selected[1]]
                                    is_selected_face = True
                                else:
                                    if props.verbose: print("Too inconsistent mesh! Aborting!")
                                    return  {"CANCELLED"}
                            else:
                                other_tri = get_adjacent_triangle(init_face)
                                if other_tri[0] != None:
                                    init_face_is_virtual = True
                                    init_virtual_face = [init_face, other_tri[0]]
                                    is_selected_face = other_tri[0].select
                                    init_virtual_face_diagonal = other_tri[1]
                                else:
                                    other_quad = get_any_adjacent_quad(init_face)
                                    if other_quad != None:
                                        init_face = other_quad
                                        init_virtual_face = [init_face]
                                        is_selected_face = init_face.select
                                        init_face_is_virtual = False
                                    else:
                                        if props.verbose: print("Too inconsistent mesh! Aborting!")
                                        return  {"CANCELLED"}
                    if len(faces_selected) == 0:
                        faces_remain.sort(key=regularity_sort)
                        init_face = faces_remain[0]
                        init_virtual_face = [init_face]
                        if len(init_face.edges) == 3:
                            other_tri = get_adjacent_triangle(init_face)
                            if other_tri[0] != None:
                                init_face_is_virtual = True
                                init_virtual_face = [init_face, other_tri[0]]
                                init_virtual_face_diagonal = other_tri[1]
                            else:
                                other_quad = get_any_adjacent_quad(init_face)
                                if other_quad != None:
                                    init_face = other_quad
                                    init_virtual_face = [init_face]
                                else:
                                    if props.verbose: print("Too inconsistent mesh! Aborting!")
                                    return  {"CANCELLED"}
                else:
                    faces_remain.sort(key=regularity_sort)
                    init_face = faces_remain[0]
                    init_virtual_face = [init_face]
                    if len(init_face.edges) == 3:
                        other_tri = get_adjacent_triangle(init_face)
                        if other_tri[0] != None:
                            init_face_is_virtual = True
                            init_virtual_face = [init_face, other_tri[0]]
                            init_virtual_face_diagonal = other_tri[1]
                        else:
                            other_quad = get_any_adjacent_quad(init_face)
                            if other_quad != None:
                                init_face = other_quad
                                init_virtual_face = [init_face]
                            else:
                                if props.verbose: print("Too inconsistent mesh! Aborting!")
                                return  {"CANCELLED"}

                current_normal = init_face.normal

                if init_face_is_virtual == True:
                    virtual_normal = mathutils.Vector((0,0,0))
                    for vf in init_virtual_face:
                        virtual_normal += vf.normal
                    current_normal = virtual_normal.normalized()

                quant_arr.sort(key=quant_sort)
                current_quant_index = 0
                current_dir_array = quant_arr[current_quant_index][1]
                init_face_dir_arr = [None, None, None, None]

                if props.verbose == True:
                    msg_start = "Most regular face for island "
                    if is_selected_face == True:
                        msg_start = "Selected face for island "
                    msg = "{:}{:}: {:} with a flatness of {:.3f} and a normal regularity of {:.3f}.\nFacing quantized to {:}\n".format(msg_start, ii, init_face.index, init_face[flatness], init_face[normal_regularity], quant_arr[0][2])
                    if init_face_is_virtual == True:
                        msg_start = "Most regular virtual face for island "
                        if is_selected_face == True:
                            msg_start = "Selected virtual face for island "
                        msg = "{:}{:}: {:}/{:} with a flatness of {:.3f} and a normal regularity of {:.3f}.\nFacing quantized to {:}\n".format(msg_start, ii, init_virtual_face[0].index, init_virtual_face[1].index, init_face[flatness], init_face[normal_regularity], quant_arr[0][2])
                    print(msg)
                
                if init_face_is_virtual == True:
                    vif_edges = []
                    for vif in init_virtual_face:
                        vif.select = True
                    for avife in init_virtual_face[0].edges:
                        if avife != init_virtual_face_diagonal:
                            vif_edges.append(avife)
                    for bvife in init_virtual_face[1].edges:
                        if bvife != init_virtual_face_diagonal:
                            vif_edges.append(bvife)
                    current_face_center = calc_edge_center(init_virtual_face_diagonal)
                    init_face_dir_arr = get_init_dir_edges(vif_edges)
                    while check_list_for_duplicate(init_face_dir_arr) == True and current_quant_index < 4:
                        current_quant_index += 1
                        try:
                            current_dir_array = quant_arr[current_quant_index][1]
                        except:
                            if props.verbose: print(f"Failed to get directions for virtual face {init_virtual_face[0].index}/{init_virtual_face[1].index}! Marking faces as trailing triangles...")
                            else:
                                pass
                    
                    if current_quant_index > 3:
                        init_face[process_sequence_attr] = process_sequence
                        if init_face in faces_remain:
                            faces_remain.remove(init_face)
                        if init_face not in faces_done:
                            faces_done.append(init_face)
                        init_face[trailing_tri] = 1
                        init_face.select = True
                        init_face[virtual_quad_attr] = 0
                        for vif in init_virtual_face:
                            vif[process_sequence_attr] = process_sequence
                            if vif in faces_remain:
                                faces_remain.remove(vif)
                            if vif not in faces_done:
                                faces_done.append(vif)
                            vif[trailing_tri] = 1
                            vif.select = True
                            vif[virtual_quad_attr] = 0
                        init_virtual_face = []
                    
                    if current_quant_index <= 3:
                        vif_verts = []
                        for vif in init_virtual_face:
                            for vifv in vif.verts:
                                if vifv not in vif_verts:
                                    vif_verts.append(vifv)
                            vif[edge_u] = init_face_dir_arr[0].index
                            vif[edge_l] = init_face_dir_arr[1].index
                            vif[edge_d] = init_face_dir_arr[2].index
                            vif[edge_r] = init_face_dir_arr[3].index
                            vif[process_sequence_attr] = process_sequence
                        these_corner_verts = get_corner_verts(init_face_dir_arr, vif_verts)
                        initial_face_found = True
                        for tcv in these_corner_verts[0]:
                            if tcv == None:
                                initial_face_found = False
                                if props.verbose == True: print(f"Failed to calculate corners for virtual face {init_virtual_face[0].index}/{init_virtual_face[1].index}. Will pick another initial face.")
                        if initial_face_found == False:
                            for vif in init_virtual_face:
                                if vif in faces_remain:
                                    faces_remain.remove(vif)
                                if vif not in faces_done:
                                    faces_done.append(vif)
                                vif[trailing_tri] = 1
                                vif.select = True
                                vif[virtual_quad_attr] = 0
                            init_virtual_face = []
                else:
                    if_verts = [ifv for ifv in init_face.verts]
                    init_face.select = True
                    ife_edges = []
                    for ife in init_face.edges:
                        ife_edges.append(ife)
                    current_face_center = init_face.calc_center_median()
                    init_face_dir_arr = get_init_dir_edges(ife_edges)
                    init_face[edge_u] = init_face_dir_arr[0].index
                    init_face[edge_l] = init_face_dir_arr[1].index
                    init_face[edge_d] = init_face_dir_arr[2].index
                    init_face[edge_r] = init_face_dir_arr[3].index
                    init_face[process_sequence_attr] = process_sequence
                    these_corner_verts = get_corner_verts(init_face_dir_arr, if_verts)
                    initial_face_found = True
                    for tcv in these_corner_verts[0]:
                        if tcv == None:
                            initial_face_found = False
                    if initial_face_found == False:
                        if init_face in faces_remain:
                            faces_remain.remove(init_face)
                        if init_face not in faces_done:
                            faces_done.append(init_face)
                        init_face[trailing_tri] = 1
                        init_face.select = True
                        init_face[virtual_quad_attr] = 0
                        init_virtual_face = []
            
            bpy.context.active_object.data.update()

            next_round = []
            init_face_corner_verts = []
            init_virtual_loops = []
            # endregion

            # region proc init face(s)
            if initial_face_found == True:
                vif_verts = []
                for vif in init_virtual_face:
                    vif[xi] = 0
                    vif[yi] = 0
                    if vif in faces_remain:
                        faces_remain.remove(vif)
                    for vifv in vif.verts:
                        if vifv not in vif_verts:
                            vif_verts.append(vifv)
                    for vifl in vif.loops:
                        if vifl not in init_virtual_loops:
                            init_virtual_loops.append(vifl)
                init_face_corner_verts = get_corner_verts(init_face_dir_arr, vif_verts)
                if props.verbose == True:
                    print(f"There are {len(init_virtual_face)} faces in the initial virtual face. Their indices:")
                    for ivf in init_virtual_face:
                        print(ivf.index)
                for vif in init_virtual_face:
                    try:
                        vif[corner_ur] = init_face_corner_verts[0][0].index
                        vif[corner_ul] = init_face_corner_verts[0][1].index
                        vif[corner_ll] = init_face_corner_verts[0][2].index
                        vif[corner_lr] = init_face_corner_verts[0][3].index
                    except:
                        if props.verbose == True: print(f"Failed to calculate corner vertices for face {vif.index}. This should not happen...")
                        if vif in init_virtual_face:
                            init_virtual_face.remove(vif)
                        pass
                    vif[xi] = 0
                    vif[yi] = 0

                avg_edge_length_x = (init_face_dir_arr[0].calc_length() + init_face_dir_arr[2].calc_length()) / 2
                avg_edge_length_y = (init_face_dir_arr[1].calc_length() + init_face_dir_arr[3].calc_length()) / 2
                for vif in init_virtual_face:
                    vif[avg_length_x] = avg_edge_length_x
                    vif[avg_length_y] = avg_edge_length_y
                
                # region init uv edit
                if self.pre_calc_edge_lengths == False:
                    avg_edge_length_x = (init_face_dir_arr[0].calc_length() + init_face_dir_arr[2].calc_length()) / 4
                    avg_edge_length_y = (init_face_dir_arr[1].calc_length() + init_face_dir_arr[3].calc_length()) / 4
                    
                    offset_co = mathutils.Vector((previous_island_end_x, previous_island_end_y))

                    init_face_corner_vert_cos = [mathutils.Vector((-avg_edge_length_x, avg_edge_length_y)) + offset_co, mathutils.Vector((-avg_edge_length_x, -avg_edge_length_y)) + offset_co, mathutils.Vector((avg_edge_length_x, -avg_edge_length_y)) + offset_co, mathutils.Vector((avg_edge_length_x, avg_edge_length_y)) + offset_co]
                    for ifvco, ifva in zip(init_face_corner_vert_cos, co_attr_arr):
                        for vif in init_virtual_face:
                            vif[ifva[0]] = ifvco.x
                            vif[ifva[1]] = ifvco.y

                    for iv, ivco in zip(init_face_corner_verts[0], init_face_corner_vert_cos):
                        for ivl in iv.link_loops:
                            if ivl.index in island_loops:
                                should_edit = True
                                if self.only_move_loops_in_face == True:
                                    should_edit = ivl in init_virtual_loops
                                if should_edit == True:
                                    ivl[uv_layer].uv = ivco
                                    new_island_max_x = max(island_max_x, ivco.x)
                                    island_max_x = new_island_max_x
                                    new_island_min_x = min(island_min_x, ivco.x)
                                    island_min_x = new_island_min_x
                                    new_island_max_y = max(island_max_y, ivco.y)
                                    island_max_y = new_island_max_y
                                    new_island_min_y = min(island_min_y, ivco.y)
                                    island_min_y = new_island_min_y
                                    ivl[uv_layer].pin_uv = True
                        if iv not in verts_done:
                            verts_done.append(iv)
                # endregion
                
                # region init pre-calc
                if self.pre_calc_edge_lengths == True:
                    this_edge_length_map_x_item = [init_face, init_face_is_virtual, init_virtual_face, [init_face_corner_verts[0][1],init_face_corner_verts[0][2]], [init_face_corner_verts[0][0],init_face_corner_verts[0][3]]]
                    this_edge_length_map_y_item = [init_face, init_face_is_virtual, init_virtual_face, [init_face_corner_verts[0][2],init_face_corner_verts[0][3]], [init_face_corner_verts[0][0],init_face_corner_verts[0][1]]]
                    map_x_index = 0
                    map_y_index = 0
                    try:
                        map_x_index = edge_length_map_x[0].index(0)
                    except:
                        map_x_index = len(edge_length_map_x[0])
                        edge_length_map_x[0].append(0)
                        edge_length_map_x[1].append([0, [], 0.0, 0])
                    edge_length_map_x[1][map_x_index][1].append(this_edge_length_map_x_item)
                    edge_length_map_x[1][map_x_index][2] += avg_edge_length_x
                    edge_length_map_x[1][map_x_index][3] += 1

                    try:
                        map_y_index = edge_length_map_y[0].index(0)
                    except:
                        map_y_index = len(edge_length_map_y[0])
                        edge_length_map_y[0].append(0)
                        edge_length_map_y[1].append([0, [], 0.0, 0])
                    edge_length_map_y[1][map_y_index][1].append(this_edge_length_map_y_item)
                    edge_length_map_y[1][map_y_index][2] += avg_edge_length_y
                    edge_length_map_y[1][map_y_index][3] += 1
                # endregion
                
                # region init ind uv attr
                for ci, corner in enumerate(init_face_corner_verts[0]):
                    for corner_loop in corner.link_loops:
                        if corner_loop in init_virtual_loops:
                            corner_loop[loop_corner_index] = ci
                            if ci == 0:
                                corner_loop[ind_face_uv][0] = 1.0
                                corner_loop[ind_face_uv][1] = 1.0
                            if ci == 1:
                                corner_loop[ind_face_uv][0] = 0.0
                                corner_loop[ind_face_uv][1] = 1.0
                            if ci == 2:
                                corner_loop[ind_face_uv][0] = 0.0
                                corner_loop[ind_face_uv][1] = 0.0
                            if ci == 3:
                                corner_loop[ind_face_uv][0] = 1.0
                                corner_loop[ind_face_uv][1] = 0.0
                # endregion
                
                # region init ind uv orig
                for ci, corner in zip(init_face_corner_verts[0], init_face_corner_verts[0]):
                    for corner_loop in corner.link_loops:
                        if corner_loop in init_virtual_loops:
                            if ci == 0:
                                corner_loop[ind_face_uv_orig][0] = 1.0
                                corner_loop[ind_face_uv_orig][1] = 1.0
                            if ci == 1:
                                corner_loop[ind_face_uv_orig][0] = 0.0
                                corner_loop[ind_face_uv_orig][1] = 1.0
                            if ci == 2:
                                corner_loop[ind_face_uv_orig][0] = 0.0
                                corner_loop[ind_face_uv_orig][1] = 0.0
                            if ci == 3:
                                corner_loop[ind_face_uv_orig][0] = 1.0
                                corner_loop[ind_face_uv_orig][1] = 0.0
                # endregion
                
                # region init face info
                init_face_winding = check_winding(init_face, init_face_corner_verts)
                if props.verbose == True:
                    if init_face_is_virtual == False:
                        print(f"Processing initial face of island {ii}: {init_face.index}. Not virtual.")
                    if init_face_is_virtual == True:
                        print(f"Processing initial face of island {ii}: {init_virtual_face[0].index}/{init_virtual_face[1].index}. Virtual quad.")
                    try:
                        print(viz_quad(init_face_dir_arr, init_face_corner_verts[0], init_face, init_virtual_face, init_face_is_virtual, init_face_winding))
                    except:
                        print(f"Failed to visualize initial face(s)! Operation will probably fail...")
                # endregion

                if init_face not in faces_done:
                    faces_done.append(init_face)

                # region fill round 0
                for di, ife in enumerate(init_face_dir_arr):
                    dirs_to_do = [0,1,2,3]
                    if self.unfold_mode == "B":
                        dirs_to_do = [1,3]
                    elif self.unfold_mode == "C":
                        dirs_to_do = [0,2]
                    if di in dirs_to_do:
                        for ifef in ife.link_faces:
                            if ifef in faces_remain and ifef not in next_round and ifef != init_face and ifef not in init_virtual_face:
                                adj_face = init_face
                                other_virtual_face = [ifef]
                                other_verts = [ov for ov in ifef.verts]
                                if init_face_is_virtual == True:
                                    for vif in init_virtual_face:
                                        if ife in vif.edges:
                                            adj_face = vif
                                other_dir_edges = get_other_dir_edges(init_face_dir_arr, init_face_corner_verts, adj_face, ifef, di)
                                is_other_face_virtual = ifef[virtual_quad_attr]
                                if is_other_face_virtual == True:
                                    other_virtual_face.append(bm.faces[ifef[vq_other_tri]])
                                    for otv in other_virtual_face[1].verts:
                                        if otv not in other_verts:
                                            other_verts.append(otv)
                                ifef[nom_dir] = di
                                ifef[nom] = adj_face.index

                                if ifef[trailing_tri] == 1:
                                    if props.verbose: print(f"Initial face of island {ii}: skipped face index {ifef.index} in direction {di} because it is a trailing triangle.")
                                    ifef[process_sequence_attr] = process_sequence
                                    ife[process_sequence_edge] = process_sequence
                                    ife[nom_dir_edge] = di
                                    for ifefl in ifef.loops:
                                        ifefl[uv_layer].pin_uv = False
                                    process_sequence += 1
                                    faces_remain.remove(ifef)
                                    ifef.select = True
                                if ifef[trailing_tri] == 0:
                                    other_corner_verts = get_corner_verts(other_dir_edges, other_verts)
                                    skip_append = False

                                    if props.verbose:
                                        print(f"Initial face of island {ii}: Added face index {ifef.index} in direction {di} to next round. Virtual face: {init_face_is_virtual}.")
                                        try:
                                            print(viz_add_quad(init_face_dir_arr, other_dir_edges, init_face_corner_verts[0], other_corner_verts[0], init_face, ifef, init_virtual_face, other_virtual_face, init_face_is_virtual, is_other_face_virtual, di))
                                        except:
                                            print(f"Face index {ifef.index} failed to generate valid corners! Face skipped.")
                                            skip_append = True
                                    
                                    if skip_append == True:
                                        ifef[trailing_tri] = 1
                                        if ifef in faces_remain:
                                            faces_remain.remove(ifef)
                                        if ifef not in faces_done:
                                            faces_done.append(ifef)
                                    if skip_append == False:
                                        next_round.append(ifef)
                                    ifef[process_sequence_attr] = process_sequence
                                    if di == 1 or di == 3:
                                        if di == 1: ifef[xi] = init_face[xi] - 1
                                        if di == 3: ifef[xi] = init_face[xi] + 1
                                        ifef[yi] = init_face[yi]
                                    if di == 0 or di == 2:
                                        ifef[xi] = init_face[xi]
                                        if di == 0: ifef[yi] = init_face[yi] + 1
                                        if di == 2: ifef[yi] = init_face[yi] - 1
                                    ife[process_sequence_edge] = process_sequence
                                    ife[nom_dir_edge] = di
                                    process_sequence += 1
                # endregion

            if props.verbose:
                if initial_face_found == False:
                    if props.verbose: print(f"No usable initial face found for island {ii}.")
            # endregion

            bpy.context.active_object.data.update()
            round = 0
            do_a = True
            retries = 0

            # region traverse mesh
            tododo = len(faces_remain)
            while tododo > 0:
                if props.verbose: print(f"{len(faces_remain)} remaining faces.")
                round += 1
                this_round = next_round.copy()
                next_round = []

                for face in faces_remain:
                    if face[trailing_tri] == 1:
                        faces_remain.remove(face)
                        face.select = True
                        for loop in face.loops:
                            loop[uv_layer].pin_uv = False
                        if props.verbose: print(f"Removed face {face.index} from remaining faces because it is a trailing triangle.")
                for ttf in this_round:
                    if ttf[trailing_tri] == 1:
                        this_round.remove(ttf)
                        ttf.select = True
                        for ttfl in ttf.loops:
                            ttfl[uv_layer].pin_uv = False
                        if props.verbose: print(f"Removed face {ttf.index} from round because it is a trailing triangle.")

                for trf in this_round:
                    # region proc t face
                    if props.verbose: print(f"Island {ii} - Round {round}: Processing face {trf.index} with direction {trf[nom_dir]}.")
                    these_edges = [bm.edges[trf[edge_u]], bm.edges[trf[edge_l]], bm.edges[trf[edge_d]], bm.edges[trf[edge_r]]]

                    is_virtual_quad = False
                    virtual_quad = [trf, None]
                    virtual_verts = []
                    virtual_loops = []
                    corners = []
                    cos = [None, None, None, None]
                    nom_cos = [
                        mathutils.Vector((bm.faces[trf[nom]][co_ur_x], bm.faces[trf[nom]][co_ur_y])),
                        mathutils.Vector((bm.faces[trf[nom]][co_ul_x], bm.faces[trf[nom]][co_ul_y])),
                        mathutils.Vector((bm.faces[trf[nom]][co_ll_x], bm.faces[trf[nom]][co_ll_y])),
                        mathutils.Vector((bm.faces[trf[nom]][co_lr_x], bm.faces[trf[nom]][co_lr_y]))
                        ]

                    if trf in faces_remain:
                        faces_remain.remove(trf)
                    if trf not in faces_done:
                        faces_done.append(trf)

                    trf.select = True

                    if trf[virtual_quad_attr] == 1:
                        trof = bm.faces[trf[vq_other_tri]]
                        trof.select = True
                        virtual_quad[0] = trf
                        virtual_quad[1] = trof
                        is_virtual_quad = True
                        if trof in faces_remain:
                            faces_remain.remove(trof)
                    
                    avg_edge_length_x = (these_edges[0].calc_length() + these_edges[2].calc_length()) / 2
                    avg_edge_length_y = (these_edges[1].calc_length() + these_edges[3].calc_length()) / 2
                    
                    for vf in virtual_quad:
                        if vf != None:
                            vf[avg_length_x] = avg_edge_length_x
                            vf[avg_length_y] = avg_edge_length_y
                            for vfv in vf.verts:
                                if vfv not in virtual_verts:
                                    virtual_verts.append(vfv)
                            for vfl in vf.loops:
                                if vfl not in virtual_loops:
                                    virtual_loops.append(vfl)
                    corners = get_corner_verts(these_edges, virtual_verts)
                    corner_missing = False

                    for corner_i, corner in enumerate(corners[0]):
                        if corner == None:
                            for vf in virtual_quad:
                                if vf != None:
                                    vf[trailing_tri] = 1
                                    vq_is = [None, None]
                                    for vqfi, vqf in enumerate(virtual_quad):
                                        try:
                                            vq_is[vqfi] = vqf.index
                                        except:
                                            vq_is.pop()
                                    if props.verbose == True and corner_missing == False:
                                        print(f"Virtual quad corner {corner_i} failed to compute! Skipping virtual quad {vq_is}...")
                                    corner_missing = True
                    
                    if corner_missing == False:
                        for ci, corner in enumerate(corners[0]):
                            for corner_loop in corner.link_loops:
                                if corner_loop in virtual_loops:
                                    corner_loop[loop_corner_index] = ci
                                    if ci == 0:
                                        corner_loop[ind_face_uv][0] = 1.0
                                        corner_loop[ind_face_uv][1] = 1.0
                                    if ci == 1:
                                        corner_loop[ind_face_uv][0] = 0.0
                                        corner_loop[ind_face_uv][1] = 1.0
                                    if ci == 2:
                                        corner_loop[ind_face_uv][0] = 0.0
                                        corner_loop[ind_face_uv][1] = 0.0
                                    if ci == 3:
                                        corner_loop[ind_face_uv][0] = 1.0
                                        corner_loop[ind_face_uv][1] = 0.0
                        
                        for ci, corner in zip(corners[1], corners[0]):
                            for corner_loop in corner.link_loops:
                                if corner_loop in virtual_loops:
                                    if ci == 0:
                                        corner_loop[ind_face_uv_orig][0] = 1.0
                                        corner_loop[ind_face_uv_orig][1] = 1.0
                                    if ci == 1:
                                        corner_loop[ind_face_uv_orig][0] = 0.0
                                        corner_loop[ind_face_uv_orig][1] = 1.0
                                    if ci == 2:
                                        corner_loop[ind_face_uv_orig][0] = 0.0
                                        corner_loop[ind_face_uv_orig][1] = 0.0
                                    if ci == 3:
                                        corner_loop[ind_face_uv_orig][0] = 1.0
                                        corner_loop[ind_face_uv_orig][1] = 0.0
                        
                        if props.verbose:
                            if is_virtual_quad == False:
                                wnd = check_winding(trf, corners)
                                print(viz_quad(these_edges, corners[0], trf, virtual_quad, is_virtual_quad, wnd))
                            else:
                                print(viz_quad(these_edges, corners[0], trf, virtual_quad, is_virtual_quad, ""))

                        if trf[nom_dir] == 0:
                            cos[0] = nom_cos[0] + mathutils.Vector((0,avg_edge_length_y))
                            cos[1] = nom_cos[1] + mathutils.Vector((0,avg_edge_length_y))
                            cos[2] = nom_cos[1]
                            cos[3] = nom_cos[0]
                        elif trf[nom_dir] == 1:
                            cos[0] = nom_cos[1]
                            cos[1] = nom_cos[1] + mathutils.Vector((-avg_edge_length_x,0))
                            cos[2] = nom_cos[2] + mathutils.Vector((-avg_edge_length_x,0))
                            cos[3] = nom_cos[2]
                        elif trf[nom_dir] == 2:
                            cos[0] = nom_cos[3]
                            cos[1] = nom_cos[2]
                            cos[2] = nom_cos[2] + mathutils.Vector((0,-avg_edge_length_y))
                            cos[3] = nom_cos[3] + mathutils.Vector((0,-avg_edge_length_y))
                        elif trf[nom_dir] == 3:
                            cos[0] = nom_cos[1] + mathutils.Vector((avg_edge_length_x,0))
                            cos[1] = nom_cos[0]
                            cos[2] = nom_cos[3]
                            cos[3] = nom_cos[3] + mathutils.Vector((avg_edge_length_x,0))
                        # endregion

                        # region t uv edit
                        if self.pre_calc_edge_lengths == False:
                            for vi, cornerv in enumerate(corners[0]):
                                for vil in cornerv.link_loops:
                                    if vil.index in island_loops:
                                        should_edit = True
                                        if self.only_move_loops_in_face == True:
                                            should_edit = vil in virtual_loops
                                        if should_edit == True:
                                            vil[uv_layer].uv = cos[vi]
                                            new_island_max_x = max(island_max_x, cos[vi].x)
                                            island_max_x = new_island_max_x
                                            new_island_min_x = min(island_min_x, cos[vi].x)
                                            island_min_x = new_island_min_x
                                            new_island_max_y = max(island_max_y, cos[vi].y)
                                            island_max_y = new_island_max_y
                                            new_island_min_y = min(island_min_y, cos[vi].y)
                                            island_min_y = new_island_min_y
                                            vil[uv_layer].pin_uv = True
                                if cornerv not in verts_done:
                                    verts_done.append(cornerv)
                            for trfvco, trfva in zip(cos, co_attr_arr):
                                for vif in virtual_quad:
                                    if vif != None:
                                        vif[trfva[0]] = trfvco.x
                                        vif[trfva[1]] = trfvco.y
                        # endregion
                        
                        # region t pre-calc prep
                        if self.pre_calc_edge_lengths == True:
                            this_edge_length_map_x_item = [trf, is_virtual_quad, virtual_quad, [corners[0][1],corners[0][2]], [corners[0][0],corners[0][3]]]
                            this_edge_length_map_y_item = [trf, is_virtual_quad, virtual_quad, [corners[0][2],corners[0][3]], [corners[0][0],corners[0][1]]]
                            map_x_index = 0
                            map_y_index = 0
                            try:
                                map_x_index = edge_length_map_x[0].index(trf[xi])
                            except:
                                map_x_index = len(edge_length_map_x[0])
                                edge_length_map_x[0].append(trf[xi])
                                edge_length_map_x[1].append([trf[xi], [], 0.0, 0])
                            edge_length_map_x[1][map_x_index][1].append(this_edge_length_map_x_item)
                            edge_length_map_x[1][map_x_index][2] += avg_edge_length_x
                            edge_length_map_x[1][map_x_index][3] += 1

                            try:
                                map_y_index = edge_length_map_y[0].index(trf[yi])
                            except:
                                map_y_index = len(edge_length_map_y[0])
                                edge_length_map_y[0].append(trf[yi])
                                edge_length_map_y[1].append([trf[yi], [], 0.0, 0])
                            edge_length_map_y[1][map_y_index][1].append(this_edge_length_map_y_item)
                            edge_length_map_y[1][map_y_index][2] += avg_edge_length_y
                            edge_length_map_y[1][map_y_index][3] += 1
                        # endregion
                        
                        # region fill next round
                        if self.unfold_mode == "A":
                            dirs_to_do =  [0,1,2,3]
                            for di, ife in enumerate(these_edges):
                                if di in dirs_to_do:
                                    for ifef in ife.link_faces:
                                        if ifef in faces_remain and ifef not in next_round and ifef != trf:
                                            adj_face = trf
                                            other_verts = [ifefv for ifefv in ifef.verts]
                                            if is_virtual_quad == True:
                                                for vif in virtual_quad:
                                                    if ife in vif.edges:
                                                        adj_face = vif
                                            other_dir_arr = get_other_dir_edges(these_edges, corners, adj_face, ifef, di)
                                            ifef[nom_dir] = di
                                            ifef[nom] = adj_face.index
                                            process_sequence += 1
                                            ifef[process_sequence_attr] = process_sequence
                                            other_virtual_face = [ifef, bm.faces[ifef[vq_other_tri]]]
                                            other_face_is_virtual = ifef[virtual_quad_attr]
                                            if other_face_is_virtual == True:
                                                for ovf in other_virtual_face:
                                                    for ovfv in ovf.verts:
                                                        if ovfv not in other_verts:
                                                            other_verts.append(ovfv)
                                            if ifef[trailing_tri] == 1:
                                                if props.verbose: print(f"Island {ii} - Round {round}: Skipped face index {ifef.index} in direction {di} because it is a trailing triangle.")
                                                ifef[process_sequence_attr] = process_sequence
                                                ife[process_sequence_edge] = process_sequence
                                                ife[nom_dir_edge] = di
                                                process_sequence += 1
                                                faces_remain.remove(ifef)
                                                ifef.select = True
                                                for ifefl in ifef.loops:
                                                    ifefl[uv_layer].pin_uv = False
                                            if ifef[trailing_tri] == 0:
                                                skip_append = False
                                                try:
                                                    other_corner_verts = get_corner_verts(other_dir_arr, other_verts)
                                                except:
                                                    skip_append = True
                                                if props.verbose:
                                                    print(f"Island {ii} - Round {round}: Added face index {ifef.index} in direction {di} to next round of island {ii}. Virtual face: {other_face_is_virtual}. Current process sequence index: {process_sequence}")
                                                    try:
                                                        print(viz_add_quad(these_edges, other_dir_arr, corners[0], other_corner_verts[0], adj_face, ifef, virtual_quad, other_virtual_face, is_virtual_quad, other_face_is_virtual, di))
                                                    except:
                                                        print(f"Failed to calculate corners for face {ifef.index}! Operation will probably fail...")
                                                if skip_append == True:
                                                    ifef[trailing_tri] = 1
                                                    if ifef in faces_remain:
                                                        faces_remain.remove(ifef)
                                                    if ifef not in faces_done:
                                                        faces_done.append(ifef)
                                                if skip_append == False:
                                                    next_round.append(ifef)
                                                ifef[process_sequence_attr] = process_sequence
                                                ife[process_sequence_edge] = process_sequence
                                                ife[nom_dir_edge] = di
                                                if di == 1 or di == 3:
                                                    if di == 1: ifef[xi] = adj_face[xi] - 1
                                                    if di == 3: ifef[xi] = adj_face[xi] + 1
                                                    ifef[yi] = adj_face[yi]
                                                if di == 0 or di == 2:
                                                    ifef[xi] = adj_face[xi]
                                                    if di == 0: ifef[yi] = adj_face[yi] + 1
                                                    if di == 2: ifef[yi] = adj_face[yi] - 1
                                                process_sequence += 1

                        if self.unfold_mode == "B" or self.unfold_mode == "C":
                            dirs_to_do = [1, 3]
                            dirs_to_do_a = [1, 3]
                            dirs_to_do_b = [0, 2]
                            if self.unfold_mode == "C":
                                dirs_to_do_a = [0, 2]
                                dirs_to_do_b = [1, 3]

                            if do_a == True:
                                dirs_to_do = dirs_to_do_a
                            if do_a == False:
                                dirs_to_do = dirs_to_do_b

                            for di, ife in enumerate(these_edges):
                                if di in dirs_to_do:
                                    for ifef in ife.link_faces:
                                        if ifef in faces_remain and ifef not in next_round and ifef != trf:
                                            adj_face = trf
                                            other_verts = [ifefv for ifefv in ifef.verts]
                                            if is_virtual_quad == True:
                                                for vif in virtual_quad:
                                                    if ife in vif.edges:
                                                        adj_face = vif
                                            try:
                                                other_dir_arr = get_other_dir_edges(these_edges, corners, adj_face, ifef, di)
                                            except:
                                                ifef[trailing_tri] = 1
                                            ifef[nom_dir] = di
                                            ifef[nom] = adj_face.index
                                            ifef[process_sequence_attr] = process_sequence
                                            other_virtual_face = [ifef, bm.faces[ifef[vq_other_tri]]]
                                            other_face_is_virtual = ifef[virtual_quad_attr]
                                            if other_face_is_virtual == True:
                                                for ovf in other_virtual_face:
                                                    for ovfv in ovf.verts:
                                                        if ovfv not in other_verts:
                                                            other_verts.append(ovfv)
                                            if ifef[trailing_tri] == 1:
                                                if props.verbose: print(f"Island {ii} - Round {round}: Skipped face index {ifef.index} in direction {di} because it is a trailing triangle.")
                                                ifef[process_sequence_attr] = process_sequence
                                                ife[process_sequence_edge] = process_sequence
                                                ife[nom_dir_edge] = di
                                                process_sequence += 1
                                                faces_remain.remove(ifef)
                                                ifef.select = True
                                                for ifefl in ifef.loops:
                                                    ifefl[uv_layer].pin_uv = False
                                            if ifef[trailing_tri] == 0:
                                                skip_append = False
                                                try:
                                                    other_corner_verts = get_corner_verts(other_dir_arr, other_verts)
                                                except:
                                                    skip_append = True
                                                if props.verbose:
                                                    print(f"Island {ii} - Round {round}: Added face index {ifef.index} in direction {di} to next round of island {ii}. Virtual face: {other_face_is_virtual}. Current process sequence index: {process_sequence}")
                                                    try:
                                                        print(viz_add_quad(these_edges, other_dir_arr, corners[0], other_corner_verts[0], adj_face, ifef, virtual_quad, other_virtual_face, is_virtual_quad, other_face_is_virtual, di))
                                                    except:
                                                        print(f"Failed to calculate corners for face {ifef.index}! Operation will probably fail...")
                                                if skip_append == True:
                                                    ifef[trailing_tri] = 1
                                                    if ifef in faces_remain:
                                                        faces_remain.remove(ifef)
                                                    if ifef not in faces_done:
                                                        faces_done.append(ifef)
                                                if skip_append == False:
                                                    next_round.append(ifef)
                                                ifef[process_sequence_attr] = process_sequence
                                                ife[process_sequence_edge] = process_sequence
                                                ife[nom_dir_edge] = di
                                                if di == 1 or di == 3:
                                                    if di == 1: ifef[xi] = adj_face[xi] - 1
                                                    if di == 3: ifef[xi] = adj_face[xi] + 1
                                                    ifef[yi] = adj_face[yi]
                                                if di == 0 or di == 2:
                                                    ifef[xi] = adj_face[xi]
                                                    if di == 0: ifef[yi] = adj_face[yi] + 1
                                                    if di == 2: ifef[yi] = adj_face[yi] - 1
                                                process_sequence += 1
                    # endregion
                
                if len(next_round) == 0:
                    if self.unfold_mode == "B" or self.unfold_mode == "C":
                        if do_a == True:
                            do_a = False
                            if props.verbose: print("No more in main axis, trying secondary axis")
                        elif do_a == False:
                            do_a = True
                            if props.verbose: print("No more in secondary axis, returning to main axis")
                        next_round = faces_done.copy()
                    else:
                        if len(faces_remain) > 0:
                            next_round.append(faces_remain[0])
                        retries += 1

                tododo = len(faces_remain)

                bpy.context.active_object.data.update()
            
            # endregion

            if self.pre_calc_edge_lengths == False:
                offset_add_x = island_max_x - island_min_x
                if props.verbose == True: print(f"Island {ii} width: {offset_add_x}")
                offset_add_y = island_max_y - island_min_y
                if props.verbose == True: print(f"Island {ii} height: {offset_add_y}")
                previous_island_end_x += (self.offset_per_island[0] * offset_add_x)
                previous_island_end_y += (self.offset_per_island[1] * offset_add_y)

            # region apply pre-calc
            if self.pre_calc_edge_lengths == True and (len(edge_length_map_x[0]) + len(edge_length_map_y[0])) > 0:
                sorted_x = edge_length_map_x[1].copy()
                sorted_x.sort(key=map_sort_by_coord)
                sorted_y = edge_length_map_y[1].copy()
                sorted_y.sort(key=map_sort_by_coord)

                list_a = sorted_x
                list_b = sorted_y

                if self.unfold_mode == "C":
                    list_a = sorted_y
                    list_b = sorted_x

                last_a_len = previous_island_end_x
                last_b_len = previous_island_end_y
                if self.unfold_mode == "C":
                    last_a_len = previous_island_end_y
                    last_b_len = previous_island_end_x
                
                this_total_a_len = 0.0
                this_total_b_len = 0.0

                a_name = "column"
                b_name = "row"
                if self.unfold_mode == "C":
                    a_name = "row"
                    b_name = "column"
                
                if props.verbose: print(f"Island {ii} has {len(list_a)} {a_name}s and {len(list_b)} {b_name}s ")

                for list_a_item in list_a:
                    a_len = list_a_item[2]/list_a_item[3]
                    old_a_len = 0.0 + last_a_len
                    new_a_len = 0.0 + last_a_len + a_len
                    if props.verbose: print(f"Average edge length for {a_name} {list_a_item[0]}: {a_len}\nTotal previous length: {last_a_len}")
                    for face_data in list_a_item[1]:
                        for alc in face_data[3]:
                            for alcl in alc.link_loops:
                                if alcl.index in island_loops:
                                    should_edit = True
                                    if self.only_move_loops_in_face == True:
                                        should_edit = alcl in virtual_loops
                                    if should_edit == True:
                                        if self.unfold_mode == "C":
                                            alcl[uv_layer].uv.y = old_a_len
                                        else:
                                            alcl[uv_layer].uv.x = old_a_len
                        for amc in face_data[4]:
                            for amcl in amc.link_loops:
                                if amcl.index in island_loops:
                                    should_edit = True
                                    if self.only_move_loops_in_face == True:
                                        should_edit = amcl in virtual_loops
                                    if should_edit == True:
                                        if self.unfold_mode == "C":
                                            amcl[uv_layer].uv.y = new_a_len
                                        else:
                                            amcl[uv_layer].uv.x = new_a_len
                    last_a_len += a_len
                    this_total_a_len += a_len
                
                bpy.context.active_object.data.update()

                for list_b_item in list_b:
                    b_len = list_b_item[2]/list_b_item[3]
                    old_b_len = 0.0 + last_b_len
                    new_b_len = 0.0 + last_b_len + b_len
                    if props.verbose: print(f"Average edge length for {b_name} {list_b_item[0]}: {b_len}\nTotal previous length: {last_b_len}")
                    for face_data in list_b_item[1]:
                        for blc in face_data[3]:
                            for blcl in blc.link_loops:
                                if blcl.index in island_loops:
                                    should_edit = True
                                    if self.only_move_loops_in_face == True:
                                        should_edit = blcl in virtual_loops
                                    if should_edit == True:
                                        if self.unfold_mode == "C":
                                            blcl[uv_layer].uv.x = old_b_len
                                        else:
                                            blcl[uv_layer].uv.y = old_b_len
                                        blcl[uv_layer].pin_uv = True
                        for bmc in face_data[4]:
                            for bmcl in bmc.link_loops:
                                if bmcl.index in island_loops:
                                    should_edit = True
                                    if self.only_move_loops_in_face == True:
                                        should_edit = bmcl in virtual_loops
                                    if should_edit == True:
                                        if self.unfold_mode == "C":
                                            bmcl[uv_layer].uv.x = new_b_len
                                        else:
                                            bmcl[uv_layer].uv.y = new_b_len
                                        bmcl[uv_layer].pin_uv = True
                    last_b_len += b_len
                    this_total_b_len += b_len
                
                if self.unfold_mode == "C":
                    previous_island_end_x += (self.offset_per_island[0] * this_total_b_len)
                    previous_island_end_y += (self.offset_per_island[1] * this_total_a_len)
                else:
                    previous_island_end_x += (self.offset_per_island[0] * this_total_a_len)
                    previous_island_end_y += (self.offset_per_island[1] * this_total_b_len)
                bpy.context.active_object.data.update()
            # endregion

            for face in island:
                face[process_sequence_max] = process_sequence
        
        bpy.context.active_object.data.update()

        try:
            bpy.ops.uv.unwrap(method='MINIMUM_STRETCH', fill_holes=True, correct_aspect=True, use_subsurf_data=False, margin=0, no_flip=False, iterations=10, use_weights=False, weight_group="uv_importance", weight_factor=1)
        except:
            bpy.ops.uv.unwrap(method='ANGLE_BASED', fill_holes=True, correct_aspect=True, use_subsurf_data=False, margin=0, no_flip=False, iterations=10, use_weights=False, weight_group="uv_importance", weight_factor=1)

        bpy.context.active_object.data.update()

        for face in bm.faces:
            for loop in face.loops:
                loop[uv_layer].pin_uv = False
        return {"FINISHED"}

class SloppyUVToMesh(bpy.types.Operator):
    bl_idname = "operator.sloppy_uv_to_mesh"
    bl_label = "Sloppy UV to Mesh"
    bl_description = "Cut mesh by the seams and flatten out the pieces in the shape of their UVs"
    bl_options = {'REGISTER', 'UNDO'}

    iP = bpy.props.IntProperty
    fP = bpy.props.FloatProperty
    fvP = bpy.props.FloatVectorProperty
    bP = bpy.props.BoolProperty
    eP = bpy.props.EnumProperty
    bvP = bpy.props.BoolVectorProperty
    sP = bpy.props.StringProperty

    # region prop def
    
    #endregion

    def execute(self, context):
        props = context.scene.sloppy_props
        bm = bmesh.from_edit_mesh(bpy.context.active_object.data)
        uv_layer = bm.loops.layers.uv.verify()
        islands = bmesh_utils.bmesh_linked_uv_islands(bm, uv_layer)
        
        # region init vars
        attr_dict = [
            {"name": "loop_orig_co", "type": "FLOAT_VECTOR", "domain": "CORNER", "layer": None},
            {"name": "v_orig_co", "type": "FLOAT_VECTOR", "domain": "POINT", "layer": None},
            {"name": "loop_uv", "type": "FLOAT_VECTOR", "domain": "CORNER", "layer": None},
            {"name": "v_uv", "type": "FLOAT_VECTOR", "domain": "POINT", "layer": None},
            {"name": "loop_co_diff", "type": "FLOAT_VECTOR", "domain": "CORNER", "layer": None},
            {"name": "v_co_diff", "type": "FLOAT_VECTOR", "domain": "POINT", "layer": None}
            ]
        
        for nam in attr_dict:
            attr = props.find_or_add_attribute(nam["name"], nam["type"], nam["domain"])
            nam["layer"] = props.get_attribute_layer(nam["name"], nam["type"], nam["domain"], bm)
        
        loop_orig_co = props.get_dict_layer("loop_orig_co", attr_dict)
        v_orig_co = props.get_dict_layer("v_orig_co", attr_dict)
        loop_uv = props.get_dict_layer("loop_uv", attr_dict)
        v_uv = props.get_dict_layer("v_uv", attr_dict)
        loop_co_diff = props.get_dict_layer("loop_co_diff", attr_dict)
        v_co_diff = props.get_dict_layer("v_co_diff", attr_dict)

        # endregion

        bm.verts.ensure_lookup_table()
        bm.edges.ensure_lookup_table()
        bm.faces.ensure_lookup_table()

        selected_edges = []

        for vert in bm.verts:
            vert.select = False
        for edge in bm.edges:
            if edge.seam == True:
                edge.select = edge.seam
                selected_edges.append(edge)
            else:
                edge.select = False
        for face in bm.faces:
            face.select = False

        bpy.context.active_object.data.update()

        bmesh.ops.split_edges(bm, edges=selected_edges)

        bpy.context.active_object.data.update()

        bm.verts.ensure_lookup_table()
        bm.edges.ensure_lookup_table()
        bm.faces.ensure_lookup_table()

        for face in bm.faces:
            for loop in face.loops:
                uvco = mathutils.Vector((loop[uv_layer].uv.x, loop[uv_layer].uv.y, 0))
                uvdiff = uvco - loop.vert.co
                loop.vert[v_uv] = uvco
                loop.vert[v_co_diff] = uvdiff

        bpy.context.active_object.data.update()

        return {"FINISHED"}

class RedoUVEdgeLength(bpy.types.Operator):
    bl_idname = "operator.sloppy_rebuild_uv_edge_length"
    bl_label = "Rebuild UV Edge Length"
    bl_description = "Rebuild edge lengths of selected (preferably straight) UV edges"
    bl_options = {'REGISTER', 'UNDO'}

    iP = bpy.props.IntProperty
    fP = bpy.props.FloatProperty
    fvP = bpy.props.FloatVectorProperty
    bP = bpy.props.BoolProperty
    eP = bpy.props.EnumProperty
    bvP = bpy.props.BoolVectorProperty
    sP = bpy.props.StringProperty

    # region prop def
    use_world_scale : bP(
        name = "Use World Scale",
        description = "Rebuild edge lengths in World-UV 1:1 scale",
        default = False
        ) # type: ignore
    
    transform_pivot : eP(
        name = "Pivot Point",
        description = "Pivot point for rebuilding UV edge",
        items = [
            ("A", "Center", "Relative to averaged center of UV selection"),
            ("B", "Cursor", "Using UV cursor")
            ],
        ) # type: ignore

    shift_pivot : fP(
        name = "Shift Pivot Point",
        description = "Shift pivot point along selected edges",
        default = 0.5,
        min = 0.0,
        max = 1.0
        ) # type: ignore

    #endregion

    def execute(self, context):
        props = context.scene.sloppy_props
        bm = bmesh.from_edit_mesh(bpy.context.active_object.data)
        uv_layer = bm.loops.layers.uv.verify()
        islands = bmesh_utils.bmesh_linked_uv_islands(bm, uv_layer)
        remapped_shift = (self.shift_pivot * 2) - 1
        print('Remapped shift:', remapped_shift)
        
        if len(islands) < 1:
            print('No islands, using all faces.')
            islands = [[face for face in bm.faces]]

        bm.verts.ensure_lookup_table()
        bm.edges.ensure_lookup_table()
        bm.faces.ensure_lookup_table()

        curpos = None
        for area in bpy.context.screen.areas:
            if area.type == 'IMAGE_EDITOR':   #find the UVeditor
                curpos = area.spaces.active.cursor_location

        for ii, island in enumerate(islands):
            island_loops = []
            for face in island:
                for loop in face.loops:
                    if loop not in island_loops:
                        island_loops.append(loop)
            print('Loops in island', ii, ':', len(island_loops))
            selected_uv_loops = []
            min_x = 9999.0
            min_y = 9999.0
            max_x = -9990.0
            max_y = -9999.0
            for loop in island_loops:
                if loop[uv_layer].select == True:
                    loop_bundle = [loop, loop[uv_layer].uv]
                    numin_x = min(min_x, loop[uv_layer].uv.x)
                    min_x = numin_x
                    numin_y = min(min_y, loop[uv_layer].uv.y)
                    min_y = numin_y
                    numax_x = max(max_x, loop[uv_layer].uv.x)
                    max_x = numax_x
                    numax_y = max(max_y, loop[uv_layer].uv.y)
                    max_y = numax_y
                    selected_uv_loops.append(loop_bundle)
            x_diff = max_x - min_x
            y_diff = max_y - min_y
            print('Number of selected loops in island', ii, ':', len(selected_uv_loops))
            selected_uv_loops_only = [loop for loop in island_loops if loop[uv_layer].select == True]
            selected_uv_loops_edge = [loop for loop in island_loops if loop[uv_layer].select_edge == True]
            selected_verts = [loop[0].vert for loop in selected_uv_loops]
            selected_edges = [edge for edge in bm.edges if edge.verts[0] in selected_verts and edge.verts[1] in selected_verts]
            if x_diff > y_diff:
                selected_uv_loops.sort(key=props.uv_sort_x)
            if x_diff < y_diff:
                selected_uv_loops.sort(key=props.uv_sort_y)
            if len(selected_uv_loops) > 0:
                min_vert = selected_uv_loops[0][0].vert
                max_vert = selected_uv_loops[-1][0].vert
                min_uv = selected_uv_loops[0][1]
                max_uv = selected_uv_loops[-1][1]
                center_uv = min_uv.lerp(max_uv, 0.5)
                print('UV Center:', center_uv.x, center_uv.y)
                min_to_max_vec = max_uv - min_uv
                min_to_max_dir = min_to_max_vec.normalized()
                print('Direction:', min_to_max_dir)
                total_uv_length = (max_uv - min_uv).length
                print('Total UV length:', total_uv_length)
                total_edge_length = 0.0
                for edge in selected_edges:
                    total_edge_length += edge.calc_length()
                print('Total world scale edge length:', total_edge_length)
                world_to_uv_ratio = total_edge_length/total_uv_length
                print('World to UV scale ratio: 1:', world_to_uv_ratio)
                last_vert = min_vert
                length_to_here = 0.0
                min_to_max_lerp = 0.0
                for uvl in selected_uv_loops:
                    this_vert = uvl[0].vert
                    if this_vert != last_vert:
                        this_edge = None
                        for edge in this_vert.link_edges:
                            if edge in last_vert.link_edges:
                                this_edge = edge
                        try:
                            length_to_here += this_edge.calc_length()
                        except:
                            pass
                    min_to_max_lerp = length_to_here/total_edge_length
                    last_vert = this_vert
                    new_uv_calc = uvl[1]
                    if self.use_world_scale == False:
                        if self.transform_pivot == "A":
                            nuv_min = center_uv - (min_to_max_dir * (total_uv_length / 2))
                            nuv_max = center_uv + (min_to_max_dir * (total_uv_length / 2))
                            new_uv_calc = nuv_min.lerp(nuv_max, min_to_max_lerp) + (min_to_max_dir * ((total_uv_length / 2) * remapped_shift))
                        if self.transform_pivot == "B":
                            nuv_min = curpos - (min_to_max_dir * (total_uv_length / 2))
                            nuv_max = curpos + (min_to_max_dir * (total_uv_length / 2))
                            new_uv_calc = nuv_min.lerp(nuv_max, min_to_max_lerp) + (min_to_max_dir * ((total_uv_length / 2) * remapped_shift))
                    if self.use_world_scale == True:
                        if self.transform_pivot == "A":
                            wuv_min = center_uv - (min_to_max_dir * (total_edge_length / 2))
                            wuv_max = center_uv + (min_to_max_dir * (total_edge_length / 2))
                            new_uv_calc = wuv_min.lerp(wuv_max, min_to_max_lerp) + (min_to_max_dir * ((total_edge_length / 2) * remapped_shift))
                        if self.transform_pivot == "B":
                            wuv_min = curpos - (min_to_max_dir * (total_edge_length / 2))
                            wuv_max = curpos + (min_to_max_dir * (total_edge_length / 2))
                            new_uv_calc = wuv_min.lerp(wuv_max, min_to_max_lerp) + (min_to_max_dir * ((total_edge_length / 2) * remapped_shift))
                    for loop in this_vert.link_loops:
                        if loop in selected_uv_loops_only and loop in island_loops:
                            loop[uv_layer].uv = new_uv_calc


        bpy.context.active_object.data.update()

        return {"FINISHED"}

class SloppyBasicUVUnfold(bpy.types.Operator):
    bl_idname = "operator.sloppy_basic_uv_unfold"
    bl_label = "Basic UV Unfold"
    bl_description = "Assumes that each non-boundary vertex should have a combined angle sum of 360 degrees (as though in a mesh of triangulated quads), working outwards from an initial selected or calculated non-boundary and non-seam edge"
    bl_options = {'REGISTER', 'UNDO'}

    iP = bpy.props.IntProperty
    fP = bpy.props.FloatProperty
    fvP = bpy.props.FloatVectorProperty
    bP = bpy.props.BoolProperty
    eP = bpy.props.EnumProperty
    bvP = bpy.props.BoolVectorProperty
    sP = bpy.props.StringProperty

    # region prop def
    transform_pivot : eP(
        name = "Pivot Point",
        description = "Pivot point for rebuilding UV edge",
        items = [
            ("A", "Init Edge", "Init edge placed at U:0.5 V:0.5"),
            ("B", "Init Edge to Cursor", "Init edge placed at cursor"),
            ("C", "Center", "Island center placed at U:0.5 V:0.5"),
            ("D", "Center to Cursor", "Island center placed at cursor"),
            ("E", "Positive UV Space", "Shift pivot to keep everything in positive UV space")
            ],
        default="E"
        ) # type: ignore
    
    scaling_mode : eP(
        name = "Scaling Mode",
        description = "Scaling of result",
        items = [
            ("A", "World Scale", "Scaled 1:1 World to UV"),
            ("B", "Scale to Bounds", "Scaled uniformly to fit within UV 0-1 using largest dimension"),
            ("C", "Stretch to Bounds", "Scaled non- uniformly to fit within UV 0-1")
            ],
        ) # type: ignore
    
    initial_edge_mode : eP(
        name = "Init Edge",
        description = "Selection of initial edge",
        items = [
            ("A", "Automatic", "Relative to averaged center of UV selection"),
            ("B", "Selected", "Selected edge (will fallback to automatic selection if selected edge is boundary or seam)")
            ],
        ) # type: ignore
    
    init_edge_mode_dir : eP(
        name = "Init Edge Direction",
        description = "Direction of initial edge - useful if result is 'upside down'",
        items = [
            ("A", "Default", "From vert 1 to vert 0"),
            ("B", "Reversed", "From vert 0 to 1")
            ],
        ) # type: ignore
    
    angle_projection_mode : eP(
        name = "Angle Projection Mode",
        description = "How the projected angles of the face corners of each vertex are used",
        items = [
            ("A", "Sort Criteria", "Projected angles are used to sort face corners, then face corners' calculated angles (divided by combined face corner angles and multiplied by 360) are used for actual rotation"),
            ("B", "Projected Angles", "Projected angles are used as is, divided by combined face corner angles and multiplied by 360")
            ],
        ) # type: ignore
    
    only_islands_with_selection : bP(
        name = "Selected Islands Only",
        description = "Restrict to UV islands with a selection",
        default=True
        ) # type: ignore
    
    recalc_full_circle : bP(
        name = "Racalculate Angles",
        description = "Recalculate angles to equal 360 degrees (or debug value found below) combined",
        default=True
        ) # type: ignore
    
    smoothing : iP(
        name = "Smooothing Iterations",
        description = "For each smoothing iteration, another pass is allowed for each vertex, averaged with last pass",
        default=2
        ) # type: ignore
    
    debug_fcd : fP(
        name = "Debug Circular Degrees",
        description = "Should be 360? Use this to test if it shouldn't...",
        default=360.0
        ) # type: ignore

    #endregion

    def edge_bundle_sort(self, edge_bundle):
        comp_vec = mathutils.Vector((0,1,0))
        sort_val = abs(edge_bundle[1].dot(comp_vec)) + edge_bundle[2] + edge_bundle[3]
        return sort_val
    
    def loop_angle_sort(self, loop_angle_bundle):
        return loop_angle_bundle[1]

    def get_dir(self, co_source, co_destination):
        vec = co_destination - co_source
        dir = vec.normalized()
        return dir
    
    def execute(self, context):
        props = context.scene.sloppy_props
        bm = bmesh.from_edit_mesh(bpy.context.active_object.data)
        oco = bpy.context.active_object.location
        uv_layer = bm.loops.layers.uv.verify()
        islands = bmesh_utils.bmesh_linked_uv_islands(bm, uv_layer)
        if len(islands) < 1:
            print('No islands, using all faces.')
            islands = [[face for face in bm.faces]]
        cur_region, cur_view = props.get_view3d_region_and_space()

        bm.verts.ensure_lookup_table()
        bm.edges.ensure_lookup_table()
        bm.faces.ensure_lookup_table()

        curpos = None
        for area in bpy.context.screen.areas:
            if area.type == 'IMAGE_EDITOR':   #find the UVeditor
                curpos = area.spaces.active.cursor_location
        
        init_pos = mathutils.Vector((0.5, 0.5))
        if self.transform_pivot == "B" or self.transform_pivot == "D":
            init_pos = curpos

        should_do_islands = []

        for island in islands:
            island_has_selection = False
            should_do = True
            for iif in island:
                if iif.select == True:
                    island_has_selection = True
                for iife in iif.edges:
                    if iife.select == True:
                        if iife.seam == False and len(iife.link_faces) >= 2:
                            island_has_selection = True
            
            if self.only_islands_with_selection == True:
                should_do = island_has_selection
            
            should_do_islands.append(should_do)

        island_iter = 0
        for island, should_do_island in zip(islands, should_do_islands):
            if should_do_island == True:
                island_edges = []
                island_loops = []
                selected_edge = None
                for iif in island:
                    for iifl in iif.loops:
                        if iifl not in island_loops:
                            island_loops.append(iifl)
                    for iife in iif.edges:
                        if iife.select == True:
                            if iife.seam == False and len(iife.link_faces) >= 2:
                                if selected_edge == None:
                                    selected_edge = iife
                        if iife not in island_edges:
                            island_edges.append(iife)
                
                island_edges_bundled = []
                for ie in island_edges:
                    ie_is_seam = 0
                    ie_is_boundary = 0
                    if len(ie.link_faces) < 2:
                        ie_is_boundary += 1
                    if ie.seam == True:
                        ie_is_seam += 1
                    for iev in ie.verts:
                        for ieve in iev.link_edges:
                            if len(ieve.link_faces) < 2:
                                ie_is_boundary += 1
                            if ieve.seam == True:
                                ie_is_seam += 1
                    edge_bundle = [ie, (ie.verts[1].co - ie.verts[0].co).normalized(), ie_is_seam, ie_is_boundary]
                    island_edges_bundled.append(edge_bundle)
            
                init_edge = None

                if self.initial_edge_mode == "B":
                    init_edge = selected_edge
                if init_edge == None:
                    island_edges_bundled.sort(key=self.edge_bundle_sort)
                    init_edge = island_edges_bundled[0][0]

                print('Loops in island', island_iter, ':', len(island_loops))

                
                
                init_uva = mathutils.Vector((init_pos.x + (init_edge.calc_length() / 2), init_pos.y))
                init_uvb = mathutils.Vector((init_pos.x - (init_edge.calc_length() / 2), init_pos.y))
                if self.init_edge_mode_dir == "B":
                    init_uva = mathutils.Vector((init_pos.x - (init_edge.calc_length() / 2), init_pos.y))
                    init_uvb = mathutils.Vector((init_pos.x + (init_edge.calc_length() / 2), init_pos.y))
                
                loops_done = []
                loops_done_times = []

                for iel in init_edge.verts[0].link_loops:
                    if iel in island_loops:
                        iel[uv_layer].uv = init_uva
                        loops_done.append(iel)
                        loops_done_times.append(self.smoothing)
                for iel in init_edge.verts[1].link_loops:
                    if iel in island_loops:
                        iel[uv_layer].uv = init_uvb
                        loops_done.append(iel)
                        loops_done_times.append(self.smoothing)

                verts_done = [init_edge.verts[0], init_edge.verts[1]]
                
                init_dir_a = mathutils.Vector((-1,0))
                init_dir_b = mathutils.Vector((1,0))
                if self.init_edge_mode_dir == "B":
                    init_dir_a = mathutils.Vector((1,0))
                    init_dir_b = mathutils.Vector((-1,0))
                
                geo_dir_a = self.get_dir(init_edge.verts[1].co, init_edge.verts[0].co)
                geo_dir_b = self.get_dir(init_edge.verts[0].co, init_edge.verts[1].co)

                # sub-bundle structure:
                # [0] = current node vertex to do
                # [1] = edge traveled to vertex
                # [2] = node UVs
                # [3] = UV direction traveled
                # [4] = geometric direction traveled before seam; does not update while traveling along seams (used to determine loops to affect (using dot product) when traveling along seams)
                next_bundle = [
                    [
                        init_edge.verts[0], init_edge, init_uva, init_dir_a, geo_dir_a
                    ],
                    [
                        init_edge.verts[1], init_edge, init_uvb, init_dir_b, geo_dir_b
                    ]
                ]
                
                bundle_iter = 0
                while len(next_bundle) > 0:
                    bundle_iter += 1
                    this_bundle = next_bundle.copy()
                    next_bundle.clear()

                    sub_bundle_iter = 0
                    for bi, bundle in enumerate(this_bundle):
                        sub_bundle_iter += 1
                        print('Bundle', bundle_iter, '- Sub-bundle', sub_bundle_iter)
                        bundle_valid = True
                        for bve in bundle[0].link_edges:
                            if bve.seam == True:
                                bundle_valid = False
                            if len(bve.link_faces) < 2:
                                bundle_valid = False
                        
                        if bundle_valid:
                            print('\nSub-bundle', (bi + 1), 'of', len(this_bundle), '- vert index:', bundle[0].index, '- edge index:', bundle[1].index, '- UV position:', bundle[2], '- UV direction:', bundle[3], '- geo. direction:', bundle[4])
                            loops_to_do = []
                            angles_to_do = []
                            angle_sum = 0.0

                            p_co = bundle[0].co + oco
                            p_up = bundle[0].normal.orthogonal()
                            r_mat = mathutils.Matrix.Rotation(math.radians(90), 4, bundle[0].normal)
                            p_right = p_up.copy()
                            p_right.rotate(r_mat)
                            orig_co = bundle[1].other_vert(bundle[0]).co + oco
                            orig_vec = props.project_point_on_plane_axes(p_co, p_right, p_up, orig_co)
                            orig_dir = orig_vec.normalized()

                            la_bundles = []

                            for bli, bl in enumerate(bundle[0].link_loops):
                                angle_sum += math.degrees(bl.calc_angle())
                                bl_co = bl.edge.other_vert(bundle[0]).co + oco
                                bl_vec = props.project_point_on_plane_axes(p_co, p_right, p_up, bl_co)
                                bl_dir = bl_vec.normalized()
                                blangle = orig_dir.angle_signed(bl_dir)
                                if self.angle_projection_mode == "A":
                                    if bl.edge == bundle[1]:
                                        blangle = math.degrees(0.0)
                                    else:
                                        if math.degrees(blangle) < 0.0:
                                            nu_blangle = math.radians(360) + blangle
                                            blangle = nu_blangle

                                loops_to_do.append(bl)
                                angles_to_do.append(math.degrees(blangle))
                                if self.angle_projection_mode == "A":
                                    la_bundles.append([bl, blangle])
                            
                            if self.angle_projection_mode == "A":
                                angle_sum = 0.0
                                la_bundles.sort(key=self.loop_angle_sort)
                                for lai, la_bundle in enumerate(la_bundles):
                                    loops_to_do[lai] = la_bundle[0]
                                    angles_to_do[lai] = angle_sum
                                    angle_sum += math.degrees(la_bundle[0].calc_angle())
                                print('Angles to do:', angles_to_do)

                            angle_ratio = self.debug_fcd/angle_sum
                            print('Angle sum:', angle_sum)
                            print('Angle ratio = 1 :', angle_ratio)


                            loop_iter = 0
                            for this_loop, this_angle in zip(loops_to_do, angles_to_do):
                                loop_iter += 1
                                if this_loop.edge != bundle[1]:
                                    print('\nBundle', bundle_iter, '- Sub-bundle', (bi + 1), 'of', len(this_bundle), '- Loop', (loop_iter), 'of', len(loops_to_do), ':')
                                    this_rot_deg = this_angle
                                    if self.recalc_full_circle:
                                        this_rot_deg = (this_angle/angle_sum) * self.debug_fcd
                                    dir_a = self.get_dir(bundle[0].co, bundle[1].other_vert(bundle[0]).co)
                                    dir_b = self.get_dir(bundle[0].co, this_loop.edge.other_vert(bundle[0]).co)
                                    print('Last/Current geo direction:', dir_a, '/' ,dir_b)
                                    print('Actual/Calculated angle:', this_angle, '/' ,this_rot_deg)
                                    this_dir = mathutils.Vector((bundle[3].x,bundle[3].y))
                                    this_matrix = mathutils.Matrix.Rotation(math.radians(this_rot_deg), 2, 'Z')
                                    this_dir.rotate(this_matrix)
                                    print('Last/Current UV direction:', bundle[3], '/', this_dir)
                                    this_vec = this_dir * this_loop.edge.calc_length()
                                    this_uvco = bundle[2] + this_vec
                                    print('Source/Target UV position:', bundle[2], '/', this_uvco)
                                    
                                    other_vert = this_loop.edge.other_vert(bundle[0])

                                    if other_vert not in verts_done:
                                        is_seam_vert = False
                                        for ove in other_vert.link_edges:
                                            if ove.seam == True:
                                                is_seam_vert = True

                                        for ovl in other_vert.link_loops:
                                            done_before = False
                                            done_times_before = 0
                                            ldi = 0
                                            if ovl in island_loops:
                                                should_do_loop = True
                                                if ovl in loops_done:
                                                    done_before = True
                                                    ldi = loops_done.index(ovl)
                                                    done_times_before = loops_done_times[ldi]
                                                    if done_times_before >= self.smoothing:
                                                        should_do_loop = False
                                                else:
                                                    if is_seam_vert == True:
                                                        should_do_loop = False
                                                        for ovlel in ovl.edge.link_loops:
                                                            if ovlel in bundle[1].link_loops:
                                                                should_do_loop = True
                                                if should_do_loop == True:
                                                    if not done_before:
                                                        print('Loop', ovl.index ,'not done before. Old UVs:', ovl[uv_layer].uv)
                                                        ovl[uv_layer].uv = this_uvco
                                                        print('Loop', ovl.index ,'not done before. New UVs:', ovl[uv_layer].uv)
                                                        loops_done.append(ovl)
                                                        loops_done_times.append(done_times_before)
                                                    if done_before:
                                                        loops_done_times[ldi] = done_times_before + 1
                                                        avg_uvco = ovl[uv_layer].uv.copy()
                                                        print('Loop', ovl.index ,'done', done_times_before, 'of', self.smoothing, 'times before. Old UVs:', avg_uvco)
                                                        avg_uvco += this_uvco
                                                        avg_uvco /= 2
                                                        ovl[uv_layer].uv = avg_uvco
                                                        print('Loop', ovl.index ,'done', done_times_before, 'of', self.smoothing, 'times before. Old UVs:', ovl[uv_layer].uv)
                                                    ovl[uv_layer].pin_uv = True

                                        new_bundle = [other_vert, this_loop.edge, this_uvco, -this_dir, dir_b]
                                        next_bundle.append(new_bundle)
                                        next_bundle_i = bundle_iter + 1
                                        print('Added sub-bundle to bundle', next_bundle_i, 'which now contains', len(next_bundle), 'sub-bundles.')

                                        max_num_done_times = len(other_vert.link_loops) * self.smoothing
                                        these_loops_done = 0
                                        for ovl in other_vert.link_loops:
                                            if ovl in loops_done:
                                                ovldi = loops_done.index(ovl)
                                                these_loops_done += loops_done_times[ovldi]
                                        if these_loops_done >= max_num_done_times:
                                            verts_done.append(other_vert)
                                else:
                                    print('\nBundle', bundle_iter, '- Sub-bundle', (bi + 1), 'of', len(this_bundle), '- Loop', (loop_iter), 'of', len(loops_to_do), 'has origin edge! Skipping...')
                        else:
                            print('Bundle', bundle_iter, '- Sub-bundle', sub_bundle_iter, 'invalid!')
                
                bpy.context.active_object.data.update()
                
                avg_x = 0.0
                avg_y = 0.0
                min_x = 9999.9
                min_y = 9999.9
                max_x = -9999.9
                max_y = -9999.9
                navg_x = 0.0
                navg_y = 0.0
                nmin_x = 9999.9
                nmin_y = 9999.9
                nmax_x = -9999.9
                nmax_y = -9999.9

                for il in loops_done:
                    old_uv_x = il[uv_layer].uv.x
                    old_uv_y = il[uv_layer].uv.y
                    numin_x = min(min_x, old_uv_x)
                    min_x = numin_x
                    numin_y = min(min_y, old_uv_y)
                    min_y = numin_y
                    numax_x = max(max_x, old_uv_x)
                    max_x = numax_x
                    numax_y = max(max_y, old_uv_y)
                    max_y = numax_y
                    avg_x += old_uv_x
                    avg_y += old_uv_y
                
                avg_x /= len(loops_done)
                avg_y /= len(loops_done)
                midpoint = mathutils.Vector((avg_x, avg_y))

                x_diff = max_x - min_x
                y_diff = max_y - min_y
                bounds_min = mathutils.Vector((min_x, min_y))
                bounds_max = mathutils.Vector((max_x, max_y))
                print('UV Bounds Min:', bounds_min, '\nUV Bounds Max', bounds_max)
                if self.scaling_mode != "A":
                    scaling_vector = mathutils.Vector((1.0, 1.0))
                    if self.scaling_mode == "B":
                        scaling_vector.x = 1.0/max(x_diff, y_diff)
                        scaling_vector.y = 1.0/max(x_diff, y_diff)
                    if self.scaling_mode == "C":
                        scaling_vector.x = 1.0/x_diff
                        scaling_vector.y = 1.0/y_diff
                    for il in loops_done:
                        old_uv = il[uv_layer].uv
                        new_uv = ((old_uv-bounds_min)*scaling_vector)
                        il[uv_layer].uv = new_uv
                if self.scaling_mode == "A":
                    # midpoint = bounds_min.lerp(bounds_max, 0.5)
                    print('Midpoint:', midpoint)
                    vec_to_center = init_pos - midpoint
                    print('Vector to init position', init_pos, ':', vec_to_center)
                    for il in loops_done:
                        nnew_uv = il[uv_layer].uv
                        old_uv = il[uv_layer].uv
                        if self.transform_pivot == "A":
                            nnew_uv = old_uv
                        if self.transform_pivot == "B":
                            nnew_uv = old_uv
                        if self.transform_pivot == "C":
                            nnew_uv = old_uv+vec_to_center
                            il[uv_layer].uv = nnew_uv
                        if self.transform_pivot == "D":
                            nnew_uv = old_uv+vec_to_center
                            il[uv_layer].uv = nnew_uv
                        if self.transform_pivot == "E":
                            nnew_uv = old_uv-bounds_min
                            il[uv_layer].uv = nnew_uv
                        if nnew_uv != None:
                            new_uv_x = nnew_uv.x
                            new_uv_y = nnew_uv.y
                            nnumin_x = min(nmin_x, new_uv_x)
                            nmin_x = nnumin_x
                            nnumin_y = min(nmin_y, new_uv_y)
                            nmin_y = nnumin_y
                            nnumax_x = max(nmax_x, new_uv_x)
                            nmax_x = nnumax_x
                            nnumax_y = max(nmax_y, new_uv_y)
                            nmax_y = nnumax_y
                            navg_x += new_uv_x
                            navg_y += new_uv_y
                
                navg_x /= len(loops_done)
                navg_y /= len(loops_done)
                nmidpoint = mathutils.Vector((navg_x, navg_y))
                nx_diff = nmax_x - nmin_x
                ny_diff = nmax_y - nmin_y
                nbounds_min = mathutils.Vector((nmin_x, nmin_y))
                nbounds_max = mathutils.Vector((nmax_x, nmax_y))
                print('New UV Bounds Min:', nbounds_min, '\nNew UV Bounds Max', nbounds_max)
                # nmidpoint = nbounds_min.lerp(nbounds_max, 0.5)
                print('New midpoint:', nmidpoint)
                nvec_to_center = init_pos - nmidpoint
                print('New vector to init position', init_pos, ':', nvec_to_center)
                for face in island:
                    face.select_set(True)
        bpy.context.active_object.data.update()

        bpy.ops.uv.unwrap(method='ANGLE_BASED', fill_holes=True, correct_aspect=True, use_subsurf_data=False, margin=0, no_flip=False, iterations=10, use_weights=False, weight_group="uv_importance", weight_factor=1)

        for face in bm.faces:
            for loop in face.loops:
                loop[uv_layer].pin_uv = False


        return {"FINISHED"}

class SloppyBoundaryFirstUVUnfold(bpy.types.Operator):
    bl_idname = "operator.sloppy_boundary_first_uv_unfold"
    bl_label = "Boundary First UV Unfold"
    bl_description = ""
    bl_options = {'REGISTER', 'UNDO'}

    iP = bpy.props.IntProperty
    fP = bpy.props.FloatProperty
    fvP = bpy.props.FloatVectorProperty
    bP = bpy.props.BoolProperty
    eP = bpy.props.EnumProperty
    bvP = bpy.props.BoolVectorProperty
    sP = bpy.props.StringProperty

    # region prop def
    transform_pivot : eP(
        name = "Pivot Point",
        description = "Pivot point for rebuilding UV edge",
        items = [
            ("A", "Init Edge", "Init edge placed at U:0.5 V:0.5"),
            ("B", "Init Edge to Cursor", "Init edge placed at cursor"),
            ("C", "Center", "Island center placed at U:0.5 V:0.5"),
            ("D", "Center to Cursor", "Island center placed at cursor"),
            ("E", "Positive UV Space", "Shift pivot to keep everything in positive UV space")
            ],
        default="E"
        ) # type: ignore
    
    scaling_mode : eP(
        name = "Scaling Mode",
        description = "Scaling of result",
        items = [
            ("A", "World Scale", "Scaled 1:1 World to UV"),
            ("B", "Scale to Bounds", "Scaled uniformly to fit within UV 0-1 using largest dimension"),
            ("C", "Stretch to Bounds", "Scaled non- uniformly to fit within UV 0-1")
            ],
        ) # type: ignore
    
    edge_length_mode : eP(
        name = "Edge Length Mode",
        description = "How the edge lengths are translated to UV space",
        items = [
            ("A", "As Is", "World space edge lengths are used"),
            ("B", "Projected", "Edge lengths are calculated in the same plane as the rotation angle")
            ],
        ) # type: ignore
    
    proj_norm_mode : eP(
        name = "Projection Normal Mode",
        description = "How the normal of the angle projection plane is calculated",
        items = [
            ("A", "Source Vertex", "Normal of the vertex traveled from"),
            ("B", "Destination Vertex", "Normal of the vertex traveled to"),
            ("C", "Edge Average", "Average of the normals of traveled edge's vertices"),
            ("D", "Average", "Average of the surrounding faces")
            ],
        ) # type: ignore
    
    quant_norm : eP(
        name = "Quantize Normals",
        description = "How the projected angles of the face corners of each vertex are used",
        items = [
            ("A", "No", ""),
            ("B", "6", "Normal is quantized to the faces of a cube"),
            ("C", "26", "Normal is quantized to all the elements (vertices, edges, faces) of a cube")
            ],
        ) # type: ignore
    
    only_islands_with_selection : bP(
        name = "Selected Islands Only",
        description = "Restrict to UV islands with a selection",
        default=True
        ) # type: ignore
    
    debug_fcd : fP(
        name = "Debug Circular Degrees",
        description = "Should be 360? Use this to test if it shouldn't...",
        default=360.0
        ) # type: ignore
    
    verbose : bP(
        name = "Verbose",
        description = "",
        default=True
        ) # type: ignore

    #endregion

    def edge_bundle_sort(self, edge_bundle):
        comp_vec = mathutils.Vector((0,1,0))
        sort_val = abs(edge_bundle[1].dot(comp_vec)) + edge_bundle[2] + edge_bundle[3]
        return sort_val
    
    def loop_angle_sort(self, loop_angle_bundle):
        return loop_angle_bundle[1]

    def get_dir(self, co_source, co_destination):
        vec = co_destination - co_source
        dir = vec.normalized()
        return dir
    
    def quantize_sort(self, ql):
        VdotN = ql[1].dot(ql[0])
        # print('Vector:', ql[0].x, ql[0].y, ql[0].z, '- Dot:', VdotN)
        return VdotN
    
    def make_quantize_six(self, in_normal):
        quantize_six = [
            [
                mathutils.Vector((1,0,0)),
                in_normal
            ],
            [
                mathutils.Vector((-1,0,0)),
                in_normal
            ],
            [
                mathutils.Vector((0,1,0)),
                in_normal
            ],
            [
                mathutils.Vector((0,-1,0)),
                in_normal
            ],
            [
                mathutils.Vector((0,0,1)),
                in_normal
            ],
            [
                mathutils.Vector((0,0,-1)),
                in_normal
            ],
        ]
        return quantize_six
    
    def make_quantize_twentysix(self, in_normal):
        quantize_six = [
            mathutils.Vector((1,0,0)),
            mathutils.Vector((-1,0,0)),
            mathutils.Vector((0,1,0)),
            mathutils.Vector((0,-1,0)),
            mathutils.Vector((0,0,1)),
            mathutils.Vector((0,0,-1))
        ]
        qtsx = []
        for qsxa in quantize_six:
            for qsxb in quantize_six:
                for qsxc in quantize_six:
                    qx = bl_math.clamp(qsxa.x + qsxb.x + qsxc.x, -1, 1)
                    qy = bl_math.clamp(qsxa.y + qsxb.y + qsxc.y, -1, 1)
                    qz = bl_math.clamp(qsxa.z + qsxb.z + qsxc.z, -1, 1)
                    
                    if (abs(qx)+abs(qy)+abs(qz)) > 0:
                        cmbqsx = mathutils.Vector((qx,qy,qz))
                        nqsx = cmbqsx / (abs(qx)+abs(qy)+abs(qz))
                        if nqsx not in qtsx:
                            qtsx.append(nqsx)
        quantize_twentysix = []
        for qtsxi in qtsx:
            quantize_twentysix.append([qtsxi, in_normal])
        return quantize_twentysix
    
    def execute(self, context):
        props = context.scene.sloppy_props
        bm = bmesh.from_edit_mesh(bpy.context.active_object.data)
        oco = bpy.context.active_object.location
        uv_layer = bm.loops.layers.uv.verify()
        islands = bmesh_utils.bmesh_linked_uv_islands(bm, uv_layer)
        if len(islands) < 1:
            print('No islands, using all faces.')
            islands = [[face for face in bm.faces]]

        bm.verts.ensure_lookup_table()
        bm.edges.ensure_lookup_table()
        bm.faces.ensure_lookup_table()

        curpos = None
        for area in bpy.context.screen.areas:
            if area.type == 'IMAGE_EDITOR':   #find the UVeditor
                curpos = area.spaces.active.cursor_location
        
        init_pos = mathutils.Vector((0.5, 0.5))
        if self.transform_pivot == "B" or self.transform_pivot == "D":
            init_pos = curpos

        should_do_islands = []

        for island in islands:
            island_has_selection = False
            should_do = True
            for iif in island:
                if iif.select == True:
                    island_has_selection = True
                for iife in iif.edges:
                    if iife.select == True:
                        if iife.seam == False and len(iife.link_faces) >= 2:
                            island_has_selection = True
            
            if self.only_islands_with_selection == True:
                should_do = island_has_selection
            
            should_do_islands.append(should_do)

        island_iter = 0
        for island, should_do_island in zip(islands, should_do_islands):
            should = should_do_island
            this_verbose = self.verbose and should_do_island
            loop_chains, loop_chains_eton, vert_chains, edge_chains, face_chains = props.find_island_boundaries_and_their_loops(bm, island, True, this_verbose)
            
            # for lc,lce in zip(loop_chains, loop_chains_eton):
            #     for l, e in zip(lc, lce):
            #         if e:
            #             print(l.index, e.index)
            #         else:
            #             print(l.index, e)

            if len(loop_chains) > 1:
                should = False
            if should == True:
                island_edges = []
                island_loops = []
                
                for iif in island:
                    for iifl in iif.loops:
                        if iifl not in island_loops:
                            island_loops.append(iifl)
                    for iife in iif.edges:
                        if iife not in island_edges:
                            island_edges.append(iife)
                
                lcf = loop_chains[0].copy()
                lcfes = loop_chains_eton[0].copy()
                lcf_loops = []
                lcf_pos = []
                lcr = loop_chains[0].copy()
                lcr.reverse()
                lcres = loop_chains_eton[0].copy()
                lcres.reverse()
                lcres_pop = lcres.pop(0)
                lcres.append(lcres_pop)
                lcr_loops = []
                lcr_pos = []

                full_rotation_f = 0.0
                full_rotation_r = 0.0
                rotations_f = 0
                rotations_r = 0

                print('Loop chain length:', len(lcf))
                print('Edge to next chain length:', len(lcfes))
                for nlcfi, nlcfl, nlcfe, nlcrl, nlcre in zip(range(len(lcf)), lcf, lcfes, lcr, lcres):
                    nlcfe_str = None
                    if nlcfe:
                        nlcfe_str = str(nlcfe.index)
                    nlcre_str = None
                    if nlcre:
                        nlcre_str = str(nlcre.index)
                    print(nlcfi, ':', nlcfl.index, '-', nlcfe_str, '<--->', nlcrl.index, '-', nlcre_str)

                #region Measure rotation
                lcfe_iter_mr = 0
                last_edge_f_mr = None
                while not last_edge_f_mr:
                    lcfe_iter_mr -= 1
                    last_edge_f_mr = lcfes[lcfe_iter_mr]
                last_uvdirf = mathutils.Vector((1,0))
                last_vec_f = None

                lcre_iter_mr = 0
                last_edge_r_mr = None
                while not last_edge_r_mr:
                    lcre_iter_mr -= 1
                    last_edge_r_mr = lcres[lcre_iter_mr]
                last_uvdirr = mathutils.Vector((1,0))
                last_vec_r = None

                iposf = mathutils.Vector((0,0))
                iposr = mathutils.Vector((0,0))

                last_f_was_edge = False

                
                f_iter = 0
                for lcfl, lcfe in zip(lcf, lcfes):
                    # lcf_loops.append(lcfl.index)
                    # lcf_pos.append(iposf)
                    if lcfe:
                        edge_length_f = lcfe.calc_length()
                        ivf = lcfl.vert
                        ovf = lcfe.other_vert(ivf)
                        ivfna = ivf.normal
                        if self.proj_norm_mode != "A":
                            if self.proj_norm_mode == "B":
                                ivfna = ovf.normal
                            if self.proj_norm_mode == "C":
                                ivfna += ovf.normal
                            if self.proj_norm_mode == "D":
                                for ivff in ivf.link_faces:
                                    if ivff in island:
                                        ivfna += ivff.normal
                        ivfn = ivfna.normalized()
                        if self.quant_norm != "A":
                            qlst_f = None
                            if self.quant_norm == "B":
                                qlst_f = self.make_quantize_six(ivfn)
                            if self.quant_norm == "C":
                                qlst_f = self.make_quantize_twentysix(ivfn)
                            if qlst_f:
                                qlst_f.sort(key=self.quantize_sort, reverse=True)
                                ivfn = qlst_f[0][0]
                                print('Quantized normal:', qlst_f[0][0].x, qlst_f[0][0].y, qlst_f[0][0].z)
                        plane_x_f, plane_y_f = props.plane_axes_from_normal(ivfn)
                        if not last_vec_f:
                            lovf = last_edge_f_mr.other_vert(ivf)
                            lvec_f = -props.project_point_on_plane_axes(ivf.co, plane_x_f, plane_y_f, lovf.co)
                            last_vec_f = lvec_f.normalized()
                        new_vec_f = props.project_point_on_plane_axes(ivf.co, plane_x_f, plane_y_f, ovf.co)
                        new_vec_nf = new_vec_f.normalized()
                        new_angle_f = last_vec_f.angle_signed(new_vec_nf)
                        this_matrix_f = mathutils.Matrix.Rotation(-new_angle_f, 2, 'Z')
                        full_rotation_f += new_angle_f
                        rotations_f += 1
                        last_uvdirf.rotate(this_matrix_f)
                        new_iposf = iposf + (last_uvdirf * edge_length_f)
                        iposf = new_iposf
                        last_edge_f_mr = lcfe
                        last_vec_f = None
                        last_f_was_edge = True
                    else:
                        last_f_was_edge = False
                    f_iter += 1

                # iposr = iposf.copy()
                last_r_was_edge = False

                r_iter = 0
                for lcrl, lcre in zip(lcr, lcres):
                    # lcr_loops.append(lcrl.index)
                    # lcr_pos.append(iposr)
                    if lcre:
                        edge_length_r = lcre.calc_length()
                        ivr = lcrl.vert
                        ovr = lcre.other_vert(ivr)
                        ivrna = ivr.normal
                        if self.proj_norm_mode != "A":
                            if self.proj_norm_mode == "B":
                                ivrna = ovr.normal
                            if self.proj_norm_mode == "C":
                                ivrna += ovr.normal
                            if self.proj_norm_mode == "D":
                                for ivrf in ivr.link_faces:
                                    if ivrf in island:
                                        ivrna += ivrf.normal
                        ivrn = ivrna.normalized()
                        if self.quant_norm != "A":
                            qlst_r = None
                            if self.quant_norm == "B":
                                qlst_r = self.make_quantize_six(ivrn)
                            if self.quant_norm == "C":
                                qlst_r = self.make_quantize_twentysix(ivrn)
                            if qlst_r:
                                qlst_r.sort(key=self.quantize_sort, reverse=True)
                                ivrn = qlst_r[0][0]
                                print('Quantized normal:', qlst_r[0][0].x, qlst_r[0][0].y, qlst_r[0][0].z)
                        plane_x_r, plane_y_r = props.plane_axes_from_normal(ivrn)
                        if not last_vec_r:
                            lovr = last_edge_r_mr.other_vert(ivr)
                            lvec_r = -props.project_point_on_plane_axes(ivr.co, plane_x_r, plane_y_r, lovr.co)
                            last_vec_r = lvec_r.normalized()
                        new_vec_r = props.project_point_on_plane_axes(ivr.co, plane_x_r, plane_y_r, ovr.co)
                        new_vec_nr = new_vec_r.normalized()
                        new_angle_r = last_vec_r.angle_signed(new_vec_nr)
                        this_matrix_r = mathutils.Matrix.Rotation(-new_angle_r, 2, 'Z')
                        full_rotation_r -= new_angle_r
                        rotations_r += 1
                        last_uvdirr.rotate(this_matrix_r)
                        new_iposr = iposr + (last_uvdirr * edge_length_r)
                        iposr = new_iposr
                        last_edge_r_mr = lcre
                        last_vec_r = None
                        last_r_was_edge = True
                    else:
                        last_r_was_edge = False
                    r_iter += 1
                
                
                rot_fac_f = self.debug_fcd/math.degrees(full_rotation_f)
                rot_fac_r = self.debug_fcd/math.degrees(full_rotation_r)

                print('Full rotation forward:', math.degrees(full_rotation_f), '- Rotation factor forward:', rot_fac_f)
                print('Full rotation reverse:', math.degrees(full_rotation_r), '- Rotation factor reverse:', rot_fac_r)
                
                #endregion

                #region Rotation
                lcfe_iter = 0
                last_edge_f = None
                while not last_edge_f:
                    lcfe_iter -= 1
                    last_edge_f = lcfes[lcfe_iter]
                print(last_edge_f)
                last_uvdirf = mathutils.Vector((1,0))
                last_vec_f = None

                lcre_iter = 0
                last_edge_r = None
                while not last_edge_r:
                    lcre_iter -= 1
                    last_edge_r = lcres[lcre_iter]
                print(last_edge_r)
                last_uvdirr = mathutils.Vector((-1,0))
                last_vec_r = None

                # ("A", "Init Edge", "Init edge placed at U:0.5 V:0.5"),
                # ("B", "Init Edge to Cursor", "Init edge placed at cursor"),
                # ("C", "Center", "Island center placed at U:0.5 V:0.5"),
                # ("D", "Center to Cursor", "Island center placed at cursor"),
                # ("E", "Positive UV Space", "Shift pivot to keep everything in positive UV space")
                
                iposf = mathutils.Vector((0,0))
                iposr = mathutils.Vector((0,0))
                if self.transform_pivot == "A":
                    iposf = mathutils.Vector((.5,.5))
                    iposr = mathutils.Vector((.5,.5))
                if self.transform_pivot == "B":
                    iposf = curpos
                    iposr = curpos

                last_f_was_edge = False

                f_iter = 0
                for lcfl, lcfe in zip(lcf, lcfes):
                    lcf_loops.append(lcfl.index)
                    lcf_pos.append(iposf)
                    if lcfe:
                        edge_length_f = lcfe.calc_length()
                        ivf = lcfl.vert
                        ovf = lcfe.other_vert(ivf)
                        ivfna = ivf.normal
                        if self.proj_norm_mode != "A":
                            if self.proj_norm_mode == "B":
                                ivfna = ovf.normal
                            if self.proj_norm_mode == "C":
                                ivfna += ovf.normal
                            if self.proj_norm_mode == "D":
                                for ivff in ivf.link_faces:
                                    if ivff in island:
                                        ivfna += ivff.normal
                        ivfn = ivfna.normalized()
                        if self.quant_norm != "A":
                            qlst_f = None
                            if self.quant_norm == "B":
                                qlst_f = self.make_quantize_six(ivfn)
                            if self.quant_norm == "C":
                                qlst_f = self.make_quantize_twentysix(ivfn)
                            if qlst_f:
                                qlst_f.sort(key=self.quantize_sort, reverse=True)
                                ivfn = qlst_f[0][0]
                                print('Quantized normal:', qlst_f[0][0].x, qlst_f[0][0].y, qlst_f[0][0].z)
                        plane_x_f, plane_y_f = props.plane_axes_from_normal(ivfn)
                        if not last_vec_f:
                            lovf = last_edge_f.other_vert(ivf)
                            lvec_f = -props.project_point_on_plane_axes(ivf.co, plane_x_f, plane_y_f, lovf.co)
                            if self.edge_length_mode == "B":
                                edge_length_f = lvec_f.length
                            last_vec_f = lvec_f.normalized()
                        new_vec_f = props.project_point_on_plane_axes(ivf.co, plane_x_f, plane_y_f, ovf.co)
                        new_vec_nf = new_vec_f.normalized()
                        new_angle_f = last_vec_f.angle_signed(new_vec_nf)
                        adj_angle_f = new_angle_f * rot_fac_f
                        this_matrix_f = mathutils.Matrix.Rotation(-adj_angle_f, 2, 'Z')
                        last_uvdirf.rotate(this_matrix_f)
                        new_iposf = iposf + (last_uvdirf * edge_length_f)
                        iposf = new_iposf
                        last_edge_f = lcfe
                        last_vec_f = None
                        print('( Fwd. step', f_iter, ') Vertex:', lcfl.vert.index, '- loop:', lcfl.index, '- edge to next vertex:', lcfe.index, '- next vertex:', ovf.index, '\n   - edge length:', edge_length_f, '- angle:', math.degrees(-new_angle_f), '- adjusted angle:', math.degrees(-adj_angle_f), '\n   - new UV direction:', last_uvdirf, '- next UV position:', iposf)
                        last_f_was_edge = True
                    else:
                        print('( Fwd. step', f_iter, ') Vertex', lcfl.vert.index, '- loop', lcfl.index, '\n   - current UV position:', iposf)
                        last_f_was_edge = False
                    f_iter += 1

                # iposr = iposf.copy()
                last_r_was_edge = False

                r_iter = 0
                for lcrl, lcre in zip(lcr, lcres):
                    lcr_loops.append(lcrl.index)
                    lcr_pos.append(iposr)
                    if lcre:
                        edge_length_r = lcre.calc_length()
                        ivr = lcrl.vert
                        ovr = lcre.other_vert(ivr)
                        ivrna = ivr.normal
                        if self.proj_norm_mode != "A":
                            if self.proj_norm_mode == "B":
                                ivrna = ovr.normal
                            if self.proj_norm_mode == "C":
                                ivrna += ovr.normal
                            if self.proj_norm_mode == "D":
                                for ivrf in ivr.link_faces:
                                    if ivrf in island:
                                        ivrna += ivrf.normal
                        ivrn = ivrna.normalized()
                        if self.quant_norm != "A":
                            qlst_r = None
                            if self.quant_norm == "B":
                                qlst_r = self.make_quantize_six(ivrn)
                            if self.quant_norm == "C":
                                qlst_r = self.make_quantize_twentysix(ivrn)
                            if qlst_r:
                                qlst_r.sort(key=self.quantize_sort, reverse=True)
                                ivrn = qlst_r[0][0]
                                print('Quantized normal:', qlst_r[0][0].x, qlst_r[0][0].y, qlst_r[0][0].z)
                        plane_x_r, plane_y_r = props.plane_axes_from_normal(ivrn)
                        if not last_vec_r:
                            lovr = last_edge_r.other_vert(ivr)
                            lvec_r = -props.project_point_on_plane_axes(ivr.co, plane_x_r, plane_y_r, lovr.co)
                            if self.edge_length_mode == "B":
                                edge_length_r = lvec_r.length
                            last_vec_r = lvec_r.normalized()
                        new_vec_r = props.project_point_on_plane_axes(ivr.co, plane_x_r, plane_y_r, ovr.co)
                        new_vec_nr = new_vec_r.normalized()
                        new_angle_r = last_vec_r.angle_signed(new_vec_nr)
                        adj_angle_r = new_angle_r * rot_fac_r
                        this_matrix_r = mathutils.Matrix.Rotation(-adj_angle_r, 2, 'Z')
                        last_uvdirr.rotate(this_matrix_r)
                        new_iposr = iposr + (last_uvdirr * edge_length_r)
                        iposr = new_iposr
                        last_edge_r = lcre
                        last_vec_r = None
                        print('( Rev. step', r_iter, ') Vertex', lcrl.vert.index, '- loop', lcrl.index, '- edge to next vertex:', lcre.index, '- next vertex:', ovf.index, '\n   - edge length:', edge_length_r, '- angle:', math.degrees(-new_angle_r), '- adjusted angle:', math.degrees(-adj_angle_r), '\n   - new UV direction:', last_uvdirr, '- current UV position:', iposr)
                        last_r_was_edge = True
                    else:
                        print('( Rev. step', r_iter, ') Vertex', lcrl.vert.index, '- loop', lcrl.index, '\n   - current UV position:', iposr)
                        last_r_was_edge = False
                    r_iter += 1
                #endregion


                # if last_f_was_edge == False or last_r_was_edge == False:
                #     lcf_loops.append(lcf[0].index)
                #     lcf_pos.append(iposf)
                #     lcr_loops.append(lcr[0].index)
                #     lcr_pos.append(iposr)

                print('Forward loops:', lcf_loops)
                print('Reverse loops:', lcr_loops)

                first_lcfi = lcf[0].index
                last_lcfi = lcf[-1].index
                for lf in lcf:
                    fpi = lcf_loops.index(lf.index)
                    fp = lcf_pos[fpi]
                    rpi = lcr_loops.index(lf.index)
                    rp = lcr_pos[rpi]
                    l_mid = fp.lerp(rp, 0.5)
                    l_mid = (fp + rp) / 2.0
                    if last_f_was_edge == False:
                        if lf.index == last_lcfi or lf.index == first_lcfi:
                            # f_mid = lcf_pos[0].lerp(lcf_pos[-1], 0.5)
                            f_mid = (lcf_pos[0] + lcf_pos[-1]) / 2.0
                            # r_mid = lcr_pos[0].lerp(lcr_pos[-1], 0.5)
                            r_mid = (lcr_pos[0] + lcr_pos[-1]) / 2.0
                            # l_mid = f_mid.lerp(r_mid, 0.5)
                            l_mid = (f_mid + r_mid) / 2.0
                            # fp = (lcf_pos[0] + lcf_pos[-1]) / 2.0
                            # fp = lcf_pos[0].lerp(lcf_pos[-1], 0.5)
                    print('Vertex', lf.vert.index, '> loop', lf.index, '--> new position:', l_mid)
                    lf[uv_layer].uv = l_mid
                    lf[uv_layer].pin_uv = True

                for face in bm.faces:
                    face.select_set(False)

                for face in island:
                    face.select_set(True)
        bpy.context.active_object.data.update()

        bpy.ops.uv.unwrap(method='ANGLE_BASED', fill_holes=True, correct_aspect=True, use_subsurf_data=False, margin=0, no_flip=False, iterations=10, use_weights=False, weight_group="uv_importance", weight_factor=1)

        for face in bm.faces:
            for loop in face.loops:
                loop[uv_layer].pin_uv = False


        return {"FINISHED"}