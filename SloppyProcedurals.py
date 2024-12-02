import bpy
import bpy_extras
from bpy_extras import bmesh_utils
import bmesh
import math
import mathutils
import bl_math
import sys

class ProceduralQuadUVUnfold(bpy.types.Operator):
    bl_idname = "operator.uv_quad_unfold"
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

    # region Property list
    unfold_mode : eP(
        name = "Mode",
        description = "Unfolding Mode",
        items = [
            ("A", "All Directions", "Work outwards from initial quad"),
            ("B", "X First", "Work outwards in (local) X first until no more faces are found, calculating edge lengths to use in next stage, then work outwards in Y. repeating first step when more faces are found in X"),
            ("C", "Y First", "Work outwards in (local) Y first until no more faces are found, calculating edge lengths to use in next stage, then work outwards in X. repeating first step when more faces are found in Y")
            ]
        ) # type: ignore

    initial_quad : eP(
        name = "Initial Quad(s)",
        description = "Which quad (for each UV island, if there are more than one,) to start from",
        items = [
            ("A", "Most Regular", "Start with quad (or virtual quad, formed by combining most regular face with adjacent triangle as calculated from face corner angles) that is flattest and/or that has normal most similar to island average normals"),
            ("B", "Selected", "First selected quad (or virtual quad, formed from two first triangles selected, or, if only 1 triangle is selected in island, formed by combining with adjacent triangle as calculated from face corner angles) in each island (if an island has no faces selected, initial quad falls back to Most Regular)")
            ]
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
        name = "Only Move Current Face's Loops",
        description = "Avoid moving loops of other faces when moving UV loops",
        default = False
        ) # type: ignore
    
    #endregion

    def execute(self, context):
        props = context.scene.sloppy_props
        bm = bmesh.from_edit_mesh(bpy.context.active_object.data)
        uv_layer = bm.loops.layers.uv.verify()
        islands = bmesh_utils.bmesh_linked_uv_islands(bm, uv_layer)

        # region Initial variables
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

        attr_dict = [
            {"name": "naboL", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "naboR", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "naboU", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "naboD", "type": "INT", "domain": "FACE", "layer": None},
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
        
        naboL = props.get_dict_layer("naboL", attr_dict)
        naboR = props.get_dict_layer("naboR", attr_dict)
        naboU = props.get_dict_layer("naboU", attr_dict)
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
        xi = props.get_dict_layer("xi", attr_dict)
        yi = props.get_dict_layer("yi", attr_dict)
        avg_length_x = props.get_dict_layer("avg_length_x", attr_dict)
        avg_length_y = props.get_dict_layer("avg_length_y", attr_dict)


        co_attr_arr = [[co_ur_x, co_ur_y], [co_ul_x, co_ul_y], [co_ll_x, co_ll_y], [co_lr_x, co_lr_y]]

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

        # region Operator functions
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
                                    main_tri_adjacent = True
                                    diagonal_edge = avfe
                        if main_tri_adjacent == True:
                            far_tri = avf
                            if len(far_tri.edges) > 3:
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
            ulc = str(qc[1].index).rjust(5)
            uline = str(qe[0].index).center(11,"-")
            urc = str(qc[0].index).ljust(5)
            vline = "|".center(7)
            vspace = " ".center(9)
            vu = "^".center(7)
            vd = "v".center(7)
            larrow = "<"
            rarrow = ">"
            lline = str(qe[1].index).rjust(5)
            fnam = str(qf.index).center(11)
            if virt == True:
                fnama = str(vqf[0].index) + "/" + str(vqf[1].index)
                fnam = fnama.center(12)
            rline = str(qe[3].index).ljust(5)
            llc = str(qc[2].index).rjust(5)
            dline = str(qe[2].index).center(11, "-")
            lrc = str(qc[3].index).ljust(5)
            quad_msg = "\n{0}{1}{2}\n{3}{4}{3}\n{3}{4}{3}\n{3}{4}{3}\n{5}{6}{7}\n{3}{4}{3}\n{3}{4}{3}\n{3}{4}{3}\n{8}{9}{10}\n".format(ulc,uline,urc,vline,vspace,lline,fnam,rline, llc, dline, lrc)
            if wnd == "CCW":
                quad_msg = "\n{0}{1}{13}{2}\n{12}{4}{3}\n{3}{4}{3}\n{3}{4}{3}\n{5}{6}{7}\n{3}{4}{3}\n{3}{4}{3}\n{3}{4}{11}\n{8}{14}{9}{10}\n".format(ulc,uline,urc,vline,vspace,lline,fnam,rline, llc, dline, lrc, vu, vd, larrow, rarrow)
            if wnd == "CW":
                quad_msg = "\n{0}{14}{1}{2}\n{3}{4}{11}\n{3}{4}{3}\n{3}{4}{3}\n{5}{6}{7}\n{3}{4}{3}\n{3}{4}{3}\n{12}{4}{3}\n{8}{9}{13}{10}\n".format(ulc,uline,urc,vline,vspace,lline,fnam,rline, llc, dline, lrc, vu, vd, larrow, rarrow)
            return quad_msg

        def viz_add_quad(qe, oqe, qc, oqc, qf, oqf, vqf, ovqf, virt, ovirt, add_dir):
            ulc = str(qc[1].index).rjust(5)
            oulc = str(oqc[1].index).rjust(5)
            uline = str(qe[0].index).center(12,"=")
            ouline = str(oqe[0].index).center(12,"-")
            urc = str(qc[0].index).ljust(5)
            ourc = str(oqc[0].index).ljust(5)
            vline = "||".center(7)
            ovline = "|".center(7)
            vspace = " ".center(9)
            ovspace = " ".center(19)
            ohspace = " ".center(7)
            ovu = "^".center(19)
            ovd = "v".center(19)
            olarrow = "<".center(7)
            orarrow = ">".center(7)
            lline = str(qe[1].index).rjust(5)
            olline = str(oqe[1].index).rjust(5)
            fnam = str(qf.index).center(12)
            if virt == True:
                fnama = str(vqf[0].index) + "/" + str(vqf[1].index)
                fnam = fnama.center(12)
            ofnam = str(oqf.index).center(12)
            if ovirt == True:
                ofnama = str(ovqf[0].index) + "/" + str(ovqf[1].index)
                ofnam = ofnama.center(12)
            rline = str(qe[3].index).ljust(5)
            orline = str(oqe[3].index).ljust(5)
            llc = str(qc[2].index).rjust(5)
            ollc = str(oqc[2].index).rjust(5)
            dline = str(qe[2].index).center(12, "=")
            odline = str(oqe[2].index).center(12, "-")
            lrc = str(qc[3].index).ljust(5)
            olrc = str(oqc[3].index).ljust(5)

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

        previous_island_length_x = 0.0
        previous_island_length_y = 0.0
        previous_island_end_x = 0.0
        previous_island_end_y = 0.0

        for ii, island in enumerate(islands):
            bm.verts.ensure_lookup_table()
            bm.edges.ensure_lookup_table()
            bm.faces.ensure_lookup_table()
            edge_length_map_x = [[],[]]
            edge_length_map_y = [[],[]]
            process_sequence = 0
            faces_remain = []
            faces_done = []
            verts_done = []
            faces_selected = []
            island_loops = []
            init_face = None
            init_virtual_face = []
            avg_norm = mathutils.Vector((0,0,0))
            i_offset = mathutils.Vector((self.offset_per_island[0] * ii, self.offset_per_island[1] * ii))
            print(f"Island {ii}:")

            # region Find initial face(s)
            for face in island:
                avg_norm += face.normal
                avg_angle = 0.0
                angle_count = 0
                for fe in face.edges:
                    angle_count += 1
                    if fe.is_boundary == True or len(face.edges) < 2:
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
                if props.verbose == True:
                    print(f"Average island normals quantized to {quant_arr[0][2]}")

            bpy.context.active_object.data.update()

            for face in island:
                new_flatness = props.remap_val(face[flatness], min_avg_angle, max_avg_angle, 0.0, 1.0)
                face[flatness] = new_flatness
                n_dot_avg = face.normal.dot(current_avg_normal)
                face[normal_regularity] = props.remap_val(n_dot_avg, -1.0, 1.0, 0.0, 1.0)
                face[island_index] = ii
            
            bpy.context.active_object.data.update()

            initial_face_found = False
            
            while initial_face_found == False and len(faces_remain) > 0:
                is_selected_face = False
                init_face_is_virtual = False
                init_virtual_face_diagonal = None
                if self.initial_quad == "B":
                    if len(faces_selected) > 0:
                        init_face = faces_selected[0]
                        if len(init_face.edges) < 4:
                            if len(faces_selected) > 1:
                                if len(faces_selected[1]) == 3:
                                    init_face_is_virtual = True
                                    init_virtual_face.append(init_face)
                                    init_virtual_face.append(faces_selected[1])
                                    is_selected_face = True
                                elif len(faces_selected[1]) == 4:
                                    init_face = faces_selected[1]
                                    is_selected_face = True
                                else:
                                    print("Too inconsistent mesh! Aborting!")
                                    return  {"CANCELED"}
                            else:
                                other_tri = get_adjacent_triangle(init_face)
                                if other_tri[0] != None:
                                    init_face_is_virtual = True
                                    init_virtual_face.append(init_face)
                                    init_virtual_face.append(other_tri[0])
                                    is_selected_face = True
                                    init_virtual_face_diagonal = other_tri[1]
                                else:
                                    print("Too inconsistent mesh! Aborting!")
                                    return  {"CANCELED"}
                    else:
                        faces_remain.sort(key=regularity_sort)
                        init_face = faces_remain[0]
                        if len(init_face.edges) == 3:
                            other_tri = get_adjacent_triangle(init_face)
                            if other_tri[0] != None:
                                init_face_is_virtual = True
                                init_virtual_face.append(init_face)
                                init_virtual_face.append(other_tri[0])
                                init_virtual_face_diagonal = other_tri[1]
                            else:
                                other_quad = get_any_adjacent_quad(init_face)
                                if other_quad != None:
                                    init_face = other_quad
                                else:
                                    print("Too inconsistent mesh! Aborting!")
                                    return  {"CANCELED"}
                else:
                    faces_remain.sort(key=regularity_sort)
                    init_face = faces_remain[0]
                    if len(init_face.edges) == 3:
                        other_tri = get_adjacent_triangle(init_face)
                        if other_tri[0] != None:
                            init_face_is_virtual = True
                            init_virtual_face.append(init_face)
                            init_virtual_face.append(other_tri[0])
                            init_virtual_face_diagonal = other_tri[1]
                        else:
                            other_quad = get_any_adjacent_quad(init_face)
                            if other_quad != None:
                                init_face = other_quad
                            else:
                                print("Too inconsistent mesh! Aborting!")
                                return  {"CANCELED"}

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
                            print(avife.index)
                    for bvife in init_virtual_face[1].edges:
                        if bvife != init_virtual_face_diagonal:
                            vif_edges.append(bvife)
                            print(bvife.index)
                    current_face_center = calc_edge_center(init_virtual_face_diagonal)
                    init_face_dir_arr = get_init_dir_edges(vif_edges)
                    while check_list_for_duplicate(init_face_dir_arr) == True and current_quant_index < 4:
                        current_quant_index += 1
                        print(current_quant_index)
                        try:
                            current_dir_array = quant_arr[current_quant_index][1]
                        except:
                            print(f"Failed to get directions for virtual face {init_virtual_face[0].index}/{init_virtual_face[1].index}! Marking faces as trailing triangles...")
                    
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
                    
                    if current_quant_index <= 3:
                        for vf in init_virtual_face:
                            vif[edge_u] = init_face_dir_arr[0].index
                            vif[edge_l] = init_face_dir_arr[1].index
                            vif[edge_d] = init_face_dir_arr[2].index
                            vif[edge_r] = init_face_dir_arr[3].index
                            vif[process_sequence_attr] = process_sequence
                        initial_face_found = True
                else:
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
                    initial_face_found = True
            
            bpy.context.active_object.data.update()

            next_round = []
            init_face_corner_verts = []
            # endregion

            # region Process initial face(s))
            if initial_face_found == True:
                if init_face_is_virtual == True:
                    vif_verts = []
                    for vif in init_virtual_face:
                        vif[xi] = 0
                        vif[yi] = 0
                        if vif in faces_remain:
                            faces_remain.remove(vif)
                        for vifv in vif.verts:
                            if vifv not in vif_verts:
                                vif_verts.append(vifv)
                    init_face_corner_verts = get_corner_verts(init_face_dir_arr, vif_verts)
                    print(init_virtual_face_diagonal.index)
                    print([i.index for i in init_face_dir_arr])
                    print([v.index for v in init_face_corner_verts[0] if v != None])
                    for vif in init_virtual_face:
                        vif[corner_ur] = init_face_corner_verts[0][0].index
                        vif[corner_ul] = init_face_corner_verts[0][1].index
                        vif[corner_ll] = init_face_corner_verts[0][2].index
                        vif[corner_lr] = init_face_corner_verts[0][3].index
                else:
                    if init_face in faces_remain:
                        faces_remain.remove(init_face)
                    init_face_corner_verts = get_corner_verts(init_face_dir_arr, init_face.verts)
                    init_face[corner_ur] = init_face_corner_verts[0][0].index
                    init_face[corner_ul] = init_face_corner_verts[0][1].index
                    init_face[corner_ll] = init_face_corner_verts[0][2].index
                    init_face[corner_lr] = init_face_corner_verts[0][3].index
                    init_face[xi] = 0
                    init_face[yi] = 0

                avg_edge_length_x = (init_face_dir_arr[0].calc_length() + init_face_dir_arr[2].calc_length()) / 2
                avg_edge_length_y = (init_face_dir_arr[1].calc_length() + init_face_dir_arr[3].calc_length()) / 2
                if init_face_is_virtual == True:
                    for vif in init_virtual_face:
                        vif[avg_length_x] = avg_edge_length_x
                        vif[avg_length_y] = avg_edge_length_y
                if init_face_is_virtual == False:
                    init_face[avg_length_x] = avg_edge_length_x
                    init_face[avg_length_y] = avg_edge_length_y

                if self.pre_calc_edge_lengths == False:
                    avg_edge_length_x = (init_face_dir_arr[0].calc_length() + init_face_dir_arr[2].calc_length()) / 4
                    avg_edge_length_y = (init_face_dir_arr[1].calc_length() + init_face_dir_arr[3].calc_length()) / 4
                    
                    init_face_corner_vert_cos = [mathutils.Vector((-avg_edge_length_x, avg_edge_length_y)) + i_offset, mathutils.Vector((-avg_edge_length_x, -avg_edge_length_y)) + i_offset, mathutils.Vector((avg_edge_length_x, -avg_edge_length_y)) + i_offset, mathutils.Vector((avg_edge_length_x, avg_edge_length_y)) + i_offset]
                    for ifvco, ifva in zip(init_face_corner_vert_cos, co_attr_arr):
                        init_face[ifva[0]] = ifvco.x
                        init_face[ifva[1]] = ifvco.y
                        if init_face_is_virtual == True:
                            for vif in init_virtual_face:
                                vif[ifva[0]] = ifvco.x
                                vif[ifva[1]] = ifvco.y

                    for iv, ivi, ivco in zip(init_face_corner_verts[0], init_face_corner_verts[1], init_face_corner_vert_cos):
                        for ivl in iv.link_loops:
                            if ivl.index in island_loops:
                                ivl[uv_layer].uv = ivco
                                ivl[uv_layer].pin_uv = True
                        if iv not in verts_done:
                            verts_done.append(iv)
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
                
                for ci, corner in enumerate(init_face_corner_verts[0]):
                    for corner_loop in corner.link_loops:
                        if init_face_is_virtual == True:
                            if corner_loop in init_virtual_face[0].loops or corner_loop in init_virtual_face[1].loops:
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
                        if init_face_is_virtual == False:
                            if corner_loop in init_face.loops:
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
                
                init_face_winding = check_winding(init_face, init_face_corner_verts)
                if init_face_is_virtual == False:
                    print(f"Processing initial face of island {ii}: {init_face.index}. Not virtual.")
                if init_face_is_virtual == True:
                    print(f"Processing initial face of island {ii}: {init_virtual_face[0].index}/{init_virtual_face[1].index}. Virtual quad.")
                print(viz_quad(init_face_dir_arr, init_face_corner_verts[0], init_face, init_virtual_face, init_face_is_virtual, init_face_winding))

                if init_face not in faces_done:
                    faces_done.append(init_face)

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
                                    for otv in other_virtual_face[1]:
                                        if otv not in other_verts:
                                            other_verts.append(otv)
                                ifef[nom_dir] = di
                                ifef[nom] = adj_face.index
                                if ifef[trailing_tri] == 1:
                                    print(f"Initial face of island {ii}: skipped face index {ifef.index} in direction {di} because it is a trailing triangle.")
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
                                    print(f"Initial face of island {ii}: Added face index {ifef.index} in direction {di} to next round. Virtual face: {init_face_is_virtual}.")
                                    print(viz_add_quad(init_face_dir_arr, other_dir_edges, init_face_corner_verts[0], other_corner_verts[0], init_face, ifef, init_virtual_face, other_virtual_face, init_face_is_virtual, is_other_face_virtual, di))
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
            if initial_face_found == False:
                print(f"No usable initial face found for island {ii}.")
            # endregion

            bpy.context.active_object.data.update()
            round = 0
            do_a = True
            retries = 0
            max_retries = 5

            # region Mesh traversal
            tododo = len(faces_remain)
            while tododo > 0:
                print(f"{len(faces_remain)} remaining faces.")
                round += 1
                this_round = next_round.copy()
                next_round = []

                for face in faces_remain:
                    if face[trailing_tri] == 1:
                        faces_remain.remove(face)
                        face.select = True
                        for loop in face.loops:
                            loop[uv_layer].pin_uv = False
                        print(f"Removed face {face.index} from remaining faces because it is a trailing triangle.")
                for ttf in this_round:
                    if ttf[trailing_tri] == 1:
                        this_round.remove(ttf)
                        ttf.select = True
                        for ttfl in ttf.loops:
                            ttfl[uv_layer].pin_uv = False
                        print(f"Removed face {ttf.index} from round because it is a trailing triangle.")

                

                for trf in this_round:
                    print(f"Island {ii} - Round {round}: Processing face {trf.index} with direction {trf[nom_dir]}.")
                    these_edges = [bm.edges[trf[edge_u]], bm.edges[trf[edge_l]], bm.edges[trf[edge_d]], bm.edges[trf[edge_r]]]
                    is_virtual_quad = False
                    virtual_quad = [None, None]
                    virtual_verts = []
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
                    
                    if is_virtual_quad == True:
                        for vf in virtual_quad:
                            vf[avg_length_x] = avg_edge_length_x
                            vf[avg_length_y] = avg_edge_length_y
                            for vfv in vf.verts:
                                if vfv not in virtual_verts:
                                    virtual_verts.append(vfv)
                        corners = get_corner_verts(these_edges, virtual_verts)
                    if is_virtual_quad == False:
                        trf[avg_length_x] = avg_edge_length_x
                        trf[avg_length_y] = avg_edge_length_y
                        corners = get_corner_verts(these_edges, trf.verts)
                    
                    for ci, corner in enumerate(corners[0]):
                        for corner_loop in corner.link_loops:
                            if is_virtual_quad == True:
                                if corner_loop in virtual_quad[0].loops or corner_loop in virtual_quad[1].loops:
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
                            if is_virtual_quad == False:
                                if corner_loop in trf.loops:
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
                    if self.pre_calc_edge_lengths == False:
                        for cornerv, vi in zip(corners[0], corners[1]):
                            # if cornerv not in verts_done:
                            for vil in cornerv.link_loops:
                                if vil.index in island_loops:
                                    vil[uv_layer].uv = cos[vi]
                                    vil[uv_layer].pin_uv = True
                            if cornerv not in verts_done:
                                verts_done.append(cornerv)
                    if self.pre_calc_edge_lengths == True:
                        for trfvco, trfva in zip(cos, co_attr_arr):
                            trf[trfva[0]] = trfvco.x
                            trf[trfva[1]] = trfvco.y
                            if is_virtual_quad == True:
                                for vif in virtual_quad:
                                    vif[trfva[0]] = trfvco.x
                                    vif[trfva[1]] = trfvco.y
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
                                            print(f"Round {round}: Skipped face index {ifef.index} in direction {di} because it is a trailing triangle.")
                                            ifef[process_sequence_attr] = process_sequence
                                            ife[process_sequence_edge] = process_sequence
                                            ife[nom_dir_edge] = di
                                            process_sequence += 1
                                            faces_remain.remove(ifef)
                                            ifef.select = True
                                            for ifefl in ifef.loops:
                                                ifefl[uv_layer].pin_uv = False
                                        if ifef[trailing_tri] == 0:
                                            other_corner_verts = get_corner_verts(other_dir_arr, other_verts)
                                            print(f"Round {round}: Added face index {ifef.index} in direction {di} to next round of island {ii}. Virtual face: {other_face_is_virtual}. Current process sequence index: {process_sequence}")
                                            print(viz_add_quad(these_edges, other_dir_arr, corners[0], other_corner_verts[0], adj_face, ifef, virtual_quad, other_virtual_face, is_virtual_quad, other_face_is_virtual, di))
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
                                        other_dir_arr = get_other_dir_edges(these_edges, corners, adj_face, ifef, di)
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
                                            print(f"Round {round}: Skipped face index {ifef.index} in direction {di} because it is a trailing triangle.")
                                            ifef[process_sequence_attr] = process_sequence
                                            ife[process_sequence_edge] = process_sequence
                                            ife[nom_dir_edge] = di
                                            process_sequence += 1
                                            faces_remain.remove(ifef)
                                            ifef.select = True
                                            for ifefl in ifef.loops:
                                                ifefl[uv_layer].pin_uv = False
                                        if ifef[trailing_tri] == 0:
                                            other_corner_verts = get_corner_verts(other_dir_arr, other_verts)
                                            print(f"Round {round}: Added face index {ifef.index} in direction {di} to next round of island {ii}. Virtual face: {other_face_is_virtual}. Current process sequence index: {process_sequence}")
                                            print(viz_add_quad(these_edges, other_dir_arr, corners[0], other_corner_verts[0], adj_face, ifef, virtual_quad, other_virtual_face, is_virtual_quad, other_face_is_virtual, di))
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
                
                if len(next_round) == 0:
                    if self.unfold_mode == "B" or self.unfold_mode == "C":
                        if do_a == True:
                            do_a = False
                            print("No more in main axis, trying secondary axis")
                        elif do_a == False:
                            do_a = True
                            print("No more in secondary axis, returning to main axis")
                        next_round = faces_done.copy()
                    else:
                        if len(faces_remain) > 0:
                            next_round.append(faces_remain[0])
                        retries += 1

                tododo = len(faces_remain)

                bpy.context.active_object.data.update()
            
            # endregion

            # region Apply pre-calculated lengths
            if self.pre_calc_edge_lengths == True and (len(edge_length_map_x[0]) + len(edge_length_map_y[0])) > 0:
                sorted_x = edge_length_map_x[1].copy()
                sorted_x.sort(key=map_sort_by_coord)
                sorted_y = edge_length_map_y[1].copy()
                sorted_y.sort(key=map_sort_by_coord)

                previous_island_end_x += previous_island_length_x
                previous_island_end_y += previous_island_length_y

                list_a = sorted_x
                list_b = sorted_y

                a_len_addend = previous_island_end_x
                b_len_addend = previous_island_end_y

                if self.unfold_mode == "C":
                    list_a = sorted_y
                    list_b = sorted_x
                    a_len_addend = previous_island_end_y
                    b_len_addend = previous_island_end_x

                last_a_len = a_len_addend
                last_b_len = b_len_addend

                a_name = "column"
                b_name = "row"
                if self.unfold_mode == "C":
                    a_name = "row"
                    b_name = "column"

                print(f"Island {ii} has {len(list_a)} {a_name}s and {len(list_b)} {b_name}s ")

                for list_a_item in list_a:
                    a_len = list_a_item[2]/list_a_item[3]
                    old_a_len = 0.0 + last_a_len
                    new_a_len = 0.0 + last_a_len + a_len
                    print(f"Average edge length for {a_name} {list_a_item[0]}: {a_len}\nTotal previous length: {last_a_len}")
                    for face_data in list_a_item[1]:
                        for alc in face_data[3]:
                            for alcl in alc.link_loops:
                                if alcl.index in island_loops:
                                    if self.unfold_mode == "C":
                                        alcl[uv_layer].uv.y = old_a_len
                                    else:
                                        alcl[uv_layer].uv.x = old_a_len
                        for amc in face_data[4]:
                            for amcl in amc.link_loops:
                                if amcl.index in island_loops:
                                    if self.unfold_mode == "C":
                                        amcl[uv_layer].uv.y = new_a_len
                                    else:
                                        amcl[uv_layer].uv.x = new_a_len
                    last_a_len += a_len
                
                bpy.context.active_object.data.update()

                for list_b_item in list_b:
                    b_len = list_b_item[2]/list_b_item[3]
                    old_b_len = 0.0 + last_b_len
                    new_b_len = 0.0 + last_b_len + b_len
                    print(f"Average edge length for {b_name} {list_b_item[0]}: {b_len}\nTotal previous length: {last_b_len}")
                    for face_data in list_b_item[1]:
                        for blc in face_data[3]:
                            for blcl in blc.link_loops:
                                if blcl.index in island_loops:
                                    if self.unfold_mode == "C":
                                        blcl[uv_layer].uv.x = old_b_len
                                    else:
                                        blcl[uv_layer].uv.y = old_b_len
                                    blcl[uv_layer].pin_uv = True
                        for bmc in face_data[4]:
                            for bmcl in bmc.link_loops:
                                if bmcl.index in island_loops:
                                    if self.unfold_mode == "C":
                                        bmcl[uv_layer].uv.x = new_b_len
                                    else:
                                        bmcl[uv_layer].uv.y = new_b_len
                                    bmcl[uv_layer].pin_uv = True
                    last_b_len += b_len
                
                if self.unfold_mode == "C":
                    previous_island_length_x = self.offset_per_island[0] * last_b_len
                    previous_island_length_y = self.offset_per_island[1] * last_a_len
                else:
                    previous_island_length_x = self.offset_per_island[0] * last_a_len
                    previous_island_length_y = self.offset_per_island[1] * last_b_len
                bpy.context.active_object.data.update()
            # endregion
            for face in island:
                face[process_sequence_max] = process_sequence
        
        bpy.context.active_object.data.update()

        try:
            bpy.ops.uv.unwrap(method='MINIMUM_STRETCH', fill_holes=True, correct_aspect=True, use_subsurf_data=False, margin=0, no_flip=False, iterations=10, use_weights=False, weight_group="uv_importance", weight_factor=1)
        except:
            bpy.ops.uv.unwrap(method='ANGLE_BASED', fill_holes=True, correct_aspect=True, use_subsurf_data=False, margin=0, no_flip=False, iterations=10, use_weights=False, weight_group="uv_importance", weight_factor=1)

        return {"FINISHED"}

# class SloppyDeTriangulate(bpy.types.Operator):
#     bl_idname = "operator.sloppy_detri"
#     bl_label = "Detriangulate From Tri Pair"
#     bl_description = "Detriangulate mesh using selected tri pair as template"
#     bl_options = {'REGISTER', 'UNDO'}