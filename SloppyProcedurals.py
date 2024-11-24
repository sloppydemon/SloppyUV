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

    quant_avg_norm : bP(
        name = "Quantize Average Normal",
        description = "Quantize island's averaged normal to up, down, right, left",
        default = False
        ) # type: ignore

    def execute(self, context):
        props = context.scene.sloppy_props
        bm = bmesh.from_edit_mesh(bpy.context.active_object.data)
        uv_layer = bm.loops.layers.uv.verify()
        islands = bmesh_utils.bmesh_linked_uv_islands(bm, uv_layer)

        current_dir_array = []
        max_avg_angle = 0.0
        min_avg_angle = 360.0
        current_avg_normal = mathutils.Vector((0,0,0))
        current_normal = mathutils.Vector((0,0,0))

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
            {"name": "virtual_quad", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "vq_other_tri", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "vq_diagonal", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "trailing_tri", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "process_sequence", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "process_sequence_max", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "xi", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "yi", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "island_index", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "flatness", "type": "FLOAT", "domain": "FACE", "layer": None},
            {"name": "normal_regularity", "type": "FLOAT", "domain": "FACE", "layer": None},
            {"name": "nom", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "nom_dir", "type": "FLOAT_VECTOR", "domain": "FACE", "layer": None},
            {"name": "prev_length_x", "type": "FLOAT", "domain": "FACE", "layer": None},
            {"name": "prev_length_y", "type": "FLOAT", "domain": "FACE", "layer": None},
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
        virtual_quad = props.get_dict_layer("virtual_quad", attr_dict)
        vq_other_tri = props.get_dict_layer("vq_other_tri", attr_dict)
        trailing_tri = props.get_dict_layer("trailing_tri", attr_dict)
        process_sequence = props.get_dict_layer("process_sequence", attr_dict)
        process_sequence_max = props.get_dict_layer("process_sequence_max", attr_dict)

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

        def quant_sort(e):
            quant_dot = e[0].dot(current_normal)
            return quant_dot

        def regularity_sort(e):
            flat = (1.0 - e[flatness]) * -self.reg_flat_fac
            norm_reg = e[normal_regularity] * -self.reg_norm_fac
            return flat + norm_reg
        
        def adj_vert_angle_sort(e):
            return e[1]
        
        def dot_sort_u(e):
            return 1.0 - e.dot(current_dir_array[0])
        
        def dot_sort_l(e):
            return 1.0 - e.dot(current_dir_array[1])
        
        def dot_sort_d(e):
            return 1.0 - e.dot(current_dir_array[2])
        
        def dot_sort_r(e):
            return 1.0 - e.dot(current_dir_array[3])
        
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
                        adj_vert_angle = loop.calc_angle()
                        adj_vert_arr_app.append(adj_vert_angle)
                        adj_verts_arr.append(adj_vert_arr_app)
                adj_verts_arr.sort(key=adj_vert_angle_sort)
                for avf in adj_verts_arr[0][0].link_faces:
                    if avf != main_tri and avf != source_face:
                        main_tri_adjacent = False
                        diagonal_index = 0
                        for avfe in avf.edges:
                            if avfe not in edge_arr:
                                if main_tri in avfe.link_faces:
                                    main_tri_adjacent = True
                                    diagonal_index = avfe.index
                                    diagonal_edge = avfe
                        if main_tri_adjacent == True:
                            far_tri = avf
                            if len(far_tri.edges) > 3:
                                tri_trailing = True
                if tri_trailing == False:
                    if far_tri[island_index] == main_tri[island_index]:
                        main_tri[virtual_quad], far_tri[virtual_quad] = 1
                        main_tri[trailing_tri], far_tri[trailing_tri] = 0
                        main_tri[vq_other_tri] = far_tri.index
                        far_tri[vq_other_tri] = main_tri.index
                        main_tri[diagonal_index], avf[diagonal_index] = diagonal_edge.index
                    if far_tri[island_index] != main_tri[island_index]:
                        main_tri[trailing_tri] = 1
                        main_tri[virtual_quad] = 0
                        tri_trailing = True
                if tri_trailing == True:
                    main_tri[trailing_tri] = 1
                    main_tri[virtual_quad] = 0
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
                    main_tri[edge_u], far_tri[edge_u] = virtual_edges[0]
                    main_tri[edge_l], far_tri[edge_l] = virtual_edges[1]
                    main_tri[edge_d], far_tri[edge_d] = virtual_edges[2]
                    main_tri[edge_r], far_tri[edge_r] = virtual_edges[3]

            if sides == 4:
                other_face[virtual_quad] = 0
                for ofe in other_face.edges:
                    if ofe != adj_edge:
                        if dir == 0:
                            if ofe in corner_verts[0].link_edges:
                                out_array[3] = ofe
                            elif ofe in corner_verts[1].link_edges:
                                out_array[1] = ofe
                            else:
                                out_array[0] = ofe
                        if dir == 1:
                            if ofe in corner_verts[1].link_edges:
                                out_array[0] = ofe
                            elif ofe in corner_verts[2].link_edges:
                                out_array[2] = ofe
                            else:
                                out_array[1] = ofe
                        if dir == 2:
                            if ofe in corner_verts[2].link_edges:
                                out_array[1] = ofe
                            elif ofe in corner_verts[3].link_edges:
                                out_array[3] = ofe
                            else:
                                out_array[2] = ofe
                        if dir == 3:
                            if ofe in corner_verts[3].link_edges:
                                out_array[2] = ofe
                            elif ofe in corner_verts[0].link_edges:
                                out_array[0] = ofe
                            else:
                                out_array[3] = ofe
                other_face[edge_u] = out_array[0]
                other_face[edge_l] = out_array[1]
                other_face[edge_d] = out_array[2]
                other_face[edge_r] = out_array[3]

            return out_array
        
        def check_winding(face, corner_verts_arr):
            poss = 0
            negs = 0
            out_str = ""
            last_i = corner_verts_arr[1][-1]
            for vi in corner_verts_arr[1]:
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
            out_tri = None
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
                out_tri = tri_edges[0][2]
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
                ul, ll, lr, ur = 0
                for ve in v.link_edges:
                    if ve == edge_arr[0]:
                        ul += 1
                        ur += 1
                    if ve == edge_arr[1]:
                        ul += 1
                        ll += 1
                    if ve == edge_arr[2]:
                        ll += 1
                        lr += 1
                    if ve == edge_arr[3]:
                        lr += 1
                        ur += 1
                if ur == 2:
                    out_array[0] = v
                    out_array_is[0] = v_iter
                if ul == 2:
                    out_array[1] = v
                    out_array_is[1] = v_iter
                if ll == 2:
                    out_array[2] = v
                    out_array_is[2] = v_iter
                if lr == 2:
                    out_array[3] = v
                    out_array_is[3] = v_iter
                v_iter += 1
            out_arrays.append(out_array)
            out_arrays.append(out_array_is)
            return out_arrays

        for i, island in enumerate(islands):
            process_sequence = 0
            faces_remain = []
            faces_selected = []
            init_face = None
            init_virtual_face = []
            avg_norm = mathutils.Vector((0,0,0))

            for face in island:
                avg_norm += face.normal
                avg_angle = 0.0
                angle_count = 0
                for fe in face.edges:
                    angle_count += 1
                    if fe.is_boundary == True:
                        avg_angle += 90
                    else:
                        avg_angle += math.degrees(fe.calc_face_angle())
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
                face[island_index] = i
            
            bpy.context.active_object.data.update()

            is_selected_face = False
            init_face_is_virtual = False

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
                            if other_tri != None:
                                init_face_is_virtual = True
                                init_virtual_face.append(init_face)
                                init_virtual_face.append(other_tri)
                                is_selected_face = True
                            else:
                                print("Too inconsistent mesh! Aborting!")
                                return  {"CANCELED"}
                else:
                    faces_remain.sort(key=regularity_sort)
                    init_face = faces_remain[0]
                    if len(init_face.edges) == 3:
                        other_tri = get_adjacent_triangle(init_face)
                        if other_tri != None:
                            init_face_is_virtual = True
                            init_virtual_face.append(init_face)
                            init_virtual_face.append(other_tri)
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
                    if other_tri != None:
                        init_face_is_virtual = True
                        init_virtual_face.append(init_face)
                        init_virtual_face.append(other_tri)
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
            current_dir_array = quant_arr[0][1]

            if props.verbose == True:
                msg_start = "Most regular face for island "
                if is_selected_face == True:
                    msg_start = "Selected face for island "
                msg = "{:}{:}: {:} with a flatness of {:.3f} and a normal regularity of {:.3f}.\nFacing quantized to {:}\n".format(msg_start, i, init_face.index, init_face[flatness], init_face[normal_regularity], quant_arr[0][2])
                if init_face_is_virtual == True:
                    msg_start = "Most regular virtual face for island "
                    if is_selected_face == True:
                        msg_start = "Selected virtual face for island "
                    msg = "{:}{:}: {:}/{:} with a flatness of {:.3f} and a normal regularity of {:.3f}.\nFacing quantized to {:}\n".format(msg_start, i, init_virtual_face[0].index, init_virtual_face[1].index, init_face[flatness], init_face[normal_regularity], quant_arr[0][2])
                print(msg)
            
            if init_face_is_virtual == True:
                for vf in init_virtual_face:
                    vf.select = True
            else:
                init_face.select = True

        bpy.context.active_object.data.update()

        return {"FINISHED"}

# class SloppyDeTriangulate(bpy.types.Operator):
#     bl_idname = "operator.sloppy_detri"
#     bl_label = "Detriangulate From Tri Pair"
#     bl_description = "Detriangulate mesh using selected tri pair as template"
#     bl_options = {'REGISTER', 'UNDO'}