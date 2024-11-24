import bpy
import bmesh
import mathutils
from bpy_extras import bmesh_utils

class SortVertByDist(bpy.types.Operator):
    bl_idname = "operator.sloppy_sort_vert_by_dist"
    bl_label = "Sort Vertices By Distance"
    bl_description = "Travel along edges until all vertices are visited and sort vertices according to the distance traveled"
    bl_options = {'REGISTER', 'UNDO'}

    iP = bpy.props.IntProperty
    fP = bpy.props.FloatProperty
    fvP = bpy.props.FloatVectorProperty
    bP = bpy.props.BoolProperty
    eP = bpy.props.EnumProperty
    bvP = bpy.props.BoolVectorProperty
    sP = bpy.props.StringProperty

    bottom_up_mix : bP(
        name = "Weigh From Bottom Up",
        description = "Weigh sorting of potential next element with position in Z",
        default = True
        ) # type: ignore

    z_only : bP(
        name = "Sort By Z",
        description = "Sort potential next element by position in Z",
        default = True
        ) # type: ignore

    respect_islands : bP(
        name = "Respect UV Islands",
        description = "Respect mesh UV islands when sorting",
        default = False
        ) # type: ignore
    

    def execute(self, context):
        props = context.scene.sloppy_props
        bm = bmesh.from_edit_mesh(bpy.context.active_object.data)
        current_co = mathutils.Vector((0,0,0))
        current_island = 0
        islands = []
        uv_layer = bm.loops.layers.uv.verify()

        if self.respect_islands == True:
            islands = bmesh_utils.bmesh_linked_uv_islands(bm, uv_layer)

        attr_dict = [
            {"name": "vertsortdist", "type": "FLOAT", "domain": "POINT", "layer": None}
            ]
        
        for nam in attr_dict:
            attr = props.find_or_add_attribute(nam["name"], nam["type"], nam["domain"])
            nam["layer"] = props.get_attribute_layer(nam["name"], nam["type"], nam["domain"], bm)

        vertsortdist = props.get_dict_layer("vertsortdist", attr_dict)
        max_total_distance = 0.0

        def edge_length_sort(e):
            return e.calc_length()

        def geo_vert_dist_sort(e):
            vec = e.co - current_co
            return vec.length

        def edge_length_sort_z(e):
            curr_z = e.verts[0].co.lerp(e.verts[1].co, 0.5).z
            return e.calc_length() + curr_z

        def edge_sort_only_z(e):
            curr_z = e.verts[0].co.lerp(e.verts[1].co, 0.5).z
            return curr_z

        def total_distance_sort(e):
            return e[vertsortdist]

        def island_sort(e):
            return e[1]
        
        def island_check(v):
            v_islands = []
            for i, island in enumerate(islands):
                island_arr = []
                island_arr.append(i)
                island_face_count = 0
                for face in v.link_faces:
                    if face in island:
                        island_face_count += 1
                island_arr.append(island_face_count)
                v_islands.append(island_arr)
            v_islands.sort(key=island_sort)
            return v_islands[0][0]

        def geo_vert_dist_plus_island_sort(e):
            vec = e.co - current_co
            e_island = island_check(e)
            diff = abs(current_island - e_island)
            return vec.length + (diff * 100)

        bm.verts.ensure_lookup_table()

        verts_remain = []
        edges_done = []
        verts_selected = []

        for vert in bm.verts:
            verts_remain.append(vert)
            if vert.select == True:
                verts_selected.append(vert)

        next_vert = bm.verts[0]
        if len(verts_selected) > 0:
            next_vert = verts_selected[0]

        total_distance = 0.0
        last_edge = None

        while len(verts_remain) > 0:
            if last_edge != None:
                if last_edge not in edges_done:
                    edges_done.append(last_edge)
            this_vert = next_vert
            if props.verbose == True:
                print(f"Processing vert {this_vert.index}, current total distance: {total_distance}. {len(verts_remain)} verts remaining.")
            this_vert[vertsortdist] = total_distance
            new_max_dist = max(max_total_distance, total_distance)
            max_total_distance = new_max_dist
            has_boundary = False
            
            edges = []
            boundary_edges = []
            
            for edge in this_vert.link_edges:
                if edge not in edges_done:
                    this_island = 0
                    other_island = 0
                    if self.respect_islands == True:
                        this_island = island_check(this_vert)
                        other_island = island_check(edge.other_vert(this_vert))
                    if edge.other_vert(this_vert) in verts_remain and this_island == other_island:
                        edges.append(edge)
                        if edge.is_boundary == True:
                            boundary_edges.append(edge)
                            has_boundary = True
            
            
            if self.bottom_up_mix == True:
                if self.z_only == True:
                    edges.sort(key=edge_sort_only_z)
                    boundary_edges.sort(key=edge_sort_only_z)
                else:
                    edges.sort(key=edge_length_sort_z)
                    boundary_edges.sort(key=edge_length_sort_z)
            if self.bottom_up_mix == False:
                edges.sort(key=edge_length_sort)
                boundary_edges.sort(key=edge_length_sort)

            
            if (len(boundary_edges) + len(edges)) > 0:
                if has_boundary == True:
                    if props.verbose == True:
                        print("Traversing boundary.")
                    next_vert = boundary_edges[0].other_vert(this_vert)
                    total_distance += boundary_edges[0].calc_length()
                    last_edge = boundary_edges[0]
                else:
                    if props.verbose == True:
                        print("Traversing inner mesh.")
                    next_vert = edges[0].other_vert(this_vert)
                    total_distance += edges[0].calc_length()
                    last_edge = edges[0]
            else:
                current_co = this_vert.co
                if self.respect_islands == True:
                    current_island = island_check(this_vert)
                    verts_remain.sort(key=geo_vert_dist_plus_island_sort)
                else:
                    verts_remain.sort(key=geo_vert_dist_sort)
                # if self.bottom_up_mix:
                #     verts_remain.sort(key=geo_vert_dist_sort_z)
                next_vert = verts_remain[0]
                if self.respect_islands == True:
                    curr_vi = 0
                    
                vec = verts_remain[0].co - this_vert.co
                if props.verbose == True:
                    print(f"Jumping from {this_vert.index} to {next_vert.index} at a distance of {vec.length}")
                total_distance += vec.length
            
            if this_vert in verts_remain:
                verts_remain.remove(this_vert)

        bpy.context.active_object.data.update()

        bm.verts.sort(key=total_distance_sort)

        bm.verts.index_update()
        # for face in bm.faces:
        #     face.verts.index_update()

        bpy.context.active_object.data.update()

        return {"FINISHED"}

class SortEdgeByDist(bpy.types.Operator):
    bl_idname = "operator.sloppy_sort_edge_by_dist"
    bl_label = "Sort Edges By Distance"
    bl_description = "Travel along edges until all edges are visited and sort edges according to the distance traveled"
    bl_options = {'REGISTER', 'UNDO'}

    iP = bpy.props.IntProperty
    fP = bpy.props.FloatProperty
    fvP = bpy.props.FloatVectorProperty
    bP = bpy.props.BoolProperty
    eP = bpy.props.EnumProperty
    bvP = bpy.props.BoolVectorProperty
    sP = bpy.props.StringProperty

    bottom_up_mix : bP(
        name = "Weigh From Botom Up",
        description = "Weigh sorting of potential next element with position in Z",
        default = True
        ) # type: ignore

    z_only : bP(
        name = "Sort By Z",
        description = "Sort potential next element by position in Z",
        default = True
        ) # type: ignore

    respect_islands : bP(
        name = "Respect UV Islands",
        description = "Respect mesh UV islands when sorting",
        default = False
        ) # type: ignore

    def execute(self, context):
        props = context.scene.sloppy_props
        bm = bmesh.from_edit_mesh(bpy.context.active_object.data)
        current_co = mathutils.Vector((0,0,0))
        current_island = 0
        islands = []
        uv_layer = bm.loops.layers.uv.verify()

        if self.respect_islands == True:
            islands = bmesh_utils.bmesh_linked_uv_islands(bm, uv_layer)

        attr_dict = [
            {"name": "edgesortdist", "type": "FLOAT", "domain": "EDGE", "layer": None}
            ]
        
        for nam in attr_dict:
            attr = props.find_or_add_attribute(nam["name"], nam["type"], nam["domain"])
            nam["layer"] = props.get_attribute_layer(nam["name"], nam["type"], nam["domain"], bm)

        edgesortdist = props.get_dict_layer("edgesortdist", attr_dict)

        def edge_length_sort(e):
            return e.calc_length()

        def edge_length_sort_z(e):
            curr_z = e.verts[0].co.lerp(e.verts[1].co, 0.5).z
            return e.calc_length() + curr_z

        def edge_sort_only_z(e):
            curr_z = e.verts[0].co.lerp(e.verts[1].co, 0.5).z
            return curr_z

        def geo_edge_dist_sort(e):
            vec = e.verts[0].co.lerp(e.verts[1].co, 0.5) - current_co
            return vec.length

        def total_distance_sort(e):
            return e[edgesortdist]

        def island_sort(e):
            return e[1]
        
        def island_check(v):
            v_islands = []
            for i, island in enumerate(islands):
                island_arr = []
                island_arr.append(i)
                island_face_count = 0
                for vert in v.verts:
                    for face in vert.link_faces:
                        if face in island:
                            island_face_count += 1
                island_arr.append(island_face_count)
                v_islands.append(island_arr)
            v_islands.sort(key=island_sort)
            return v_islands[0][0]

        def geo_edge_dist_plus_island_sort(e):
            vec = e.verts[0].co.lerp(e.verts[1].co, 0.5) - current_co
            e_island = island_check(e)
            diff = abs(current_island - e_island)
            return vec.length + (diff * 100)

        bm.edges.ensure_lookup_table()

        edges_remain = []
        edges_done = []
        edges_selected = []

        for edge in bm.edges:
            edges_remain.append(edge)
            if edge.select == True:
                edges_selected.append(edge)

        next_edge = bm.edges[0]
        if len(edges_selected) > 0:
            next_edge = edges_selected[0]

        total_distance = 0.0
        last_edge = None

        while len(edges_remain) > 0:
            if last_edge != None:
                if last_edge not in edges_done:
                    edges_done.append(last_edge)
                if last_edge in edges_remain:
                    edges_remain.remove(last_edge)
            this_edge = next_edge
            if props.verbose == True:
                print(f"Processing edge {this_edge.index}, current total distance: {total_distance}. {len(edges_remain)} edges remaining.")
            this_edge[edgesortdist] = total_distance
            has_boundary = False
            found_other_edge = False
            
            edges = []
            boundary_edges = []
            
            for v in this_edge.verts:
                for ve in v.link_edges:
                    this_island = 0
                    other_island = 0
                    if self.respect_islands:
                        this_island = island_check(this_edge)
                        other_island = island_check(ve)
                    if ve not in edges_done and ve != this_edge and this_island == other_island:
                        edges.append(ve)
                        if ve.is_boundary == True:
                            boundary_edges.append(ve)
                            has_boundary = True
            
            if self.bottom_up_mix == True:
                if self.z_only == True:
                    edges.sort(key=edge_sort_only_z)
                    boundary_edges.sort(key=edge_sort_only_z)
                else:
                    edges.sort(key=edge_length_sort_z)
                    boundary_edges.sort(key=edge_length_sort_z)
            if self.bottom_up_mix == False:
                edges.sort(key=edge_length_sort)
                boundary_edges.sort(key=edge_length_sort)
            
            if (len(boundary_edges) + len(edges)) > 0:
                if has_boundary == True:
                    if props.verbose == True:
                        print("Traversing boundary.")
                    next_edge = boundary_edges[0]
                    total_distance += boundary_edges[0].calc_length()
                    last_edge = boundary_edges[0]
                    found_other_edge = True
                else:
                    if props.verbose == True:
                        print("Traversing inner mesh.")
                    next_edge = edges[0]
                    total_distance += edges[0].calc_length()
                    last_edge = edges[0]
                    found_other_edge = True
            
            if found_other_edge == False:
                current_co = this_edge.verts[0].co.lerp(this_edge.verts[1].co, 0.5)
                if self.respect_islands:
                    current_island = island_check(this_edge)
                    edges_remain.sort(key=geo_edge_dist_plus_island_sort)
                else:
                    edges_remain.sort(key=geo_edge_dist_sort)
                # if self.bottom_up_mix:
                #     edges_remain.sort(key=geo_edge_dist_sort_z)
                next_edge = edges_remain[0]
                vec = edges_remain[0].verts[0].co.lerp(this_edge.verts[1].co, 0.5) - this_edge.verts[0].co.lerp(this_edge.verts[1].co, 0.5)
                if props.verbose == True:
                    print(f"Jumping from {this_edge.index} to {next_edge.index} at a distance of {vec.length}")
                total_distance += vec.length
            
            if this_edge in edges_remain:
                edges_remain.remove(this_edge)

        bpy.context.active_object.data.update()

        bm.edges.sort(key=total_distance_sort)

        bm.edges.index_update()

        bpy.context.active_object.data.update()

        return {"FINISHED"}

class SortFaceByDist(bpy.types.Operator):
    bl_idname = "operator.sloppy_sort_face_by_dist"
    bl_label = "Sort Faces By Distance"
    bl_description = "Travel along faces until all faces are visited and sort faces according to the distance traveled"
    bl_options = {'REGISTER', 'UNDO'}

    iP = bpy.props.IntProperty
    fP = bpy.props.FloatProperty
    fvP = bpy.props.FloatVectorProperty
    bP = bpy.props.BoolProperty
    eP = bpy.props.EnumProperty
    bvP = bpy.props.BoolVectorProperty
    sP = bpy.props.StringProperty

    bottom_up_mix : bP(
        name = "Weigh From Botom Up",
        description = "Weigh sorting of potential next element with position in Z",
        default = True
        ) # type: ignore

    z_only : bP(
        name = "Sort By Z",
        description = "Sort potential next element by position in Z",
        default = True
        ) # type: ignore

    respect_islands : bP(
        name = "Respect UV Islands",
        description = "Respect mesh UV islands when sorting",
        default = False
        ) # type: ignore

    def execute(self, context):
        props = context.scene.sloppy_props
        bm = bmesh.from_edit_mesh(bpy.context.active_object.data)
        current_co = mathutils.Vector((0,0,0))
        current_island = 0
        islands = []
        uv_layer = bm.loops.layers.uv.verify()

        if self.respect_islands == True:
            islands = bmesh_utils.bmesh_linked_uv_islands(bm, uv_layer)

        attr_dict = [
            {"name": "facesortdist", "type": "FLOAT", "domain": "FACE", "layer": None}
            ]
        
        for nam in attr_dict:
            attr = props.find_or_add_attribute(nam["name"], nam["type"], nam["domain"])
            nam["layer"] = props.get_attribute_layer(nam["name"], nam["type"], nam["domain"], bm)

        facesortdist = props.get_dict_layer("facesortdist", attr_dict)

        def face_avg_length_sort(e):
            avg_length = 0.0
            num_lengths = 0
            for edge in e.edges:
                avg_length += edge.calc_length()
                num_lengths += 1
            avg_length /= num_lengths
            return avg_length

        def face_avg_length_sort_z(e):
            avg_length = 0.0
            num_lengths = 0
            for edge in e.edges:
                avg_length += edge.calc_length()
                num_lengths += 1
            avg_length /= num_lengths
            curr_z = e.calc_center_median().z
            return avg_length + curr_z

        def face_sort_only_z(e):
            curr_z = e.calc_center_median().z
            return curr_z

        def geo_face_dist_sort(e):
            vec = e.calc_center_median() - current_co
            return vec.length

        def total_distance_sort(e):
            return e[facesortdist]

        def island_check(f):
            f_island = 0
            for i, island in enumerate(islands):
                if f in island:
                    f_island = i
            return f_island

        def geo_face_dist_plus_island_sort(f):
            vec = f.calc_center_median() - current_co
            f_island = island_check(f)
            diff = abs(current_island - f_island)
            return vec.length + (diff * 100)

        bm.faces.ensure_lookup_table()

        faces_remain = []
        faces_done = []
        faces_selected = []

        for face in bm.faces:
            faces_remain.append(face)
            if face.select == True:
                faces_selected.append(face)

        next_face = bm.faces[0]
        if len(faces_selected) > 0:
            next_face = faces_selected[0]

        total_distance = 0.0
        last_face = None

        while len(faces_remain) > 0:
            if last_face != None:
                if last_face not in faces_done:
                    faces_done.append(last_face)
            this_face = next_face
            if props.verbose == True:
                print(f"Processing face {this_face.index}, current total distance: {total_distance}. {len(faces_remain)} faces remaining.")
            this_face[facesortdist] = total_distance
            has_boundary = False
            
            potential_faces = []
            boundary_faces = []
            
            for e in this_face.edges:
                for f in e.link_faces:
                    this_island = 0
                    other_island = 0
                    if self.respect_islands:
                        this_island = island_check(this_face)
                        other_island = island_check(f)
                    if f not in faces_done and this_island == other_island:
                        potential_faces.append(f)
                        for fe in f.edges:
                            if fe.is_boundary == True:
                                if f not in boundary_faces:
                                    boundary_faces.append(f)
                                has_boundary = True
            
            
            if self.bottom_up_mix == True:
                if self.z_only == True:
                    potential_faces.sort(key=face_sort_only_z)
                    boundary_faces.sort(key=face_sort_only_z)
                else:
                    potential_faces.sort(key=face_avg_length_sort_z)
                    boundary_faces.sort(key=face_avg_length_sort_z)
            if self.bottom_up_mix == False:
                potential_faces.sort(key=face_avg_length_sort)
                boundary_faces.sort(key=face_avg_length_sort)
            
            if (len(boundary_faces) + len(potential_faces)) > 0:
                if has_boundary == True:
                    if props.verbose == True:
                        print("Traversing boundary.")
                    next_face = boundary_faces[0]
                    vec = boundary_faces[0].calc_center_median() - this_face.calc_center_median()
                    total_distance += vec.length
                    last_face = boundary_faces[0]
                else:
                    if props.verbose == True:
                        print("Traversing inner mesh.")
                    next_face = potential_faces[0]
                    vec = potential_faces[0].calc_center_median() - this_face.calc_center_median()
                    total_distance += vec.length
                    last_face = potential_faces[0]
            else:
                current_co = this_face.calc_center_median()
                if self.respect_islands == True:
                    current_island = island_check(this_face)
                    faces_remain.sort(key=geo_face_dist_plus_island_sort)
                else:
                    faces_remain.sort(key=geo_face_dist_sort)
                # if self.bottom_up_mix:
                #     faces_remain.sort(key=geo_face_dist_sort_z)
                next_face = faces_remain[0]
                vec = faces_remain[0].calc_center_median() - this_face.calc_center_median()
                if props.verbose == True:
                    print(f"Jumping from {this_face.index} to {next_face.index} at a distance of {vec.length}")
                total_distance += vec.length
            
            if this_face in faces_remain:
                faces_remain.remove(this_face)

        bpy.context.active_object.data.update()

        bm.faces.sort(key=total_distance_sort)

        bm.faces.index_update()

        bpy.context.active_object.data.update()

        return {"FINISHED"}

class SortVertByFacing(bpy.types.Operator):
    bl_idname = "operator.sloppy_sort_vert_by_facing"
    bl_label = "Sort Vertices By Facing"
    bl_description = "Travel along edges until all vertices are visited, selecting next element according to facing, and sort vertices according to the distance traveled"
    bl_options = {'REGISTER', 'UNDO'}

    iP = bpy.props.IntProperty
    fP = bpy.props.FloatProperty
    fvP = bpy.props.FloatVectorProperty
    bP = bpy.props.BoolProperty
    eP = bpy.props.EnumProperty
    bvP = bpy.props.BoolVectorProperty
    sP = bpy.props.StringProperty

    compare_facing : eP(
        name = "Compare With Vector",
        description = "Vector to compare facing with",
        items = [
            ("A", "View", "Compare with vector pointing towards view"),
            ("B", "Selection Average", "Compare with selection average normal"),
            ("C", "Selection Average Quantized", "Compare with selection average normal, quantized to negative or positive X, Y or Z"),
            ("D", "Island Average", "Compare with UV island average"),
            ("E", "Island Average Quantized", "Compare with UV island average, quantized to negative or positive X, Y or Z")
            ]
        ) # type: ignore

    respect_islands : bP(
        name = "Respect UV Islands",
        description = "Attempt to respect mesh UV islands when sorting, even when not using island average",
        default = False
        ) # type: ignore
    

    def execute(self, context):
        props = context.scene.sloppy_props
        bm = bmesh.from_edit_mesh(bpy.context.active_object.data)
        current_co = mathutils.Vector((0,0,0))
        current_norm = mathutils.Vector((0,0,0))
        islands = []
        uv_layer = bm.loops.layers.uv.verify()

        quant_arr = [
            [mathutils.Vector((1,0,0)), "Left"],
            [mathutils.Vector((-1,0,0)), "Right"],
            [mathutils.Vector((0,1,0)), "Back"],
            [mathutils.Vector((0,-1,0)), "Front"],
            [mathutils.Vector((0,0,1)), "Down"],
            [mathutils.Vector((0,0,-1)), "Up"],
            ]

        if self.respect_islands == True or self.compare_facing == "D" or self.compare_facing == "E":
            islands = bmesh_utils.bmesh_linked_uv_islands(bm, uv_layer)
        else:
            islands = [bm.faces]

        attr_dict = [
            {"name": "vertsortdist", "type": "FLOAT", "domain": "POINT", "layer": None}
            ]
        
        for nam in attr_dict:
            attr = props.find_or_add_attribute(nam["name"], nam["type"], nam["domain"])
            nam["layer"] = props.get_attribute_layer(nam["name"], nam["type"], nam["domain"], bm)

        vertsortdist = props.get_dict_layer("vertsortdist", attr_dict)
        max_total_distance = 0.0

        def dot_sort(e):
            norm = e.normal
            ndot = norm.dot(current_norm)
            return ndot
        
        def quant_sort(e):
            return current_norm.dot(e[0])

        def geo_vert_dist_sort(e):
            vec = e.co - current_co
            return vec.length

        def edge_sort_only(e):
            return dot_sort(e[1])

        def total_distance_sort(e):
            return e[vertsortdist]

        bm.verts.ensure_lookup_table()

        total_distance = 0.0

        for island in islands:
            total_distance += 0.01
            verts_remain = []
            edges_done = []
            verts_selected = []
            island_norm = mathutils.Vector((0,0,0))
            selected_norm = mathutils.Vector((0,0,0))
            island_min_co = mathutils.Vector((9999,9999))
            island_max_co = mathutils.Vector((-9999,-99990))
            selected_min_co = mathutils.Vector((9999,9999))
            selected_max_co = mathutils.Vector((-9999,-99990))

            for face in island:
                for fv in face.verts:
                    if fv not in verts_remain:
                        verts_remain.append(fv)
                        island_norm += fv.normal
                        if self.compare_facing == "A":
                            fvco = props.viewco(fv.co)
                            new_island_min_co_x = min(island_min_co.x, fvco.x)
                            new_island_min_co_y = min(island_min_co.y, fvco.y)
                            new_island_max_co_x = max(island_max_co.x, fvco.x)
                            new_island_max_co_y = max(island_max_co.y, fvco.y)
                            island_min_co.x = new_island_min_co_x
                            island_min_co.y = new_island_min_co_y
                            island_max_co.x = new_island_max_co_x
                            island_max_co.y = new_island_max_co_y
                    if fv.select == True and fv not in verts_selected:
                        verts_selected.append(fv)
                        selected_norm += fv.normal
                        if self.compare_facing == "A":
                            fvco = props.viewco(fv.co)
                            new_sel_min_co_x = min(selected_min_co.x, fvco.x)
                            new_sel_min_co_y = min(selected_min_co.y, fvco.y)
                            new_sel_max_co_x = max(selected_max_co.x, fvco.x)
                            new_sel_max_co_y = max(selected_max_co.y, fvco.y)
                            selected_min_co.x = new_sel_min_co_x
                            selected_min_co.y = new_sel_min_co_y
                            selected_max_co.x = new_sel_max_co_x
                            selected_max_co.y = new_sel_max_co_y
            
            if self.compare_facing == "A":
                island_center = island_min_co.lerp(island_max_co, 0.5)
                if len(verts_selected) > 0:
                    island_center = selected_min_co.lerp(selected_max_co, 0.5)
                current_norm = props.viewvec(props.viewco(island_center))
            if self.compare_facing == "B":
                current_norm = island_norm.normalized() * -1
                if len(verts_selected) > 0:
                    current_norm = selected_norm.normalized() * -1
            if self.compare_facing == "C":
                current_norm = island_norm.normalized() * -1
                if len(verts_selected) > 0:
                    current_norm = selected_norm.normalized() * -1
                quant_arr.sort(key=quant_sort)
                current_norm = quant_arr[0]
            if self.compare_facing == "D":
                current_norm = island_norm.normalized() * -1
            if self.compare_facing == "E":
                current_norm = island_norm.normalized() * -1
                quant_arr.sort(key=quant_sort)
                current_norm = quant_arr[0]
            
            verts_remain.sort(key=dot_sort)
            verts_selected.sort(key=dot_sort)
            next_vert = verts_remain[0]
            if self.compare_facing == "A" or self.compare_facing == "B" or self.compare_facing == "C":
                next_vert = verts_selected[0]
            
            last_edge = None

            while len(verts_remain) > 0:
                if last_edge != None:
                    if last_edge not in edges_done:
                        edges_done.append(last_edge)
                this_vert = next_vert
                if props.verbose == True:
                    print(f"Processing vert {this_vert.index}, current total distance: {total_distance}. {len(verts_remain)} verts remaining.")
                this_vert[vertsortdist] = total_distance
                new_max_dist = max(max_total_distance, total_distance)
                max_total_distance = new_max_dist

                edges = []
            
                for edge in this_vert.link_edges:
                    if edge not in edges_done:
                        if edge.other_vert(this_vert) in verts_remain:
                            edge_arr = []
                            edge_arr.append(edge)
                            edge_arr.append(edge.other_vert(this_vert))
                            edges.append(edge_arr)
                edges.sort(key=edge_sort_only)

                if len(edges) > 0:
                    next_vert = edges[0][1]
                    total_distance += edges[0][0].calc_length()
                    last_edge = edges[0][0]
                else:
                    current_co = this_vert.co
                    verts_remain.sort(key=geo_vert_dist_sort)
                    next_vert = verts_remain[0]
                    vec = verts_remain[0].co - this_vert.co
                    if props.verbose == True:
                        print(f"Jumping from {this_vert.index} to {next_vert.index} at a distance of {vec.length}")
                    total_distance += vec.length
                
                if this_vert in verts_remain:
                    verts_remain.remove(this_vert)
            
            bpy.context.active_object.data.update()

        bm.verts.sort(key=total_distance_sort)

        bm.verts.index_update()

        bpy.context.active_object.data.update()

        return {"FINISHED"}

class SortEdgeByFacing(bpy.types.Operator):
    bl_idname = "operator.sloppy_sort_edge_by_facing"
    bl_label = "Sort Edges By Facing"
    bl_description = "Travel along edges until all edges are visited, selecting next element according to facing, and sort edges according to the distance traveled"
    bl_options = {'REGISTER', 'UNDO'}

    iP = bpy.props.IntProperty
    fP = bpy.props.FloatProperty
    fvP = bpy.props.FloatVectorProperty
    bP = bpy.props.BoolProperty
    eP = bpy.props.EnumProperty
    bvP = bpy.props.BoolVectorProperty
    sP = bpy.props.StringProperty

    compare_facing : eP(
        name = "Compare With Vector",
        description = "Vector to compare facing with",
        items = [
            ("A", "View", "Compare with vector pointing towards view"),
            ("B", "Selection Average", "Compare with selection average normal"),
            ("C", "Selection Average Quantized", "Compare with selection average normal, quantized to negative or positive X, Y or Z"),
            ("D", "Island Average", "Compare with UV island average"),
            ("E", "Island Average Quantized", "Compare with UV island average, quantized to negative or positive X, Y or Z")
            ]
        ) # type: ignore

    respect_islands : bP(
        name = "Respect UV Islands",
        description = "Attempt to respect mesh UV islands when sorting, even when not using island average",
        default = False
        ) # type: ignore
    

    def execute(self, context):
        props = context.scene.sloppy_props
        bm = bmesh.from_edit_mesh(bpy.context.active_object.data)
        current_co = mathutils.Vector((0,0,0))
        current_norm = mathutils.Vector((0,0,0))
        islands = []
        uv_layer = bm.loops.layers.uv.verify()

        quant_arr = [
            [mathutils.Vector((1,0,0)), "Left"],
            [mathutils.Vector((-1,0,0)), "Right"],
            [mathutils.Vector((0,1,0)), "Back"],
            [mathutils.Vector((0,-1,0)), "Front"],
            [mathutils.Vector((0,0,1)), "Down"],
            [mathutils.Vector((0,0,-1)), "Up"],
            ]

        if self.respect_islands == True or self.compare_facing == "D" or self.compare_facing == "E":
            islands = bmesh_utils.bmesh_linked_uv_islands(bm, uv_layer)
        else:
            islands = [bm.faces]

        attr_dict = [
            {"name": "edgesortdist", "type": "FLOAT", "domain": "POINT", "layer": None}
            ]
        
        for nam in attr_dict:
            attr = props.find_or_add_attribute(nam["name"], nam["type"], nam["domain"])
            nam["layer"] = props.get_attribute_layer(nam["name"], nam["type"], nam["domain"], bm)

        edgesortdist = props.get_dict_layer("edgesortdist", attr_dict)
        max_total_distance = 0.0

        def calc_edge_normal(e):
            norms = e.verts[0].normal + e.verts[1].normal
            return norms.normalized()
        
        def calc_edge_co(e):
            return e.verts[0].co.lerp(e.verts[1].co, 0.5)

        def dot_sort(e):
            norm = calc_edge_normal(e)
            ndot = norm.dot(current_norm)
            return ndot
        
        def quant_sort(e):
            return current_norm.dot(e[0])

        def geo_edge_dist_sort(e):
            vec = calc_edge_co(e) - current_co
            return vec.length

        def total_distance_sort(e):
            return e[edgesortdist]

        bm.edges.ensure_lookup_table()

        total_distance = 0.0

        for island in islands:
            total_distance += 0.01
            edges_remain = []
            edges_done = []
            edges_selected = []
            island_norm = mathutils.Vector((0,0,0))
            selected_norm = mathutils.Vector((0,0,0))
            island_min_co = mathutils.Vector((9999,9999))
            island_max_co = mathutils.Vector((-9999,-99990))
            selected_min_co = mathutils.Vector((9999,9999))
            selected_max_co = mathutils.Vector((-9999,-99990))

            for face in island:
                for fe in face.edges:
                    if fe not in edges_remain:
                        edges_remain.append(fe)
                        island_norm += calc_edge_normal(fe)
                        if self.compare_facing == "A":
                            feco = props.viewco(calc_edge_co(fe))
                            new_island_min_co_x = min(island_min_co.x, feco.x)
                            new_island_min_co_y = min(island_min_co.y, feco.y)
                            new_island_max_co_x = max(island_max_co.x, feco.x)
                            new_island_max_co_y = max(island_max_co.y, feco.y)
                            island_min_co.x = new_island_min_co_x
                            island_min_co.y = new_island_min_co_y
                            island_max_co.x = new_island_max_co_x
                            island_max_co.y = new_island_max_co_y
                    if fe.select == True and fe not in edges_selected:
                        edges_selected.append(fe)
                        selected_norm += calc_edge_normal(fe)
                        if self.compare_facing == "A":
                            feco = props.viewco(calc_edge_co(fe))
                            new_sel_min_co_x = min(selected_min_co.x, feco.x)
                            new_sel_min_co_y = min(selected_min_co.y, feco.y)
                            new_sel_max_co_x = max(selected_max_co.x, feco.x)
                            new_sel_max_co_y = max(selected_max_co.y, feco.y)
                            selected_min_co.x = new_sel_min_co_x
                            selected_min_co.y = new_sel_min_co_y
                            selected_max_co.x = new_sel_max_co_x
                            selected_max_co.y = new_sel_max_co_y
            
            if self.compare_facing == "A":
                island_center = island_min_co.lerp(island_max_co, 0.5)
                if len(edges_selected) > 0:
                    island_center = selected_min_co.lerp(selected_max_co, 0.5)
                current_norm = props.viewvec(props.viewco(island_center))
            if self.compare_facing == "B":
                current_norm = island_norm.normalized() * -1
                if len(edges_selected) > 0:
                    current_norm = selected_norm.normalized() * -1
            if self.compare_facing == "C":
                current_norm = island_norm.normalized() * -1
                if len(edges_selected) > 0:
                    current_norm = selected_norm.normalized() * -1
                quant_arr.sort(key=quant_sort)
                current_norm = quant_arr[0]
            if self.compare_facing == "D":
                current_norm = island_norm.normalized() * -1
            if self.compare_facing == "E":
                current_norm = island_norm.normalized() * -1
                quant_arr.sort(key=quant_sort)
                current_norm = quant_arr[0]
            
            edges_remain.sort(key=dot_sort)
            edges_selected.sort(key=dot_sort)
            next_edge = edges_remain[0]
            if self.compare_facing == "A" or self.compare_facing == "B" or self.compare_facing == "C":
                next_vert = edges_selected[0]
            
            last_edge = None

            while len(edges_remain) > 0:
                if last_edge != None:
                    if last_edge not in edges_done:
                        edges_done.append(last_edge)
                    if last_edge in edges_remain:
                        edges_remain.remove(last_edge)
                this_edge = next_edge
                if props.verbose == True:
                    print(f"Processing edge {this_edge.index}, current total distance: {total_distance}. {len(edges_remain)} edges remaining.")
                this_edge[edgesortdist] = total_distance
                new_max_dist = max(max_total_distance, total_distance)
                max_total_distance = new_max_dist
                found_other_edge = False

                edges = []
            
                for vert in this_edge.verts:
                    for ev in vert.link_edges:
                        if ev not in edges_done:
                            if ev in edges_remain:
                                edges.append(ev)
                edges.sort(key=dot_sort)

                if len(edges) > 0:
                    next_edge = edges[0]
                    total_distance += edges[0].calc_length()
                    last_edge = edges[0]
                    found_other_edge = True

                if found_other_edge == False:
                    current_co = calc_edge_co(this_edge)
                    edges_remain.sort(key=geo_edge_dist_sort)
                    next_edge = edges_remain[0]
                    new_length = geo_edge_dist_sort(next_edge)
                    if props.verbose == True:
                        print(f"Jumping from {this_edge.index} to {next_edge.index} at a distance of {new_length}")
                    total_distance += new_length
                
                if this_edge in edges_remain:
                    edges_remain.remove(this_edge)
            
            bpy.context.active_object.data.update()

        bm.edges.sort(key=total_distance_sort)

        bm.edges.index_update()

        bpy.context.active_object.data.update()

        return {"FINISHED"}

class SortFaceByFacing(bpy.types.Operator):
    bl_idname = "operator.sloppy_sort_face_by_facing"
    bl_label = "Sort Faces By Facing"
    bl_description = "Travel along faces until all faces are visited, selecting next element according to facing, and sort faces according to the distance traveled"
    bl_options = {'REGISTER', 'UNDO'}

    iP = bpy.props.IntProperty
    fP = bpy.props.FloatProperty
    fvP = bpy.props.FloatVectorProperty
    bP = bpy.props.BoolProperty
    eP = bpy.props.EnumProperty
    bvP = bpy.props.BoolVectorProperty
    sP = bpy.props.StringProperty

    compare_facing : eP(
        name = "Compare With Vector",
        description = "Vector to compare facing with",
        items = [
            ("A", "View", "Compare with vector pointing towards view"),
            ("B", "Selection Average", "Compare with selection average normal"),
            ("C", "Selection Average Quantized", "Compare with selection average normal, quantized to negative or positive X, Y or Z"),
            ("D", "Island Average", "Compare with UV island average"),
            ("E", "Island Average Quantized", "Compare with UV island average, quantized to negative or positive X, Y or Z")
            ]
        ) # type: ignore

    respect_islands : bP(
        name = "Respect UV Islands",
        description = "Attempt to respect mesh UV islands when sorting, even when not using island average",
        default = False
        ) # type: ignore
    

    def execute(self, context):
        props = context.scene.sloppy_props
        bm = bmesh.from_edit_mesh(bpy.context.active_object.data)
        current_co = mathutils.Vector((0,0,0))
        current_norm = mathutils.Vector((0,0,0))
        islands = []
        uv_layer = bm.loops.layers.uv.verify()

        quant_arr = [
            [mathutils.Vector((1,0,0)), "Left"],
            [mathutils.Vector((-1,0,0)), "Right"],
            [mathutils.Vector((0,1,0)), "Back"],
            [mathutils.Vector((0,-1,0)), "Front"],
            [mathutils.Vector((0,0,1)), "Down"],
            [mathutils.Vector((0,0,-1)), "Up"],
            ]

        if self.respect_islands == True or self.compare_facing == "D" or self.compare_facing == "E":
            islands = bmesh_utils.bmesh_linked_uv_islands(bm, uv_layer)
        else:
            islands = [bm.faces]

        attr_dict = [
            {"name": "facesortdist", "type": "FLOAT", "domain": "POINT", "layer": None}
            ]
        
        for nam in attr_dict:
            attr = props.find_or_add_attribute(nam["name"], nam["type"], nam["domain"])
            nam["layer"] = props.get_attribute_layer(nam["name"], nam["type"], nam["domain"], bm)

        facesortdist = props.get_dict_layer("facesortdist", attr_dict)
        max_total_distance = 0.0
        
        def calc_face_dist(fa, fb):
            vec = fb.calc_center_median() - fa.calc_center_median()
            return vec.length

        def dot_sort(e):
            norm = e.normal
            ndot = norm.dot(current_norm)
            return ndot
        
        def quant_sort(e):
            return current_norm.dot(e[0])

        def geo_face_dist_sort(e):
            vec = e.calc_center_median() - current_co
            return vec.length

        def total_distance_sort(e):
            return e[facesortdist]

        bm.faces.ensure_lookup_table()

        total_distance = 0.0

        for island in islands:
            total_distance += 0.01
            faces_remain = []
            faces_done = []
            faces_selected = []
            island_norm = mathutils.Vector((0,0,0))
            selected_norm = mathutils.Vector((0,0,0))
            island_min_co = mathutils.Vector((9999,9999))
            island_max_co = mathutils.Vector((-9999,-99990))
            selected_min_co = mathutils.Vector((9999,9999))
            selected_max_co = mathutils.Vector((-9999,-99990))

            for face in island:
                if face not in faces_remain:
                    faces_remain.append(face)
                    island_norm += face.normal
                    if self.compare_facing == "A":
                        fco = face.calc_center_median()
                        new_island_min_co_x = min(island_min_co.x, fco.x)
                        new_island_min_co_y = min(island_min_co.y, fco.y)
                        new_island_max_co_x = max(island_max_co.x, fco.x)
                        new_island_max_co_y = max(island_max_co.y, fco.y)
                        island_min_co.x = new_island_min_co_x
                        island_min_co.y = new_island_min_co_y
                        island_max_co.x = new_island_max_co_x
                        island_max_co.y = new_island_max_co_y
                    if face.select == True and face not in faces_selected:
                        faces_selected.append(face)
                        selected_norm += face.normal
                        if self.compare_facing == "A":
                            fco = face.calc_center_median()
                            new_sel_min_co_x = min(selected_min_co.x, fco.x)
                            new_sel_min_co_y = min(selected_min_co.y, fco.y)
                            new_sel_max_co_x = max(selected_max_co.x, fco.x)
                            new_sel_max_co_y = max(selected_max_co.y, fco.y)
                            selected_min_co.x = new_sel_min_co_x
                            selected_min_co.y = new_sel_min_co_y
                            selected_max_co.x = new_sel_max_co_x
                            selected_max_co.y = new_sel_max_co_y
            
            if self.compare_facing == "A":
                island_center = island_min_co.lerp(island_max_co, 0.5)
                if len(faces_selected) > 0:
                    island_center = selected_min_co.lerp(selected_max_co, 0.5)
                current_norm = props.viewvec(props.viewco(island_center))
            if self.compare_facing == "B":
                current_norm = island_norm.normalized() * -1
                if len(faces_selected) > 0:
                    current_norm = selected_norm.normalized() * -1
            if self.compare_facing == "C":
                current_norm = island_norm.normalized() * -1
                if len(faces_selected) > 0:
                    current_norm = selected_norm.normalized() * -1
                quant_arr.sort(key=quant_sort)
                current_norm = quant_arr[0]
            if self.compare_facing == "D":
                current_norm = island_norm.normalized() * -1
            if self.compare_facing == "E":
                current_norm = island_norm.normalized() * -1
                quant_arr.sort(key=quant_sort)
                current_norm = quant_arr[0]
            
            faces_remain.sort(key=dot_sort)
            faces_selected.sort(key=dot_sort)
            next_face = faces_remain[0]
            if self.compare_facing == "A" or self.compare_facing == "B" or self.compare_facing == "C":
                next_face = faces_selected[0]
            
            last_face = None

            while len(faces_remain) > 0:
                if last_face != None:
                    if last_face not in faces_done:
                        faces_done.append(last_face)
                    if last_face in faces_remain:
                        faces_remain.remove(last_face)
                this_face = next_face
                if props.verbose == True:
                    print(f"Processing edge {this_face.index}, current total distance: {total_distance}. {len(faces_remain)} edges remaining.")
                this_face[facesortdist] = total_distance
                new_max_dist = max(max_total_distance, total_distance)
                max_total_distance = new_max_dist
                found_other_face = False

                potential_faces = []
            
                for edge in this_face.edges:
                    for ef in edge.link_faces:
                        if ef not in faces_done:
                            if ef in faces_remain:
                                potential_faces.append(ef)
                potential_faces.sort(key=dot_sort)

                if len(potential_faces) > 0:
                    next_face = potential_faces[0]
                    total_distance += calc_face_dist(this_face, next_face)
                    last_face = potential_faces[0]
                    found_other_face = True

                if found_other_face == False:
                    current_co = this_face.calc_center_median()
                    faces_remain.sort(key=geo_face_dist_sort)
                    next_face = faces_remain[0]
                    new_length = calc_face_dist(this_face, next_face)
                    if props.verbose == True:
                        print(f"Jumping from {this_face.index} to {next_face.index} at a distance of {new_length}")
                    total_distance += new_length
                
                if this_face in faces_remain:
                    faces_remain.remove(this_face)
            
            bpy.context.active_object.data.update()

        bm.faces.sort(key=total_distance_sort)

        bm.faces.index_update()

        bpy.context.active_object.data.update()

        return {"FINISHED"}