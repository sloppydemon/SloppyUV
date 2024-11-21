import bpy
import bmesh
import mathutils
from bpy_extras import bmesh_utils

class SortVertByDist(bpy.types.Operator):
    bl_idname = "operator.sloppy_sort_vert_by_dist"
    bl_label = "Sort Vertices By Distance"
    bl_description = "Travel along edges until all vertices are visited and sort vertices according to the distance traveled."
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
            return vec.length + diff

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
    bl_description = "Travel along edges until all edges are visited and sort edges according to the distance traveled."
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
            return vec.length + diff

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
    bl_description = "Travel along faces until all faces are visited and sort faces according to the distance traveled."
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
            return vec.length + diff

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