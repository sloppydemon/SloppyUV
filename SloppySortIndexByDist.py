import bpy
import bmesh
import mathutils


class SortVertByDist(bpy.types.Operator):
    bl_idname = "operator.sloppy_sort_vert_by_dist"
    bl_label = "Sort Vertices By Distance"
    bl_description = "Travel along edges until all vertices are visited and sort vertices according to the distance traveled."
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        props = context.scene.sloppy_props
        bm = bmesh.from_edit_mesh(bpy.context.active_object.data)
        current_co = mathutils.Vector((0,0,0))

        attr_dict = [
            {"name": "vertsortdist", "type": "FLOAT", "domain": "POINT", "layer": None}
            ]
        
        for nam in attr_dict:
            attr = props.find_or_add_attribute(nam["name"], nam["type"], nam["domain"])
            nam["layer"] = props.get_attribute_layer(nam["name"], nam["type"], nam["domain"], bm)

        def edge_length_sort(e):
            return e.calc_length()

        def geo_vert_dist_sort(e):
            vec = e.co - current_co
            return vec.length

        def total_distance_sort(e):
            return e[props.get_dict_layer("vertsortdist", attr_dict)]

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
            this_vert[props.get_dict_layer("vertsortdist", attr_dict)] = total_distance
            has_boundary = False
            
            edges = []
            boundary_edges = []
            
            for edge in this_vert.link_edges:
                if edge not in edges_done:
                    if edge.other_vert(this_vert) in verts_remain:
                        edges.append(edge)
                        if edge.is_boundary == True:
                            boundary_edges.append(edge)
                            has_boundary = True
            
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
        # for face in bm.faces:
        #     face.verts.index_update()

        bpy.context.active_object.data.update()

        return {"FINISHED"}

class SortEdgeByDist(bpy.types.Operator):
    bl_idname = "operator.sloppy_sort_edge_by_dist"
    bl_label = "Sort Edges By Distance"
    bl_description = "Travel along edges until all edges are visited and sort edges according to the distance traveled."
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        props = context.scene.sloppy_props
        bm = bmesh.from_edit_mesh(bpy.context.active_object.data)
        current_co = mathutils.Vector((0,0,0))

        attr_dict = [
            {"name": "edgesortdist", "type": "FLOAT", "domain": "EDGE", "layer": None}
            ]
        
        for nam in attr_dict:
            attr = props.find_or_add_attribute(nam["name"], nam["type"], nam["domain"])
            nam["layer"] = props.get_attribute_layer(nam["name"], nam["type"], nam["domain"], bm)

        def edge_length_sort(e):
            return e.calc_length()

        def geo_edge_dist_sort(e):
            vec = e.verts[0].co.lerp(this_edge.verts[1].co, 0.5) - current_co
            return vec.length

        def total_distance_sort(e):
            return e[props.get_dict_layer("edgesortdist", attr_dict)]

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
            this_edge[props.get_dict_layer("edgesortdist", attr_dict)] = total_distance
            has_boundary = False
            found_other_edge = False
            
            edges = []
            boundary_edges = []
            
            for v in this_edge.verts:
                for ve in v.link_edges:
                    if ve not in edges_done and ve != this_edge:
                        edges.append(ve)
                        if ve.is_boundary == True:
                            boundary_edges.append(ve)
                            has_boundary = True
            
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
                edges_remain.sort(key=geo_edge_dist_sort)
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

    def execute(self, context):
        props = context.scene.sloppy_props
        bm = bmesh.from_edit_mesh(bpy.context.active_object.data)
        current_co = mathutils.Vector((0,0,0))

        attr_dict = [
            {"name": "facesortdist", "type": "FLOAT", "domain": "FACE", "layer": None}
            ]
        
        for nam in attr_dict:
            attr = props.find_or_add_attribute(nam["name"], nam["type"], nam["domain"])
            nam["layer"] = props.get_attribute_layer(nam["name"], nam["type"], nam["domain"], bm)

        def face_avg_length_sort(e):
            avg_length = 0.0
            num_lengths = 0
            for edge in e.edges:
                avg_length += edge.calc_length()
                num_lengths += 1
            avg_length /= num_lengths
            return avg_length

        def geo_face_dist_sort(e):
            vec = e.calc_center_median() - current_co
            return vec.length

        def total_distance_sort(e):
            return e[props.get_dict_layer("facesortdist", attr_dict)]

        bm.verts.ensure_lookup_table()

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
            this_face[props.get_dict_layer("facesortdist", attr_dict)] = total_distance
            has_boundary = False
            
            potential_faces = []
            boundary_faces = []
            
            for e in this_face.edges:
                for f in e.link_faces:
                    if f not in faces_done:
                        potential_faces.append(f)
                        for fe in f.edges:
                            if fe.is_boundary == True:
                                if f not in boundary_faces:
                                    boundary_faces.append(f)
                                has_boundary = True
            
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
                faces_remain.sort(key=geo_face_dist_sort)
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