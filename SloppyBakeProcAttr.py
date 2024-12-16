import bpy
import bpy_extras
import bmesh
import math
import mathutils

# region ProcAttrBke Op
class SloppyProcAttrBake(bpy.types.Operator):
    bl_idname = "operator.sloppy_proc_attr_bake"
    bl_label = "Bake Procedural Attributes"
    bl_description = "Bake procedural values (distance from seam, distance along seam, island index, curvature and thickness) to mesh element attributes"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        props = context.scene.sloppy_props
        bm = bmesh.from_edit_mesh(bpy.context.active_object.data)

        attr_dict = [
            {"name" : "SeamDistance", "type":"FLOAT", "domain":"POINT", "layer": None},
            {"name" : "DistanceAlongSeam", "type":"FLOAT", "domain":"CORNER", "layer": None},
            {"name" : "NumPaths", "type":"INT", "domain":"POINT", "layer": None},
            {"name" : "NumHits", "type":"INT", "domain":"POINT", "layer": None},
            {"name" : "IslandIndex", "type":"INT", "domain":"FACE", "layer": None},
            {"name" : "FaceAlongSeamIndex", "type":"INT", "domain":"FACE", "layer": None},
            {"name" : "NearSeamTangent", "type":"FLOAT_VECTOR", "domain":"CORNER", "layer": None},
            {"name" : "Curvature", "type":"FLOAT", "domain":"POINT", "layer": None},
            {"name" : "Thickness", "type":"FLOAT", "domain":"POINT", "layer": None},
            {"name" : "FaceThickness", "type":"FLOAT", "domain":"POINT", "layer": None},
            ]
        
        for nam in attr_dict:
            attr = props.find_or_add_attribute(nam["name"], nam["type"], nam["domain"])
            nam["layer"] = props.get_attribute_layer(nam["name"], nam["type"], nam["domain"], bm)
            
        uv_layer = bm.loops.layers.uv.verify()

        seam_dist = props.get_dict_layer("SeamDistance", attr_dict)
        seam_dist_along = props.get_dict_layer("DistanceAlongSeam", attr_dict)
        near_seam_tangent = props.get_dict_layer("NearSeamTangent", attr_dict)
        curvature = props.get_dict_layer("Curvature", attr_dict)
        num_paths = props.get_dict_layer("NumPaths", attr_dict)
        thickness = props.get_dict_layer("Thickness", attr_dict)
        face_thickness = props.get_dict_layer("FaceThickness", attr_dict)
        num_hits = props.get_dict_layer("NumHits", attr_dict)
        island_index = props.get_dict_layer("IslandIndex", attr_dict)
        face_index = props.get_dict_layer("FaceAlongSeamIndex", attr_dict)

        bmi = bpy_extras.bmesh_utils.bmesh_linked_uv_islands(bm, uv_layer)
        if props.verbose:
            print(f"Number of UV islands: {len(bmi)}")
        seam_loops = []

        depsgraph = bpy.context.evaluated_depsgraph_get()

        done_loops = []
        done_verts = []
        seed_verts = []

        max_curv = 0.0
        min_curv = 360.0

        max_thick = 0
        min_thick = 100
        max_fthick = 0
        min_fthick = 100

        # DISTANCE ALONG SEAM (+ 

        along_edges = []
        along_faces = []
        along_faces_remaining = []
        along_loops = []
        along_verts = []
        all_along_loops = []
        all_along_loops_remaining = []
        along_loops_sorted = []
        along_loops_remaining = []
        distances_along_seams = []

        iter = 0

        for island in bmi:
            island_edges = []
            island_faces = []
            island_verts = []
            island_faces_remaining = []
            island_loops = []
            all_island_loops = []
            all_island_loops_remaining = []
            island_loops_remaining = []
            for f in island:
                f[island_index] = iter
                has_seam = 0
                for v in f.verts:
                    for edge in v.link_edges:
                        if edge.seam == True:
                            has_seam = 1
                            if edge not in island_edges:
                                island_edges.append(edge)
                            for ev in edge.verts:
                                if ev not in island_verts:
                                    island_verts.append(ev)
                for loop in f.loops:
                    if loop not in all_island_loops:
                        all_island_loops.append(loop)
                        all_island_loops_remaining.append(loop)
                    if has_seam == 1:
                        if loop not in island_loops:
                            island_loops.append(loop)
                            island_loops_remaining.append(loop)
                if has_seam == 1:
                    if f not in island_faces:
                        island_faces.append(f)
                        island_faces_remaining.append(f)
            along_edges.append(island_edges)
            along_loops.append(island_loops)
            along_verts.append(island_verts)
            along_faces_remaining.append(island_faces_remaining)
            along_faces.append(island_faces)
            all_along_loops.append(all_island_loops)
            all_along_loops_remaining.append(all_island_loops_remaining)
            along_loops_remaining.append(island_loops_remaining)
            iter += 1

        latest_face_pos = mathutils.Vector((-100,-100,-100))

        def sort_remaining(e):
            vec = e.calc_center_median() - latest_face_pos
            return vec.length

        max_face_indices = []
        iter = 0
        for aes, afs, avs, afrs, als, alrs, ifs, avs in zip(along_edges, along_faces, along_verts, along_faces_remaining, along_loops, along_loops_remaining, bmi, along_verts):
            latest_face_pos = bm.verts[0].co
            afs.sort(key=sort_remaining)
            ifrs = ifs.copy()
            sorted_loops = []
            sorted_faces = []
            sorted_edges = []
            sorted_faces_tangent = []
            sorted_faces_start_distance = []
            sorted_faces_end_distance = []
            sorted_faces_next = []
            next_face = afs[0]
            this_tangent = mathutils.Vector((1,0,0))
            for f in afs:
                is_seamed = False
                for edge in f.edges:
                    if edge.seam == True:
                        is_seamed = True
                        for loop in edge.link_loops:
                            if loop in f.loops:
                                if loop.edge == edge:
                                    this_tangent = edge.calc_tangent(loop)
                if is_seamed == True:
                    next_face = f
            max_face_index = 0
            seam_edge = aes[0]
            dist_along = 0
            print(len(afrs))
            
            face_iter = 0
            while len(afrs) > 0:
                this_face = next_face
                start_distance = dist_along
                    
                
                this_face[face_index] = face_iter
                max_face_index = face_iter
                for edge in this_face.edges:
                    if edge.seam == True:
                        dist_along += edge.calc_length()
                        for loop in edge.link_loops:
                            if loop in f.loops:
                                this_tangent = edge.calc_tangent(loop)
                if this_face not in sorted_faces:
                    sorted_faces.append(this_face)
                    sorted_faces_start_distance.append(start_distance)
                    sorted_faces_end_distance.append(dist_along)
                    sorted_faces_tangent.append(this_tangent)
                    for edge in this_face.edges:
                        if edge not in sorted_edges:
                            sorted_edges.append(edge)
                if this_face in afrs:
                    afrs.remove(this_face)
                adjacent_faces = []
                
                
                for edge in this_face.edges:
                    if edge.seam == False:
                        for face in edge.link_faces:
                            if face in afs:
                                if face not in sorted_faces:
                                    adjacent_faces.append(face)
                if len(adjacent_faces) == 0:
                    for v in this_face.verts:
                        for face in v.link_faces:
                            if face in afs:
                                if face not in sorted_faces:
                                    adjacent_faces.append(face)
                face_iter += 1
                if len(adjacent_faces) > 0:
                    face_seams = []
                    for face in adjacent_faces:
                        seams = 0
                        for edge in face.edges:
                            if edge.seam == True:
                                seams += 1
                        face_seams.append(seams)
                    max_seams = max(face_seams)
                    get_i = face_seams.index(max_seams)
                    next_face = adjacent_faces[get_i]
                    sorted_faces_next.append(next_face)
                else:
        #            face_iter = 0
                    if len(afrs) > 0:
                        if len(afrs) > 0:
                            latest_face_pos = this_face.calc_center_median()
                            afrs.sort(key=sort_remaining)
                            next_face = afrs[0]
                            sorted_faces_next.append(next_face)
                        else:
                            sorted_faces_next.append(sorted_faces[0])
                            break
            
            for face, start, end, tan, next in zip(sorted_faces, sorted_faces_start_distance, sorted_faces_end_distance, sorted_faces_tangent, sorted_faces_next):
                for loop in face.loops:
                    if loop in als:
                        loop[near_seam_tangent] = tan
                        if loop.vert in next.verts:
                            loop[seam_dist_along] = end
                        else:
                            loop[seam_dist_along] = start
                        if loop not in sorted_loops:
                            sorted_loops.append(loop)
            
            bpy.context.active_object.data.update() 
            
            next_faces = sorted_faces.copy()
            referring_faces = sorted_faces.copy()
            
            face_iter = 0
            
            while len(ifrs) > 0:
                print(f"Island {iter + 1} of {len(bmi)}. Faces remaining: {len(ifrs)} Processing {len(next_faces)} faces. {len(sorted_faces)} faces sorted. Retry iteration: {face_iter}")
                these_faces = next_faces.copy()
                these_referring_faces = referring_faces.copy()
                next_faces = []
                referring_faces = []
                
                for f, rf in zip(these_faces, these_referring_faces):
                    if f not in sorted_faces:
                        f[face_index] = rf[face_index]
                        off_edges = []
                        
                        for edge in f.edges:
                            off = False
                            for v in edge.verts:
                                if v in rf.verts:
                                    if edge not in rf.edges:
                                        off = True
                            if off == True:
                                if edge not in off_edges:
                                    off_edges.append(edge)
                        
                        for edge in off_edges:
                            dist_to_use = 0
                            tan = edge.link_loops[0][near_seam_tangent]
                            for v in edge.verts:
                                if v in rf.verts:
                                    for loop in v.link_loops:
                                        if loop.face == rf:
                                            dist_to_use = loop[seam_dist_along]
                                            tan = loop[near_seam_tangent]
                            for v in edge.verts:
                                if v in f.verts:
                                    for loop in v.link_loops:
                                        if loop.face == f:
                                            if loop not in sorted_loops:
                                                loop[seam_dist_along] = dist_to_use
                                                loop[near_seam_tangent] = tan
                                                sorted_loops.append(loop)
                            bpy.context.active_object.data.update() 
                            print(dist_to_use)
                        for loop in f.loops:
                            if loop not in sorted_loops:
                                sorted_loops.append(loop)
                        sorted_faces.append(f)
                        for edge in f.edges:
                            if edge not in sorted_edges:
                                sorted_edges.append(edge)
                        if f in ifrs:
                            ifrs.remove(f)
                                        
                    for edge in f.edges:
                        for face in edge.link_faces:
                            if face in ifs:
                                if face not in sorted_faces:
                                    if face not in next_faces:
                                        next_faces.append(face)
                                        referring_faces.append(f)
                    if f in sorted_faces:
                        if f in ifrs:
                            ifrs.remove(f)
                
                if len(next_faces) == 0:
                    if len(ifrs) > 0:
                        face_iter += 1
                        next_faces.append(ifrs[face_iter])
            
            for f in ifs:
                for v in f.verts:
                    if v not in avs:
                        has_zero = False
                        num_vloops = 0
                        sum_length = 0
                        for loop in v.link_loops:
                            num_vloops += 1
                            sum_length += loop[seam_dist_along]
                            if loop[seam_dist_along] == 0.0:
                                has_zero = True
                        if has_zero == False:
                            sum_length /= num_vloops
                            for loop in v.link_loops:
                                loop[seam_dist_along] = sum_length
            
            bpy.context.active_object.data.update() 
            
            distances_along_seams.append(dist_along)
                
            
            max_face_indices.append(max_face_index)
            iter += 1

        max_face_index_sum = max(max_face_indices)
        max_distance_along_seam = max(distances_along_seams)

        # FACE THICKNESS

        bpy.context.active_object.data.update() 

        for v in bm.verts:
            v[thickness] = 0.0
            v[face_thickness] = 0.0
            v[num_hits] = 0
            
        bpy.context.active_object.data.update() 

        for f in bm.faces:
            thick = 0
            dirdot = 0
            hits = 0
            nordir = -f.normal
            dirs = [nordir]
            dirs.append(mathutils.Vector((0,0,nordir.z)).normalized())
            dirs.append(mathutils.Vector((nordir.x,nordir.y,0)).normalized())
            fco = f.calc_center_median() + (nordir * 0.0001)
            for d in dirs:
                hit, loc, norm, index, obj, ab = bpy.context.scene.ray_cast(depsgraph, fco, d)
                if hit:
                    if index != f.index:
                        vec = loc - fco
                        thick += vec.length
                        hits += 1
                        dir = vec.normalized()
                        nordir = f.normal.normalized()
                        dirdot = dir.dot(nordir)
                        if props.verbose:
                            print(f"Hit: {hit}, at {loc}. Thickness: {thick}. Dot: {dirdot}. Hit face {index}.")
            if hits > 0:
                thick /= hits
            for v in f.verts:
                v[face_thickness] += thick

        bpy.context.active_object.data.update()

        # VERT THICKNESS + CURVATURE

        for v in bm.verts:
            curv = 0
            for edge in v.link_edges:
                curv += math.degrees(edge.calc_face_angle_signed(0))
            newthick = 0
            d = -v.normal
            dirs = [d]
            hits = 0
            dirs.append(mathutils.Vector((0,0,d.z)).normalized())
            dirs.append(mathutils.Vector((d.x,d.y,0)).normalized())
            fco = v.co + (d * 0.0001)
            for dl in dirs:
                hit, loc, norm, index, obj, ab = bpy.context.scene.ray_cast(depsgraph, fco, dl)
                if hit:
                    if bm.faces[index] not in v.link_faces:
                        vec = loc - fco
                        newthick += vec.length
                        hits += 1
                        if props.verbose:
                            print(f"Hit: {hit}, at {loc}. Thickness: {newthick}. Hit face {index}.")
            if hits > 0:
                newthick /= hits
            newthick /= v.calc_shell_factor()
            newmax_thick = max(max_thick, newthick)
            newmin_thick = min(min_thick, newthick)
            max_thick = newmax_thick
            min_thick = newmin_thick
            v[thickness] = newthick
            v[face_thickness] /= len(v.link_faces)
            newmax_fthick = max(max_fthick, v[face_thickness])
            newmin_fthick = min(min_fthick, v[face_thickness])
            max_fthick = newmax_fthick
            min_fthick = newmin_fthick
            curv /= len(v.link_edges)
            v[curvature] = curv
            newmax_curv = max(max_curv, curv)
            newmin_curv = min(min_curv, curv)
            max_curv = newmax_curv
            min_curv = newmin_curv

        bpy.context.active_object.data.update() 

        # DISTANCE FROM SEAM

        for v in bm.verts:
            v[seam_dist] = 0.0
            v[num_paths] = 0

        for f in bm.faces:
            for loop in f.loops:
                uva = loop[uv_layer]
                for edge in loop.vert.link_edges:
                    if edge.seam == True:
                        for v in edge.verts:
                            v[seam_dist] = 0.0
                            v[num_paths] = 0
                            if v not in done_verts:
                                done_verts.append(v)
                            if v not in seed_verts:
                                seed_verts.append(v)
                        done_loops.append(loop)
                        seam_loops.append(loop)


                        
        bpy.context.active_object.data.update()

        remain = len(bm.verts) - len(done_verts)

        max_dist = 0
        max_paths = 0

        while remain > 0:
            print(f"{remain} verts not yet done.")
            
            next_verts = []

            for v in seed_verts:
                if v not in done_verts:
                    if v[num_paths] > 0:
                        v[seam_dist] /= v[num_paths]
                        newmax = max(max_dist, v[seam_dist])
                        newmaxpaths = max(max_paths, v[num_paths])
                        max_dist = newmax
                        max_paths = newmaxpaths
                    done_verts.append    
                for edge in v.link_edges:
                    vo = edge.other_vert(v)
                    if vo not in done_verts and vo not in seed_verts:
                        last_dist = v[seam_dist]
                        dist = last_dist + edge.calc_length()
                        vo[seam_dist] += dist
                        vo[num_paths] += 1
                        if vo not in next_verts:
                            next_verts.append(vo)
            
            for v in seed_verts:
                if v not in done_verts:
                    done_verts.append(v)
                
            seed_verts = []
            for v in next_verts:
                if v not in seed_verts:
                    seed_verts.append(v)
            
            remain = len(bm.verts) - len(done_verts)
            
            bpy.context.active_object.data.update()

        print(f"Maximum distance from seam: {max_dist}. Maximum number of paths: {max_paths}")
        bpy.context.object["MaxDistanceFromSeam"] = max_dist
        bpy.context.object["MaxDistanceAlongSeam"] = max_distance_along_seam
        bpy.context.object["MaxCurvature"] = max_curv
        bpy.context.object["MinCurvature"] = min_curv
        bpy.context.object["MaxThickness"] = max_thick
        bpy.context.object["MinThickness"] = min_thick
        bpy.context.object["MaxFaceThickness"] = max_fthick
        bpy.context.object["MinFaceThickness"] = min_fthick
        bpy.context.object["NumberOfUVIslands"] = len(bmi)
        bpy.context.object["MaxAlongFaceIndex"] = max_face_index_sum

        bpy.context.active_object.data.update()

        return {"FINISHED"}
# endregion

# region Blur Attribute Op
class SloppyBlurAttribute(bpy.types.Operator):
    bl_idname = "operator.sloppy_attr_blur"
    bl_label = "Blur Attribute"
    bl_description = "Blur active attribute values"
    bl_options = {'REGISTER', 'UNDO'}

    iP = bpy.props.IntProperty
    fP = bpy.props.FloatProperty
    fvP = bpy.props.FloatVectorProperty
    bP = bpy.props.BoolProperty
    eP = bpy.props.EnumProperty
    bvP = bpy.props.BoolVectorProperty
    sP = bpy.props.StringProperty

    # region AttrBlur Properties
    respect_islands : bP(
        name = "Respect UV Islands",
        description = "If true, only blur with elemeents within same UV island",
        default = False
        ) # type: ignore
    # endregion
    
    def execute(self, context):
        props = context.scene.sloppy_props
        bm = bmesh.from_edit_mesh(bpy.context.active_object.data)

        active_attribute = bpy.context.active_object.data.attributes.active
        new_attr_name = active_attribute.name + '_blur'
        if '_blur' in active_attribute.name:
            prev_blur_iter = active_attribute.name[-1]
            if prev_blur_iter.isdigit() == True:
                new_blur_iter = int(prev_blur_iter) + 1
                new_attr_name = active_attribute.name + str(new_blur_iter)
            else:
                new_attr_name = active_attribute.name + '1'
        new_attr_domain = active_attribute.domain
        new_attr_type = active_attribute.data_type

        if new_attr_domain == 'CORNER':
            return {"CANCELLED"}

        attr_dict = [
            {"name" : active_attribute.name, "type":active_attribute.data_type, "domain":active_attribute.domain, "layer": None},
            {"name" : new_attr_name, "type":new_attr_domain, "domain":new_attr_type, "layer": None},
            ]
        
        for nam in attr_dict:
            attr = props.find_or_add_attribute(nam["name"], nam["type"], nam["domain"])
            nam["layer"] = props.get_attribute_layer(nam["name"], nam["type"], nam["domain"], bm)
            
        uv_layer = bm.loops.layers.uv.verify()

        clean_attr = props.get_dict_layer(active_attribute.name, attr_dict)
        blur_attr = props.get_dict_layer(new_attr_name, attr_dict)

        bmi = []
        if self.respect_islands == True:
            bmi = bpy_extras.bmesh_utils.bmesh_linked_uv_islands(bm, uv_layer)
        else:
            bmi.append([i for i in bm.faces])
        if props.verbose:
            print(f"Number of islands: {len(bmi)}")

        for island in bmi:
            if new_attr_domain == 'FACE':
                for face in island:
                    done_i = [face.index]
                    attr_base = face[clean_attr]
                    for fv in face.verts:
                        for fvf in fv.link_faces:
                            if fvf.index not in done_i and fvf in island:
                                done_i.append(fvf.index)
                                attr_base += fvf[clean_attr]
                    new_attr_calc = attr_base / len(done_i)
                    face[blur_attr] = new_attr_calc
            if new_attr_domain == 'EDGE' or new_attr_domain == 'POINT':
                island_edges = []
                for face in island:
                    for fe in face.edges:
                        if fe not in island_edges:
                            island_edges.append(fe)
                if new_attr_domain == 'EDGE':
                    for edge in island_edges:
                        done_i = [edge.index]
                        attr_base = edge[clean_attr]
                        for ev in edge.verts:
                            for eve in ev.link_edges:
                                if eve.index not in done_i and eve in island_edges:
                                    done_i.append(eve.index)
                                    attr_base += eve[clean_attr]
                        new_attr_calc = attr_base / len(done_i)
                        edge[blur_attr] = new_attr_calc
                if new_attr_domain == 'POINT':
                    island_verts = []
                    for edge in island_edges:
                        for ev in edge.verts:
                            if ev not in island_verts:
                                island_verts.append(ev)
                    for vert in island_verts:
                        done_i = [vert.index]
                        attr_base = vert[clean_attr]
                        for ve in vert.link_edges:
                            if ve in island_edges:
                                veov = ve.other_vert(vert)
                                if veov.index not in done_i:
                                    done_i.append(veov.index)
                                    attr_base += veov[clean_attr]
                        new_attr_calc = attr_base / len(done_i)
                        vert[blur_attr] = new_attr_calc

        bpy.context.active_object.data.update()

        return {"FINISHED"}
# endregion