import bpy
import bpy_extras
from bpy_extras import bmesh_utils
import bmesh
import math
import mathutils
import bl_math
import sys

bm = bmesh.from_edit_mesh(bpy.context.active_object.data)

uv_layer = bm.loops.layers.uv.verify()

bmi = bpy_extras.bmesh_utils.bmesh_linked_uv_islands(bm, uv_layer)
print(f"Number of UV islands: {len(bmi)}")
seam_loops = []


try:
    bpy.context.object.data.attributes["SeamDistance"]
except:
    bpy.context.object.data.attributes.new(name="SeamDistance", type="FLOAT", domain="POINT")
    
try:
    bpy.context.object.data.attributes["DistanceAlongSeam"]
except:
    bpy.context.object.data.attributes.new(name="DistanceAlongSeam", type="FLOAT", domain="CORNER")
    
try:
    bpy.context.object.data.attributes["NumPaths"]
except:
    bpy.context.object.data.attributes.new(name="NumPaths", type="INT", domain="POINT") 
   
try:
    bpy.context.object.data.attributes["NumHits"]
except:
    bpy.context.object.data.attributes.new(name="NumHits", type="INT", domain="POINT")
   
try:
    bpy.context.object.data.attributes["IslandIndex"]
except:
    bpy.context.object.data.attributes.new(name="IslandIndex", type="INT", domain="FACE")
   
try:
    bpy.context.object.data.attributes["FaceAlongSeamIndex"]
except:
    bpy.context.object.data.attributes.new(name="FaceAlongSeamIndex", type="INT", domain="FACE")
    
try:
    bpy.context.object.data.attributes["NearSeamTangent"]
except:
    bpy.context.object.data.attributes.new(name="NearSeamTangent", type="FLOAT_VECTOR", domain="CORNER")

try:
    bpy.context.object.data.attributes["Curvature"]
except:
    bpy.context.object.data.attributes.new(name="Curvature", type="FLOAT", domain="POINT")
    
try:
    bpy.context.object.data.attributes["Thickness"]
except:
    bpy.context.object.data.attributes.new(name="Thickness", type="FLOAT", domain="POINT")
    
try:
    bpy.context.object.data.attributes["FaceThickness"]
except:
    bpy.context.object.data.attributes.new(name="FaceThickness", type="FLOAT", domain="POINT")

seam_tang = bm.verts.layers.float_vector.get("NearSeamTangent")

seam_dist = bm.verts.layers.float.get("SeamDistance")
seam_dist_along = bm.loops.layers.float.get("DistanceAlongSeam")
near_seam_tangent = bm.loops.layers.float_vector.get("NearSeamTangent")
curvature = bm.verts.layers.float.get("Curvature")
num_paths = bm.verts.layers.int.get("NumPaths")
thickness = bm.verts.layers.float.get("Thickness")
face_thickness = bm.verts.layers.float.get("FaceThickness")
num_hits = bm.verts.layers.int.get("NumHits")
island_index = bm.faces.layers.int.get("IslandIndex")
face_index = bm.faces.layers.int.get("FaceAlongSeamIndex")

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

#for als, aals in zip(along_loops, all_along_loops):
#    sorted = []
#    sorted_faces = []
#    sorted_verts = []
#    sorted_verts_is = []
#    sorted_verts_lengths = []
#    sorted_verts_tangents = []
#    next_loop = als[0]
#    for loop in als:
#        if loop.edge.seam == True:
#            next_loop = loop
#            break
#    max_dist_along = 0.0
#    dist_along = 0.0
#    alrs_index = 0
#    iter_this = 0
#    alrs = along_loops_remaining[iter]
#    aalrs = all_along_loops_remaining[iter]
#    
#    iter_while = len(alrs)
#    chain = 0
    
#    while len(alrs) > 0:
#        this_loop = next_loop
#        this_edge = this_loop.edge
#        this_face = this_loop.face
#        this_island_index = this_face[island_index]
#        aes = along_edges[this_island_index]
#        this_vert = this_loop.vert
#        next_vert = this_edge.other_vert(this_vert)
#        this_distance = dist_along
#        this_tangent = this_edge.calc_tangent(this_loop)
#        dist_along += this_edge.calc_length()
#        print(f"Current distance: {dist_along}")
#        new_max_along = max(max_dist_along, dist_along)
#        max_dist_along = new_max_along
#        if this_edge in aes:
#            print("Edge is seam!")
#            chain += 1
#            sorted_faces.append(this_face)
#            if this_vert not in sorted_verts:
#                sorted_verts.append(this_vert)
#                temp_arr = []
#                temp_arr.append(this_distance)
#                temp_tan_arr = []
#                temp_tan_arr.append(this_tangent)
#                sorted_verts_is.append(iter_this)
#                sorted_verts_lengths.append(temp_arr)
#                sorted_verts_tangents.append(temp_tan_arr)
#                iter_this += 1
#            if this_vert in sorted_verts:
#                get_i = sorted_verts.index(this_vert)
#                sorted_verts_lengths[get_i].append(this_distance)
#                sorted_verts_tangents[get_i].append(this_tangent)
#                
#            if this_vert not in sorted_verts:
#                sorted_verts.append(this_vert)
#                temp_arr = []
#                temp_arr.append(dist_along)
#                temp_tan_arr = []
#                temp_tan_arr.append(this_tangent)
#                sorted_verts_is.append(iter_this)
#                sorted_verts_lengths.append(temp_arr)
#                sorted_verts_tangents.append(temp_tan_arr)
#                iter_this += 1
#            if next_vert in sorted_verts:
#                get_i = sorted_verts.index(next_vert)
#                sorted_verts_lengths[get_i].append(dist_along)
#                sorted_verts_tangents[get_i].append(this_tangent)
#                
#            this_loop[seam_dist_along] = this_distance
#            this_loop[near_seam_tangent] = this_tangent
#            if this_loop not in sorted:
#                sorted.append(this_loop)
#            if this_loop in alrs:
#                alrs.remove(this_loop)
#                aalrs.remove(this_loop)
#            
#            for loop in next_vert.link_loops:
#                if loop in this_face.loops:
#                    loop[seam_dist_along] = dist_along
#                    loop[near_seam_tangent] = this_tangent
#                    if loop not in sorted:
#                        sorted.append(loop)
#                    if loop in alrs:
#                        alrs.remove(loop)
#                        aalrs.remove(loop)
#                else:
#                    if loop in als:
#                        if loop.edge in aes:
#                            if loop.edge != this_edge:
#                                if loop not in sorted:
#                                    next_loop = loop
#                                    break
#                                if loop in sorted:
#                                    alrs_index = 0
#                                    next_loop = alrs[alrs_index]
#                                    break
#                                
#        
#        if this_edge not in aes:
#            print("Edge is not seam!")
#            chain = 0
#            dist_along = 0.0
#            if this_vert in sorted_verts:
#                print("Vert sorted!")
#                before = False
#                first_face = sorted_faces[0]
#                last_face = sorted_faces[-1]
#                loop_to_copy = this_loop
#                if this_loop == sorted[0]:
#                    for edge in this_face.edges:
#                        if last_face in edge.link_faces:
#                            before = True
#                    get_i = sorted_verts.index(this_vert)
#                    if before == True:
#                        for loop in this_vert.link_loops:
#                            if loop in als:
#                                if loop in sorted:
#                                    if loop.edge in aes:
#                                        if loop.edge in last_face.edges:
#                                            if loop != this_loop:
#                                                loop_to_copy = loop
#                        print("Taking previous loops length.")
#                        this_loop[seam_dist_along] = loop_to_copy[seam_dist_along]
#                        this_loop[near_seam_tangent] = loop_to_copy[near_seam_tangent]
#                    if before == False:
#                        for loop in this_vert.link_loops:
#                            if loop in als:
#                                if loop in sorted:
#                                    if loop.edge in aes:
#                                        if loop.edge in first_face.edges:
#                                            if loop != this_loop:
#                                                loop_to_copy = loop
#                        print("Taking next loops length.")
#                        this_loop[seam_dist_along] = loop_to_copy[seam_dist_along]
#                        this_loop[near_seam_tangent] = loop_to_copy[near_seam_tangent]
#                else:
#                    for loop in this_vert.link_loops:
#                        if loop in als:
#                            if loop in sorted:
#                                if loop.edge in aes:
#                                    if loop.edge in this_face.edges:
#                                        if loop != this_loop:
#                                            loop_to_copy = loop
#                    this_loop[seam_dist_along] = loop_to_copy[seam_dist_along]
#                    this_loop[near_seam_tangent] = loop_to_copy[near_seam_tangent]
#                sorted.append(this_loop)
#                print(loop[seam_dist_along])
#                if this_loop in alrs:
#                    alrs.remove(this_loop)
#                    aalrs.remove(this_loop)
#                if this_loop in sorted:
#                    alrs_index = 0
#                    next_loop = alrs[alrs_index]
#                    break
#            if this_vert not in sorted_verts:
#                print("Vert not sorted!")
#                alrs_index += 1
#                next_loop = alrs[alrs_index]
#        bpy.context.active_object.data.update() 
#        print(f"Island: {iter + 1} of {len(bmi)} Loops remaining: {len(alrs)} of {len(als)} Next loop: {next_loop.index} ALRS Index: {alrs_index} Current chain links: {chain}")
    
#    distances_along_seams.append(max_dist_along)
#    iter += 1

#max_distance_along_seam = max(distances_along_seams)

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
#    dirs.append(mathutils.Vector((0,nordir.y,nordir.z)).normalized())
#    dirs.append(mathutils.Vector((nordir.x,0,nordir.z)).normalized())
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
#                print(f"Hit: {hit}, at {loc}. Thickness: {thick}. Dot: {dirdot}. Hit face {index}.")
#            if index == f.index:
#                print("Hit face is this face!")
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
#    newthick = v[thickness] / (len(v.link_faces) + (len(v.link_faces) * v[num_hits]))
#    newthick = v[thickness] / len(v.link_faces)
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
#                print(f"Hit: {hit}, at {loc}. Thickness: {newthick}. Hit face {index}.")
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