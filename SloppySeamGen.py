import bpy
import bpy_extras
from bpy_extras import bmesh_utils
import bmesh
import math
import mathutils
import bl_math
import sys

class SloppySeamGen(bpy.types.Operator):
    bl_idname = "operator.sloppy_seam_gen"
    bl_label = "Generate Seans"
    bl_description = "Generate seams according to edge concavity, among other things"
    
    def execute(self, context):
        props = context.scene.sloppy_props
        bm = bmesh.from_edit_mesh(bpy.context.active_object.data)
        bpy.ops.mesh.select_linked(delimit={'SEAM'})
        uv_layer = bm.loops.layers.uv.verify()
        
        # circle to square function, maps circular coordinates to square coordinates
        def circle_to_square(u, v):
            x = 0
            y = 0
            u2 = u * u
            v2 = v * v
            twosqrt2 = 2.0 * math.sqrt(2.0)
            subtermx = 2.0 + u2 - v2
            subtermy = 2.0 - u2 + v2
            termx1 = subtermx + u * twosqrt2
            termx2 = subtermx - u * twosqrt2
            termy1 = subtermy + v * twosqrt2
            termy2 = subtermy - v * twosqrt2
            if termx1 < 0 or termx2 < 0:
                x = 0
            else:
                x = (0.5 * math.sqrt(termx1)) - (0.5 * math.sqrt(termx2))
            if termy1 < 0 or termy2 < 0:
                y = 0
            else:
                y = (0.5 * math.sqrt(termy1)) - (0.5 * math.sqrt(termy2))
            result = mathutils.Vector((x,y))
            return result

        # initialize console commands
        CURSOR_UP = '\033[F'
        ERASE_LINE = '\033[K'
        
        # initialize variables
        island_center = mathutils.Vector((0,0))
        directions = []
        initial_boundary_co_x = []
        initial_boundary_co_y = []
        boundary_co = []
        distances = []
        real_vert_co = []
        boundary_loops = []
        inner_loops = []
        linked_verts = []
        unique_verts = []
        linked_loops = []
        edge_threshold = 0.001

        # calculate island_center and append real vert coords and linked verts for each loop
        iter = 0
        seam_iter = 0
        if props.verbose == True:
            print("\nPelt UV Generation started\n")
            print("Step 1: Calculate UV island center and initialize island boundary")
            print("=================================================================")
        for v in bm.faces:
            iter += 1
            if v.select == True:
                for loop in v.loops:
                    uva = loop[uv_layer]
                    has_seam = 0
                    for edge in loop.vert.link_edges:
                        if edge.seam == True:
                            has_seam = 1
                    seam_iter += has_seam
                    if has_seam == 0:
                        if loop not in inner_loops:
                            inner_loops.append(loop)
                    if has_seam == 1:
#                    if uva.select == True:
                        uvco = uva.uv
                        initial_boundary_co_x.append(uvco.x)
                        initial_boundary_co_y.append(uvco.y)
                        
                        if v.index not in unique_verts:
                            unique_verts.append(v.index)
                            linked_loops.append([])
                        
                        if loop not in boundary_loops:
                            boundary_loops.append(loop)
                            real_vert_co.append(loop.vert.co)
                            linked_verts.append(v.index)

                        get_i = unique_verts.index(v.index)
                        linked_loops[get_i].append(loop)
                        if props.verbose == True:
                            done_pct = (iter / len(bm.verts)) * 100
                            if len(boundary_loops) > 1:
                                print(CURSOR_UP + CURSOR_UP + CURSOR_UP + CURSOR_UP + CURSOR_UP)
                            msg = "Added loop {:} to boundary list. \nCurrent number of boundary loops: {:}. \nVertex {:} of {:} checked. \nStep 1 {:.1f}% done.".format(loop.index, len(boundary_loops), iter, len(bm.verts), done_pct)
                            print(msg, flush=True)
                            sys.stdout.flush()
        
        # cancel if no geometry selected
        if len(boundary_loops) < 1:
            props.sloppy_error_msg_heading = "Can't generate pelt!"
            props.sloppy_error_msg = "    No geometry selected."
            bpy.ops.operator.sloppy_dialog('INVOKE_DEFAULT')
            if props.verbose:
                print("Pelt Generation cancelled...\n")
            return {'CANCELLED'}
        
        max_x = max(initial_boundary_co_x)
        min_x = min(initial_boundary_co_x)
        max_y = max(initial_boundary_co_y)
        min_y = min(initial_boundary_co_y)

        island_center.x = bl_math.lerp(min_x, max_x, 0.5)
        island_center.y = bl_math.lerp(min_y, max_y, 0.5)
        
        if props.verbose == True:
            print(CURSOR_UP + ERASE_LINE + CURSOR_UP)
            print("Step 1 complete!\n")
            print("Step 1 results:")
            print("---------------")
            print(f"Center of UV island: {island_center.x}, {island_center.y}")
            print(f"Number of selected verts: {len(unique_verts)}")
            print(f"{seam_iter} of boundary loops were linked to seam edges through their vertex.")
            print(f"Number of boundary UV loops: {len(boundary_loops)}")
            print(f"Number of non-boundary UV loops: {len(inner_loops)}\n")

        # initialize corners and corner distance arrays
        corners = [mathutils.Vector((0,1)), mathutils.Vector((1,1)), mathutils.Vector((0,0)), mathutils.Vector((1,0))]
        corner_dists = [[], [], [], []]
        
        if props.verbose == True:
            print("Step 2: Square boundary and calculate distance from center for each boundary point")
            print("==================================================================================")

        # calculate distances and directions from center for boundary loop points,
        # move points to 1 UV unit radius circle and then square that circle
        iter = 0
        for loop in boundary_loops:
            iter += 1
            uva = loop[uv_layer]
            uvco = uva.uv
            uvco_old = uvco
            vec = uvco - island_center
            dir = vec.normalized()
            dist = vec.length
            directions.append(dir)
            distances.append(dist)
            vec_sqr = circle_to_square(dir.x, dir.y)
            uvco.x = 0.5 + (vec_sqr.x * 0.5)
            uvco.y = 0.5 + (vec_sqr.y * 0.5)
            
            for c,d in zip(corners, corner_dists):
                vec_cx = c.x - uvco.x
                vec_cy = c.y - uvco.y
                dist_c = mathutils.Vector((vec_cx,vec_cy))
                d.append(dist_c)
            
            boundary_co.append(uvco)
            
            if props.verbose == True:
                done_pct = (iter / len(boundary_loops)) * 100
                if iter > 1:
                    print(CURSOR_UP + CURSOR_UP + CURSOR_UP + CURSOR_UP + ERASE_LINE + CURSOR_UP)
                msg = "Moved point at {:.5f}, {:.5f} to new position at {:.5f}, {:.5f}. \nInitial distance from island center: {:.5f}. \nBoundary point direction: {:.5f}, {:.5f} \nStep 2 {:.1f}% done.".format(uvco_old.x, uvco_old.y, uvco.x, uvco.y, dist, dir.x, dir.y, done_pct)
                print(msg, flush=True)
                sys.stdout.flush()
            bpy.context.active_object.data.update()
            
        if props.verbose == True:
            print(CURSOR_UP + ERASE_LINE + CURSOR_UP)
            print("Step 2 complete!\n")
            print("Step 2 results:")
            print("---------------")

        max_distance = max(distances)
        get_max_i = distances.index(max_distance)
        max_dist_dir = directions[get_max_i]
        min_distance = min(distances)
        get_min_i = distances.index(min_distance)
        min_dist_dir = directions[get_min_i]
        if props.verbose == True:
            print("Max. distance from center is {:.5f} in direction: {:.5f}, {:.5f}".format(max_distance, max_dist_dir.x, max_dist_dir.y))
            print("Min. distance from center is {:.5f} in direction: {:.5f}, {:.5f}\n".format(min_distance, min_dist_dir.x, min_dist_dir.y))
            print("Step 3: Align boundary points along UV edges")
            print("============================================")
        
        # move points nearest to corners to corners
        for c,d in zip(corners, corner_dists):
            near_c = min(d)
            cis = [i for i,x in enumerate(d) if x == near_c]
            for ci in cis:
                boundary_loops[ci][uv_layer].uv.x = c.x
                boundary_loops[ci][uv_layer].uv.y = c.y

        bpy.context.active_object.data.update()

        # align boundary points within edge_threshold distance to edges to edges
        ud = 1.0 - edge_threshold
        ld = 0.0 + edge_threshold

        iter = 0
        for loop in boundary_loops:
            if loop[uv_layer].uv.y >= ud:
                iter += 1
                loop[uv_layer].uv.y = 1.0
            if loop[uv_layer].uv.y <= ld:
                iter += 1
                loop[uv_layer].uv.y = 0.0
            if loop[uv_layer].uv.x >= ud:
                iter += 1
                loop[uv_layer].uv.x = 1.0
            if loop[uv_layer].uv.x <= ld:
                iter += 1
                loop[uv_layer].uv.x = 0.0

        bpy.context.active_object.data.update()
        
        if props.verbose == True:
            print(f"{iter} boundary points within edge threshold {edge_threshold} of edge aligned.")
        
        if props.space_edges == True:
            # space boundary points equally along their edge
            def sortloopx(e):
                return e[uv_layer].uv.x

            def sortloopy(e):
                return e[uv_layer].uv.y

            boundary_sorted_x = []
            boundary_sorted_y = []
            boundary_unique_lx = []
            boundary_unique_rx = []
            boundary_unique_uy = []
            boundary_unique_ly = []

            sorting_boundary_loops_x = boundary_loops
            sorting_boundary_loops_x.sort(key=sortloopx)

            sorting_boundary_loops_y = boundary_loops
            sorting_boundary_loops_y.sort(key=sortloopy)

            for loop in sorting_boundary_loops_x:
                if loop[uv_layer].uv.y == 1.0:
                    boundary_sorted_x.append(loop)
                    if loop[uv_layer].uv.x not in boundary_unique_uy:
                        boundary_unique_uy.append(loop[uv_layer].uv.x)
                if loop[uv_layer].uv.y == 0.0:
                    boundary_sorted_x.append(loop)
                    if loop[uv_layer].uv.x not in boundary_unique_ly:
                        boundary_unique_ly.append(loop[uv_layer].uv.x)
                        
            for loop in sorting_boundary_loops_y:
                if loop[uv_layer].uv.x == 1.0:
                    boundary_sorted_y.append(loop)
                    if loop[uv_layer].uv.y not in boundary_unique_rx:
                        boundary_unique_rx.append(loop[uv_layer].uv.y)
                if loop[uv_layer].uv.x == 0.0:
                    boundary_sorted_y.append(loop)
                    if loop[uv_layer].uv.y not in boundary_unique_lx:
                        boundary_unique_lx.append(loop[uv_layer].uv.y)
                        
            max_iter = len(boundary_sorted_x) + len(boundary_sorted_y)
            iter = 0

            for loop in boundary_sorted_x:
                iter += 1
                if loop[uv_layer].uv.y == 1.0 and loop[uv_layer].uv.x != 1.0 and loop[uv_layer].uv.x != 0.0:
                    get_i = boundary_unique_uy.index(loop[uv_layer].uv.x)
                    newx = (1/(len(boundary_unique_uy)-1)) * get_i
                    loop[uv_layer].uv.x = newx
                if loop[uv_layer].uv.y == 0.0 and loop[uv_layer].uv.x != 1.0 and loop[uv_layer].uv.x != 0.0:
                    get_i = boundary_unique_ly.index(loop[uv_layer].uv.x)
                    newx = (1/(len(boundary_unique_ly)-1)) * get_i
                    loop[uv_layer].uv.x = newx
                if props.verbose == True:
                    done_pct = (iter / max_iter) * 100
                    if iter > 1:
                        print(CURSOR_UP + CURSOR_UP + CURSOR_UP)
                    msg = "Processing boundary loop {:} of {:} \nStep 3 {:.1f}% done.".format(iter, max_iter, done_pct)
                    print(msg, flush=True)
                    sys.stdout.flush()

            for loop in boundary_sorted_y:
                iter += 1
                if loop[uv_layer].uv.x == 1.0 and loop[uv_layer].uv.y != 1.0 and loop[uv_layer].uv.y != 0.0:
                    get_i = boundary_unique_rx.index(loop[uv_layer].uv.y)
                    newy = (1/(len(boundary_unique_rx)-1)) * get_i
                    loop[uv_layer].uv.y = newy
                if loop[uv_layer].uv.x == 0.0 and loop[uv_layer].uv.y != 1.0 and loop[uv_layer].uv.y != 0.0:
                    get_i = boundary_unique_lx.index(loop[uv_layer].uv.y)
                    newy = (1/(len(boundary_unique_lx)-1)) * get_i
                    loop[uv_layer].uv.y = newy
                if props.verbose == True:
                    done_pct = (iter / max_iter) * 100
                    if iter > 1:
                        print(CURSOR_UP + CURSOR_UP + CURSOR_UP)
                    msg = "Processing boundary loop {:} of {:} \nStep 3 {:.1f}% done.".format(iter, max_iter, done_pct)
                    print(msg, flush=True)
                    sys.stdout.flush()

            bpy.context.active_object.data.update()
        
        if props.space_edges == False:
            if props.verbose == True:
                print("Skipped even spacing of boundary loops.\n")
        
        if props.verbose == True:
            print(CURSOR_UP + ERASE_LINE + CURSOR_UP)
            print("Step 3 complete!\n")
            if props.space_edges == True:
                print("Step 3 results:")
                print("---------------")
                print(f"Points along left edge: {len(boundary_unique_lx)}")
                print(f"Points along right edge: {len(boundary_unique_rx)}")
                print(f"Points along bottom edge: {len(boundary_unique_ly)}")
                print(f"Points along top edge: {len(boundary_unique_uy)}\n")
            print("Step 4: Calculate direction and distance from island center for remaining loops")
            print("===============================================================================")
            
        # get remaining loops and calculate distance multiplier from distances of boundary loop points
        rest_loop_dist_factors = []
        rest_boundary_co = []

        iter = 0
        iter_max = len(inner_loops)

        for loop in inner_loops:
            iter += 1
            uva = loop[uv_layer]
            uvco = uva.uv
            vec = uvco - island_center
            dir = vec.normalized()
            dirdots = []
            for d in directions:
                dirdot = 1 - dir.dot(d)
                dirdots.append(dirdot)
            min_dirdot = min(dirdots)
            min_dirdot_i = dirdots.index(min_dirdot)
            dist_fac = vec.length / distances[min_dirdot_i]
            rest_loop_dist_factors.append(dist_fac)
            rest_boundary_co.append(boundary_co[min_dirdot_i])
            if props.verbose == True:
                done_pct = (iter / iter_max) * 100
                if iter > 1:
                    print(CURSOR_UP + CURSOR_UP + CURSOR_UP + CURSOR_UP)
                else:
                    print(CURSOR_UP)
                msg = "Calculated distance multiplier of {:.5f} for loop {:}. \nProcessed loop {:} of {:}. \nStep 4 {:.1f}% done.".format(dist_fac, loop.index, iter, iter_max, done_pct)
                print(msg, flush=True)
                sys.stdout.flush()
        
        if props.verbose == True:
            print(CURSOR_UP + ERASE_LINE + CURSOR_UP)
            print("Step 4 complete!\n")
            print("Step 5: Adjust positions of non-boundary loops")
            print("==============================================")
        
        # move non-boundary points
        iter = 0

        for loop, fac, bco in zip(inner_loops, rest_loop_dist_factors, rest_boundary_co):
            iter += 1
            uva = loop[uv_layer]
            uvco = uva.uv
            uvco_old = uvco
            mid = mathutils.Vector((0.5, 0.5))
            newx = bl_math.lerp(0.5, bco.x, fac)
            newy = bl_math.lerp(0.5, bco.y, fac)
            uvco.x = newx
            uvco.y = newy
            if props.verbose == True:
                done_pct = round((iter/len(inner_loops)) * 1000) * 0.1
                if iter > 1:
                    print(CURSOR_UP + CURSOR_UP + ERASE_LINE + CURSOR_UP + CURSOR_UP)
                msg = "Processing loop {:} of {:}. \nMoved point at {:.5}, {:.5} to new position at {:.5}, {:.5}. \nStep 5 {:.1f}% done.".format(iter, len(inner_loops), uvco_old.x, uvco_old.y, uvco.x, uvco.y, done_pct)
                print(msg, flush=True)
                sys.stdout.flush()
            
        bpy.context.active_object.data.update()
        
        if props.verbose == True:
            print(CURSOR_UP + ERASE_LINE + CURSOR_UP)
            print("Step 5 complete!\n")
        
        if props.relax_inner == True:
            for loop in inner_loops:
                uva = loop[uv_layer]
                uva.select = True
            try:
                bpy.ops.uv.univ_relax()
            except AttributeError:
                if props.verbose:
                    print("Could not find UniV add-on, skipping relax.\n")
            else:
                if props.verbose:
                    print("UniV add-on found, relaxing.\n")
            bpy.ops.uv.select_all(action='DESELECT')
            
        if props.verbose:
            print("Pelt Generation complete!\n")
        
        return {"FINISHED"}