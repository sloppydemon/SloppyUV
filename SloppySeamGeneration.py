import bpy
import bmesh
import math
import mathutils
import bl_math
from bpy_extras import bmesh_utils

class SloppySeamGen(bpy.types.Operator):
    bl_idname = "operator.sloppy_seam_gen"
    bl_label = "Generate Seams"
    bl_description = "Generate seams according to edge concavity, among other things"
    bl_options = {'REGISTER', 'UNDO'}

    iP = bpy.props.IntProperty
    fP = bpy.props.FloatProperty
    fvP = bpy.props.FloatVectorProperty
    bP = bpy.props.BoolProperty
    eP = bpy.props.EnumProperty
    bvP = bpy.props.BoolVectorProperty
    sP = bpy.props.StringProperty

    # region SeamGen Properties
    angle_factor : fP(
        name = "Concavity Influence",
        description = "Influence of concavity on seam generation",
        default = 1,
        min = -1,
        max = 1
        ) # type: ignore

    avg_angle_factor : fP(
        name = "Average Concavity Influence",
        description = "Influence of averaged area concavity on seam generation",
        default = 0,
        min = -1,
        max = 1
        ) # type: ignore

    ao_factor : fP(
        name = "AO Influence",
        description = "Influence of ambient occlusion on seam generation",
        default = 0,
        min = -1,
        max = 1
        ) # type: ignore

    ed_factor : fP(
        name = "Edge Density Influence",
        description = "Influence of edge density on seam generation",
        default = 0,
        min = -1,
        max = 1
        ) # type: ignore

    depth_factor : fP(
        name = "Depth Influence",
        description = "Influence of edge depth (compared to average Z of surrounding faces) on seam generation",
        default = 0,
        min = -1,
        max = 1
        ) # type: ignore
    
    depth_mode : eP(
        name = "Depth Calculation Mode",
        description = "Mode of calculating average depth",
        items = [
            ("A", "Radius", "Calculate average depth from vertices within radius"),
            ("B", "Connected", "Calculate average depth from vertices within number of connected edges")
            ],
        default="B"
        ) # type: ignore
    
    depth_avg_rad : fP(
        name = "Depth Average Radius",
        description = "Calculate depth average within this radius",
        default = 0.1,
        min = 0.001,
        max = 10000.0
        ) # type: ignore

    depth_avg_iters : iP(
        name = "Depth Average Iterations",
        description = "Number of edges out from current edge to calculate depth average",
        default = 2,
        min = 1,
        max = 1000
        ) # type: ignore
    
    depth_weight_by_dist : bP(
        name = "Weight by Distance",
        description = "Weigh depth average by distance from target",
        default = False
        ) # type: ignore

    no_rounds : iP(
        name = "Iterations",
        description = "Number of initial iterations",
        default = 20,
        min = 1,
        max = 1000
        ) # type: ignore

    no_retries : iP(
        name = "Max Retries",
        description = "Maximum number of retry iterations",
        default = 100,
        min = 1,
        max = 1000
        ) # type: ignore

    angle_thresh_start : fP(
        name = "Min Angle Threshold",
        description = "Initial angular threshold",
        default = -50,
        min = -180,
        max = 180
        ) # type: ignore

    angle_thresh_end : fP(
        name = "Max Angle Threshold",
        description = "Maximum angular threshold",
        default = -20,
        min = -180,
        max = 180
        ) # type: ignore

    add_normalized_attr : bP(
        name = "Add Normalized Attributes",
        description = "Generate additional normalized attributes - can be usual for debugging",
        default = True
        ) # type: ignore

    clear_seam : bP(
        name = "Clear Seams",
        description = "Clear any existing seams",
        default = True
        ) # type: ignore

    unwrap : bP(
        name = "Unwrap",
        description = "Run auto-unwrap. (Mainly to check if islands turn out as desired.)",
        default = True
        ) # type: ignore
    # endregion

    def find_depth(self, in_bm, edg, coo, nor, radius, num_out, mode, weight_by_dist):
        evs = []
        for ev in edg.verts:
            evs.append(ev)
        avg_depth = coo.copy()
        avg_norm = nor.copy()
        num_contributors = 1.0
        if mode == "A":
            max_dist = 0.0
            min_dist = 9999999.0
            dist_dat = []
            for vert in in_bm.verts:
                if vert not in evs:
                    vec = coo - vert.co
                    dist = vec.length
                    if dist <= radius:
                        numax = max(max_dist, dist)
                        numin = min(min_dist, dist)
                        max_dist = numax
                        min_dist = numin
                        dist_dat.append([dist, vert.co, vert.normal])
            dist_span = max_dist - min_dist
            for d in dist_dat:
                d_weight = 1.0
                if weight_by_dist == True:
                    if dist_span > 0.0:
                        d_weight = 1.0 - (d[0]-min_dist)/dist_span
                co_contribution = d[1] * d_weight
                no_contribution = d[2] * d_weight
                avg_depth += co_contribution
                avg_norm += no_contribution
                num_contributors += d_weight
        if mode == "B":
            done_verts = evs.copy()
            dist_dat = []
            out_iter = 0
            next_batch = evs.copy()
            weight_increment = 1.0/num_out
            while out_iter < num_out:
                this_batch = next_batch.copy()
                next_batch.clear()
                for bv in this_batch:
                    done_verts.append(bv)
                    for bve in bv.link_edges:
                        for bvev in bve.verts:
                            if bvev not in done_verts:
                                dist_dat.append([out_iter, bvev.co, bvev.normal])
                out_iter += 1
            for d in dist_dat:
                d_weight = 1.0
                if weight_by_dist == True:
                    d_weight = ((num_out + 1) - d[0]) * weight_increment
                co_contribution = d[1] * d_weight
                no_contribution = d[2] * d_weight
                avg_depth += co_contribution
                avg_norm += no_contribution
                num_contributors += d_weight
        avg_depth /= num_contributors
        avg_normal = avg_norm.normalized()
        covec = avg_depth - coo
        codir = covec.normalized()
        nor_dot = avg_normal.dot(codir)
        dots = covec.length
        dots *= (nor_dot * -1.0)
        print('Average ontributors:', num_contributors, ' - dots combined:', dots, '- average normal:', avg_normal, '- Source coordinates:', coo, '- avg. coordinates:', avg_depth)
        return(dots)



    def execute(self, context):
        props = context.scene.sloppy_props
        bm = bmesh.from_edit_mesh(bpy.context.active_object.data)

        verbose = props.verbose
        
        attr_dict = [
            {"name": "AO", "type": "FLOAT_COLOR", "domain": "POINT", "layer": None},
            {"name": "eAO", "type": "FLOAT", "domain": "EDGE", "layer": None},
            {"name": "depth", "type": "FLOAT", "domain": "EDGE", "layer": None},
            {"name": "edge_density", "type": "FLOAT", "domain": "EDGE", "layer": None},
            {"name": "avg_angle", "type": "FLOAT", "domain": "EDGE", "layer": None},
            {"name": "cost_to_boundary", "type": "FLOAT", "domain": "EDGE", "layer": None},
            {"name": "dist_to_boundary", "type": "FLOAT", "domain": "EDGE", "layer": None},
            ]
        
        if self.add_normalized_attr == True:
            attr_dict.append({"name": "normalized_depth", "type": "FLOAT", "domain": "EDGE", "layer": None})
            attr_dict.append({"name": "normalized_edge_density", "type": "FLOAT", "domain": "EDGE", "layer": None})
            attr_dict.append({"name": "normalized_avg_angle", "type": "FLOAT", "domain": "EDGE", "layer": None})

        for nam in attr_dict:
            attr = props.find_or_add_attribute(nam["name"], nam["type"], nam["domain"])
            nam["layer"] = props.get_attribute_layer(nam["name"], nam["type"], nam["domain"], bm)

        vAO = props.get_dict_layer("AO", attr_dict)
        eAO = props.get_dict_layer("eAO", attr_dict)
        depth = props.get_dict_layer("depth", attr_dict)
        edge_density = props.get_dict_layer("edge_density", attr_dict)
        avg_angle_lyr = props.get_dict_layer("avg_angle", attr_dict)
        cost_to_boundary = props.get_dict_layer("cost_to_boundary", attr_dict)
        dist_to_boundary = props.get_dict_layer("dist_to_boundary", attr_dict)

        min_edge_density = 9999.0
        max_edge_density = 0.0
        min_avg_angle = 0.0
        max_avg_angle = 0.0
        min_depth = 0.0
        max_depth = 0.0
        ao_fac = self.ao_factor
        ed_fac = self.ed_factor
        dpth_fac = self.depth_factor
        angle_fac = self.angle_factor
        avg_angle_fac = self.avg_angle_factor

        def angle_sort(e):
            result = 0.0
            ao = props.remap_val(e[eAO], 0, 1, -180, 180)
            ed = props.remap_val(e[edge_density], min_edge_density, max_edge_density, -180, 180)
            dpth = props.remap_val(e[depth], min_depth, max_depth, -180, 180)
            if len(e.link_faces) > 1:
                result += (ao * ao_fac)
                result += (ed * ed_fac)
                result += (dpth * dpth_fac)
                result += (math.degrees(e.calc_face_angle_signed()) * angle_fac)
                result += (e[avg_angle_lyr] * avg_angle_fac)
            return result

        def cost_sort(e):
            return e[cost_to_boundary]
        def dist_sort(e):
            return e[dist_to_boundary]

        def calc_cost_to_boundary(this_edge, curr_seam, other_seam, has_max_cost, max_cost):
            done = []
            boundary_found = False
            cost = 0
            next_round = [this_edge]
            while boundary_found == False:
                this_round = next_round.copy()
                next_round = []
                for edge in this_round:
                    for v in edge.verts:
                        for ve in v.link_edges:
                            if ve not in done and ve.index != edge.index:
                                next_round.append(ve)
                                if other_seam == True:
                                    if ve.seam == True:
                                        if ve not in curr_seam:
                                            boundary_found = True
                                            return cost
                                else:
                                    if ve.seam == True or len(ve.link_faces) == 1 or ve.is_boundary:
                                        boundary_found = True
                                        return cost
                            if boundary_found == True:
                                break
                        if boundary_found == True:
                            break
                    if boundary_found == True:
                        break
                    done.append(edge)
                cost += 1
                if has_max_cost == True:
                    if cost >= max_cost:
                        return cost
                if verbose == True:
                    print(f"Current cost for edge {this_edge.index}: {cost}")
            return cost

        def calc_dist_to_boundary(this_edge, curr_seam, other_seam, has_max_dist, max_dist):
            done = []
            boundary_found = False
            dist = 0.0
            next_round = [this_edge]
            while boundary_found == False:
                this_round = next_round.copy()
                next_round = []
                this_dist = 0.0
                num_dists = 0
                for edge in this_round:
                    this_dist += edge.calc_length()
                    num_dists += 1
                    for v in edge.verts:
                        for ve in v.link_edges:
                            if ve not in done and ve.index != edge.index:
                                next_round.append(ve)
                                if other_seam == True:
                                    if ve.seam == True:
                                        if ve not in curr_seam:
                                            boundary_found = True
                                            return dist
                                else:
                                    if ve.seam == True or len(ve.link_faces) == 1 or ve.is_boundary:
                                        boundary_found = True
                                        return dist
                            if boundary_found == True:
                                break
                        if boundary_found == True:
                            break
                    if boundary_found == True:
                        break
                    done.append(edge)
                
                dist += this_dist / num_dists
                if has_max_dist == True:
                    if dist >= max_dist:
                        return dist
                if verbose == True:
                    print(f"Current distance for edge {this_edge.index}: {dist}")
            return dist

        all_edges = []

        if self.clear_seam == True:
            for edge in bm.edges:
                edge.seam = False

        e_iter = 0
        for edge in bm.edges:
            if edge.seam == False:
                all_edges.append(edge)
                ao = 0.0
                avg_edge_length = 0.0
                avg_angle = 0.0
                num_edges = 0
                num_angles = 0
                for v in edge.verts:
                    ao += v[vAO].x
                    for e in v.link_edges:
                        if e.index != edge.index:
                            num_edges += 1
                            avg_edge_length += e.calc_length()
                near_pars = props.find_near_parallels(edge)
                if len(near_pars) > 1:
                    for e in near_pars:
                        if len(e.link_faces) == 2:
                            num_angles += 1
                            avg_angle += math.degrees(e.calc_face_angle_signed())
                    avg_edge_length /= num_edges
                    avg_angle /= num_angles
                this_depth = 0.0
                if self.depth_factor > 0.0:
                    this_depth = self.find_depth(bm, edge, props.calc_edge_center(edge), props.calc_edge_avg_normal(edge), self.depth_avg_rad, self.depth_avg_iters, self.depth_mode, self.depth_weight_by_dist)
                new_max_density = max(avg_edge_length, max_edge_density)
                new_min_density = min(avg_edge_length, min_edge_density)
                new_max_angle = max(avg_angle, max_avg_angle)
                new_min_angle = min(avg_angle, min_avg_angle)
                new_max_depth = max(max_depth, this_depth)
                new_min_depth = min(min_depth, this_depth)
                max_depth = new_max_depth
                min_depth = new_min_depth
                max_edge_density = new_max_density
                min_edge_density = new_min_density
                max_avg_angle = new_max_angle
                min_avg_angle = new_min_angle
                ao /= 2
                edge[edge_density] = avg_edge_length
                edge[avg_angle_lyr] = avg_angle
                edge[eAO] = ao
                edge[depth] = this_depth
            e_iter += 1
            remain_num = len(bm.edges) - e_iter
            if verbose == True:
                print(f"{e_iter} edges done, {remain_num} remaining.")

        all_edges.sort(key=angle_sort)

        remain_edges = all_edges.copy()

        if self.add_normalized_attr == True:
            ndpth = props.get_dict_layer("normalized_depth", attr_dict)
            ned = props.get_dict_layer("normalized_edge_density", attr_dict)
            naa = props.get_dict_layer("normalized_avg_angle", attr_dict)
            for edg in all_edges:
                if self.depth_factor > 0.0:
                    if (max_depth - min_depth) > 0.0:
                        edg[ndpth] = (edg[depth] - min_depth) / (max_depth - min_depth)
                if (max_edge_density - min_edge_density) > 0.0:
                    edg[ned] = (edg[edge_density] - min_edge_density) / (max_edge_density - min_edge_density)
                if (max_avg_angle - min_avg_angle) > 0.0:
                    edg[naa] = (edg[avg_angle_lyr] - min_avg_angle) / (max_avg_angle - min_avg_angle)

        seams = []

        rounds = self.no_rounds
        max_retries = self.no_retries

        near_parallel_edges = []

        angle_threshold_start = self.angle_thresh_start
        angle_threshold_end = self.angle_thresh_end
        angle_threshold_interval = angle_threshold_end - angle_threshold_start
        angle_threshold_point = 0
        if max_retries > 0:
            angle_threshold_point = angle_threshold_interval / max_retries

        for i in range(rounds):
            remain_edges.sort(key=angle_sort)
            
            seam_ended = False
            next_round = [remain_edges[0]]
            verts_done = []
            current_adjacent = []
            current_adjacent_faces = []
            current_seam = []
            retries = 0

            while seam_ended == False:
                this_round = next_round.copy()
                next_round = []
                angle_threshold = angle_threshold_start
                
                for edge in this_round:
                    edge.seam = True
                    if verbose == True:
                        print(f"Added {edge.index} to seam.")
                    current_seam.append(edge)
                    near_para = props.find_near_parallels(edge)
                    for i in near_para:
                        near_parallel_edges.append(i)
                    if edge in remain_edges:
                        remain_edges.remove(edge)
                    for f in edge.link_faces:
                        current_adjacent_faces.append(f)
                    for v in edge.verts:
                        if v not in verts_done:
                            vert_edges = []
                            for eb in v.link_edges:
                                if eb.link_faces not in current_adjacent_faces:
                                    current_adjacent_faces.append(f)
                                    if len(eb.link_faces) > 1:
                                        if eb.index != edge.index and eb.seam == False and eb not in current_seam and eb not in current_adjacent and eb not in near_parallel_edges and angle_sort(eb) < angle_threshold:
                                            vert_edges.append(eb)
                                            current_adjacent.append(eb)
                            if len(vert_edges) > 0:
                                vert_edges.sort(key=angle_sort)
                                next_round.append(vert_edges[0])
                            verts_done.append(v)
                
                if len(next_round) == 0:
                    if retries < max_retries:
                        retries += 1
                        retry_round = []
                        retry_round.append(current_seam[-1])
                        angle_threshold = angle_threshold_start + (angle_threshold_point * retries)
                        if props.verbose == True:
                            print(f"Retrying edge {current_seam[-1].index} with angle threshold {angle_threshold}. Retry #{retries}.")
                        for edge in retry_round:
                            for v in edge.verts:
                                vert_edges = []
                                for eb in v.link_edges:
                                    if eb.link_faces not in current_adjacent_faces:
                                        current_adjacent_faces.append(f)
                                        if len(eb.link_faces) > 1:
                                            if eb.index != edge.index and eb.seam == False and eb not in current_seam and eb not in current_adjacent and eb not in near_parallel_edges and angle_sort(eb) < angle_threshold:
                                                vert_edges.append(eb)
                                                current_adjacent.append(eb)
                                if len(vert_edges) > 0:
                                    vert_edges.sort(key=angle_sort)
                                    next_round.append(vert_edges[0])
                    else:
                        seam_ended = True

            if len(current_seam) < 2:
                for edge in current_seam:
                    edge.seam = False
            else:
                seams.append(current_seam)

        bpy.context.active_object.data.update()

        if len(seams) > 1:
            for seam in seams:
                for edge in seam:
                    this_dist = calc_dist_to_boundary(edge, current_seam, True, True, 0.5)
                    edge[dist_to_boundary] = this_dist
                seam.sort(key=dist_sort)
                
                next_round = [seam[0]]
                seam_ended = False
                verts_done = []
                current_adjacent = []
                current_adjacent_faces = []
                current_seam = []

                while seam_ended == False:
                    this_round = next_round.copy()
                    next_round = []
                    
                    for edge in this_round:
                        edge.seam = True
                        if edge not in seam:
                            seam.append(edge)
                        if edge not in current_seam:
                            current_seam.append(edge)
                        for v in edge.verts:
                            vert_edges = []
                            for eb in v.link_edges:
                                if eb.seam == True and eb not in seam:
                                    seam_ended = True
                                    break
                                if eb.link_faces not in current_adjacent_faces:
                                    current_adjacent_faces.append(f)
                                    if len(eb.link_faces) > 1:
                                        if eb.index != edge.index and eb.seam == False and eb not in current_seam and eb not in current_adjacent and eb not in near_parallel_edges:
                                            if props.verbose == True:
                                                print(f"Edge index: {edge.index} Other edge index: {eb.index}")
                                            vert_edges.append(eb)
                                            current_adjacent.append(eb)
                                            this_dist = calc_dist_to_boundary(eb, seam, False, True, 0.5)
                                            eb[dist_to_boundary] = this_dist
                            if len(vert_edges) > 0:
                                vert_edges.sort(key=dist_sort)
                                next_round.append(vert_edges[0])
                    if len(next_round) == 0:
                        seam_ended = True
            
                    if len(current_seam) > 5:
                        for edge in current_seam:
                            edge.seam = False
        bpy.context.active_object.data.update()

        if self.unwrap == True:
            selected_before = []
            loops_selected_before = []
            for face in bm.faces:
                if face.select == True:
                    selected_before.append(face)
                face.select = True
            bpy.ops.uv.unwrap(method='ANGLE_BASED', margin=0)
            for face in bm.faces:
                if face in selected_before:
                    pass
                else:
                    face.select = False

        
        return {"FINISHED"}

class SloppyIslandBoundaryConnect(bpy.types.Operator):
    bl_idname = "operator.sloppy_boundary_connect"
    bl_label = "Connect Island Boundaries"
    bl_description = "Try to connect island boundaries if there are more than one (an example would be a cylindrical shape with no seam to split it into a strip, instead having a ring-shaped UV-map)"
    bl_options = {'REGISTER', 'UNDO'}

    iP = bpy.props.IntProperty
    fP = bpy.props.FloatProperty
    fvP = bpy.props.FloatVectorProperty
    ivP = bpy.props.IntVectorProperty
    bP = bpy.props.BoolProperty
    eP = bpy.props.EnumProperty
    bvP = bpy.props.BoolVectorProperty
    sP = bpy.props.StringProperty

    island_selection_mode: eP(
        name = "Island Selection Mode",
        description = "Mode of selecting which island to do",
        items = [
            ("A", "Selection As Island", ""),
            ("B", "From Selected Face", ""),
            ("C", "Index", "")
            ],
        default="A"
        ) # type: ignore

    island_index : iP(
        name = "Island Index",
        description = "",
        default = 0
        ) # type: ignore

    z_factor : fP(
        name = "Z Factor",
        description = "How much influence Z axis has on connection calculation",
        default = 1.0,
        max=100.0,
        min=-100.0
        ) # type: ignore

    concavity_factor : fP(
        name = "Concavity Factor",
        description = "How much influence concavity has on connection calculation",
        default = 1.0,
        max=100.0,
        min=-100.0
        ) # type: ignore

    proximity_factor : fP(
        name = "Proximity Factor",
        description = "How much influence proximity to other boundary has on connection calculation",
        default = 1.0,
        max=100.0,
        min=-100.0
        ) # type: ignore

    verbose : bP(
        name = "Verbose",
        description = "Print debug information to system console",
        default = True
        ) # type: ignore
    
    def sort_by_z_angle_cost(self, dat):
        z_val = ((dat[5]-dat[6])+(dat[2]-dat[5]))*self.z_factor
        a_val = ((dat[3]+180.0)/180.0)*self.concavity_factor
        c_val = float(dat[4])*self.proximity_factor
        combined_val = z_val + a_val + c_val
        print('Edge', dat[1].index, 'to vertex', dat[0].index, '- Sorting by values', z_val, '+', a_val, '+', c_val, '=', combined_val)
        return combined_val

    def execute(self, context):
        props = context.scene.sloppy_props
        bm = bmesh.from_edit_mesh(bpy.context.active_object.data)
        uv_layer = bm.loops.layers.uv.verify()
        islands = bmesh_utils.bmesh_linked_uv_islands(bm, uv_layer)
        do_connect = False

        island_index_clamped = int(bl_math.clamp(float(self.island_index), 0, len(islands) - 1))

        island = islands[island_index_clamped]
        if self.island_selection_mode == "A":
            island = []
            for face in bm.faces:
                if face.select == True:
                    island.append(face)
        if self.island_selection_mode == "B":
            for ii, isl in enumerate(islands):
                for face in isl:
                    if face.select == True:
                        island = islands[ii]
                        break
        
        loop_chains, loop_chains_eton, vert_chains, edge_chains, face_chains = props.find_island_boundaries_and_their_loops(bm, island, True, self.verbose)

        if len(loop_chains) > 1:
            do_connect = True
        else:
            if self.verbose:
                print('\nConnecting of boundaries unnecessary since there is only one boundary.')
        
        island_edges = []
        new_seam = []
        

        if do_connect == True:
            this_loop_chain = []
            this_vert_chain = []
            this_edge_chain = []
            next_chain = len(loop_chains) - 1
            other_verts = []
            other_edges = []
            other_loops = []
            this_loop_chain = loop_chains[1].copy()
            print(this_loop_chain)
            this_vert_chain = vert_chains[1].copy()
            print(this_vert_chain)
            this_edge_chain = edge_chains[1].copy()
            print(this_edge_chain)
            for isf in island:
                for isfe in isf.edges:
                    if isfe not in island_edges:
                        island_edges.append(isfe)
            for ec in edge_chains:
                for ece in ec:
                    if ece not in this_edge_chain:
                        other_edges.append(ece)
                        for ecev in ece.verts:
                            if ecev not in other_verts:
                                other_verts.append(ecev)
                    if ece in this_edge_chain:
                        for ecev in ece.verts:
                            if ecev not in this_vert_chain:
                                this_vert_chain.append(ecev)
            for lv in loop_chains:
                for lvl in lv:
                    if lvl not in this_loop_chain:
                        other_loops.append(lvl)
                    else:
                        if lvl.vert not in this_vert_chain:
                            this_vert_chain.append(lvl.vert)

            print('\nThere are still', len(loop_chains), 'loop chains in island.')
            
            other_boundary_found = False
            next_verts = []
            for vcv in this_vert_chain:
                if vcv not in next_verts:
                    next_verts.append(vcv)

            first_z = None
            max_z = None

            b_iter = 0
            while other_boundary_found == False:
                if b_iter > 0:
                    print('Iteration', b_iter, '- Other boundary still not found.\n')
                b_iter += 1
                next_indices = []
                pot_max_z = -999999.9
                for nv in next_verts:
                    next_indices.append(nv.index)
                    if not max_z:
                        new_pot_max = max(pot_max_z, nv.co.z)
                        pot_max_z = new_pot_max
                if not max_z:
                    max_z = pot_max_z
                        
                print('Vertices to check:', next_indices)
                these_verts = next_verts.copy()
                next_verts.clear()
                tve_dats = []
                for tv in these_verts:
                    for tve in tv.link_edges:
                        if tve not in new_seam:
                            if tve in island_edges:
                                print('Vertex', tv.index, '> Edge', tve.index, '>')
                                other_vert = tve.other_vert(tv)
                                other_vert_loops = props.edge_corners_per_vert(tve, other_vert)
                                edge_viable = True
                                for ovl in other_vert_loops:
                                    if ovl in this_loop_chain:
                                        if edge_viable == True:
                                            edge_viable = False
                                            print('Edge', tve.index, 'not viable.')
                                for ovl in other_vert.link_loops:
                                    if ovl in other_loops:
                                        if other_boundary_found == False:
                                            other_boundary_found = True
                                            print('Other boundary found!\n')
                                            tve_dat = []
                                            tve_dat.append(other_vert)
                                            tve_dat.append(tve)
                                            tve_dat.append(other_vert.co.z)
                                            try:
                                                tve_dat.append(tve.calc_face_angle_signed())
                                            except:
                                                tve_dat.append(0.0)
                                            tve_dat.append(int(-2.0*self.proximity_factor))
                                            if first_z:
                                                tve_dat.append(first_z)
                                            else:
                                                tve_dat.append(tv.co.z)
                                            tve_dat.append(max_z)
                                            tve_dats.append(tve_dat)
                                if edge_viable == True:
                                    if other_boundary_found == False:
                                        tve_dat = []
                                        tve_dat.append(other_vert)
                                        tve_dat.append(tve)
                                        tve_dat.append(other_vert.co.z)
                                        try:
                                            tve_dat.append(tve.calc_face_angle_signed())
                                        except:
                                            tve_dat.append(0.0)
                                        found_cost = False
                                        next_vert_batch = [other_vert]
                                        visited_edges = [tve]
                                        current_cost = 0
                                        while found_cost == False:
                                            print('Total cost for edge', tve.index, 'still not found. Current cost:', current_cost)
                                            this_vert_batch = next_vert_batch.copy()
                                            next_vert_batch.clear()
                                            current_cost += 1
                                            for tvbv in this_vert_batch:
                                                for tvbve in tvbv.link_edges:
                                                    if tvbve in island_edges:
                                                        if tvbve not in visited_edges:
                                                            tvbveov = tvbve.other_vert(tvbv)
                                                            for tvbveovl in tvbveov.link_loops:
                                                                if tvbveovl in other_loops:
                                                                    if found_cost == False:
                                                                        found_cost = True
                                                                        print('Found cost, it was', current_cost, 'at vert', tvbveov.index)
                                                            if found_cost == False:
                                                                next_vert_batch.append(tvbveov)
                                                            visited_edges.append(tvbve)
                                        tve_dat.append(current_cost)
                                        if first_z:
                                            tve_dat.append(first_z)
                                        else:
                                            tve_dat.append(tv.co.z)
                                        tve_dat.append(max_z)
                                        tve_dats.append(tve_dat)
                tve_dats.sort(key=self.sort_by_z_angle_cost,reverse=True)
                print('\nSorted:')
                tve_dats.sort(key=self.sort_by_z_angle_cost)
                next_verts = [tve_dats[0][0]]
                first_z = tve_dats[0][5]
                print('Edge', tve_dats[0][1].index, 'set as seam.\n')
                tve_dats[0][1].seam = True
                if tve_dats[0][1] not in new_seam:
                    new_seam.append(tve_dats[0][1])
            
            bpy.context.active_object.data.update()
            
            for bmf in bm.faces:
                bmf.select = False

            for isf in island:
                isf.select = True
            bpy.ops.uv.unwrap(method='ANGLE_BASED', margin=0)
            
        bpy.context.active_object.data.update()
        if do_connect == True:
            bpy.ops.uv.unwrap(method='ANGLE_BASED', margin=0)

        return {"FINISHED"}

class SloppySeamGenVis(bpy.types.Operator):
    bl_idname = "operator.sloppy_seam_gen_vis"
    bl_label = "Generate Visibility Seams"
    bl_description = "Generate seams according to visibility"
    bl_options = {'REGISTER', 'UNDO'}

    iP = bpy.props.IntProperty
    fP = bpy.props.FloatProperty
    fvP = bpy.props.FloatVectorProperty
    ivP = bpy.props.IntVectorProperty
    bP = bpy.props.BoolProperty
    eP = bpy.props.EnumProperty
    bvP = bpy.props.BoolVectorProperty
    sP = bpy.props.StringProperty

    # region SeamGen Properties
    angular_resolution : ivP(
        name = "Angular Resolution",
        description = "Visibility resolution in number of raycasting points in full rotation",
        default=(1, 4, 4),
        min=2,
        max=512,
        size = 3
        ) # type: ignore

    dome : bP(
        name = "Dome",
        description = "No raycast from underneath object",
        default = True
        ) # type: ignore
    
    init_height_factor : fP(
        name = "Height Factor",
        description = "Height as factor times object's height",
        default = 0.5,
        min = 0,
        max = 1
        ) # type: ignore
    
    vis_threshold : fP(
        name = "Visibility Threshold",
        description = "Edges below this threshold will be marked as seams",
        default = 0.5,
        min = 0,
        max = 1
        ) # type: ignore
    
    angle_factor : fP(
        name = "Concavity Influence",
        description = "Influence of concavity on seam generation",
        default = 1,
        min = 0,
        max = 1
        ) # type: ignore

    avg_angle_factor : fP(
        name = "Average Concavity Influence",
        description = "Influence of averaged area concavity on seam generation",
        default = 0,
        min = 0,
        max = 1
        ) # type: ignore

    ao_factor : fP(
        name = "AO Influence",
        description = "Influence of ambient occlusion on seam generation",
        default = 0,
        min = 0,
        max = 1
        ) # type: ignore

    ed_factor : fP(
        name = "Edge Density Influence",
        description = "Influence of edge density on seam generation",
        default = 0,
        min = 0,
        max = 1
        ) # type: ignore

    no_rounds : iP(
        name = "Iterations",
        description = "Number of initial iterations",
        default = 20,
        min = 1,
        max = 1000
        ) # type: ignore

    no_retries : iP(
        name = "Max Retries",
        description = "Maximum number of retry iterations",
        default = 100,
        min = 1,
        max = 1000
        ) # type: ignore

    angle_thresh_start : fP(
        name = "Min Angle Threshold",
        description = "Initial angular threshold",
        default = -50,
        min = -180,
        max = 180
        ) # type: ignore

    angle_thresh_end : fP(
        name = "Max Angle Threshold",
        description = "Maximum angular threshold",
        default = -20,
        min = -180,
        max = 180
        ) # type: ignore

    clear_seam : bP(
        name = "Clear Seams",
        description = "Clear any existing seams",
        default = True
        ) # type: ignore

    unwrap : bP(
        name = "Unwrap",
        description = "Run auto-unwrap. (Mainly to check if islands turn out as desired.)",
        default = False
        ) # type: ignore
    # endregion

    def execute(self, context):
        props = context.scene.sloppy_props
        bm = bmesh.from_edit_mesh(bpy.context.active_object.data)
        depsgraph = bpy.context.evaluated_depsgraph_get()
        verbose = props.verbose

        init_dist = max(bpy.context.active_object.dimensions) * 1.5
        init_height = bpy.context.active_object.dimensions.z * self.init_height_factor
        init_pov = mathutils.Vector((init_dist, 0.0, init_height))
        curr_pov = init_pov
        
        attr_dict = [
            {"name": "edge_vis", "type": "FLOAT", "domain": "EDGE", "layer": None},
            {"name": "vert_vis", "type": "FLOAT", "domain": "POINT", "layer": None},
            ]

        for nam in attr_dict:
            attr = props.find_or_add_attribute(nam["name"], nam["type"], nam["domain"])
            nam["layer"] = props.get_attribute_layer(nam["name"], nam["type"], nam["domain"], bm)

        edge_vis = props.get_dict_layer("edge_vis", attr_dict)
        vert_vis = props.get_dict_layer("vert_vis", attr_dict)

        y_rot_inc = 360 / (self.angular_resolution[1] - 1)
        if self.dome == True:
            y_rot_inc = 180 / (self.angular_resolution[1] - 1)
        z_rot_inc = 360 / (self.angular_resolution[2])

        max_iter = self.angular_resolution[1] * self.angular_resolution[2] * len(bm.verts)

        iter = 0

        for yi in range(self.angular_resolution[1]):
            y_rot = mathutils.Euler((0,y_rot_inc*yi,0), 'XYZ')
            y_rotated = init_pov.rotate(y_rot)
            if y_rotated == None:
                y_rotated = init_pov
            for zi in range(self.angular_resolution[2]):
                z_rot = mathutils.Euler((0,0,z_rot_inc*zi))
                z_rotated = y_rotated.rotate(z_rot)
                if z_rotated == None:
                    z_rotated = y_rotated
                curr_pov = z_rotated

                for vert in bm.verts:
                    percent_done = (iter/max_iter) * 100
                    print(f"{percent_done} percent done")
                    iter += 1
                    vert[vert_vis] += 1
                    vec = vert.co - curr_pov
                    dir = vec.normalized()
                    hit, hit_loc, hit_norm, hit_i, hit_o, hit_ab = bpy.context.scene.ray_cast(depsgraph, curr_pov, dir)
                    if hit == True:
                        if bm.faces[hit_i] not in vert.link_faces:
                            vert[vert_vis] -= 1
        
        vis_div = self.angular_resolution[1] * self.angular_resolution[2]

        bpy.context.active_object.data.update()

        for vert in bm.verts:
            vert[vert_vis] /= vis_div
            print(f"Visibility for vert {vert.index}: {vert[vert_vis]}")

        bpy.context.active_object.data.update()

        for edge in bm.edges:
            edge[edge_vis] = (edge.verts[0][vert_vis] + edge.verts[1][vert_vis]) / 2
            print(f"Visibility for edge {edge.index}: {edge[edge_vis]}")

        bpy.context.active_object.data.update()

        if self.unwrap == True:
            selected_before = []
            for face in bm.faces:
                if face.select == True:
                    selected_before.append(face)
                face.select = True
            bpy.ops.uv.unwrap(method='ANGLE_BASED', margin=0)
            for face in bm.faces:
                if face in selected_before:
                    pass
                else:
                    face.select = False

        
        return {"FINISHED"}