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
        dat = bpy.context.object.data
        bm = bmesh.from_edit_mesh(bpy.context.active_object.data)
        uv_layer = bm.loops.layers.uv.verify()
        islands = bmesh_utils.bmesh_linked_uv_islands(bm, uv_layer)
        
        attribute_dict = [
            {"name": "AO", "type": "FLOAT_COLOR", "domain": "POINT", "layer": None},
            {"name": "eAO", "type": "FLOAT", "domain": "EDGE", "layer": None},
            {"name": "edge_density", "type": "FLOAT", "domain": "EDGE", "layer": None},
            {"name": "avg_angle", "type": "FLOAT", "domain": "EDGE", "layer": None},
            {"name": "cost_to_boundary", "type": "FLOAT", "domain": "EDGE", "layer": None},
            {"name": "dist_to_boundary", "type": "FLOAT", "domain": "EDGE", "layer": None},
            ]

        def find_or_add_attribute(name, attr_type, attr_domain):
            attribute = dat.attributes[0]
            get_i = dat.attributes.find(name)
            if get_i == -1:
                attribute = dat.attributes.new(name=name, type=attr_type, domain=attr_domain)
            else:
                attribute = dat.attributes[get_i]
            return attribute

        def get_attribute_layer(name, attr_type, attr_domain):
            layer = None
            if attr_domain == "FACE":
                if attr_type == "INT":
                    layer = bm.faces.layers.int.get(name)
                if attr_type == "FLOAT_VECTOR":
                    layer = bm.faces.layers.float_vector.get(name)
                if attr_type == "FLOAT":
                    layer = bm.faces.layers.float.get(name)
            if attr_domain == "POINT":
                if attr_type == "FLOAT_COLOR":
                    layer = bm.verts.layers.float_color.get(name)
            if attr_domain == "EDGE":
                if attr_type == "FLOAT":
                    layer = bm.edges.layers.float.get(name)
            return layer

        def get_dict_layer(name):
            layer = None
            for dict in attribute_dict:
                if dict["name"] == name:
                    layer = dict["layer"]
                    return layer

        def remap_val(val, in_min, in_max, out_min, out_max):
            in_interval = in_max - in_min
            out_interval = out_max - out_min
            in_val = val - in_min
            in_fac = in_val / in_interval
            out_val = out_min + (in_fac * out_interval)
            if out_val > out_max:
                out_val = out_max
            if out_val < out_min:
                out_val = out_min
            return out_val

        for nam in attribute_dict:
            attr = find_or_add_attribute(nam["name"], nam["type"], nam["domain"])
            nam["layer"] = get_attribute_layer(nam["name"], nam["type"], nam["domain"])

        min_edge_density = 9999.0
        max_edge_density = 0.0
        min_avg_angle = 0.0
        max_avg_angle = 0.0
        ao_fac = props.seamgen_ao_fac
        ed_fac = props.seamgen_ed_fac
        angle_fac = props.seamgen_angle_fac
        avg_angle_fac = props.seamgen_avg_angle_fac

        def angle_sort(e):
            result = 0.0
            ao = remap_val(e[get_dict_layer("eAO")], 0, 1, -180, 180)
            ed = remap_val(e[get_dict_layer("edge_density")], min_edge_density, max_edge_density, -180, 180)
            if len(e.link_faces) > 1:
                result += (ao * ao_fac)
                result += (ed * ed_fac)
                result += (math.degrees(e.calc_face_angle_signed()) * angle_fac)
                result += (e[get_dict_layer("avg_angle")] * avg_angle_fac)
            return result

        def cost_sort(e):
            return e[get_dict_layer("cost_to_boundary")]
        def dist_sort(e):
            return e[get_dict_layer("dist_to_boundary")]

        def find_near_parallels(e):
            near_parallels = []
            for f in e.link_faces:
                for edge in f.edges:
                    across = False
                    for v in edge.verts:
                        if e in v.link_edges:
                            across = True
                            break
                    if across == False:
                        near_parallels.append(edge)
            return near_parallels

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
                print(f"Current distance for edge {this_edge.index}: {dist}")
            return dist

        all_edges = []
        max_cost_to_boundary = 0

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
                    ao += v[bm.verts.layers.float_color['AO']].x
                    for e in v.link_edges:
                        if e.index != edge.index:
                            num_edges += 1
                            avg_edge_length += e.calc_length()
        #                    if len(e.link_faces) == 2:
        #                        avg_angle += math.degrees(e.calc_face_angle_signed())
        #                        num_angles += 1
                near_pars = find_near_parallels(edge)
                if len(near_pars) > 1:
                    for e in near_pars:
                        if len(e.link_faces) == 2:
                            num_angles += 1
                            avg_angle += math.degrees(e.calc_face_angle_signed())
                    avg_edge_length /= num_edges
                    avg_angle /= num_angles
                new_max_density = max(avg_edge_length, max_edge_density)
                new_min_density = min(avg_edge_length, min_edge_density)
                new_max_angle = max(avg_angle, max_avg_angle)
                new_min_angle = min(avg_angle, min_avg_angle)
                max_edge_density = new_max_density
                min_edge_density = new_min_density
                max_avg_angle = new_max_angle
                min_avg_angle = new_min_angle
                ao /= 2
                edge[get_dict_layer("edge_density")] = avg_edge_length
                edge[get_dict_layer("avg_angle")] = avg_angle
                edge[get_dict_layer("eAO")] = ao
            e_iter += 1
            remain_num = len(bm.edges) - e_iter
            print(f"{e_iter} edges done, {remain_num} remaining.")
                

        print(min_edge_density)
        print(max_edge_density)

        all_edges.sort(key=angle_sort)

        #all_edges_angle = []
        #for edge in all_edges:
        #    if len(edge.link_faces) > 1:
        #        str = "{:.3}".format(math.degrees(angle_sort))
        #        all_edges_angle.append(str)
        #print(f"Edges: {all_edges_angle}")

        remain_edges = all_edges.copy()

        seams = []

        rounds = props.seamgen_rounds
        max_retries = props.seamgen_retries

        near_parallel_edges = []

        angle_threshold_start = props.seamgen_angle_threshold_start
        angle_threshold_end = props.seamgen_angle_threshold_end
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
            max_cost_to_boundary = 0
            retries = 0

            while seam_ended == False:
                this_round = next_round.copy()
                next_round = []
                angle_threshold = angle_threshold_start
                
                for edge in this_round:
                    print(angle_sort(edge))
                    edge.seam = True
                    print(f"Added {edge.index} to seam.")
                    current_seam.append(edge)
                    near_para = find_near_parallels(edge)
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
        #                                    this_cost = calc_cost_to_boundary(eb, current_seam)
        #                                    new_max_cost_to_boundary = max(max_cost_to_boundary, this_cost)
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
        #                                        this_cost = calc_cost_to_boundary(eb, current_seam, False, True, 8)
        #                                        eb[get_dict_layer("cost_to_boundary")] = this_cost
        #                                        this_dist = calc_dist_to_boundary(eb, current_seam, False, True, 0.5)
        #                                        eb[get_dict_layer("dist_to_boundary")] = this_dist
                                if len(vert_edges) > 0:
        #                            vert_edges.sort(key=cost_sort)
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
        #            this_cost = calc_cost_to_boundary(eb, current_seam, True, True, 8)
        #            eb[get_dict_layer("cost_to_boundary")] = this_cost
        #        seam.sort(key=cost_sort)
                    this_dist = calc_dist_to_boundary(edge, current_seam, True, True, 0.5)
                    edge[get_dict_layer("dist_to_boundary")] = this_dist
                seam.sort(key=dist_sort)
                
                next_round = [seam[0]]
                seam_ended = False
                verts_done = []
                current_adjacent = []
                current_adjacent_faces = []
                current_seam = []
                max_cost_to_boundary = 0

                while seam_ended == False:
                    this_round = next_round.copy()
                    next_round = []
                    
                    for edge in this_round:
                        other_seam_found = False
                        edge.seam = True
                        if edge not in seam:
                            seam.append(edge)
                        if edge not in current_seam:
                            current_seam.append(edge)
                        for v in edge.verts:
                            vert_edges = []
                            for eb in v.link_edges:
                                if eb.seam == True and eb not in seam:
                                    other_seam_found = True
                                    seam_ended = True
                                    break
                                if eb.link_faces not in current_adjacent_faces:
                                    current_adjacent_faces.append(f)
                                    if len(eb.link_faces) > 1:
                                        if eb.index != edge.index and eb.seam == False and eb not in current_seam and eb not in current_adjacent and eb not in near_parallel_edges:
                                            print(f"Edge index: {edge.index} Other edge index: {eb.index}")
                                            vert_edges.append(eb)
                                            current_adjacent.append(eb)
        #                                    this_cost = calc_cost_to_boundary(eb, seam, False, True, 8)
        #                                    eb[get_dict_layer("cost_to_boundary")] = this_cost
                                            this_dist = calc_dist_to_boundary(eb, seam, False, True, 0.5)
                                            eb[get_dict_layer("dist_to_boundary")] = this_dist
                            if len(vert_edges) > 0:
        #                        vert_edges.sort(key=cost_sort)
                                vert_edges.sort(key=dist_sort)
                                next_round.append(vert_edges[0])
                    if len(next_round) == 0:
                        seam_ended = True
            
                    if len(current_seam) > 5:
                        for edge in current_seam:
                            edge.seam = False
        bpy.context.active_object.data.update()
        
        return {"FINISHED"}