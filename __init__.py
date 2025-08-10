import bpy
import bpy_extras
import bmesh
from bpy_extras import bmesh_utils
from importlib import reload
import math
import mathutils
import bl_math
import sys
import subprocess

bl_info = {
    "name": "SloppyUV",
    "description": "An assortment of simplish UV tools",
    "author": "Kim Falck",
    "version": (0, 0, 1),
    "blender": (4, 2, 0),
    "location": "UV > Sidebar",
    "wiki_url": "...",
    "category": "User Interface",
}

from SloppyUV.SloppySeamGeneration import SloppySeamGen, SloppySeamGenVis, SloppyIslandBoundaryConnect # type: ignore
from SloppyUV.SloppySortIndexByDist import SortVertByDist, SortEdgeByDist, SortFaceByDist # type: ignore
from SloppyUV.SloppyBakeProcAttr import SloppyProcAttrBake, SloppyBlurAttribute # type: ignore
from SloppyUV.SloppyProcedurals import SloppyQuadUVUnfold, SloppyUVToMesh, RedoUVEdgeLength, SloppyBasicUVUnfold,  SloppyBoundaryFirstUVUnfold # type: ignore


class SloppyProperties(bpy.types.PropertyGroup):
    
    iP = bpy.props.IntProperty
    fP = bpy.props.FloatProperty
    fvP = bpy.props.FloatVectorProperty
    bP = bpy.props.BoolProperty
    eP = bpy.props.EnumProperty
    bvP = bpy.props.BoolVectorProperty
    sP = bpy.props.StringProperty
    
    # region Common Functions
    def viewco(self, vec):
        '''Fetches (TODO: active or) any View3D region and space and outputs the input vector mapped
        to that 3D view.'''
        view3d_region = None
        view3d = None
        for ar in bpy.context.screen.areas:
            if ar.type == 'VIEW_3D':
                view3d_region = ar.regions[0]
                for spc in ar.spaces:
                    if spc.type == 'VIEW_3D':
                        view3d = spc.region_3d
                        break
        outco = bpy_extras.view3d_utils.location_3d_to_region_2d(view3d_region, view3d, vec)
        return outco
    
    def get_view3d_region_and_space(self):
        '''Fetches (TODO: active or) any View3D region and space.'''
        view3d_region = None
        view3d = None
        for ar in bpy.context.screen.areas:
            if ar.type == 'VIEW_3D':
                view3d_region = ar.regions[0]
                for spc in ar.spaces:
                    if spc.type == 'VIEW_3D':
                        view3d = spc.region_3d
                        break
        return view3d_region, view3d
    
    def set_view3d_lookat(self, position):
        '''Fetches (TODO: active or) any View3D region and space.'''
        view3d_region = None
        view3d = None
        for ar in bpy.context.screen.areas:
            if ar.type == 'VIEW_3D':
                view3d_region = ar.regions[0]
                for spc in ar.spaces:
                    if spc.type == 'VIEW_3D':
                        view3d = spc.region_3d
                        view3d.view_location = position
    
    def view3d_position(self, matrix):
        """ From 4x4 matrix, calculate view3d location """
        t = (matrix[0][3], matrix[1][3], matrix[2][3])
        r = (
        (matrix[0][0], matrix[0][1], matrix[0][2]),
        (matrix[1][0], matrix[1][1], matrix[1][2]),
        (matrix[2][0], matrix[2][1], matrix[2][2])
        )
        rp = (
        (-r[0][0], -r[1][0], -r[2][0]),
        (-r[0][1], -r[1][1], -r[2][1]),
        (-r[0][2], -r[1][2], -r[2][2])
        )
        output = (
        rp[0][0] * t[0] + rp[0][1] * t[1] + rp[0][2] * t[2],
        rp[1][0] * t[0] + rp[1][1] * t[1] + rp[1][2] * t[2],
        rp[2][0] * t[0] + rp[2][1] * t[1] + rp[2][2] * t[2],
        )
        return output

    def view3d_details(self, view):
        """ Return position, rotation data about a given view for the first space attached to it """
        look_at = view.view_location
        matrix = view.view_matrix
        view_pos = self.view3d_position(matrix)
        rotation = view.view_rotation
        return look_at, matrix, view_pos, rotation
    
    def project_point_on_plane_normal(self, plane_center, plane_normal, point_co):
        """ Returns a 2D vector of a point projected onto a plane given its center and normal """
        plane_up = plane_normal.orthogonal()
        right_rot_mat = mathutils.Matrix.Rotation(math.radians(90), 4, plane_normal)
        plane_right = plane_up.copy()
        plane_right.rotate(right_rot_mat)
        V = point_co - plane_center
        pt_x = V.dot(plane_right)
        pt_y = V.dot(plane_up)
        out_vec = mathutils.Vector((pt_x, pt_y))
        return out_vec
    
    def plane_axes_from_normal(self, plane_normal):
        """ Returns a 2D vector of a point projected onto a plane given its center and normal """
        plane_up = plane_normal.orthogonal()
        right_rot_mat = mathutils.Matrix.Rotation(math.radians(90), 4, plane_normal)
        plane_right = plane_up.copy()
        plane_right.rotate(right_rot_mat)
        return plane_right, plane_up
    
    def project_point_on_plane_axes(self, plane_center, plane_x, plane_y, point_co):
        """ Returns a 2D vector of a point projected onto a plane given its two axes """
        V = point_co - plane_center
        pt_x = V.dot(plane_x)
        pt_y = V.dot(plane_y)
        out_vec = mathutils.Vector((pt_x, pt_y))
        return out_vec
    
    def get_dir(self, co_source, co_destination):
        vec = co_destination - co_source
        dir = vec.normalized()
        return dir
    
    def sort_key_second_item_in_list_of_lists_item(self, list_of_lists_item):
        return list_of_lists_item[1]
    
    def sort_list_of_lists_by_list_length(self, list_in_a_list_of_lists):
        return len(list_in_a_list_of_lists)
    
    def sort_list_of_lists_by_first_item_length(self, list_in_a_list_of_lists):
        return len(list_in_a_list_of_lists[0])
    
    def sort_verts_by_nearest_most_opposite(self, two_vert_lst):
        distance = math.dist(two_vert_lst[0].co, two_vert_lst[1].co)
        pco_dot_vco = -two_vert_lst[0].co.dot(two_vert_lst[1].co)
        pno_dot_vno = two_vert_lst[0].normal.dot(two_vert_lst[1].normal)
        # return(pco_dot_vco + (pno_dot_vno * two_vert_lst[2]))
        return(pco_dot_vco)
    
    def get_verts_filtered_loops_in_island_using_perp_edge(self, v, e, lf, isl):
        out_loops = []
        fs = [f for f in e.link_faces if f in isl]
        fsls = []
        for fsf in fs:
            for fsfl in fsf.loops:
                if fsfl in lf:
                    fsls.append(fsfl)
        for vl in v.link_loops:
            if vl in fsls:
                out_loops.append(vl)
        return(out_loops)
    
    def generate_uv_loop_connection_map(self, input_loops, uv_layer):
        connected_loops = []
        for il in input_loops:
            conls = []
            for ill in input_loops:
                ld = math.dist(il[uv_layer].uv, ill[uv_layer].uv)
                if ld == 0.0:
                    conls.append(ill)
            connected_loops.append(conls)
        return(connected_loops)
    
    def find_contiguous_seams(self, bm, select_seams_by_index = False, select_index_start = 0, select_index_end = 0, sort_by_length = True, verbose = False):
        """ Returns lists (seams, face corners to the right, face corners to the left, boolean that is True if seam is open-ended (that is if one end does not end in a seam or boundary edge), boolean that is True if seam is isolated (that is if it never touches another seam)) of contiguous seams (which way is contiguous at crossroad seams will be decided according to which next seam is most well-aligned with the last) and two lists of face corners belonging on either side of the seam """
        seams = []
        loops_l = []
        loops_r = []
        seams_open_ended = []
        seams_isolated = []
        seams_meet_boundary = []
        seam_edges_total = []
        seam_edges_todo = []
        seam_edges_done = []

        for edge in bm.edges:
            if edge.seam == True:
                seam_edges_total.append(edge)
                seam_edges_todo.append(edge)
        
        while (len(seam_edges_total) - len(seam_edges_done)) > 0:
            this_seam = []
            these_loops_l = []
            these_loops_r = []
            open_end = False
            isolated = False
            ever_met_other_seam = False
            ever_met_boundary_edge = False
            open_end_count = 0
            init_edge = seam_edges_todo[0]
            init_edge_normal = self.calc_edge_avg_normal(init_edge, True)
            r_rot_mat = mathutils.Matrix.Rotation(math.radians(90), 4, init_edge_normal)
            init_edge_dir = self.get_dir(init_edge.verts[0].co, init_edge.verts[1].co)
            this_r = init_edge_dir.copy()
            this_r.rotate(r_rot_mat)
            
            for iev in init_edge.verts:
                for ievl in iev.link_loops:
                    ievl_tan = ievl.calc_tangent()
                    # ievl_tan += self.get_dir(self.calc_edge_center(init_edge), ievl.face.calc_center_bounds())
                    ievl_tan_nor = ievl_tan.normalized()
                    ie_tan_dot = ievl_tan.dot(this_r)
                    if ie_tan_dot > 0.0:
                        if ievl not in these_loops_r and ievl not in these_loops_l:
                            these_loops_r.append(ievl)
                    if ie_tan_dot <= 0.0:
                        if ievl not in these_loops_l and ievl not in these_loops_r:
                            these_loops_l.append(ievl)

            init_sub_bundle_a = [init_edge.verts[0], 0, -1, init_edge, -init_edge_dir]
            init_sub_bundle_b = [init_edge.verts[1], 0, 1, init_edge, init_edge_dir]
            next_bundle = [init_sub_bundle_a, init_sub_bundle_b]

            this_seam.append([init_edge, 0])

            if init_edge not in seam_edges_done:
                seam_edges_done.append(init_edge)
            if init_edge in seam_edges_todo:
                seam_edges_todo.remove(init_edge)

            while len(next_bundle) > 0:
                this_bundle = next_bundle.copy()
                next_bundle.clear()
                for sb in this_bundle:
                    next_edge = None
                    seam_cand = []
                    met_other_seam = False
                    met_boundary_edge = False

                    for sbe in sb[0].link_edges:
                        if sbe != sb[3]:
                            if len(sbe.link_faces) < 2:
                                met_boundary_edge = True
                                ever_met_boundary_edge = True
                            if sbe.seam == True:
                                met_other_seam = True
                                if sbe not in seam_edges_done:
                                    sbe_dir = self.get_dir(sb[0].co, sbe.other_vert(sb[0]).co)
                                    sbe_dot = sbe_dir.dot(sb[4])
                                    seam_cand.append([sbe, sbe_dot])
                                else:
                                    ever_met_other_seam = True
                    if len(seam_cand) > 0:
                        if len(seam_cand) > 1:
                            ever_met_other_seam = True
                            seam_cand.sort(key=self.sort_key_second_item_in_list_of_lists_item, reverse=True)
                        next_edge = seam_cand[0][0]
                    
                    if next_edge:
                        sb_ov = next_edge.other_vert(sb[0])
                        sb_dir = self.get_dir(sb[0].co, sb_ov.co)
                        if sb[2] == -1:
                            sb_dir *= -1
                        sb_normal = self.calc_edge_avg_normal(next_edge, True)
                        sb_r_rot_mat = mathutils.Matrix.Rotation(math.radians(90), 4, sb_normal)
                        sb_r = sb_dir.copy()
                        sb_r.rotate(sb_r_rot_mat)

                        for nev in next_edge.verts:
                            for nevl in nev.link_loops:
                                nevl_tan = nevl.calc_tangent()
                                # nevl_tan += self.get_dir(self.calc_edge_center(next_edge), nevl.face.calc_center_bounds())
                                nevl_tan_nor = nevl_tan.normalized()
                                ne_tan_dot = nevl_tan.dot(sb_r)
                                if ne_tan_dot > 0.0:
                                    if nevl not in these_loops_r and nevl not in these_loops_l:
                                        these_loops_r.append(nevl)
                                if ne_tan_dot <= 0.0:
                                    if nevl not in these_loops_l and nevl not in these_loops_r:
                                        these_loops_l.append(nevl)

                        sb_sort_val = sb[1] + sb[2]

                        ne_bundle = [sb_ov, sb_sort_val, sb[2], next_edge, sb_dir]
                        next_bundle.append(ne_bundle)

                        this_seam.append([next_edge, sb_sort_val])

                        if next_edge not in seam_edges_done:
                            seam_edges_done.append(next_edge)
                        if next_edge in seam_edges_todo:
                            seam_edges_todo.remove(next_edge)
                    else:
                        if met_other_seam == False and met_boundary_edge == False:
                            open_end = True
                            open_end_count += 1
            this_seam.sort(key=self.sort_key_second_item_in_list_of_lists_item)

            this_seam_clean = [i[0] for i in this_seam]
            
            if open_end_count >= 2:
                if ever_met_other_seam == False:
                    isolated = True

            seams.append(this_seam_clean)
            loops_l.append(these_loops_l)
            loops_r.append(these_loops_r)
            seams_open_ended.append(open_end)
            seams_isolated.append(isolated)
            seams_meet_boundary.append(ever_met_boundary_edge)
            

            num_remain = len(seam_edges_total) - len(seam_edges_done)
            if verbose:
                print('Number of seams total:', len(seam_edges_total), '\nNumber of seams done:', len(seam_edges_done), '\nNumber of seams remaining:', num_remain)

        sorting_array = []
        for smi, sm in enumerate(seams):
            sorting_sub = [seams[smi], loops_r[smi], loops_l[smi], seams_open_ended[smi], seams_isolated[smi], seams_meet_boundary[smi]]
            sorting_array.append(sorting_sub)

        if sort_by_length == True:
            sorting_array.sort(key=self.sort_list_of_lists_by_first_item_length, reverse=True)
            for sai, sa in enumerate(sorting_array):
                seams[sai] = sa[0]
                loops_r[sai] = sa[1]
                loops_l[sai] = sa[2]
                seams_open_ended[sai] = sa[3]
                seams_isolated[sai] = sa[4]
                seams_meet_boundary[sai] = sa[5]
        
        if select_seams_by_index:
            i_s = int(bl_math.clamp(float(select_index_start), 0, float(len(seams)-1)))
            i_e = int(bl_math.clamp(float(select_index_end), 0, float(len(seams)-1)))

            if i_s == i_e:
                for sme in seams[i_s]:
                    sme.select_set(True)
            elif i_e < i_s:
                for i in range(i_e, i_s + 1):
                    for sme in seams[i]:
                        sme.select_set(True)
            else:
                for i in range(i_s, i_e + 1):
                    for sme in seams[i]:
                        sme.select_set(True)
            bpy.context.active_object.data.update()

        if verbose:
            sorted_by_length_string = ":"
            if sort_by_length == True:
                sorted_by_length_string = "- Sorted by length:"
            print('Number of contiguous seams:', len(seams), sorted_by_length_string)
            for smi, sm in enumerate(seams):
                smis = []
                for sme in sm:
                    smis.append(sme.index)
                open_ended = "(closed"
                if seams_open_ended[smi] == True:
                    open_ended = "(open-ended"
                if seams_isolated[smi] == True:
                    open_ended = "(isolated"
                meets_boundary_string = open_ended + ")"
                if seams_meet_boundary[smi] == True:
                    meets_boundary_string = open_ended + " - reaches boundary)"
                print('Seam',(smi), meets_boundary_string,'- length',len(smis),'=',smis)

        return seams, loops_r, loops_l, seams_open_ended, seams_isolated, seams_meet_boundary
    
    def find_island_boundaries_and_their_loops(self, bm, island, sort_by_length = True, verbose = False):
        """ Returns lists (loop chains, vertex chains, edge chains, face chains) of contiguous boundaries of input island/collection of faces """
        island_loops = []
        boundary_faces = []
        boundary_edges = []
        boundary_verts = []
        boundary_loops = []
        loops_todo = []
        loops_done = []

        for iface in island:
            for ifl in iface.loops:
                if ifl not in island_loops:
                    island_loops.append(ifl)

        for iface in island:
            for ife in iface.edges:
                if ife not in boundary_edges:
                    is_boundary = False
                    if len(ife.link_faces) < 2:
                        is_boundary = True
                    else:
                        faces_in_island = 0
                        for ifef in ife.link_faces:
                            if ifef in island:
                                faces_in_island += 1
                        if faces_in_island <= 1:
                            is_boundary = True
                    if ife.seam == True:
                        is_boundary = True
                    if is_boundary:
                        if iface not in boundary_faces:
                            boundary_faces.append(iface)
                        if ife not in boundary_edges:
                            boundary_edges.append(ife)
                            for ifev in ife.verts:
                                if ifev not in boundary_verts:
                                    boundary_verts.append(ifev)
                                for ifevl in ifev.link_loops:
                                    if ifevl in island_loops:
                                        if ifevl not in boundary_loops:
                                            boundary_loops.append(ifevl)
                                            loops_todo.append(ifevl)
                                    if ifevl.face in island:
                                        if ifevl.face not in boundary_faces:
                                            boundary_faces.append(ifevl.face)
        if verbose == True:
            print('Boundary loops:', len(boundary_loops))
            print('Boundary vertices:', len(boundary_verts))
            print('Boundary edges:', len(boundary_edges))
            print('Boundary faces:', len(boundary_faces))
            print('Island loops:', len(island_loops))

        loop_chains = []
        loop_chains_eton = []
        vert_chains = []
        edge_chains = []
        face_chains = []

        while len(loops_todo) > 0:
            # print('Verts done:', len(verts_done), 'Verts to do:', len(verts_todo), 'Vert chains:', len(vert_chains))
            loop_chain = []
            loop_chain_eton = []
            vert_chain = []
            edge_chain = []
            face_chain = []

            next_loop = loops_todo[0]
            better_loop_found = False
            for ltd in loops_todo:
                if better_loop_found == False:
                    if len(ltd.vert.link_faces) <= 2:
                        next_loop = ltd
                        better_loop_found = True
                    else:
                        faces_in_island = 0
                        for ltdf in ltd.vert.link_faces:
                            if ltdf in island:
                                faces_in_island += 1
                        if faces_in_island <= 1:
                            next_loop = ltd
                            better_loop_found = True
            first_loop = next_loop
            vert_chain.append(next_loop.vert)
            face_chain.append(next_loop.face)

            edge_between = None
            
            lc_iter = 0
            while next_loop != None:
                full_circle = None
                this_loop = next_loop
                next_loop = None
                loop_chain.append(this_loop)
                if this_loop not in loops_done:
                    loops_done.append(this_loop)
                if this_loop in loops_todo:
                    loops_todo.remove(this_loop)

                if not next_loop:
                    for tle in this_loop.vert.link_edges:
                        if next_loop == None:
                            if tle in boundary_edges:
                                for tleovl in tle.other_vert(this_loop.vert).link_loops:
                                    if tleovl not in loops_done:
                                        if tleovl.face == this_loop.face:
                                            next_loop = tleovl
                                            edge_between = tle
                                    else:
                                        if lc_iter > 1:
                                            if tleovl.face == this_loop.face:
                                                if tleovl == first_loop:
                                                    full_circle = True
                                                    edge_between = tle
                    if next_loop:
                        edge_chain.append(edge_between)
                        loop_chain_eton.append(edge_between)

                if not next_loop:
                    for tlvl in this_loop.vert.link_loops:
                        if not next_loop:
                            share_face = False
                            for tlvlfe in tlvl.face.edges:
                                if tlvlfe not in boundary_edges:
                                    for tlvlfef in tlvlfe.link_faces:
                                        if tlvlfef == this_loop.face:
                                            share_face = True
                            if share_face == True:
                                if tlvl not in loops_done:
                                    next_loop = tlvl
                                else:
                                    if lc_iter > 1:
                                        if tlvl == first_loop:
                                            full_circle = True
                                            edge_between = None
                    if next_loop:
                        face_chain.append(this_loop.face)
                        loop_chain_eton.append(None)

                if full_circle == True:
                    loop_chain_eton.append(edge_between)
                
                lc_iter +=1
            loop_chains.append(loop_chain)
            loop_chains_eton.append(loop_chain_eton)
            vert_chains.append(vert_chain)
            edge_chains.append(edge_chain)
            face_chains.append(face_chain)

            if sort_by_length == True:
                sorting_array = []
                for lci, lc in enumerate(loop_chains):
                    sorting_sub = [lc, loop_chains_eton[lci], vert_chains[lci], edge_chains[lci], face_chains[lci]]
                    sorting_array.append(sorting_sub)

                if sort_by_length == True:
                    sorting_array.sort(key=self.sort_list_of_lists_by_first_item_length, reverse=True)
                    for sai, sa in enumerate(sorting_array):
                        loop_chains[sai] = sa[0]
                        loop_chains_eton[sai] = sa[1]
                        vert_chains[sai] = sa[2]
                        edge_chains[sai] = sa[3]
                        face_chains[sai] = sa[4]
        
        if verbose == True:
            print('Number of loop chains:', len(loop_chains))

        return loop_chains, loop_chains_eton, vert_chains, edge_chains, face_chains
    
    def align_view3d_against_normal(self, normal, view3d):
        rot_quat = -normal.to_track_quat('Z', 'Y')
        view3d_region = None
        view3d = None
        for ar in bpy.context.screen.areas:
            if ar.type == 'VIEW_3D':
                view3d_region = ar.regions[0]
                for spc in ar.spaces:
                    if spc.type == 'VIEW_3D':
                        view3d = spc.region_3d
                        view3d.view_rotation = rot_quat

    def line_up_exe(times = 1):
        str = ''
        CURSOR_UP = '\033[F'
        for i in range(times):
            str += CURSOR_UP
        print(str)

    def line_del_exe(times = 1):
        str = ''
        ERASE_LINE = '\033[K'
        for i in range(times):
            str += ERASE_LINE
        print(str)

    def line_up(times = 1):
        str = ''
        CURSOR_UP = '\033[F'
        for i in range(times):
            str += CURSOR_UP
        return(str)

    def line_del(times = 1):
        str = ''
        ERASE_LINE = '\033[K'
        for i in range(times):
            str += ERASE_LINE
        return(str)

    def viewvec(self, loc):
        '''Fetch vectpr pointing from 3D View origin to point on screen'''
        view3d_region = None
        view3d = None
        for ar in bpy.context.screen.areas:
            if ar.type == 'VIEW_3D':
                view3d_region = ar.regions[0]
                for spc in ar.spaces:
                    if spc.type == 'VIEW_3D':
                        view3d = spc.region_3d
                        break
        outvec = bpy_extras.view3d_utils.region_2d_to_vector_3d(view3d_region, view3d, loc)
        
        return outvec.normalized()
    
    def find_or_add_attribute(self, name, attr_type, attr_domain):
        dat = bpy.context.object.data
        attribute = dat.attributes[0]
        get_i = dat.attributes.find(name)
        if get_i == -1:
            attribute = dat.attributes.new(name=name, type=attr_type, domain=attr_domain)
        else:
            try:
                attribute = dat.attributes[get_i]
            except:
                bpy.context.active_object.data.update()
                attribute = dat.attributes.new(name=name, type=attr_type, domain=attr_domain)
        return attribute
    
    def find_or_add_attribute_other_obj(self, name, attr_type, attr_domain, other_obj):
        dat = other_obj.data
        attribute = dat.attributes[0]
        get_i = dat.attributes.find(name)
        if get_i == -1:
            attribute = dat.attributes.new(name=name, type=attr_type, domain=attr_domain)
        else:
            try:
                attribute = dat.attributes[get_i]
            except:
                dat.update()
                attribute = dat.attributes.new(name=name, type=attr_type, domain=attr_domain)
        return attribute
    
    def get_attribute_layer(self, name, attr_type, attr_domain, bmi):
            layer = None
            if attr_domain == "FACE":
                if attr_type == "INT":
                    layer = bmi.faces.layers.int.get(name)
                if attr_type == "FLOAT_VECTOR":
                    layer = bmi.faces.layers.float_vector.get(name)
                if attr_type == "FLOAT":
                    layer = bmi.faces.layers.float.get(name)
            if attr_domain == "CORNER":
                if attr_type == "INT":
                    layer = bmi.loops.layers.int.get(name)
                if attr_type == "FLOAT":
                    layer = bmi.loops.layers.float.get(name)
                if attr_type == "FLOAT_COLOR":
                    layer = bmi.loops.layers.float_color.get(name)
                if attr_type == "FLOAT_VECTOR":
                    layer = bmi.loops.layers.float_vector.get(name)
            if attr_domain == "POINT":
                if attr_type == "FLOAT_COLOR":
                    layer = bmi.verts.layers.float_color.get(name)
                if attr_type == "FLOAT_VECTOR":
                    layer = bmi.verts.layers.float_vector.get(name)
                if attr_type == "FLOAT":
                    layer = bmi.verts.layers.float.get(name)
                if attr_type == "INT":
                    layer = bmi.verts.layers.int.get(name)
            if attr_domain == "EDGE":
                if attr_type == "FLOAT":
                    layer = bmi.edges.layers.float.get(name)
                if attr_type == "INT":
                    layer = bmi.edges.layers.int.get(name)
            return layer
    
    def get_dict_layer(self, name, attribute_dict):
        layer = None
        for dict in attribute_dict:
            if dict["name"] == name:
                layer = dict["layer"]
                return layer

    def find_near_parallels(self, e):
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
    
    def select_by_index_bm(self, domain, choose_domain, uvs, i_start, i_end, sel_range, extend, output_only, in_bm):
        '''Select mesh element by index or range of indices. 
        (Can also be used to only output a list of elements.)'''
        uv_layer = in_bm.loops.layers.uv.verify()
        out_arr = []
        indices = [i_start]
        if sel_range:
            indices = range(i_start, i_end + 1)
        dom = None
        dom_i = 0
        modes = bpy.context.tool_settings.mesh_select_mode[:]
        if choose_domain:
            dom_i = domain
            if domain == 0:
                dom = in_bm.verts
            if domain == 1:
                dom = in_bm.edges
            if domain == 2:
                dom = in_bm.faces
        else:
            if modes[0]:
                dom = in_bm.verts
                dom_i = 0
            if modes[1]:
                dom = in_bm.edges
                dom_i = 1
            if modes[2]:
                dom = in_bm.faces
                dom_i = 2
        if output_only:
            for el in dom:
                if el.index in indices:
                    out_arr.append(el)
        else:
            for el in dom:
                if el.index in indices:
                    out_arr.append(el)
                if uvs:
                    if dom_i == 2:
                        for ell in el.loops:
                            if el.index in indices:
                                ell[uv_layer].select = True
                            else:
                                if extend == False:
                                    ell[uv_layer].select = False
                    else:
                        for ell in el.link_loops:
                            if el.index in indices:
                                ell[uv_layer].select = True
                            else:
                                if extend == False:
                                    ell[uv_layer].select = False
                else:
                    if el.index in indices:
                        el.select = True
                    else:
                        if extend == False:
                            el.select = False
        bpy.context.active_object.data.update()
        return out_arr
    
    def shift_select_by_index_bm(self, domain, choose_domain, shift_amount, shift_by_len, shift_by_len_neg, extend, both_directions, output_only, in_bm):
        '''Shift mesh element selection by their indices. 
        (Can also be used to only output a list of elements.)'''
        
        out_arr = []
        dom = None
        modes = bpy.context.tool_settings.mesh_select_mode[:]
        if choose_domain:
            if domain == 0:
                dom = in_bm.verts
            if domain == 1:
                dom = in_bm.edges
            if domain == 2:
                dom = in_bm.faces
        else:
            if modes[0]:
                dom = in_bm.verts
            if modes[1]:
                dom = in_bm.edges
            if modes[2]:
                dom = in_bm.faces

        dom.ensure_lookup_table()

        indices = [i.index for i in dom if i.select == True]

        if extend == False:
            for i in indices:
                dom[i].select = False

        if len(indices) == 0:
            return out_arr
        else:
            not_indices = indices.copy()
            for i, index in enumerate(not_indices):
                mod_i = index
                mod_i_n = index
                if shift_by_len:
                    arr_len = len(indices)
                    if shift_by_len_neg:
                        arr_len *= -1
                    mod_i += arr_len
                    mod_i_n -= arr_len
                else:
                    mod_i += shift_amount
                    mod_i_n -= shift_amount
                new_i = mod_i%len(dom)
                if both_directions:
                    new_i_n = mod_i_n%len(dom)
                    indices.append(new_i_n)
                indices[i] = new_i
        for i in indices:
            out_arr.append(dom[i])
            if output_only == False:
                dom[i].select = True
        bpy.context.active_object.data.update()
        return out_arr
    
    def calc_average_normal(self, faces):
        '''Calculate average normal of a list of faces.'''
        norm = mathutils.Vector((0,0,0))
        for face in faces:
            norm += face.normal
        return norm.normalized()
    
    def uv_sort_x(self, loop_uvco_list):
        return(loop_uvco_list[1].x)
    
    def uv_sort_y(self, loop_uvco_list):
        return(loop_uvco_list[1].y)

    def calc_edge_center(self, edge):
        return(edge.verts[0].co.lerp(edge.verts[1].co, 0.5))
    
    def edge_corners_per_vert(self, edge, vert):
        corners = []
        all_corners = []
        for ef in edge.link_faces:
            for efl in ef.loops:
                all_corners.append(efl)
        for v in edge.verts:
            if v == vert:
                for vl in v.link_loops:
                    if vl in all_corners:
                        corners.append(vl)
        return corners
    
    def face_corner_per_vert(self, face, vert):
        for loop in face.loops:
            if loop.vert == vert:
                return loop
        return None
    
    def calc_edge_avg_normal(self, edge, use_faces = False):
        if use_faces:
            add_nor = mathutils.Vector((0,0,0))
            for face in edge.link_faces:
                add_nor += face.normal
            avg_nor = add_nor.normalized()
            return(avg_nor)
        else:
            add_nor = edge.verts[0].normal + edge.verts[1].normal
            avg_nor = add_nor.normalized()
            return(avg_nor)
    
    def calc_selected_edges_total_length(self):
        totlen = 0
        for obj in bpy.context.objects_in_mode:
            for edge in obj.edges:
                if edge.select == True:
                    parent_mesh = edge.id_data
                    v1 = parent_mesh.vertices[edge.vertices[0]]
                    v2 = parent_mesh.vertices[edge.vertices[1]]
                    totlen += math.dist(v1.co, v2.co)
        return(totlen)
    
    def calc_total_and_selected_surface_area(self, objects):
        total_area = 0
        selected_area = 0
        for obj in objects:
            for pol in obj.data.polygons:
                total_area += pol.area
                if pol.select == True:
                    selected_area += pol.area
        return total_area, selected_area
    
    def copy_float_to_clipboard(self, val_to_copy):
        cmd=f'echo {val_to_copy}|clip'
        return subprocess.check_call(cmd, shell=True)
    
    def remap_val(self, val, in_min, in_max, out_min, out_max):
        in_interval = in_max - in_min
        out_interval = out_max - out_min
        in_val = val - in_min
        in_fac = 0
        if in_interval > 0:
            in_fac = in_val / in_interval
        out_val = out_min + (in_fac * out_interval)
        if out_val > out_max:
            out_val = out_max
        if out_val < out_min:
            out_val = out_min
        return out_val

    def update_uvs_geo_axis(self, context):
        '''Update function for UV alignment 3D view axis, to ensure only one is selected.'''
        if self.update_axis_switch == False:
            self.update_axis_switch = True
            unchanged = []
            change_i = 0
            
            for i in range(2):
                if self.align_uv_geo_axis[i] == True:
                        if self.align_uv_geo_axis_prev[i] == True:
                            unchanged.append(i)
            
            for j in unchanged:
                self.align_uv_geo_axis[j] = False
            
            for k in range(2):
                if self.align_uv_geo_axis[k] == True:
                    self.align_uv_geo_axis_prev[k] = True
                if self.align_uv_geo_axis[k] == False:
                    self.align_uv_geo_axis_prev[k] = False
            self.update_axis_switch = False
        return None
    
    def update_uv_axis(self, context):
        '''Update function for UV alignment axis, to ensure only one is selected.'''
        if self.update_axis_switch == False:
            self.update_axis_switch = True
            unchanged = []
            change_i = 0
            
            for i in range(2):
                if self.align_uv_axis[i] == True:
                        if self.align_uv_axis_prev[i] == True:
                            unchanged.append(i)
            
            for j in unchanged:
                self.align_uv_axis[j] = False
            
            for k in range(2):
                if self.align_uv_axis[k] == True:
                    self.align_uv_axis_prev[k] = True
                if self.align_uv_axis[k] == False:
                    self.align_uv_axis_prev[k] = False
            self.update_axis_switch = False
        return None
    
    def update_uv_axis_scaling(self, context):
        '''Update function for UV scaling axis, to ensure only one is selected.'''
        if self.update_axis_switch == False:
            self.update_axis_switch = True
            unchanged = []
            change_i = 0
            
            for i in range(2):
                if self.scale_calc_uv_axis[i] == True:
                        if self.scale_calc_uv_axis_prev[i] == True:
                            unchanged.append(i)
            
            for j in unchanged:
                self.scale_calc_uv_axis[j] = False
            
            for k in range(2):
                if self.scale_calc_uv_axis[k] == True:
                    self.scale_calc_uv_axis_prev[k] = True
                if self.scale_calc_uv_axis[k] == False:
                    self.scale_calc_uv_axis_prev[k] = False
            self.update_axis_switch = False
            self.update_scalings(context)
        return None
    
    def update_swizzle_x_axis(self, context):
        '''Update function for 3D axis swizzle for X axis, to ensure only one is selected.'''
        if self.update_axis_switch == False:
            self.update_axis_switch = True
            unchanged = []
            change_i = 0
            
            for i in range(3):
                if self.scale_uv_x_geo_axis[i] == True:
                        if self.scale_uv_x_geo_axis_prev[i] == True:
                            unchanged.append(i)
            
            for j in unchanged:
                self.scale_uv_x_geo_axis[j] = False
            
            for k in range(3):
                if self.scale_uv_x_geo_axis[k] == True:
                    self.scale_uv_x_geo_axis_prev[k] = True
                if self.scale_uv_x_geo_axis[k] == False:
                    self.scale_uv_x_geo_axis_prev[k] = False
            self.update_axis_switch = False
            self.update_scalings(context)
        return None
    
    def update_swizzle_y_axis(self, context):
        '''Update function for 3D axis swizzle for Y axis, to ensure only one is selected.'''
        if self.update_axis_switch == False:
            self.update_axis_switch = True
            unchanged = []
            change_i = 0
            
            for i in range(3):
                if self.scale_uv_y_geo_axis[i] == True:
                        if self.scale_uv_y_geo_axis_prev[i] == True:
                            unchanged.append(i)
            
            for j in unchanged:
                self.scale_uv_y_geo_axis[j] = False
            
            for k in range(3):
                if self.scale_uv_y_geo_axis[k] == True:
                    self.scale_uv_y_geo_axis_prev[k] = True
                if self.scale_uv_y_geo_axis[k] == False:
                    self.scale_uv_y_geo_axis_prev[k] = False
            self.update_axis_switch = False
            self.update_scalings(context)
        return None
    
    def update_scale_island(self, context):
        '''Update function for seleecting island mode, to ensure only one mode is active.'''
        if self.update_scale_switch == False:
            self.update_scale_switch = True
            if self.scale_calc_island == True:
                self.scale_calc_selected = False
                self.scale_calc_edge_length = False
        self.update_scale_switch = False
        self.update_scalings(context)
        return None
    
    def update_scale_selected(self, context):
        '''Update function for seleecting selected UVs mode, to ensure only one mode is active.'''
        if self.update_scale_switch == False:
            self.update_scale_switch = True
            if self.scale_calc_selected == True:
                self.scale_calc_island = False
                self.scale_calc_edge_length = False
        self.update_scale_switch = False
        self.update_scalings(context)
        return None
    
    def update_scale_length(self, context):
        '''Update function for seleecting edge length mode, to ensure only one mode is active.'''
        if self.update_scale_switch == False:
            self.update_scale_switch = True
            if self.scale_calc_edge_length == True:
                self.scale_calc_island = False
                self.scale_calc_selected = False
                self.update_axis_switch = True
                self.scale_uv_x_geo_axis[0] = True
                self.scale_uv_x_geo_axis[1] = False
                self.scale_uv_x_geo_axis[2] = False
                self.scale_uv_y_geo_axis[0] = False
                self.scale_uv_y_geo_axis[1] = True
                self.scale_uv_y_geo_axis[2] = False
                self.update_axis_switch = False
        self.update_scale_switch = False
        self.update_scalings(context)
        return None
    
    def update_scale_view(self, context):
        '''update function for seleecting viewport transform space, to ensure only one
        transform space is active.'''
        if self.update_scale_switch == False:
            self.update_scale_switch = True
            if self.scale_geo_space_view == True:
                self.scale_geo_space_global = False
                self.scale_geo_space_object = False
                self.update_axis_switch = True
                self.scale_uv_x_geo_axis[0] = True
                self.scale_uv_x_geo_axis[1] = False
                self.scale_uv_x_geo_axis[2] = False
                self.scale_uv_y_geo_axis[0] = False
                self.scale_uv_y_geo_axis[1] = True
                self.scale_uv_y_geo_axis[2] = False
                self.update_axis_switch = False
        self.update_scale_switch = False
        self.update_scalings()
        return None
    
    def update_scale_global(self, context):
        '''Update function for seleecting global transform space, to ensure only one
        transform space is active.'''
        if self.update_scale_switch == False:
            self.update_scale_switch = True
            if self.scale_geo_space_global  == True:
                self.scale_geo_space_view = False
                self.scale_geo_space_object = False
        self.update_scale_switch = False
        self.update_scalings(context)
        return None
    
    def update_scale_object(self, context):
        '''Update function for seleecting object transform space, to ensure only one
        transform space is active.'''
        if self.update_scale_switch == False:
            self.update_scale_switch = True
            if self.scale_geo_space_object  == True:
                self.scale_geo_space_view = False
                self.scale_geo_space_global = False
        self.update_scale_switch = False
        self.update_scalings(context)
        return None
    
    def update_scale_uv_to_bounds(self, context):
        self.update_scalings(context)
        return None
    
    def update_scale_uv_real_world(self, context):
        self.update_scalings(context)
        return None
    
    def update_scalings(self, context):
        '''Update calculation results.'''
        some_as = []
        some_bs = []
        for a in self.scale_uv_x_geo_axis:
            some_as.append(a)
        for b in self.scale_uv_y_geo_axis:
            some_bs.append(b)
        xi = some_as.index(True)
        yi = some_bs.index(True)
        if self.scale_calc_edge_length == True:
            self.scale_geo_size[0] = self.scale_geo_edge_lengths_real[xi]
            self.scale_geo_size[1] = self.scale_geo_edge_lengths_real[yi]
        if self.scale_calc_selected == True:
            self.scale_geo_size[0] = self.scale_geo_selected_real[xi]
            self.scale_geo_size[1] = self.scale_geo_selected_real[yi]
        if self.scale_calc_island == True:
            self.scale_geo_size[0] = self.scale_geo_size_real[xi]
            self.scale_geo_size[1] = self.scale_geo_size_real[yi]
        geo_size_max = max(self.scale_geo_size[0], self.scale_geo_size[1])
        geo_size_min = min(self.scale_geo_size[0], self.scale_geo_size[1])
        
        if self.scale_calc_edge_length == True:
            uv_size_trans_x = (self.scale_uv_size_current[0] / self.scale_uv_edge_lengths[0]) * self.scale_geo_size[0]
            uv_size_trans_y = (self.scale_uv_size_current[1] / self.scale_uv_edge_lengths[1]) * self.scale_geo_size[1]
            if self.scale_uv_real_world == True:
                self.scale_uv_size[0] = uv_size_trans_x
                self.scale_uv_size[1] = uv_size_trans_y
            else:
                self.scale_uv_size[0] = uv_size_trans_x / geo_size_max
                self.scale_uv_size[1] = uv_size_trans_y / geo_size_max
        uv_size_max = max(self.scale_uv_size[0], self.scale_uv_size[1])
        self.scale_uv_aspect[0] = self.scale_uv_size[0] / uv_size_max
        self.scale_uv_aspect[1] = self.scale_uv_size[1] / uv_size_max
        if self.scale_uv_to_bounds == True:
            self.scale_uv_size[0] /= uv_size_max
            self.scale_uv_size[1] /= uv_size_max
        if self.scale_uv_real_world == True:
            self.scale_uv_geo_relative[0] = 1.0
            self.scale_uv_geo_relative[1] = 1.0
        else:
            self.scale_uv_geo_relative[0] = self.scale_uv_size[0] / self.scale_geo_size[0]
            self.scale_uv_geo_relative[1] = self.scale_uv_size[1] / self.scale_geo_size[1]
        

        return None
    
    def update_distsort_vars(self, context):
        SortVertByDist.bottom_up_mix = self.distsort_bottom_up
        SortVertByDist.z_only = self.distsort_z_only
        SortEdgeByDist.bottom_up_mix = self.distsort_bottom_up
        SortEdgeByDist.z_only = self.distsort_z_only
        SortFaceByDist.bottom_up_mix = self.distsort_bottom_up
        SortFaceByDist.z_only = self.distsort_z_only
        return None
    # endregion

    # region Properties Def
    distsort_bottom_up : bP (
        name = "From Bottom Up",
        description = "Weigh sorting of distance sort from bottom up",
        default = False,
        update = update_distsort_vars
        ) # type: ignore
    
    distsort_z_only : bP (
        name = "Only Sort From Z",
        description = "Only sort potential elements from Z",
        default = False,
        update = update_distsort_vars
        ) # type: ignore
    
    space_edges : bP(
        name = "Space Seam",
        description = "Space UV points along edges evenly",
        default = True
        ) # type: ignore
    
    relax_inner : bP(
        name = "Relax Inner",
        description = "Relax inner (non-boundary) UVs (UniV add-on required!)",
        default = True
        ) # type: ignore
    
    verbose : bP(
        name = "Verbose",
        description = "Prints progress info to System Console",
        default = False
        ) # type: ignore
    
    align_uv_geo_axis_prev : bvP(
        name = "3D Axis (Previous)",
        description = "Just for coding purposes",
        subtype = 'XYZ'
        ) # type: ignore
    
    align_uv_geo_axis : bvP(
        name = "View Axis",
        description = "Align UVs with 3D view in this axis",
        default = (True, False),
        subtype = 'XYZ',
        size = 2,
        update=update_uvs_geo_axis
        ) # type: ignore
    
    align_uv_axis_prev : bvP(
        name = "UV Axis (Previous)",
        description = "Just for coding purposes",
        default = (False, True),
        subtype = 'XYZ',
        size = 2,
        ) # type: ignore
    
    align_uv_axis : bvP(
        name = "UV Axis",
        description = "Align UVs with UV points on this axis",
        default = (False, True),
        subtype = 'XYZ',
        size = 2,
        update = update_uv_axis
        ) # type: ignore
    
    align_uvs_use_geo : bP(
        name = "Use geometry to align UVs",
        default = False,
        ) # type: ignore
    
    scale_geo_size : fvP(
        name = "3D Size",
        description = "Dimensions of geometry",
        default = (1.0, 1.0),
        precision = 4,
        subtype = 'XYZ',
        size = 2,
        ) # type: ignore
    
    scale_geo_selected_real : fvP(
        name = "Actual 3D Size",
        description = "Actual dimensions of geometry",
        default = (1.0, 1.0, 1.0),
        precision = 4,
        subtype = 'XYZ',
        size = 3,
        ) # type: ignore
    
    scale_geo_edge_lengths_real : fvP(
        name = "Actual 3D Size",
        description = "Actual dimensions of geometry",
        default = (1.0, 1.0),
        precision = 4,
        subtype = 'XYZ',
        size = 2,
        ) # type: ignore
    
    scale_geo_size_real : fvP(
        name = "Actual 3D Size",
        description = "Actual dimensions of geometry",
        default = (1.0, 1.0, 1.0),
        precision = 4,
        subtype = 'XYZ',
        size = 3,
        ) # type: ignore
    
    scale_uv_edge_lengths : fvP(
        name = "Selected UV Edges Length",
        default = (1.0, 1.0),
        subtype = 'XYZ',
        size = 2,
        ) # type: ignore
    
    scale_uv_selected_lengths : fvP(
        name = "Selected UVs Length",
        default = (1.0, 1.0),
        subtype = 'XYZ',
        size = 2,
        ) # type: ignore
    
    scale_uv_size : fvP(
        name = "UV Size",
        description = "UV dimensions in UV units",
        default = (1.0, 1.0),
        precision = 4,
        subtype = 'XYZ',
        size = 2,
        ) # type: ignore
    
    scale_uv_size_current : fvP(
        name = "Current UV Size",
        description = "UV dimensions in UV units",
        default = (1.0, 1.0),
        precision = 4,
        subtype = 'XYZ',
        size = 2,
        ) # type: ignore
    
    scale_uv_aspect : fvP(
        name = "UV Aspect Ratio",
        description = "Aspect ratio of UV island",
        default = (1.0, 1.0),
        min = 0.0,
        max = 1.0,
        soft_min = 0.0,
        soft_max = 1.0,
        precision = 4,
        subtype = 'XYZ',
        size = 2,
        ) # type: ignore
    
    scale_uv_geo_relative : fvP(
        name = "UV Relative to 3D",
        description = "UV size relative to geometry in 3D space",
        default = (1.0, 1.0),
        precision = 4,
        subtype = 'XYZ',
        size = 2,
        ) # type: ignore
    
    scale_calc_uv_axis : bvP(
        name = "UV Axis",
        description = "Axis to scale UVs in according to length of linked edges",
        default = (True, False),
        subtype = 'XYZ',
        size = 2,
        update = update_uv_axis_scaling
        ) # type: ignore
    
    scale_calc_uv_axis_prev : bvP(
        name = "UV Axis",
        description = "Axis to scale UVs in according to length of linked edges",
        default = (True, False),
        subtype = 'XYZ',
        size = 2,
        ) # type: ignore
    
    scale_calc_uv_corner : fvP(
        name = "Least Counds Corner",
        description = "Corner of UV bounds to use when applying scale",
        default = (True, False),
        subtype = 'XYZ',
        size = 2,
        ) # type: ignore
    
    scale_geo_space_view : bP(
        name = "View",
        description = "Calculate scaling in 3D View (automatically maps view XY to UV)", 
        default = True,
        update = update_scale_view
        ) # type: ignore
    
    scale_geo_space_global : bP(
        name = "Global",
        description = "Calculate scaling in global space", 
        default = False,
        update = update_scale_global
        ) # type: ignore
    
    scale_geo_space_object : bP(
        name = "Object",
        description = "Calculate scaling in object space", 
        default = False,
        update = update_scale_object
        ) # type: ignore
    
    scale_calc_island : bP(
        name = "Island",
        description = "Calculate X and Y from vertices linked to UV island boundary points", 
        default = False,
        update = update_scale_island
        ) # type: ignore
    
    scale_calc_selected : bP(
        name = "Select",
        description = "Calculate X and Y from vertices linked to boundary points of selected UVs", 
        default = False,
        update = update_scale_selected
        ) # type: ignore
    
    scale_calc_edge_length : bP(
        name = "Edges Length",
        description = "Calculate X or Y from sum of lengths of edges connected to selected UVs", 
        default = True,
        update = update_scale_length
        ) # type: ignore
    
    scale_uv_x_geo_axis : bvP(
        name = "X",
        description = "3D axis to map to U (UV X) (Ignored if Transform Space is 3D View)",
        default = (True, False, False),
        subtype = 'XYZ',
        size = 3,
        update = update_swizzle_x_axis
        ) # type: ignore
    
    scale_uv_x_geo_axis_prev : bvP(
        name = "X",
        description = "3D axis to map to U (UV X) (Ignored if Transform Space is 3D View)",
        default = (True, False, False),
        subtype = 'XYZ',
        size = 3,
        ) # type: ignore
    
    scale_uv_y_geo_axis : bvP(
        name = "Y",
        description = "3D axis to map to V (UV Y) (Ignored if Transform Space is 3D View)",
        default = (False, True, False),
        subtype = 'XYZ',
        size = 3,
        update = update_swizzle_y_axis
        ) # type: ignore
    
    scale_uv_y_geo_axis_prev : bvP(
        name = "Y",
        description = "3D axis to map to V (UV Y) (Ignored if Transform Space is 3D View)",
        default = (False, True, False),
        subtype = 'XYZ',
        size = 3,
        ) # type: ignore
    
    scale_uv_real_world : bP(
        name = "Real World Scale",
        description = "Scale Uvs to real world scale", 
        default = False,
        update = update_scale_uv_real_world
        ) # type: ignore
    
    scale_uv_to_bounds : bP(
        name = "Scale to Bounds",
        description = "Scale result to bounds", 
        default = True,
        update = update_scale_uv_to_bounds
        ) # type: ignore
    
    align_uv_scale_to_bounds : bP(
        name = "Scale to Bounds",
        description = "Scale result to bounds", 
        default = True
        ) # type: ignore
    
    align_uv_scale_to_bounds_by_greatest : bP(
        name = "Scale by Greatest Bound",
        description = "Scale uniformly to greatest bound", 
        default = True
        ) # type: ignore
    
    scale_uv_edge_direction : bP(
        name = "Scale Along Vector of Edges",
        description = "When in Edge Length Mode, scale UVs in directions of edges used for calculation", 
        default = False,
        ) # type: ignore
    
    scale_uv_edge_direction_x : fvP(
        name = "Edge Direction X",
        description = "Actual dimensions of geometry",
        default = (1.0, 0.0),
        precision = 4,
        subtype = 'XYZ',
        size = 2
        ) # type: ignore
    
    scale_uv_edge_direction_y : fvP(
        name = "Edge Direction Y",
        description = "Actual dimensions of geometry",
        default = (0.0, 1.0),
        precision = 4,
        subtype = 'XYZ',
        size = 2
        ) # type: ignore
    
    update_axis_switch : bP(
        name = "Axis Updating",
        default = False,
        ) # type: ignore
    
    update_scale_switch : bP(
        name = "Scale Updating",
        default = False,
        ) # type: ignore
    
    sloppy_error_msg_heading : sP(
        name = "Error Message Heading",
        default = "Error!",
        ) # type: ignore
    
    sloppy_error_msg : sP(
        name = "Error Message",
        default = "Unknonwn error occured.",
        ) # type: ignore
    
    delta_rotation_mode: eP(
        name = "Delta Rotation Mode",
        description = "",
        items = [
            ("A", "Euler", ""),
            ("B", "Quaternion", "")
            ],
        default="A"
        ) # type: ignore

    current_sort_axis : sP(
        name = "Current Sort Axis",
        default = "",
        ) # type: ignore

    # region SeamGen Properties
    seamgen_angle_factor : fP(
        name = "Concavity Influence",
        description = "Influence of concavity on seam generation",
        default = 1,
        min = 0,
        max = 1
        ) # type: ignore

    seamgen_avg_angle_factor : fP(
        name = "Average Concavity Influence",
        description = "Influence of averaged area concavity on seam generation",
        default = 0,
        min = 0,
        max = 1
        ) # type: ignore

    seamgen_ao_factor : fP(
        name = "AO Influence",
        description = "Influence of ambient occlusion on seam generation",
        default = 0,
        min = 0,
        max = 1
        ) # type: ignore

    seamgen_ed_factor : fP(
        name = "Edge Density Influence",
        description = "Influence of edge density on seam generation",
        default = 0,
        min = 0,
        max = 1
        ) # type: ignore

    seamgen_no_rounds : iP(
        name = "Iterations",
        description = "Number of initial iterations",
        default = 20,
        min = 1,
        max = 1000
        ) # type: ignore

    seamgen_no_retries : iP(
        name = "Max Retries",
        description = "Maximum number of retry iterations",
        default = 100,
        min = 1,
        max = 1000
        ) # type: ignore

    seamgen_angle_thresh_start : fP(
        name = "Min Angle Threshold",
        description = "Initial angular threshold",
        default = -50,
        min = -180,
        max = 180
        ) # type: ignore

    seamgen_angle_thresh_end : fP(
        name = "Max Angle Threshold",
        description = "Maximum angular threshold",
        default = -20,
        min = -180,
        max = 180
        ) # type: ignore

    seamgen_clear_seam : bP(
        name = "Clear Seams",
        description = "Clear any existing seams",
        default = True
        ) # type: ignore

    seamgen_unwrap : bP(
        name = "Unwrap",
        description = "Run auto-unwrap. (Mainly to check if islands turn out as desired.)",
        default = True
        ) # type: ignore
    # endregion

    # region Batch Render Properties
    br_output_path : sP(
        name = "Output Path",
        description = "Base path to batch render to",
        default = "",
        subtype="DIR_PATH"
        ) # type: ignore
    
    br_use_viewport : bP(
        name="Viewport Render",
        default=False
        ) # type: ignore
    #endregion

    # region CopyPaste Properties
    cpt_location : fvP(
        name = "Location",
        description = "",
        default = (0,0,0),
        subtype='XYZ_LENGTH',
        precision=10
        ) # type: ignore
    
    cpt_rotation : fvP(
        name = "Rotation",
        description = "",
        default = (1,0,0,0),
        subtype='QUATERNION',
        size=4,
        precision=10
        ) # type: ignore
    
    cpt_scale : fvP(
        name = "Scale",
        description = "",
        default = (1,1,1),
        subtype='XYZ',
        precision=10
        ) # type: ignore
    
    cpt_dimensions : fvP(
        name = "Dimensions",
        description = "",
        default = (0,0,0),
        subtype='XYZ_LENGTH',
        precision=10
        ) # type: ignore
    
    cpt_spline_length : fP(
        name = "Clipboard Spline Length",
        description = "",
        default = 0
        ) # type: ignore
    
    cpt_median : fvP(
        name = "Clipboard Median",
        description = "",
        default = (0,0,0),
        subtype='XYZ'
        ) # type: ignore
    
    ep_total_area : fP(
        name = "Total Surface Area",
        description = "",
        default = 0
        ) # type: ignore
    
    ep_selected_area : fP(
        name = "Selected Surface Area",
        description = "",
        default = 0
        ) # type: ignore
    
    ep_total_selected_edge_length : fP(
        name = "Selected Edges Length",
        description = "",
        default = 0
        ) # type: ignore
    #endregion


    # region ProcUV Properties
    procuv_unfold_mode : eP(
        name = "Mode",
        description = "Unfolding Mode",
        items = [
            ("A", "All Directions", "Work outwards from initial quad"),
            ("B", "X First", "Work outwards in (local) X first until no more faces are found, calculating edge lengths to use in next stage, then work outwards in Y. repeating first step when more faces are found in X"),
            ("C", "Y First", "Work outwards in (local) Y first until no more faces are found, calculating edge lengths to use in next stage, then work outwards in X. repeating first step when more faces are found in Y")
            ]
        ) # type: ignore

    procuv_initial_quad : eP(
        name = "Initial Quad(s)",
        description = "Which quad (for each UV island, if there are more than one,) to start from",
        items = [
            ("A", "Most Regular", "Start with quad (or virtual quad, formed by combining most regular face with adjacent triangle as calculated from face corner angles) that is flattest and/or that has normal most similar to island average normals"),
            ("B", "Selected", "First selected quad (or virtual quad, formed from two first triangles selected, or, if only 1 triangle is selected in island, formed by combining with adjacent triangle as calculated from face corner angles) in each island (if an island has no faces selected, initial quad falls back to Most Regular)")
            ]
        ) # type: ignore

    procuv_reg_flat_fac : fP(
        name = "Flatness Factor",
        description = "How much flatness influences initial quad selection",
        default = 1.0,
        min = 0.0,
        max = 2.0
        ) # type: ignore

    procuv_reg_norm_fac : fP(
        name = "Normal Factor",
        description = "How much island average normal influences initial quad selection",
        default = 1.0,
        min = 0.0,
        max = 2.0
        ) # type: ignore

    procuv_pre_calc_edge_lengths : bP(
        name = "Pre-Calculate Edge Lengths",
        description = "Calculate average edge lengths of columns and rows before setting UV coordinates",
        default = False
        ) # type: ignore

    procuv_offset_per_island : fvP(
        name = "Offset/Island",
        description = "Offset UVs per island",
        size = 2
        ) # type: ignore

    procuv_quant_avg_norm : bP(
        name = "Quantize Average Normal",
        description = "Quantize island's averaged normal to up, down, right, left",
        default = False
        ) # type: ignore

    procuv_only_move_loops_in_face : bP(
        name = "Only Edit Current Loops",
        description = "Avoid moving loops of other faces than current face when editing UV loops",
        default = False
        ) # type: ignore
    #endregion
    # endregion

#region Panel Classes
class SloppyUVPanel(bpy.types.Panel):
    bl_idname = "UV_PT_SloppyUVPanel"
    bl_label = "Sloppy UV"
    bl_space_type = "IMAGE_EDITOR"
    bl_region_type = "UI"
    bl_category = "SloppyUV"
    
    def draw(self, context):
        layout = self.layout
        props = context.scene.sloppy_props

class SloppyUVAlignPanel(bpy.types.Panel):
    bl_idname = "UV_PT_SloppyUVAlignPanel"
    bl_parent_id = "UV_PT_SloppyUVPanel"
    bl_label = "Align UVs"
    bl_space_type = "IMAGE_EDITOR"
    bl_region_type = "UI"
    bl_category = "SloppyUV"
    
    def draw(self, context):
        layout = self.layout
        align_box = layout.box()
        label = align_box.label(text="Align with:")
        align_grid = align_box.grid_flow(row_major = True, columns=2)
        props = context.scene.sloppy_props
        
        if props.align_uv_geo_axis[0] or props.align_uv_geo_axis[1] == True:
            align_grid.operator("operator.sloppy_align_uvs_geo")
        else:
            align_grid.label(text="Please Select Axis!")
        align_grid.prop(props, "align_uv_geo_axis")
        if props.align_uv_axis[0] or props.align_uv_axis[1] == True:
            align_grid.operator("operator.sloppy_align_uvs_uv")
        else:
            align_grid.label(text="Please Select Axis!")
        align_grid.prop(props, "align_uv_axis")
        
        align_box.prop(props, "align_uv_scale_to_bounds")
        if props.align_uv_scale_to_bounds == True:
            align_box.prop(props, "align_uv_scale_to_bounds_by_greatest")

class SloppyUVScalePanel(bpy.types.Panel):
    bl_idname = "UV_PT_SloppyUVScalePanel"
    bl_parent_id = "UV_PT_SloppyUVPanel"
    bl_label = "Scale UVs"
    bl_space_type = "IMAGE_EDITOR"
    bl_region_type = "UI"
    bl_category = "SloppyUV"
    
    def draw(self, context):
        props = context.scene.sloppy_props
        layout = self.layout
        scale_box = layout.box()
        
        calc_grid = scale_box.grid_flow(columns=2, even_columns=True)
        if any([props.scale_calc_island, props.scale_calc_selected, props.scale_calc_edge_length]) == True and any([props.scale_geo_space_view, props.scale_geo_space_object, props.scale_geo_space_global]) == True and any(props.scale_calc_uv_axis) == True and any(props.scale_uv_x_geo_axis) == True and any(props.scale_uv_y_geo_axis) == True:
            if all([props.scale_uv_x_geo_axis[0], props.scale_uv_y_geo_axis[0]]) == True or  all([props.scale_uv_x_geo_axis[1], props.scale_uv_y_geo_axis[1]]) == True or  all([props.scale_uv_x_geo_axis[2], props.scale_uv_y_geo_axis[2]]) == True:
                col_err1 = calc_grid.box()
                col_err1.label(text="Error:")
                col_err2 = calc_grid.box()
                col_err1.label(text="Same Swizzle")
                col_err1.label(text="in x & Y")
                col_err2.label(text="not allowed!")
            else:
                calc_grid.operator("operator.sloppy_scale_calc")
#                calc_grid.operator("operator.sloppy_scale_apply")
        else:
            col_err1 = calc_grid.box()
            col_err1.label(text="Error:")
            col_err2 = calc_grid.box()
            col_err2.label(text="Error:")
            if any([props.scale_calc_island, props.scale_calc_selected, props.scale_calc_edge_length]) == True:
                pass
            else:
                col_err1.label(text="Select Mode!")
                col_err2.label(text="Select Mode!")
            if any([props.scale_geo_space_view, props.scale_geo_space_object, props.scale_geo_space_global]) == True:
                pass
            else:
                col_err1.label(text="Select Space!")
                col_err2.label(text="Select Space!")
            if props.scale_calc_edge_length == True:
                if any(props.scale_calc_uv_axis) == True:
                    pass
                else:
                    col_err1.label(text="Select UV Axis!")
                    col_err2.label(text="Select UV Axis!")
            else:
                if any(props.scale_uv_x_geo_axis) == True:
                    
                    pass
                else:
                    col_err1.label(text="Select Swizzle X!")
                    col_err2.label(text="Select Swizzle X!")
                if any(props.scale_uv_y_geo_axis) == True:
                    pass
                else:
                    col_err1.label(text="Select Swizzle Y!")
                    col_err2.label(text="Select Swizzle Y!")
# started to get too compliocated                    
#        col_mode = calc_grid.column(heading="Mode:")
#        col_mode.prop(props, "scale_calc_island")
#        col_mode.prop(props, "scale_calc_selected")
#        col_mode.prop(props, "scale_calc_edge_length")
#        col_space = calc_grid.column(heading="Space:")
#        col_space.prop(props, "scale_geo_space_view")
#        col_space.prop(props, "scale_geo_space_global")
        
        scale_box.label(text="Calculation Results:")
        axes_grid = scale_box.grid_flow(columns=3, even_columns=True)
        if props.scale_calc_edge_length == True:
            col_axis = axes_grid.column()
            col_axis.prop(props, "scale_calc_uv_axis")
        else:
            col_swizzle = axes_grid.column()
            lbl_swizzle = col_swizzle.label(text="SwizzleXYZ:")
            row_swizzle1 = col_swizzle.row()
            row_swizzle1.prop(props, "scale_uv_x_geo_axis")
            row_swizzle2 = col_swizzle.row()
            row_swizzle2.prop(props, "scale_uv_y_geo_axis")
        col_geo = axes_grid.column()
        prop_geo = props.scale_geo_size
        prop_geo_x = "{:.3f}".format(prop_geo[0])
        prop_geo_y = "{:.3f}".format(prop_geo[1])
        col_geo.label(text="3D Size:")
        col_geo.label(text=prop_geo_x)
        col_geo.label(text=prop_geo_y)
        col_uv = axes_grid.column()
        prop_uv = props.scale_uv_size
        prop_uv_x = "{:.3f}".format(prop_uv[0])
        prop_uv_y = "{:.3f}".format(prop_uv[1])
        col_uv.label(text="UV Size:")
        col_uv.label(text=prop_uv_x)
        col_uv.label(text=prop_uv_y)
        
        rel_grid = scale_box.grid_flow(columns=2, even_columns=True)
        col_georel = rel_grid.column()
        prop_georel = props.scale_uv_geo_relative
        col_georel.label(text="UVs Relative to 3D:")
        prop_georel_x = "{:.3f}".format(prop_georel[0])
        prop_georel_y = "{:.3f}".format(prop_georel[1])
        col_georel.label(text=prop_georel_x)
        col_georel.label(text=prop_georel_y)
        
        col_aspect = rel_grid.column()
        prop_aspect = props.scale_uv_aspect
        col_aspect.label(text="UVs Aspect Ratio:")
        prop_aspect_x = "{:.3f}".format(prop_aspect[0])
        prop_aspect_y = "{:.3f}".format(prop_aspect[1])
        col_aspect.label(text=prop_aspect_x)
        col_aspect.label(text=prop_aspect_y)
        
        scale_box.prop(props, "scale_uv_to_bounds")
        
        if props.scale_uv_to_bounds == False:
            scale_box.prop(props, "scale_uv_real_world")
        if props.scale_calc_edge_length == True:
            scale_box.prop(props, "scale_uv_edge_direction")
        
        scale_box.operator("operator.sloppy_scale_apply")

class SloppyPeltPanel(bpy.types.Panel):
    bl_idname = "UV_PT_SloppyPeltPanel"
    bl_parent_id = "UV_PT_SloppyUVPanel"
    bl_region_type = "UI"
    bl_space_type = "IMAGE_EDITOR"
    bl_category = "SloppyUV"
    bl_label = "Square Pelt UVs"
    
    def draw(self, context):
        layout = self.layout
        pelt_box = layout.box()
        pelt_grid = pelt_box.grid_flow(columns=2, even_columns=True)
        props = context.scene.sloppy_props
        
        pelt_grid.operator("operator.pelt_uvs")
        col_opt = pelt_grid.column(heading="Pelt Options:")
        col_opt.prop(props, "space_edges")
        col_opt.prop(props, "relax_inner")

class SloppySeamGenPanel(bpy.types.Panel):
    bl_idname = "UV_PT_SloppySeamGenPanel"
    bl_parent_id = "UV_PT_SloppyUVPanel"
    bl_region_type = "UI"
    bl_space_type = "IMAGE_EDITOR"
    bl_category = "SloppyUV"
    bl_label = "Sloppy Seam Generation"

    def draw(self, context):
        layout = self.layout
        seamgen_box = layout.box()
        props = context.scene.sloppy_props
        
        seamgen_props = seamgen_box.operator("operator.sloppy_seam_gen")
        sg_col_opt = seamgen_box.column(heading="Seam Generation Options:")
        sg_col_opt.prop(props, "seamgen_angle_factor")
        sg_col_opt.prop(props, "seamgen_avg_angle_factor")
        sg_col_opt.prop(props, "seamgen_ao_factor")
        sg_col_opt.prop(props, "seamgen_ed_factor")
        sg_col_opt.prop(props, "seamgen_no_rounds")
        sg_col_opt.prop(props, "seamgen_no_retries")
        sg_col_opt.prop(props, "seamgen_angle_thresh_start")
        sg_col_opt.prop(props, "seamgen_angle_thresh_end")
        sg_col_opt.prop(props, "seamgen_clear_seam")
        sg_col_opt.prop(props, "seamgen_unwrap")

        seamgen_props["angle_factor"] = props.seamgen_angle_factor
        seamgen_props["avg_angle_factor"] = props.seamgen_avg_angle_factor
        seamgen_props["ao_factor"] = props.seamgen_ao_factor
        seamgen_props["ed_factor"] = props.seamgen_ed_factor
        seamgen_props["no_rounds"] = props.seamgen_no_rounds
        seamgen_props["no_retries"] = props.seamgen_no_retries
        seamgen_props["angle_thresh_start"] = props.seamgen_angle_thresh_start
        seamgen_props["angle_thresh_end"] = props.seamgen_angle_thresh_end
        seamgen_props["clear_seam"] = props.seamgen_clear_seam
        seamgen_props["unwrap"] = props.seamgen_unwrap

class SloppyFlatQuadPanel(bpy.types.Panel):
    bl_idname = "UV_PT_SloppyFlatQuadPanel"
    bl_parent_id = "UV_PT_SloppyUVPanel"
    bl_region_type = "UI"
    bl_space_type = "IMAGE_EDITOR"
    bl_category = "SloppyUV"
    bl_label = "Flat Quad UV Unfold"

    def draw(self, context):
        layout = self.layout
        fq_box = layout.box()
        props = context.scene.sloppy_props

        fq_props = fq_box.operator("operator.sloppy_quad_unfold")
        fq_col_opt = fq_box.column(heading="Virtual Quad Unwrap Options:")
        fq_col_opt.prop(props, "procuv_unfold_mode")
        fq_col_opt.prop(props, "procuv_initial_quad")
        fq_col_opt.prop(props, "procuv_reg_flat_fac")
        fq_col_opt.prop(props, "procuv_reg_norm_fac")
        fq_col_opt.prop(props, "procuv_pre_calc_edge_lengths")
        fq_col_opt.prop(props, "procuv_offset_per_island")
        fq_col_opt.prop(props, "procuv_quant_avg_norm")
        fq_col_opt.prop(props, "procuv_only_move_loops_in_face")
        fq_props["unfold_mode"] = props.procuv_unfold_mode
        fq_props["initial_quad"] = props.procuv_initial_quad
        fq_props["reg_flat_fac"] = props.procuv_reg_flat_fac
        fq_props["reg_norm_fac"] = props.procuv_reg_norm_fac
        fq_props["pre_calc_edge_lengths"] = props.procuv_pre_calc_edge_lengths
        fq_props["offset_per_island"] = props.procuv_offset_per_island
        fq_props["quant_avg_norm"] = props.procuv_quant_avg_norm
        fq_props["only_move_loops_in_face"] = props.procuv_only_move_loops_in_face

class SloppySortDistPanel(bpy.types.Panel):
    bl_idname = "UV_PT_SloppySortDistPanel"
    bl_parent_id = "UV_PT_SloppyUVPanel"
    bl_region_type = "UI"
    bl_space_type = "IMAGE_EDITOR"
    bl_category = "SloppyUV"
    bl_label = "Sloppy Sort Index by Distance Traveled"

    def draw(self, context):
        layout = self.layout
        sortdist_box = layout.box()
        props = context.scene.sloppy_props

        sortdist_box.operator("operator.sloppy_sort_vert_by_dist")
        sortdist_box.operator("operator.sloppy_sort_edge_by_dist")
        sortdist_box.operator("operator.sloppy_sort_face_by_dist")
        sortdist_opt_row = sortdist_box.row()
        sortdist_opt_row.prop(props, "distsort_bottom_up")
        if props.distsort_bottom_up == True:
            sortdist_opt_row.prop(props, "distsort_z_only")

class SloppyDebugPanel(bpy.types.Panel):
    bl_idname = "UV_PT_SloppyDebugPanel"
    bl_region_type = "UI"
    bl_space_type = "IMAGE_EDITOR"
    bl_category = "SloppyUV"
    bl_label = "Debug"
    
    def draw(self, context):
        layout = self.layout
        row = layout.row()
        props = context.scene.sloppy_props
        
        layout.prop(props, "verbose")

class SloppyBRPanel(bpy.types.Panel):
    bl_idname = "UV_PT_SloppyBRPanel"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "SloppyUV"
    bl_label = "Batch Render"

    def draw(self, context):
        props = context.scene.sloppy_props
        layout = self.layout
        layout.prop(props, "br_output_path")
        layout.prop(props, "br_use_viewport")
        layout.operator("operator.sloppy_batch_render")

class SloppyExtraObjPropsPanel(bpy.types.Panel):
    bl_idname = "UV_PT_SloppyExtraObjPropsPanel"
    bl_label = "Extra Properties"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Item"
    
    def draw(self, context):
        props = context.scene.sloppy_props
        layout = self.layout
        context = bpy.context
        obj = context.object
        parent = obj.parent
        row = layout.column(align=True)
        row.scale_y = 0.8

        if parent:
            row.label(text="Parent Inverse Location:")
            row.prop(obj, "matrix_parent_inverse", text="X", index=12)
            row.prop(obj, "matrix_parent_inverse", text="Y", index=13)
            row.prop(obj, "matrix_parent_inverse", text="Z", index=14)
            row.separator()
            row.label(text="World Location:")
            row.prop(obj, "matrix_world", text="X", index=12)
            row.prop(obj, "matrix_world", text="Y", index=13)
            row.prop(obj, "matrix_world", text="Z", index=14)
            row.separator()
        
        if obj.mode == "OBJECT":
            copy_grid = row.grid_flow(columns=2)

            copy_grid.operator("operator.sloppy_copy_location")
            copy_grid.operator("operator.sloppy_copy_rotation")
            copy_grid.operator("operator.sloppy_copy_scale")
            copy_grid.operator("operator.sloppy_copy_dimensions")

            row.separator()
            row.operator("operator.sloppy_copy_transform")
            row.separator()

            paste_grid = row.grid_flow(columns=2)

            paste_grid.operator("operator.sloppy_paste_location")
            paste_grid.operator("operator.sloppy_paste_rotation")
            paste_grid.operator("operator.sloppy_paste_scale")
            paste_grid.operator("operator.sloppy_paste_dimensions")
            row.separator()
            row.operator("operator.sloppy_paste_transform")
            row.separator()

            if props.delta_rotation_mode == "A":
                row.prop(obj, "delta_rotation_euler", text="Delta Rotation")
            else:
                row.prop(obj, "delta_rotation_quaternion", text="Delta Rotation")
            row.prop(props, "delta_rotation_mode", text="")

        if obj.type == "CURVE":
            total_curve_length = 0
            row.separator()
            row.label(text="Spline Lengths:")
            roro = row.row(align=True)
            rororo = roro.column(align=True)
            roroco = roro.column(align=True)
            roroop = roro.column(align=False)
            for spli, spl in enumerate(obj.data.splines):
                this_length = spl.calc_length()
                total_curve_length += this_length
                rororo.label(text=f"Spline {spli}:")
                roroco.label(text=f"{this_length:.3f}")
                this_op_props = roroop.operator("operator.sloppy_copy_to_clipboard")
                this_op_props.copied_value = this_length
            rororo.label(text=f"Total:")
            roroco.label(text=f"{total_curve_length:.3f}")
            tot_op_props = roroop.operator("operator.sloppy_copy_to_clipboard")
            tot_op_props.copied_value = total_curve_length

        if obj.mode == "EDIT":
            total_area = 0
            selected_area = 0
            total_edge_length = 0
            rowr = row.row(align=True)
            rowrow = rowr.column(align=True)
            rowcow = rowr.column(align=True)
            rowops = rowr.column(align=True)
            for ob in bpy.context.objects_in_mode:
                ob.update_from_editmode()
                for pol in ob.data.polygons:
                    total_area += pol.area
                    if pol.select == True:
                        selected_area += pol.area
                for e in ob.data.edges:
                    if e.select == True:
                        v1 = e.id_data.vertices[e.vertices[0]]
                        v2 = e.id_data.vertices[e.vertices[1]]
                        total_edge_length += math.dist(v1.co, v2.co)
            
            rowrow.label(text=f"Total Area:")
            rowcow.label(text=f"{total_area:.3f}")
            totar_props = rowops.operator("operator.sloppy_copy_to_clipboard")
            totar_props.copied_value = total_area
            
            rowrow.label(text=f"Selected Area:")
            rowcow.label(text=f"{selected_area:.3f}")
            selar_props = rowops.operator("operator.sloppy_copy_to_clipboard")
            selar_props.copied_value = selected_area
            
            rowrow.label(text=f"Selected Edges Length:")
            rowcow.label(text=f"{total_edge_length:.3f}")
            totel_props = rowops.operator("operator.sloppy_copy_to_clipboard")
            totel_props.copied_value = total_edge_length

class SloppyCPTPanel(bpy.types.Panel):
    bl_idname = "UV_PT_SloppyCPTPanel"
    bl_parent_id = "UV_PT_SloppyExtraObjPropsPanel"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Item"
    bl_label = "Copied Transform"

    def draw(self, context):
        props = context.scene.sloppy_props
        layout = self.layout
        row = layout.column(align=True)
        row.prop(props, "cpt_location")
        row.prop(props, "cpt_rotation")
        row.prop(props, "cpt_scale")
        row.prop(props, "cpt_dimensions")


# endregion

class SloppyErrorDialog(bpy.types.Operator):
    bl_idname = "operator.sloppy_dialog"
    bl_label = "Error!"
    
    def execute(self, context):
        props = context.scene.sloppy_props
        heading = props.sloppy_error_msg_heading
        message = props.sloppy_error_msg
        
        self.report({"INFO"}, heading)
        self.report({"INFO"}, message)
        
        return {"FINISHED"}
    
    def invoke(self, context, event):
        wm = context.window_manager
        return wm.invoke_props_dialog(self)
    
    def draw(self,context):
        props = context.scene.sloppy_props
        msg_heading = props.sloppy_error_msg_heading
        msg = props.sloppy_error_msg
        layout = self.layout
        box_err = layout.box()
        box_err.label(text = msg_heading)
        box_err.label(text = msg)

# Scaling section very unfinished
class SloppyCalcScale(bpy.types.Operator):
    bl_idname = "operator.sloppy_scale_calc"
    bl_label = "Calculate"
    bl_description = "Calculate UV scaling"
    
    def execute(self, context):
        props = context.scene.sloppy_props
        bpy.ops.mesh.select_linked(delimit={'SEAM'})
        bm = bmesh.from_edit_mesh(bpy.context.active_object.data)
        uv_layer = bm.loops.layers.uv.verify()
        
        # find view3d region and space
        def sort_loops_geo(e):
            props = context.scene.sloppy_props
            if props.scale_geo_space_view == True:
                if props.current_sort_axis == "X":
                    return props.viewco(e.vert.co).x
                if props.current_sort_axis == "Y":
                    return props.viewco(e.vert.co).y
            if props.scale_geo_space_global == True:
                if props.current_sort_axis == "X":
                    return e.vert.co.x
                if props.current_sort_axis == "Y":
                    return e.vert.co.y
                if props.current_sort_axis == "Z":
                    return e.vert.co.z
                
        def sort_loops_uv(e):
            props = context.scene.sloppy_props
            if props.current_sort_axis == "X":
                return e[uv_layer].uv.x
            if props.current_sort_axis == "Y":
                    return e[uv_layer].uv.y
        
        island_loops = []
        selected_faces = []
        selected_uv_loops = []
        selected_uv_edges = []
        
        for f in bm.faces:
            if f.select == True:
                selected_faces.append(f)
                for loop in f.loops:
                    uva = loop[uv_layer]
                    if loop not in island_loops:
                        island_loops.append(loop)
                    if uva.select == True:
                        if loop not in selected_uv_loops:
                            selected_uv_loops.append(loop)
                    if uva.select_edge == True:
                        if loop not in selected_uv_edges:
                            selected_uv_edges.append(loop)
        
        if len(selected_faces) < 1:
            props.sloppy_error_msg_heading = "Can't calculate scaling!"
            props.sloppy_error_msg = "    No faces selected."
            bpy.ops.operator.sloppy_dialog('INVOKE_DEFAULT')
            if props.verbose:
                print("UV Scale Calculation cancelled...\n")
            return {'CANCELLED'}
        
        if len(selected_uv_edges) < 1:
            props.sloppy_error_msg_heading = "Can't calculate scaling!"
            props.sloppy_error_msg = "    No UVs selected."
            bpy.ops.operator.sloppy_dialog('INVOKE_DEFAULT')
            if props.verbose:
                print("UV Scale Calculation cancelled...\n")
            return {'CANCELLED'}
        
        island_bounds_loops = []
        selected_bounds_loops = []
        edges_bounds_loops = []
        island_bounds_geo = []
        selected_bounds_geo = []
           
        props.current_sort_axis = "X"
        island_loops.sort(key=sort_loops_uv)
        selected_uv_loops.sort(key=sort_loops_uv)
        island_bounds_loops.append(island_loops[0])
        island_bounds_loops.append(island_loops[-1])
        selected_bounds_loops.append(selected_uv_loops[0])
        selected_bounds_loops.append(selected_uv_loops[-1])
        selected_bounds_geo.append(island_loops[0])
        selected_bounds_geo.append(island_loops[-1])
        props.current_sort_axis = "Y"
        island_loops.sort(key=sort_loops_uv)
        selected_uv_loops.sort(key=sort_loops_uv)
        island_bounds_loops.append(island_loops[0])
        island_bounds_loops.append(island_loops[-1])
        selected_bounds_loops.append(selected_uv_loops[0])
        selected_bounds_loops.append(selected_uv_loops[-1])
        selected_bounds_geo.append(island_loops[0])
        selected_bounds_geo.append(island_loops[-1])
        
        props.current_sort_axis = "X"
        island_loops.sort(key=sort_loops_geo)
        island_bounds_geo.append(island_loops[0])
        island_bounds_geo.append(island_loops[-1])
        props.current_sort_axis = "Y"
        island_loops.sort(key=sort_loops_geo)
        island_bounds_geo.append(island_loops[0])
        island_bounds_geo.append(island_loops[-1])
#        props.current_sort_axis = "Z"
#        island_loops.sort(key=sort_loops_geo)
#        island_bounds_geo.append(island_loops[0])
#        island_bounds_geo.append(island_loops[-1])
        
        props.scale_uv_size_current[0] = island_bounds_loops[1][uv_layer].uv.x - island_bounds_loops[0][uv_layer].uv.x
        props.scale_uv_size_current[1] = island_bounds_loops[3][uv_layer].uv.y - island_bounds_loops[2][uv_layer].uv.y
        props.scale_calc_uv_corner[0] = island_bounds_loops[0][uv_layer].uv.x
        props.scale_calc_uv_corner[1] = island_bounds_loops[2][uv_layer].uv.y
        
        axis_int = 0
        if props.scale_calc_uv_axis[0] == True:
            props.current_sort_axis = "X"
            axis_int = 0
        if props.scale_calc_uv_axis[1] == True:
            props.current_sort_axis = "Y"
            axis_int = 1
        selected_uv_edges.sort(key=sort_loops_uv)
        edges_bounds_loops.append(selected_uv_edges[0])
        edges_bounds_loops.append(selected_uv_edges[-1])
        
        edge_length_total = 0.0
        
        vec = edges_bounds_loops[1][uv_layer].uv - edges_bounds_loops[0][uv_layer].uv
        vec_len = vec.length
        dir = vec.normalized()
        for loop in selected_uv_edges:
            edge_length_total += loop.edge.calc_length()
        props.scale_geo_edge_lengths_real[axis_int] = edge_length_total
        props.scale_uv_edge_lengths[axis_int] = vec_len
        if axis_int == 0:
            print(f"Direction of X: {dir}")
            props.scale_uv_edge_direction_x = dir
        if axis_int == 1:
            print(f"Direction of Y: {dir}")
            props.scale_uv_edge_direction_y = dir
                
        props.update_scalings(context)
        
        return {"FINISHED"}

class SloppyApplyScale(bpy.types.Operator):
    bl_idname = "operator.sloppy_scale_apply"
    bl_label = "Apply"
    bl_description = "Apply calculated scaling to UV map"
    
    def execute(self, context):
        props = context.scene.sloppy_props
        bpy.ops.mesh.select_linked(delimit={'SEAM'})
        bm = bmesh.from_edit_mesh(bpy.context.active_object.data)
        uv_layer = bm.loops.layers.uv.verify()
        
        scale_factor_x = props.scale_uv_size.x / props.scale_uv_size_current.x
        scale_factor_y = props.scale_uv_size.y / props.scale_uv_size_current.y
        minx = props.scale_calc_uv_corner.x
        miny = props.scale_calc_uv_corner.y
        
        selected_faces = []
        xs = []
        ys = []
        
        for f in bm.faces:
            if f.select == True:
                selected_faces.append(f)
                for loop in f.loops:
                    uva = loop[uv_layer]
                    uvco = uva.uv
                    trans_x = uvco.x
                    trans_y = uvco.y
                    if props.scale_uv_real_world == True or props.scale_uv_to_bounds == True:
                        trans_x = uvco.x - minx
                        trans_y = uvco.y - miny
                    if props.scale_uv_edge_direction == True:
                        trans_co = mathutils.Vector((trans_x, trans_y))
                        print(trans_co)
                        newco1 = trans_co * props.scale_uv_edge_direction_x * scale_factor_x
                        print(newco1)
                        newco2 = newco1 * props.scale_uv_edge_direction_y * scale_factor_y
                        print(newco2)
                        uvco.x = newco2.x
                        xs.append(newco2.x)
                        uvco.y = newco2.y
                        ys.append(newco2.y)
                    else:
                        new_x = trans_x * scale_factor_x
                        new_y = trans_y * scale_factor_y
                        uvco.x = new_x
                        xs.append(new_x)
                        uvco.y = new_y
                        ys.append(new_y)
        
        if len(selected_faces) < 1:
            props.sloppy_error_msg_heading = "Can't apply scaling!"
            props.sloppy_error_msg = "    No faces selected."
            bpy.ops.operator.sloppy_dialog('INVOKE_DEFAULT')
            if props.verbose:
                print("UV Scaling cancelled...\n")
            return {'CANCELLED'}
        
        bpy.context.active_object.data.update()
        
        if props.scale_uv_to_bounds == True:
            minx = min(xs)
            maxx = max(xs)
            x_div = maxx - minx
            miny = min(ys)
            maxy = max(ys)
            y_div = maxy - miny
            div_max = max(x_div, y_div)
            
            iter = 0
            iter_max = len(selected_faces)
            l_iter = 0
            
            for f in selected_faces:
                iter += 1
                for loop in f.loops:
                    l_iter += 1
                    uva = loop[uv_layer]
                    uvco = uva.uv
                    uvco_x_old = uva.uv.x
                    uvco_y_old = uva.uv.y
                    trans_x = uvco.x - minx
                    new_x = trans_x / x_div
                    trans_y = uvco.y - miny
                    new_y = trans_y / y_div
                    if props.align_uv_scale_to_bounds_by_greatest == True:
                        new_x = trans_x / div_max
                        new_y = trans_y / div_max
                    uvco.x = new_x
                    uvco.y = new_y
        
        bpy.context.active_object.data.update()
        
        props.scale_uv_size_current.x = props.scale_uv_size.x
        props.scale_uv_size_current.y = props.scale_uv_size.y
        if props.scale_uv_real_world == True or props.scale_uv_to_bounds == True:
            props.scale_calc_uv_corner.x = 0.0
            props.scale_calc_uv_corner.x = 0.0
        
        props.update_scalings(context)
        
        return {"FINISHED"}

class SloppyAlignUVsGeo(bpy.types.Operator):
    bl_idname = "operator.sloppy_align_uvs_geo"
    bl_label = "3D View"
    bl_description = "Rotate UVs to match geometry in 3D view using axis chosen below"
    
    def execute(self, context):
        props = context.scene.sloppy_props
        props.align_uvs_use_geo = True
        
        bpy.ops.operator.sloppy_align_uvs()
        
        return {"FINISHED"}

class SloppyAlignUVsUV(bpy.types.Operator):
    bl_idname = "operator.sloppy_align_uvs_uv"
    bl_label = "Selected UVs"
    bl_description = "Rotate UVs to align selected UV points along axis chosen below"
    
    def execute(self, context):
        props = context.scene.sloppy_props
        props.align_uvs_use_geo = False
        
        bpy.ops.operator.sloppy_align_uvs()
        
        return {"FINISHED"}

class SloppyAlignUVs(bpy.types.Operator):
    bl_idname = "operator.sloppy_align_uvs"
    bl_label = "Align UVs"
    bl_description = "Rotate UVs to make them match either UV axis or 3D axis"
    
    def execute(self, context):
        props = context.scene.sloppy_props
        bpy.ops.mesh.select_linked(delimit={'SEAM'})
        bm = bmesh.from_edit_mesh(bpy.context.active_object.data)
        uv_layer = bm.loops.layers.uv.verify()
        
        # initialize console commands
        CURSOR_UP = '\033[F'
        ERASE_LINE = '\033[K'
        
        def sort_loops(e):
            props = context.scene.sloppy_props
            if props.align_uvs_use_geo == True:
                if props.align_uv_geo_axis[0] == True:
                    return props.viewco(e.vert.co).x
                if props.align_uv_geo_axis[1] == True:
                    return props.viewco(e.vert.co).y
            if props.align_uvs_use_geo == False:
                if props.align_uv_axis[0] == True:
                    return e[uv_layer].uv.x
                if props.align_uv_axis[1] == True:
                    return e[uv_layer].uv.y

        loops = []
        
        vec = mathutils.Vector((0,0))
        dir = mathutils.Vector((0,0))
        
        if props.verbose == True:
            print("\nUV Alignment started.\n")
            print("Step 1: Collect pertinent loops")
            print("===============================")
        
        iter = 0
        iter_max = len(bm.faces)
        selected_faces = []
        
        if props.align_uvs_use_geo == True:
            for f in bm.faces:
                iter += 1
                if f.select == True:
                    selected_faces.append(f)
                    for loop in f.loops:
                        if loop not in loops:
                            loops.append(loop)
                if props.verbose == True:
                    done_pct = (iter / iter_max) * 100
                    if iter > 1:
                        print(CURSOR_UP + CURSOR_UP + CURSOR_UP + CURSOR_UP)
                    msg = "Processing face {:} of {:}.\nNumber of loops collected: {:}. \nStep 1 {:.1f}% done.".format(iter, iter_max, len(loops), done_pct)
                    print(msg, flush=True)
                    sys.stdout.flush()
            
        if props.align_uvs_use_geo == False:
            for f in bm.faces:
                iter += 1
                if f.select == True:
                    selected_faces.append(f)
                    for loop in f.loops:
                        uva = loop[uv_layer]
                        if uva.select == True:
                            if loop not in loops:
                                loops.append(loop)
                if props.verbose == True:
                    done_pct = (iter / iter_max) * 100
                    if iter > 1:
                        print(CURSOR_UP + CURSOR_UP + CURSOR_UP + CURSOR_UP)
                    msg = "Processing face {:} of {:}.\nNumber of loops collected: {:}. \nStep 1 {:.1f}% done.".format(iter, iter_max, len(loops), done_pct)
                    print(msg, flush=True)
                    sys.stdout.flush()
        
        if len(loops) < 1:
            props.sloppy_error_msg_heading = "Can't align UVs!"
            if props.align_uvs_use_geo == True:
                props.sloppy_error_msg = "    No faces selected."
                bpy.ops.operator.sloppy_dialog('INVOKE_DEFAULT')
            if props.align_uvs_use_geo == False:
                props.sloppy_error_msg = "    No UVs selected."
                bpy.ops.operator.sloppy_dialog('INVOKE_DEFAULT')
            if props.verbose:
                print("UV Alignment cancelled...\n")
            return {'CANCELLED'}
        
        if props.verbose == True:
            print(CURSOR_UP + ERASE_LINE + CURSOR_UP)
            print("Step 1 complete!\n")
            print("Step 2: Calculate rotational difference")
            print("=======================================")
        
        if props.align_uvs_use_geo == True:
            loops.sort(key=sort_loops)
            vec = props.viewco(loops[-1].vert.co) - props.viewco(loops[0].vert.co)
            dir = vec.normalized()
        if props.align_uvs_use_geo == False:
            loops.sort(key=sort_loops)
            if props.align_uv_axis[0] == True:
                dir = mathutils.Vector((1,0))
            if props.align_uv_axis[1] == True:
                dir = mathutils.Vector((0,1))
        
        uvec = loops[-1][uv_layer].uv - loops[0][uv_layer].uv
        udir = uvec.normalized()
        uangle = udir.angle_signed(dir)
        
        if props.verbose == True:
            print("{:} loops sorted.\nCalculated rotational difference: {:.1f}°".format(len(loops), math.degrees(uangle)))
            print("Step 2 complete!\n")
            print("Step 3: Rotate UVs")
            print("==================")
        
        mid_x = bl_math.lerp(loops[-1][uv_layer].uv.x, loops[0][uv_layer].uv.x, 0.5)
        mid_y = bl_math.lerp(loops[-1][uv_layer].uv.y, loops[0][uv_layer].uv.y, 0.5)
        mid = mathutils.Vector((mid_x, mid_y))
        
        xs = []
        ys = []
        
        iter = 0
        iter_max = len(selected_faces)
        l_iter = 0
        
        for f in selected_faces:
            iter += 1
            for loop in f.loops:
                l_iter += 1
                uva = loop[uv_layer]
                uvco = uva.uv
                uvco_x_old = uva.uv.x
                uvco_y_old = uva.uv.y
                vec = uvco - mid
                rot = mathutils.Matrix.Rotation(-uangle, 2, 'Z')
                vec.rotate(rot)
                if props.align_uvs_use_geo == True:
                    mid = mathutils.Vector((0.5, 0.5))
                uvco.x = mid.x + vec.x
                uvco.y = mid.y + vec.y
                xs.append(uvco.x)
                ys.append(uvco.y)
                if props.verbose == True:
                    done_pct = (iter / iter_max) * 100
                    if l_iter > 1:
                        print(CURSOR_UP + CURSOR_UP + CURSOR_UP + CURSOR_UP)
                    msg = "Processing face {:} of {:}.\nMoved loop {:} from {:.3f}, {:.3f} to {:.3f}, {:.3f}. \nStep 3 {:.1f}% done.".format(iter, iter_max, loop.index, uvco_x_old, uvco_y_old, uvco.x, uvco.y, done_pct)
                    print(msg, flush=True)
                    sys.stdout.flush()
        bpy.context.active_object.data.update()
        
        if props.verbose == True:
            print(CURSOR_UP + ERASE_LINE + CURSOR_UP)
            print("Step 3 complete!\n")
        
        if props.align_uv_scale_to_bounds == True:
            if props.verbose == True:
                print("Step 4: Scale To Bounds")
                print("==================")
            minx = min(xs)
            maxx = max(xs)
            x_div = maxx - minx
            miny = min(ys)
            maxy = max(ys)
            y_div = maxy - miny
            div_max = max(x_div, y_div)
            
            iter = 0
            iter_max = len(selected_faces)
            l_iter = 0
            
            for f in selected_faces:
                iter += 1
                for loop in f.loops:
                    l_iter += 1
                    uva = loop[uv_layer]
                    uvco = uva.uv
                    uvco_x_old = uva.uv.x
                    uvco_y_old = uva.uv.y
                    trans_x = uvco.x - minx
                    new_x = trans_x / x_div
                    trans_y = uvco.y - miny
                    new_y = trans_y / y_div
                    if props.align_uv_scale_to_bounds_by_greatest == True:
                        new_x = trans_x / div_max
                        new_y = trans_y / div_max
                    uvco.x = new_x
                    uvco.y = new_y
                    if props.verbose == True:
                        done_pct = (iter / iter_max) * 100
                        if l_iter > 1:
                            print(CURSOR_UP + CURSOR_UP + CURSOR_UP + CURSOR_UP + CURSOR_UP)
                        msg = ""
                        if props.align_uv_scale_to_bounds_by_greatest == True:
                            scale_factor = 1.0 / div_max
                            msg = "Processing face {:} of {:}.\nUniform scale factor: {:.3f}\nMoved loop {:} from {:.3f}, {:.3f} to {:.3f}, {:.3f}. \nStep 4 {:.1f}% done.".format(iter, iter_max, scale_factor, loop.index, uvco_x_old, uvco_y_old, uvco.x, uvco.y, done_pct)
                        else:
                            scale_factor_x = 1.0 / x_div
                            scale_factor_y = 1.0 / y_div
                            msg = "Processing face {:} of {:}.\nScale factor: {:.3f}, {:.3f}\nMoved loop {:} from {:.3f}, {:.3f} to {:.3f}, {:.3f}. \nStep 4 {:.1f}% done.".format(iter, iter_max, scale_factor_x, scale_factor_y, loop.index, uvco_x_old, uvco_y_old, uvco.x, uvco.y, done_pct)
                        print(msg, flush=True)
                        sys.stdout.flush()
        if props.verbose == True:
            print(CURSOR_UP + ERASE_LINE + CURSOR_UP)
            print("Step 4 complete!\n")  
        bpy.context.active_object.data.update()
        
        if props.verbose:
            print("UV Alignment complete!\n")
        
        return {"FINISHED"}

class SloppyUVSmartAlign(bpy.types.Operator):
    bl_idname = "operator.sloppy_uv_smart_align"
    bl_label = "Smart Align UVs"
    bl_description = "Align line of UV edges to axis and auto-unwrap result"
    bl_options = {'REGISTER', 'UNDO'}
    
    # region UVSmartAlign Properties
    iP = bpy.props.IntProperty
    fP = bpy.props.FloatProperty
    fvP = bpy.props.FloatVectorProperty
    bP = bpy.props.BoolProperty
    eP = bpy.props.EnumProperty
    bvP = bpy.props.BoolVectorProperty
    sP = bpy.props.StringProperty

    align_axis : eP(
        name = "Axis",
        description = "Axis to align UVs with",
        items = [
            ("A", "Average", "Align with average axis of selection"),
            ("B", "Endpoints", "Align with vector between endpoints of selection"),
            ("C", "U", "Align selected UVs with U/X"),
            ("D", "V", "Align selected UVs with V/Y")
            ],
        default="A"
        ) # type: ignore

    align_pivot : eP(
        name = "Pivot",
        description = "Pivot to align UVs with (note that any intersecting pinned UVs will override this)",
        items = [
            ("A", "Middle", "Middle of length is used as pivot"),
            ("B", "Endpoints", "Midpoint between endpoints is used as pivot"),
            ("C", "Endpoint A", "Endpoint A is used as pivot"),
            ("D", "Endpoint B", "Endpoint B is used as pivot")
            ],
        default="A"
        ) # type: ignore

    unwrap_mode : eP(
        name = "Auto-Unwrap Mode",
        description = "",
        items = [
            ("A", "Angle-Based", ""),
            ("B", "Conformal", "")
            ],
        default="A"
        ) # type: ignore

    keep_pins : bP(
        name = "Keep Pins",
        description = "Keep UV pins",
        default = True
        ) # type: ignore

    respect_pins : bP(
        name = "Respect Pins",
        description = "Respect existing pins",
        default = True
        ) # type: ignore
    #endregion


    def execute(self, context):
        props = context.scene.sloppy_props
        bms = []
        for ed_ob in bpy.context.objects_in_mode:
            nbm = bmesh.from_edit_mesh(ed_ob.data)
            bms.append(nbm)
        
        for bm in bms:
            uv_layer = bm.loops.layers.uv.verify()
            islands = bmesh_utils.bmesh_linked_uv_islands(bm, uv_layer)
            
            # initialize console commands
            CURSOR_UP = '\033[F'
            ERASE_LINE = '\033[K'
            
            def sort_isu_map(e):
                return e[0]
            
            def sort_by_u(e):
                return e[2].x
            
            def sort_by_v(e):
                return e[2].y

            if props.verbose == True:
                print("\nUV Alignment started.\n")
                print("Step 1: Collect pertinent edges and loops")
                print("=========================================")

            island_loop_map = []
            selected_edges = []
            selected_uvs = []
            pinned_uvs = []
            
            for island in islands:
                island_loops = []
                island_selected_edges = []
                island_selected_uvs = []
                island_pinned_uvs = []
                for i_f in island:
                    for il in i_f.loops:
                        if il not in island_loops:
                            island_loops.append(il)
                            if il[uv_layer].select == True:
                                if il not in island_selected_uvs:
                                    island_selected_uvs.append(il)
                                if il[uv_layer].select_edge == True:
                                    if il.edge not in island_selected_edges:
                                        island_selected_edges.append(il.edge)
                            if il[uv_layer].pin_uv == True:
                                if il not in island_pinned_uvs:
                                    island_pinned_uvs.append(il)
                island_loop_map.append(island_loops)
                selected_edges.append(island_selected_edges)
                selected_uvs.append(island_selected_uvs)
                pinned_uvs.append(island_pinned_uvs)
            
            ii = 0
            for ilm, ise, isu, ipu, isls in zip(island_loop_map, selected_edges, selected_uvs, pinned_uvs, islands):
                
                abort_island = False
                
                if len(isu) < 1:
                    abort_island = True

                if abort_island == False:
                    print("Island", ii, "- selected UV edges:", len(ise))
                    uvlcm = props.generate_uv_loop_connection_map(isu, uv_layer)

                    pin_intersect = None
                    for su in isu:
                        if su in ipu:
                            if pin_intersect:
                                if pin_intersect[-1][uv_layer].uv != su[uv_layer].uv:
                                    pin_intersect.append(su)
                            else:
                                pin_intersect = [su]
                    
                    length_relation = 0.0
                    uv_length = 0.0
                    uv_length_mid = 0.0
                    geo_length = 0.0
                    uv_pivot = None
                    uv_pivot_offset = mathutils.Vector((0,0))
                    too_many_intersections = False
                    if pin_intersect:
                        print("Number of pinned intersections:", len(pin_intersect))
                        if self.respect_pins == True:
                            if len(pin_intersect) > 1:
                                uvlerpa = pin_intersect[0][uv_layer].uv
                                uvlerpb = pin_intersect[1][uv_layer].uv
                                uv_pivot = uvlerpa.lerp(uvlerpb, 0.5)
                                if len(pin_intersect) > 2:
                                    too_many_intersections = True
                            else:
                                if len(pin_intersect) == 1:
                                    uv_pivot = pin_intersect[0][uv_layer].uv
                    
                        if too_many_intersections == True:
                            if props.verbose:
                                print("Too many pin intersections (", len(pin_intersect), ")! Island", ii, "marked for abort...")
                            abort_island = True
                    
                    if abort_island == False:
                        for se in ise:
                            geo_length += se.calc_length()
                            lenuvs = [None, None]
                            for sevi, sev in enumerate(se.verts):
                                for sevl in sev.link_loops:
                                    if sevl in isu:
                                        if not lenuvs[sevi]:
                                            lenuvs[sevi] = sevl[uv_layer].uv
                            uv_length += math.dist(lenuvs[0], lenuvs[1])


                        length_relation = uv_length/geo_length
                        print("Length relationship:", length_relation)
                        uv_length_mid = uv_length / 2

                    isu_map = []
                    found_ends = 2
                    loops_todo = isu.copy()
                    edges_todo = ise.copy()
                    next_vdns = []
                    first_loop = loops_todo[0]
                    first_vert = first_loop.vert
                    
                    next_edges = []
                    for lve in loops_todo[0].vert.link_edges:
                        if lve in edges_todo:
                            next_edges.append(lve)
                            edges_todo.remove(lve)
                            found_ends -= 1

                    if len(next_edges) > 2:
                        print("Intersecting selection! Island", ii, "marked for abort...")
                        abort_island = True

                    for ni, ne in enumerate(next_edges):
                        current_dir = (0 + 1) - (ni * 2)
                        ov = ne.other_vert(first_vert)
                        arrival_loops = props.get_verts_filtered_loops_in_island_using_perp_edge(ov, ne, isu, isls)
                        next_vdns.append([ov, current_dir, current_dir, arrival_loops, ne.calc_length() * current_dir])
                    
                    isu_index = isu.index(first_loop)

                    isu_map.append([0, uvlcm[isu_index].copy(), first_loop[uv_layer].uv, 0.0])

                    for uvl in uvlcm[isu_index]:
                        loops_todo.remove(uvl)

                    while len(edges_todo) > 0:
                        these_vdns = next_vdns.copy()
                        next_vdns = []
                        for vdn in these_vdns:
                            tvedges = []
                            tv = vdn[0]
                            td = vdn[1]
                            tn = vdn[2]
                            al = vdn[3][0]
                            glen = vdn[4]

                            isu_index = isu.index(al)
                            isu_map.append([tn, uvlcm[isu_index].copy(), al[uv_layer].uv, glen])

                            for tve in tv.link_edges:
                                if tve in edges_todo:
                                    tvedges.append(tve)
                                    edges_todo.remove(tve)
                            
                            if len(tvedges) > 1:
                                print("Intersecting selection! Island", ii, "marked for abort...")
                                abort_island = True
                            
                            for tvi, tve in enumerate(tvedges):
                                tvov = tve.other_vert(tv)
                                tval = props.get_verts_filtered_loops_in_island_using_perp_edge(tvov, tve, isu, isls)
                                next_vdns.append([tvov, td, tn + td, tval, glen + (tve.calc_length() * td)])

                            if len(tvedges) == 0:
                                found_ends += 1
                        
                        # if found_ends >= 2:
                        #     break

                    isu_map.sort(key=sort_isu_map)

                    for isui in isu_map:
                        print(isui[0], "- loops:", [iu.index for iu in isui[1]], "- U:", isui[2].x, "V:", isui[2].y, "- Distance from 0:", isui[3])

                    axis_vec = mathutils.Vector((0,0))
                    avg_axis_vec = mathutils.Vector((0,0))
                    curr_uv_length = 0.0
                    next_uv_length = 0.0
                    curr_guv_length = 0.0
                    next_guv_length = 0.0
                    next_guv_part = 0.0

                    isu_guv_parts = [0.0]

                    midpoint_element_lerp = None
                    pivot_element_lerp = None
                    pivot_offset_dist = 0.0
                    mid_guv_length = 0.0
                    mid_guv_part = 0.0
                    part_offset = 0.0
                    
                    if uv_pivot:
                        for ppvi in range(len(isu_map)-1):
                            vec = isu_map[ppvi][2]
                            nexvec = isu_map[ppvi+1][2]
                            if uv_pivot == vec:
                                pivot_element_lerp = [isu_map[ppvi][0], isu_map[ppvi+1][0], 0.0, uv_pivot]
                            elif uv_pivot == nexvec:
                                pivot_element_lerp = [isu_map[ppvi][0], isu_map[ppvi+1][0], 1.0, uv_pivot]
                            else:
                                pvec = uv_pivot - vec
                                pnvec = uv_pivot - nexvec
                                pdir = pvec.normalized()
                                pndir = pnvec.normalized()
                                pdotn = pdir.dot(pndir)
                                if pdotn < 0:
                                    povershot_distance = math.dist(vec, uv_pivot)
                                    puv_section_length = math.dist(vec, nexvec)
                                    plerp_alpha = povershot_distance / puv_section_length
                                    pivot_element_lerp = [isu_map[ppvi][0], isu_map[ppvi+1][0], plerp_alpha, uv_pivot]

                    for pvi in range(len(isu_map)-1):
                        curr_uv_length = next_uv_length
                        curr_guv_length = next_guv_length
                        curr_guv_part = next_guv_part
                        vec = isu_map[pvi][2]
                        nexvec = isu_map[pvi+1][2]
                        uv_section_length = math.dist(vec, nexvec)
                        guv_section_length = isu_map[pvi+1][3] - isu_map[pvi][3]
                        next_uv_length = curr_uv_length + uv_section_length
                        next_guv_length =  curr_guv_length + guv_section_length
                        next_guv_part = next_guv_length / geo_length
                        isu_guv_parts.append(next_guv_part)
                        avg_axis_vec += nexvec - vec

                        if uv_length_mid > curr_uv_length and uv_length_mid < next_uv_length:
                            overshot_distance = uv_length_mid - curr_uv_length
                            lerp_alpha = overshot_distance / uv_section_length
                            lerped_pos = vec.lerp(nexvec, lerp_alpha)
                            mid_guv_length = bl_math.lerp(curr_guv_length, next_guv_length, lerp_alpha)
                            mid_guv_part = mid_guv_length / uv_length
                            midpoint_element_lerp = [isu_map[pvi][0], isu_map[pvi+1][0], lerp_alpha, lerped_pos]
                    
                    if self.align_axis == "A":
                        axis_vec = avg_axis_vec.normalized()
                    if self.align_axis == "B":
                        vec = isu_map[0][2]
                        lastvec = isu_map[-1][2]
                        endpoints_axis_vec = lastvec - vec
                        axis_vec = endpoints_axis_vec.normalized()
                    if self.align_axis == "C":
                        xa = isu_map[0][2].x
                        xb = isu_map[-1][2].x
                        axis_vec = mathutils.Vector((1,0))
                        if xa > xb:
                            axis_vec *= -1
                    if self.align_axis == "D":
                        ya = isu_map[0][2].y
                        yb = isu_map[-1][2].y
                        axis_vec = mathutils.Vector((0,1))
                        if ya > yb:
                            axis_vec *= -1
                    
                    if pivot_element_lerp:
                        # uv_pivot_offset = uv_pivot - midpoint_element_lerp[3]
                        parta = isu_guv_parts[pivot_element_lerp[0]]
                        partb = isu_guv_parts[pivot_element_lerp[1]]
                        parts_lerped = bl_math.lerp(parta, partb, pivot_element_lerp[2])
                        print("\nPart lerped:", parts_lerped)
                        # mid_guv_part = 0.5
                        # part_offset = -(1.0 - parts_lerped) + (mid_guv_part - 0.5)
                        # part_offset = -(1.0 - parts_lerped) - (0.5 - mid_guv_part)
                        if parts_lerped > mid_guv_part:
                            part_offset = -(1.0 - parts_lerped)
                        else:
                            mulfac = max(mid_guv_part, parts_lerped) - min(mid_guv_part, parts_lerped)
                            part_offset = -parts_lerped*mulfac
                            # pass
                        # part_offset = -parts_lerped
                        uv_pivot_offset = mathutils.Vector((0,0))
                    else:
                        if self.align_pivot == "A":
                            uv_pivot = midpoint_element_lerp[3]
                            uv_pivot_offset = mathutils.Vector((0,0))
                        if self.align_pivot == "B":
                            vec = isu_map[0][2]
                            nexvec = isu_map[-1][2]
                            vmid = vec.lerp(nexvec, 0.5)
                            uv_pivot = vmid
                            uv_pivot_offset = mathutils.Vector((0,0))
                        if self.align_pivot == "C":
                            vec = isu_map[0][2]
                            nexvec = isu_map[-1][2]
                            vmid = vec.lerp(nexvec, 0.5)
                            uv_pivot = vec
                            part_offset = -mid_guv_part
                            uv_pivot_offset = mathutils.Vector((0,0))
                        if self.align_pivot == "D":
                            vec = isu_map[0][2]
                            nexvec = isu_map[-1][2]
                            vmid = vec.lerp(nexvec, 0.5)
                            uv_pivot = nexvec
                            part_offset = mid_guv_part
                            uv_pivot_offset = mathutils.Vector((0,0))

                    print("Current mid geo-UV part:", mid_guv_part, "- Current part offset:", part_offset, "- Sum of part modification:", -(mid_guv_part-part_offset))

                    for isum, isump in zip(isu_map, isu_guv_parts):
                        nupart = isump - mid_guv_part - part_offset
                        new_uv = uv_pivot + uv_pivot_offset + (axis_vec * nupart * uv_length)
                        print("New UVs for element", isum[0], "(at", isump, "/", nupart, "of length): U:", new_uv.x, "V:", new_uv.y)
                        for isul in isum[1]:
                            isul[uv_layer].uv = new_uv
                            isul[uv_layer].pin_uv = True
                    

                if abort_island == False:
                    pass
                else:
                    pass

                ii += 1

        for ed_ob in bpy.context.objects_in_mode:
            ed_ob.data.update()
        
        uvw_method = 'ANGLE_BASED'
        if self.unwrap_mode == "B":
            uvw_method = 'CONFORMAL'
        
        bpy.ops.uv.unwrap(method=uvw_method, fill_holes=True, correct_aspect=True, use_subsurf_data=False, margin=0, no_flip=False, iterations=10, use_weights=False, weight_group="uv_importance", weight_factor=1)

        return {"FINISHED"}

class PeltUVs(bpy.types.Operator):
    bl_idname = "operator.pelt_uvs"
    bl_label = "Generate Pelt"
    bl_description = "Transform UVs of currently selected faces into a square"
    
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

class SloppySelectByIndex(bpy.types.Operator):
    bl_idname = "operator.select_by_index"
    bl_label = "Select By Index"
    bl_description = "Select mesh element by index or range of indices. (Can also be used to only output a list of elements.)"
    bl_options = {'REGISTER', 'UNDO'}

    iP = bpy.props.IntProperty
    fP = bpy.props.FloatProperty
    fvP = bpy.props.FloatVectorProperty
    bP = bpy.props.BoolProperty
    eP = bpy.props.EnumProperty
    bvP = bpy.props.BoolVectorProperty
    sP = bpy.props.StringProperty

    max_i = 999999999

    domain_enum : eP(
        name = "Domain",
        description = "Element domain to select",
        items = [
            ("A", "Vertex", ""),
            ("B", "Edge", ""),
            ("C", "Face", "")
            ]
        ) # type: ignore

    choose_domain : bP(
        name = "Force Domain",
        description = "Force the element domain to select",
        default = False
        ) # type: ignore

    uvs : bP(
        name = "UV",
        description = "Select UVs",
        default = False
        ) # type: ignore

    i_start : iP(
        name = "First Index",
        description = "First index to select",
        default = 0,
        min = 0,
        max = max_i
        ) # type: ignore

    i_end : iP(
        name = "Last Index",
        description = "Last index to select, if range",
        default = 0,
        min = 0,
        max = max_i
        ) # type: ignore

    sel_range : bP(
        name = "Select Range",
        description = "Select a range of indices",
        default = False
        ) # type: ignore

    extend : bP(
        name = "Extend",
        description = "Extend existing selection",
        default = False
        ) # type: ignore

    def execute(self, context):
        props = context.scene.sloppy_props
        in_bm = bmesh.from_edit_mesh(bpy.context.active_object.data)
        domains = ["A", "B", "C"]
        domain = domains.index(self.domain_enum)
        self.max_i = max(len(in_bm.verts), len(in_bm.edges), len(in_bm.faces))
        props.select_by_index_bm(domain, self.choose_domain, self.uvs, self.i_start, self.i_end, self.sel_range, self.extend, False, in_bm)
        
        return {"FINISHED"}

class SloppyShiftSelectByIndex(bpy.types.Operator):
    bl_idname = "operator.shift_select_by_index"
    bl_label = "Shift Selection By Index"
    bl_description = "Shift mesh element selection by index. (Can also be used to only output a list of elements.)"
    bl_options = {'REGISTER', 'UNDO'}

    iP = bpy.props.IntProperty
    fP = bpy.props.FloatProperty
    fvP = bpy.props.FloatVectorProperty
    bP = bpy.props.BoolProperty
    eP = bpy.props.EnumProperty
    bvP = bpy.props.BoolVectorProperty
    sP = bpy.props.StringProperty

    max_i = 999999999

    domain_enum : eP(
        name = "Domain",
        description = "Element domain to select",
        items = [
            ("A", "Vertex", ""),
            ("B", "Edge", ""),
            ("C", "Face", "")
            ]
        ) # type: ignore

    choose_domain : bP(
        name = "Force Domain",
        description = "Force the element domain to select",
        default = False
        ) # type: ignore

    shift_amount : iP(
        name = "Shift Amount",
        description = "Shift indices by this amount",
        default = 1,
        min = -max_i,
        max = max_i
        ) # type: ignore

    shift_by_len : bP(
        name = "Shift By Number of Selected",
        description = "Shift by number of previously selected elements",
        default = False
        ) # type: ignore

    shift_by_len_neg : bP(
        name = "Negative Selection Amount",
        description = "Multiply number of previously selected elements by -1",
        default = False
        ) # type: ignore

    extend : bP(
        name = "Extend",
        description = "Extend existing selection",
        default = False
        ) # type: ignore

    both_directions : bP(
        name = "Up and Down",
        description = "Shift selection in both positive and negative directions",
        default = False
        ) # type: ignore

    def execute(self, context):
        props = context.scene.sloppy_props
        in_bm = bmesh.from_edit_mesh(bpy.context.active_object.data)
        domains = ["A", "B", "C"]
        domain = domains.index(self.domain_enum)
        self.max_i = max(len(in_bm.verts), len(in_bm.edges), len(in_bm.faces))
        props.shift_select_by_index_bm(domain, self.choose_domain, self.shift_amount, self.shift_by_len, self.shift_by_len_neg, self.extend, self.both_directions, False, in_bm)
        
        return {"FINISHED"}

class SloppyFindContiguousSeams(bpy.types.Operator):
    bl_idname = "operator.find_contiguous_seams"
    bl_label = "Find Contiguous Seams"
    bl_description = "(Note! Operator just for debugging.) Makes lists of contiguous seams (which way is contiguous at crossroad seams will be decided according to which next seam is most well-aligned with the last) and two lists of face corners belonging on either side of the seam"
    bl_options = {'REGISTER', 'UNDO'}

    iP = bpy.props.IntProperty
    fP = bpy.props.FloatProperty
    fvP = bpy.props.FloatVectorProperty
    bP = bpy.props.BoolProperty
    eP = bpy.props.EnumProperty
    bvP = bpy.props.BoolVectorProperty
    sP = bpy.props.StringProperty

    max_i = 999999999

    sel_seams_by_i : bP(
        name = "Select Seams by Index",
        description = "",
        default = False
        ) # type: ignore

    sel_i_start : iP(
        name = "Start Index",
        description = "",
        default = 0
        ) # type: ignore

    sel_i_end : iP(
        name = "End Index",
        description = "",
        default = 0
        ) # type: ignore

    add_debug_attributes : bP(
        name = "Generate Debug Attribute",
        description = "Generate color attributes to debug seam sides",
        default = False
        ) # type: ignore

    sort_by_length : bP(
        name = "Sort by Length",
        description = "Sort seams by their length (longest first)",
        default = True
        ) # type: ignore

    verbose : bP(
        name = "Verbose",
        description = "Print debug information to system console",
        default = False
        ) # type: ignore

    def execute(self, context):
        props = context.scene.sloppy_props
        in_bm = bmesh.from_edit_mesh(bpy.context.active_object.data)
        sms, lsr, lsl, smsoe, smsi, smsmb = props.find_contiguous_seams(in_bm, self.sel_seams_by_i, self.sel_i_start, self.sel_i_end, self.sort_by_length, self.verbose)

        # if self.verbose:
        #     print('Number of contiguous seams:', len(sms))
        #     for smi, sm in enumerate(sms):
        #         smis = []
        #         for sme in sm:
        #             smis.append(sme.index)
        #         open_ended = "(closed)"
        #         if smsoe[smi] == True:
        #             open_ended = "(open-ended)"
        #         if smsi[smi] == True:
        #             open_ended = "(isolated)"
        #         print('Seam',(smi), open_ended,'- length',len(smis),'=',smis)
        
        if self.add_debug_attributes:
            props.find_or_add_attribute("seam_side", "FLOAT_COLOR", "CORNER")
            props.find_or_add_attribute("normalized_seam_index", "FLOAT", "EDGE")
            props.find_or_add_attribute("normalized_seam_sequence", "FLOAT", "EDGE")
            seam_side = props.get_attribute_layer("seam_side", "FLOAT_COLOR", "CORNER", in_bm)
            normalized_seam_index = props.get_attribute_layer("normalized_seam_index", "FLOAT", "EDGE", in_bm)
            normalized_seam_sequence = props.get_attribute_layer("normalized_seam_sequence", "FLOAT", "EDGE", in_bm)

            colr = mathutils.Color((0,1,0))
            coll = mathutils.Color((1,0,0))

            for smsi, sm in enumerate(sms):
                nor_sm_i = 0
                if len(sms) > 1:
                    nor_sm_i = smsi/(len(sms) - 1)
                for smi, sme in enumerate(sm):
                    sme[normalized_seam_index] = nor_sm_i
                    nor_sm_seq = 1.0
                    if len(sm) > 1:
                        nor_sm_seq = smi/(len(sm) - 1)
                    sme[normalized_seam_sequence] = nor_sm_seq

            for lrs, lls in zip(lsr, lsl):
                for lrl, lll in zip(lrs, lls):
                    lrl[seam_side].x = colr.r
                    lrl[seam_side].y = colr.g
                    lrl[seam_side].z = colr.b
                    lll[seam_side].x = coll.r
                    lll[seam_side].y = coll.g
                    lll[seam_side].z = coll.b

            bpy.context.active_object.data.update()

        return {"FINISHED"}

class SloppyFindIslandBoundaries(bpy.types.Operator):
    bl_idname = "operator.find_island_boundaries"
    bl_label = "Find Island Boundaries"
    bl_description = "Returns lists (vertex chains, face corner subchains, edge to next vertex in vertex chains, complete face corner chains, complete edge chains, complete vert chains, complete face chains) of contiguous boundaries of input island/collection of faces"
    bl_options = {'REGISTER', 'UNDO'}

    iP = bpy.props.IntProperty
    fP = bpy.props.FloatProperty
    fvP = bpy.props.FloatVectorProperty
    bP = bpy.props.BoolProperty
    eP = bpy.props.EnumProperty
    bvP = bpy.props.BoolVectorProperty
    sP = bpy.props.StringProperty

    max_i = 999999999

    island_selection_mode: eP(
        name = "Island Selection Mode",
        description = "Mode of selecting which island to do",
        items = [
            ("A", "Selection As Island", ""),
            ("B", "From Selected Face", ""),
            ("C", "Index", "")
            ],
        default="B"
        ) # type: ignore

    island_index : iP(
        name = "Island Index",
        description = "",
        default = 0
        ) # type: ignore

    add_debug_attributes : bP(
        name = "Generate Debug Attribute",
        description = "Generate color attributes to debug seam sides",
        default = True
        ) # type: ignore

    verbose : bP(
        name = "Verbose",
        description = "Print debug information to system console",
        default = True
        ) # type: ignore

    def execute(self, context):
        props = context.scene.sloppy_props
        in_bm = bmesh.from_edit_mesh(bpy.context.active_object.data)
        uv_layer = in_bm.loops.layers.uv.verify()
        islands = bmesh_utils.bmesh_linked_uv_islands(in_bm, uv_layer)

        island_index_clamped = int(bl_math.clamp(float(self.island_index), 0, len(islands) - 1))

        island = islands[island_index_clamped]
        if self.island_selection_mode == "A":
            island = []
            for face in in_bm.faces:
                if face.select == True:
                    island.append(face)
        if self.island_selection_mode == "B":
            for ii, isl in enumerate(islands):
                for face in isl:
                    if face.select == True:
                        island = islands[ii]
                        break
        # vert_chains, vert_loop_subchains, vert_chains_eton, complete_loop_chains, complete_edge_chains, complete_face_chains
        lcs, lcetons, vcs, ecs, fcs = props.find_island_boundaries_and_their_loops(in_bm, island, True, self.verbose)
        
        if self.add_debug_attributes:
            props.find_or_add_attribute("normalized_loop_sequence", "FLOAT_COLOR", "CORNER")
            normalized_loop_sequence = props.get_attribute_layer("normalized_loop_sequence", "FLOAT_COLOR", "CORNER", in_bm)

            color_init = mathutils.Color((1,0.33,0.33))
            for iface in island:
                for ifl in iface.loops:
                    ifl[normalized_loop_sequence].x = 0.0
                    ifl[normalized_loop_sequence].y = 0.0
                    ifl[normalized_loop_sequence].z = 0.0

            for lci,lc in enumerate(lcs):
                this_base_color = color_init.copy()
                this_base_color.h = lci/len(lcs)
                for cli, cl in enumerate(lc):
                    this_cl_color = this_base_color.copy()
                    # this_cl_color.v = (cli/len(lc) * 0.5) + 0.5
                    cl[normalized_loop_sequence].x = this_cl_color.r
                    cl[normalized_loop_sequence].y = this_cl_color.g
                    cl[normalized_loop_sequence].z = this_cl_color.b

        bpy.context.active_object.data.update()

        return {"FINISHED"}

class SloppyAlignViewToSelected(bpy.types.Operator):
    bl_idname = "operator.align_view_to_selection"
    bl_label = "Align View to Selected"
    bl_description = "Points view at averaged selection, aligned against averaged normals"
    bl_options = {'REGISTER', 'UNDO'}

    iP = bpy.props.IntProperty
    fP = bpy.props.FloatProperty
    fvP = bpy.props.FloatVectorProperty
    bP = bpy.props.BoolProperty
    eP = bpy.props.EnumProperty
    bvP = bpy.props.BoolVectorProperty
    sP = bpy.props.StringProperty

    def execute(self, context):
        props = context.scene.sloppy_props
        mode = bpy.context.mode
        cur_reg, cur_view = props.get_view3d_region_and_space()
        
        if mode == 'OBJECT':
            objects = [obj for obj in bpy.context.scene.objects if obj.select_get() == True]
            avg_loc = mathutils.Vector()
            for obj in objects:
                avg_loc += obj.location
            avg_loc /= len(objects)
            props.set_view3d_lookat(avg_loc)
        if mode == 'EDIT_MESH':
            objects = []
            bms = []
            for obj in bpy.context.objects_in_mode:
                objects.append(obj)
                bm = bmesh.from_edit_mesh(obj.data)
                bms.append(bm)
            avg_obj_loc = mathutils.Vector()
            avg_loc = mathutils.Vector()
            num_elements = 0
            num_objects = 0
            avg_nor = mathutils.Vector()
            for ob, obm in zip(objects, bms):
                has_selection = False
                for v in obm.verts:
                    if v.select == True:
                        has_selection = True
                        num_elements += 1
                        avg_loc += v.co
                        avg_nor += v.normal
                for e in obm.edges:
                    if e.select == True:
                        has_selection = True
                        num_elements += 1
                        avg_loc += props.calc_edge_center(e)
                        avg_nor += props.calc_edge_avg_normal(e)
                for f in obm.faces:
                    if f.select == True:
                        has_selection = True
                        num_elements += 1
                        avg_loc += f.calc_center_median_weighted()
                        avg_nor += f.normal
                if has_selection:
                    avg_obj_loc += ob.location
                    num_objects += 1
            avg_obj_loc /= num_objects
            avg_loc /= num_elements
            nor = avg_nor.normalized()
            props.align_view3d_against_normal(nor, cur_view)
            props.set_view3d_lookat(avg_loc + avg_obj_loc)
        return {"FINISHED"}
    
class SloppyFindMidpoints(bpy.types.Operator):
    bl_idname = "operator.find_midpoints"
    bl_label = "Find Mesh Midpoints"
    bl_description = ""
    bl_options = {'REGISTER', 'UNDO'}

    iP = bpy.props.IntProperty
    fP = bpy.props.FloatProperty
    fvP = bpy.props.FloatVectorProperty
    bP = bpy.props.BoolProperty
    eP = bpy.props.EnumProperty
    bvP = bpy.props.BoolVectorProperty
    sP = bpy.props.StringProperty

    def fmp_dist_sort(self, two_verts_lst):
        sort_dist = math.dist(two_verts_lst[0].co, two_verts_lst[1].co)
        return sort_dist

    def execute(self, context):
        props = context.scene.sloppy_props
        bm = bmesh.from_edit_mesh(bpy.context.active_object.data)
        max_dimension = max(bpy.context.object.dimensions)
        
        verts_to_do = []

        for vert in bm.verts:
            verts_to_do.append(vert)

        # for i, tv in enumerate(verts_to_do):
        #     progress_pct = 100*(i/len(verts_to_do))
        #     progress_str = str(progress_pct) + '%'
        #     print('Processing vertex', tv.index, 'of', len(verts_to_do), '-', progress_str, 'done.')
        #     vits = 0
        #     totvec = mathutils.Vector()
        #     for ov in verts_to_do:
        #         if ov != tv:
        #            odott = -ov.normal.dot(-tv.normal)
        #            if odott < -0.25:
        #                 vits += 1
        #                 midco = tv.co.lerp(ov.co, 0.5)
        #                 totvec += midco
        #     if vits >= 1:
        #         totvec /= vits
        #         bmesh.ops.create_vert(bm, co=totvec)

        for i, tv in enumerate(verts_to_do):
            progress_pct = 100*(i/len(verts_to_do))
            progress_str = str(progress_pct) + '%'
            print('Processing vertex', tv.index, 'of', len(verts_to_do), '-', progress_str, 'done.')
            overts = [] 
            for ov in verts_to_do:
                if ov != tv:
                    overts.append([tv, ov, max_dimension])
            overts.sort(key=props.sort_verts_by_nearest_most_opposite)
            midco = tv.co.lerp(overts[0][1].co, 0.5)
            bmesh.ops.create_vert(bm, co=midco)

        bpy.context.active_object.data.update()

        new_verts = []
        done_verts = []

        for av in bm.verts:
            if av not in verts_to_do:
                new_verts.append(av)
        
        for nv in new_verts:
            noverts = []
            for nov in new_verts:
                if nov != nv:
                    if nov not in done_verts:
                        noverts.append([nv, nov])
            noverts.sort(key=self.fmp_dist_sort)
            bmesh.ops.contextual_create(bm, geom=[noverts[0][0], noverts[0][1]])
            done_verts.append(noverts[0][1])
        
        bpy.context.active_object.data.update()

        return {"FINISHED"}

class SloppyBatchRender(bpy.types.Operator):
    bl_idname = "operator.sloppy_batch_render"
    bl_label = "Batch Render"
    bl_description = ""
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        props = context.scene.sloppy_props

        prev_render_path = bpy.context.scene.render.filepath

        viewport_render = props.br_use_viewport
        
        coll_parent = bpy.context.collection
        name_parent = coll_parent.name
        lyrcoll_parent = bpy.context.scene.view_layers['ViewLayer'].layer_collection.children[name_parent]

        for sub_coll in coll_parent.children:
            print(sub_coll.name)
            if 'Collection' in str(sub_coll):
                name_sub = sub_coll.name
                lyrcoll_sub = None
                for lc in lyrcoll_parent.children:
                    if lc.name == name_sub:
                        lc.exclude = False
                        lyrcoll_sub = lc
                        print('Subcollection:', lc.name)
                    else:
                        lc.exclude = True
                for model_coll in sub_coll.children:
                    if 'Collection' in str(model_coll):
                        name_model = model_coll.name
                        print('Subcollection:', name_sub + ' > ' + 'Model:', name_model)
                        lyrcoll_model = None
                        for lcm in lyrcoll_sub.children:
                            if lcm.name == name_model:
                                lcm.exclude = False
                                lyrcoll_model = lcm
                            else:
                                lcm.exclude = True
                        render_path = props.br_output_path + '\\' + sub_coll.name + '\\' + model_coll.name + '_'
                        if viewport_render == True:
                            render_path = props.br_output_path + '\\' + sub_coll.name + '\\' + model_coll.name + '_VP_'
                        bpy.context.scene.render.filepath = render_path
                        if viewport_render == True:
                            bpy.ops.render.opengl(animation=True)
                        else:
                            bpy.ops.render.render(animation=True, use_viewport=True)
                        
            
        print('\nAll renders done!')
        bpy.context.scene.render.filepath = prev_render_path

        return {"FINISHED"}

class SloppyCopyFloatToClipboard(bpy.types.Operator):
    bl_idname = "operator.sloppy_copy_to_clipboard"
    bl_label = "Copy Value"
    bl_description = ""
    bl_options = {'REGISTER', 'UNDO'}
    
    fP = bpy.props.FloatProperty

    copied_value : fP(
        name="Copied Value",
        precision=10
        ) # type: ignore

    def execute(self, context):
        props = context.scene.sloppy_props

        props.copy_float_to_clipboard(self.copied_value)

        return {"FINISHED"}

class SloppyCopyLocation(bpy.types.Operator):
    bl_idname = "operator.sloppy_copy_location"
    bl_label = "Copy Location"
    bl_description = ""
    bl_options = {'REGISTER', 'UNDO'}
    
    def execute(self, context):
        props = context.scene.sloppy_props

        props.cpt_location = bpy.context.object.location
        print(props.cpt_location)

        return {"FINISHED"}

class SloppyCopyRotation(bpy.types.Operator):
    bl_idname = "operator.sloppy_copy_rotation"
    bl_label = "Copy Rotation"
    bl_description = ""
    bl_options = {'REGISTER', 'UNDO'}
    
    def execute(self, context):
        props = context.scene.sloppy_props

        if bpy.context.object.rotation_mode == 'QUATERNION':
            quat = bpy.context.object.rotation_quaternion
            props.cpt_rotation = mathutils.Vector((quat[0], quat[1], quat[2], quat[3]))
            print(quat)
        elif bpy.context.object.rotation_mode == 'AXIS_ANGLE':
            axisangle = bpy.context.object.rotation_axis_angle
            rot = mathutils.Matrix.Rotation(axisangle[0], 4, (axisangle[1], axisangle[2], axisangle[3]))
            quat = rot.to_quaternion()
            props.cpt_rotation = mathutils.Vector((quat[0], quat[1], quat[2], quat[3]))
            print('Axis: X:', axisangle[1], ', Y:', axisangle[2], ', Z:', axisangle[3], '| Angle:', math.degrees(axisangle[0]),'->', quat)
        else:
            quat = bpy.context.object.rotation_euler.to_quaternion()
            props.cpt_rotation = mathutils.Vector((quat[0], quat[1], quat[2], quat[3]))
            print(bpy.context.object.rotation_euler, '->', quat)

        return {"FINISHED"}

class SloppyCopyScale(bpy.types.Operator):
    bl_idname = "operator.sloppy_copy_scale"
    bl_label = "Copy Scale"
    bl_description = ""
    bl_options = {'REGISTER', 'UNDO'}
    
    def execute(self, context):
        props = context.scene.sloppy_props

        props.cpt_scale = bpy.context.object.scale
        print(props.cpt_scale)

        return {"FINISHED"}

class SloppyCopyDimensions(bpy.types.Operator):
    bl_idname = "operator.sloppy_copy_dimensions"
    bl_label = "Copy Dimensions"
    bl_description = ""
    bl_options = {'REGISTER', 'UNDO'}
    
    def execute(self, context):
        props = context.scene.sloppy_props

        props.cpt_dimensions = bpy.context.object.dimensions
        print(props.cpt_dimensions)

        return {"FINISHED"}

class SloppyCopyTransform(bpy.types.Operator):
    bl_idname = "operator.sloppy_copy_transform"
    bl_label = "Copy Transform"
    bl_description = ""
    bl_options = {'REGISTER', 'UNDO'}
    
    def execute(self, context):
        props = context.scene.sloppy_props

        props.cpt_location = bpy.context.object.location
        print(props.cpt_location)

        if bpy.context.object.rotation_mode == 'QUATERNION':
            quat = bpy.context.object.rotation_quaternion
            props.cpt_rotation = mathutils.Vector((quat[0], quat[1], quat[2], quat[3]))
            print(quat)
        elif bpy.context.object.rotation_mode == 'AXIS_ANGLE':
            axisangle = bpy.context.object.rotation_axis_angle
            rot = mathutils.Matrix.Rotation(axisangle[0], 4, (axisangle[1], axisangle[2], axisangle[3]))
            quat = rot.to_quaternion()
            props.cpt_rotation = mathutils.Vector((quat[0], quat[1], quat[2], quat[3]))
            print('Axis: X:', axisangle[1], ', Y:', axisangle[2], ', Z:', axisangle[3], '| Angle:', math.degrees(axisangle[0]),'->', quat)
        else:
            quat = bpy.context.object.rotation_euler.to_quaternion()
            props.cpt_rotation = mathutils.Vector((quat[0], quat[1], quat[2], quat[3]))
            print(bpy.context.object.rotation_euler, '->', quat)

        props.cpt_scale = bpy.context.object.scale
        print(props.cpt_scale)

        return {"FINISHED"}
    
class SloppyPasteLocation(bpy.types.Operator):
    bl_idname = "operator.sloppy_paste_location"
    bl_label = "Paste Location"
    bl_description = ""
    bl_options = {'REGISTER', 'UNDO'}
    
    def execute(self, context):
        props = context.scene.sloppy_props

        bpy.context.object.location = props.cpt_location

        return {"FINISHED"}

class SloppyPasteRotation(bpy.types.Operator):
    bl_idname = "operator.sloppy_paste_rotation"
    bl_label = "Paste Rotation"
    bl_description = ""
    bl_options = {'REGISTER', 'UNDO'}
    
    def execute(self, context):
        props = context.scene.sloppy_props

        quat = mathutils.Quaternion((props.cpt_rotation[0], props.cpt_rotation[1], props.cpt_rotation[2], props.cpt_rotation[3]))

        if bpy.context.object.rotation_mode == 'QUATERNION':
            bpy.context.object.rotation_quaternion = quat
            print(quat)
        elif bpy.context.object.rotation_mode == 'AXIS_ANGLE':
            axisangle = quat.to_axis_angle()
            print(quat, '-->', axisangle)
            bpy.context.object.rotation_axis_angle[0] = axisangle[1]
            bpy.context.object.rotation_axis_angle[1] = axisangle[0][0]
            bpy.context.object.rotation_axis_angle[1] = axisangle[0][1]
            bpy.context.object.rotation_axis_angle[1] = axisangle[0][2]
        else:
            bpy.context.object.rotation_euler = quat.to_euler()
            print(quat, '-->', quat.to_euler())

        return {"FINISHED"}

class SloppyPasteScale(bpy.types.Operator):
    bl_idname = "operator.sloppy_paste_scale"
    bl_label = "Paste Scale"
    bl_description = ""
    bl_options = {'REGISTER', 'UNDO'}
    
    def execute(self, context):
        props = context.scene.sloppy_props

        bpy.context.object.scale = props.cpt_scale

        return {"FINISHED"}

class SloppyPasteDimensions(bpy.types.Operator):
    bl_idname = "operator.sloppy_paste_dimensions"
    bl_label = "Paste Dimensions"
    bl_description = ""
    bl_options = {'REGISTER', 'UNDO'}
    
    def execute(self, context):
        props = context.scene.sloppy_props

        bpy.context.object.dimensions = props.cpt_dimensions

        return {"FINISHED"}

class SloppyPasteTransform(bpy.types.Operator):
    bl_idname = "operator.sloppy_paste_transform"
    bl_label = "Paste Transform"
    bl_description = ""
    bl_options = {'REGISTER', 'UNDO'}
    
    def execute(self, context):
        props = context.scene.sloppy_props

        bpy.context.object.location = props.cpt_location

        quat = mathutils.Quaternion((props.cpt_rotation[0], props.cpt_rotation[1], props.cpt_rotation[2], props.cpt_rotation[3]))

        if bpy.context.object.rotation_mode == 'QUATERNION':
            bpy.context.object.rotation_quaternion = quat
            print(quat)
        elif bpy.context.object.rotation_mode == 'AXIS_ANGLE':
            axisangle = quat.to_axis_angle()
            print(quat, '-->', axisangle)
            bpy.context.object.rotation_axis_angle[0] = axisangle[1]
            bpy.context.object.rotation_axis_angle[1] = axisangle[0][0]
            bpy.context.object.rotation_axis_angle[1] = axisangle[0][1]
            bpy.context.object.rotation_axis_angle[1] = axisangle[0][2]
        else:
            bpy.context.object.rotation_euler = quat.to_euler()
            print(quat, '-->', quat.to_euler())

        bpy.context.object.scale = props.cpt_scale

        return {"FINISHED"}

#region Initialization
classes = [SloppyProperties,
           SloppyErrorDialog,
           SloppyUVPanel,
           SloppyUVAlignPanel,
           SloppyUVScalePanel,
           SloppyPeltPanel,
           SloppyDebugPanel,
           PeltUVs,
           SloppyCalcScale,
           SloppyApplyScale,
           SloppyAlignUVsGeo,
           SloppyAlignUVsUV,
           SloppyAlignUVs,
           SloppySeamGen,
           SloppySeamGenVis,
           SloppySeamGenPanel,
           SortVertByDist,
           SortEdgeByDist,
           SortFaceByDist,
           SloppySortDistPanel,
           SloppyProcAttrBake,
           SloppyBlurAttribute,
           SloppySelectByIndex,
           SloppyShiftSelectByIndex,
           SloppyQuadUVUnfold,
           SloppyFlatQuadPanel,
           SloppyUVToMesh,
           RedoUVEdgeLength,
           SloppyAlignViewToSelected,
           SloppyFindContiguousSeams,
           SloppyFindIslandBoundaries,
           SloppyFindMidpoints,
           SloppyIslandBoundaryConnect,
           SloppyBasicUVUnfold,
           SloppyBoundaryFirstUVUnfold,
           SloppyUVSmartAlign,
           SloppyBatchRender,
           SloppyBRPanel,
           SloppyExtraObjPropsPanel,
           SloppyCopyFloatToClipboard,
           SloppyCopyLocation,
           SloppyCopyRotation,
           SloppyCopyScale,
           SloppyCopyDimensions,
           SloppyCopyTransform,
           SloppyPasteLocation,
           SloppyPasteRotation,
           SloppyPasteScale,
           SloppyPasteDimensions,
           SloppyPasteTransform,
           SloppyCPTPanel
           ]

def register():
    for cls in classes:
        try: 
            bpy.utils.register_class(cls)
        except:
            bpy.utils.unregister_class(cls)
            bpy.utils.register_class(cls)
    bpy.types.Scene.sloppy_props = bpy.props.PointerProperty(type = SloppyProperties)

        
def unregister():
    for cls in classes:
        bpy.utils.unregister_class(cls)
    del bpy.types.Scene.sloppy_props
        
if __name__ == "__main__":
    register()
#endregion