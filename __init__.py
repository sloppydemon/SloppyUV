import bpy
import bpy_extras
import bmesh
import math
import mathutils
import bl_math
import sys

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

from SloppyUV.SloppySeamGeneration import SloppySeamGen # type: ignore
from SloppyUV.SloppySortIndexByDist import SortVertByDist # type: ignore
from SloppyUV.SloppySortIndexByDist import SortEdgeByDist # type: ignore
from SloppyUV.SloppySortIndexByDist import SortFaceByDist # type: ignore
from SloppyUV.SloppyBakeProcAttr import SloppyProcAttrBake # type: ignore
from SloppyUV.SloppyProcedurals import SloppyQuadUVUnfold # type: ignore

class SloppyProperties(bpy.types.PropertyGroup):
    
    iP = bpy.props.IntProperty
    fP = bpy.props.FloatProperty
    fvP = bpy.props.FloatVectorProperty
    bP = bpy.props.BoolProperty
    eP = bpy.props.EnumProperty
    bvP = bpy.props.BoolVectorProperty
    sP = bpy.props.StringProperty
    
    def viewco(self, vec):
        '''Fetches (TODO: active or) any View3D region and space and outputs the input vector mapped
        to that 3D view.'''
        view3d_region = None
        view3d = None
        for scr in bpy.data.screens:
            for ar in scr.areas:
                if ar.type == 'VIEW_3D':
                    view3d_region = ar.regions[0]
                    for spc in ar.spaces:
                        if spc.type == 'VIEW_3D':
                            view3d = spc.region_3d
                            break
        outco = bpy_extras.view3d_utils.location_3d_to_region_2d(view3d_region, view3d, vec)
        return outco

    def viewvec(self, loc):
        '''Fetch vectpr pointing from 3D View origin to point on screen'''
        view3d_region = None
        view3d = None
        for scr in bpy.data.screens:
            for ar in scr.areas:
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

    def remap_val(self, val, in_min, in_max, out_min, out_max):
        in_interval = in_max - in_min
        out_interval = out_max - out_min
        in_val = val - in_min
        in_fac = 0
        if in_interval < 0:
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
    
    current_sort_axis : sP(
        name = "Current Sort Axis",
        default = "",
        ) # type: ignore

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
        seamgen_grid = seamgen_box.grid_flow(columns=2, even_columns=True)

        seamgen_props = seamgen_grid.operator("operator.sloppy_seam_gen")
        sg_col_opt = seamgen_grid.column(heading="Seam Generation Options:")
        sg_col_opt.prop(seamgen_props, "angle_fac")
        sg_col_opt.prop(seamgen_props, "avg_angle_fac")
        sg_col_opt.prop(seamgen_props, "ao_factor")
        sg_col_opt.prop(seamgen_props, "ed_factor")
        sg_col_opt.prop(seamgen_props, "no_rounds")
        sg_col_opt.prop(seamgen_props, "no_retries")
        sg_col_opt.prop(seamgen_props, "angle_thresh_start")
        sg_col_opt.prop(seamgen_props, "angle_thresh_end")
        sg_col_opt.prop(seamgen_props, "clear_seam")
        sg_col_opt.prop(seamgen_props, "unwrap")

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

        fq_props = fq_box.operator("operator.sloppy_quad_unfold")
        fq_col_opt = fq_box.column(heading="Virtual Quad Unwrap Options:")
        fq_col_opt.prop(fq_props, "unfold_mode")
        fq_col_opt.prop(fq_props, "initial_quad")
        fq_col_opt.prop(fq_props, "reg_flat_fac")
        fq_col_opt.prop(fq_props, "reg_norm_fac")
        fq_col_opt.prop(fq_props, "pre_calc_edge_lengths")
        fq_col_opt.prop(fq_props, "offset_per_island")
        fq_col_opt.prop(fq_props, "quant_avg_norm")
        fq_col_opt.prop(fq_props, "only_move_loops_in_face")

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
           SloppySeamGenPanel,
           SortVertByDist,
           SortEdgeByDist,
           SortFaceByDist,
           SloppySortDistPanel,
           SloppyProcAttrBake,
           SloppySelectByIndex,
           SloppyShiftSelectByIndex,
           SloppyQuadUVUnfold,
           SloppyFlatQuadPanel
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