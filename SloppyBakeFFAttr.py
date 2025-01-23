import bpy
import bmesh
import mathutils
import math
import bl_math
from bpy_extras import bmesh_utils

props = bpy.context.scene.sloppy_props

depsgraph = bpy.context.evaluated_depsgraph_get()

per_vert = False
use_ray = False

x_res = 8
y_res = 8
z_res = 18

min_x = bpy.context.object.bound_box[3][0]
max_x = bpy.context.object.bound_box[4][0]
dim_x = max_x - min_x
min_y = bpy.context.object.bound_box[5][1]
max_y = bpy.context.object.bound_box[6][1]
dim_y = max_y - min_y
min_z = bpy.context.object.bound_box[0][2]
max_z = bpy.context.object.bound_box[1][2]
dim_z = max_z - min_z

corners = [
    [
        mathutils.Vector((0, min_y, min_z)),
        mathutils.Vector((0, max_y, min_z)),
        mathutils.Vector((0, max_y, max_z)),
        mathutils.Vector((0, min_y, max_z))
    ],
    [
        mathutils.Vector((min_x, 0, min_z)),
        mathutils.Vector((max_x, 0, min_z)),
        mathutils.Vector((max_x, 0, max_z)),
        mathutils.Vector((min_x, 0, max_z))
    ],
    [
        mathutils.Vector((min_x, min_y, 0)),
        mathutils.Vector((max_x, min_y, 0)),
        mathutils.Vector((max_x, max_y, 0)),
        mathutils.Vector((min_x, max_y, 0))
    ]
]

uv_corners = [
    [
        0.0,0.0
    ],
    [
        1.0,0.0
    ],
    [
        1.0,1.0
    ],
    [
        0.0,1.0
    ]
]

mins = [
    min_x,
    min_y,
    min_z
]
maxs = [
    max_x,
    max_y,
    max_z
]

x_inc = dim_x / x_res
y_inc = dim_y / y_res
z_inc = dim_z / z_res

incs = [
    x_inc,
    y_inc,
    z_inc
]

lvls_x = [(min_x + i*x_inc) for i in range(x_res + 1)]
lvls_y = [(min_y + i*y_inc) for i in range(y_res + 1)]
lvls_z = [(min_z + i*z_inc) for i in range(z_res + 1)]

lvls_norm_x = [(1/x_res)*i for i in range(x_res + 1)]
lvls_norm_y = [(1/y_res)*i for i in range(y_res + 1)]
lvls_norm_z = [(1/z_res)*i for i in range(z_res + 1)]

lvls = [
        lvls_x,
        lvls_y,
        lvls_z
    ]

lvlns = [
        lvls_norm_x,
        lvls_norm_y,
        lvls_norm_z
    ]


x_div = x_res - 1
y_div = y_res - 1
z_div = z_res - 1

mo = bpy.context.object
bake_plane = bpy.context.scene.objects['Plane']

olax = [
    'x',
    'y',
    'z'
]

axos = [None, None, None]
vaxos = [None, None, None]
bake_planes = [None, None, None]

for axi, axis in enumerate(olax):
    this_axis_o = mo.copy()
    this_axis_o.data = mo.data.copy()
    this_axis_o.name = mo.name + '_lvls.' + axis
    bpy.context.collection.objects.link(this_axis_o)
    this_axis_o.select_set(True)
    axos[axi] = this_axis_o

    this_vaxis_o = mo.copy()
    this_vaxis_o.data = mo.data.copy()
    this_vaxis_o.name = mo.name + '_vol.' + axis
    bpy.context.collection.objects.link(this_vaxis_o)
    this_vaxis_o.select_set(True)
    vaxos[axi] = this_vaxis_o

    this_bake_plane = bake_plane.copy()
    this_bake_plane.data = bake_plane.data.copy()
    this_bake_plane.name = mo.name + '_bakePlane.' + axis
    bpy.context.collection.objects.link(this_bake_plane)
    this_bake_plane.select_set(True)
    bake_planes[axi] = this_bake_plane

lvl_cuts = [[],[],[]]

for lvli, lvlsa in enumerate(lvls):
    for lvlvi, lvlv in enumerate(lvlsa):
        this_cut_o = mo.copy()
        this_cut_o.data = mo.data.copy()
        this_cut_o.name = mo.name + '_cut.' + olax[lvli] + str(lvlvi)
        bpy.context.collection.objects.link(this_cut_o)
        props.find_or_add_attribute_other_obj('xmid', 'FLOAT_VECTOR', 'POINT', this_cut_o)
        props.find_or_add_attribute_other_obj('ymid', 'FLOAT_VECTOR', 'POINT', this_cut_o)
        props.find_or_add_attribute_other_obj('zmid', 'FLOAT_VECTOR', 'POINT', this_cut_o)
        props.find_or_add_attribute_other_obj('norm_xmid', 'FLOAT_VECTOR', 'POINT', this_cut_o)
        props.find_or_add_attribute_other_obj('norm_ymid', 'FLOAT_VECTOR', 'POINT', this_cut_o)
        props.find_or_add_attribute_other_obj('norm_zmid', 'FLOAT_VECTOR', 'POINT', this_cut_o)
        props.find_or_add_attribute_other_obj('dist_x', 'FLOAT', 'POINT', this_cut_o)
        props.find_or_add_attribute_other_obj('dist_y', 'FLOAT', 'POINT', this_cut_o)
        props.find_or_add_attribute_other_obj('dist_z', 'FLOAT', 'POINT', this_cut_o)
        props.find_or_add_attribute_other_obj('norm_pos_x', 'FLOAT', 'POINT', this_cut_o)
        props.find_or_add_attribute_other_obj('norm_pos_y', 'FLOAT', 'POINT', this_cut_o)
        props.find_or_add_attribute_other_obj('norm_pos_z', 'FLOAT', 'POINT', this_cut_o)
        props.find_or_add_attribute_other_obj('norm_dist_x', 'FLOAT', 'POINT', this_cut_o)
        props.find_or_add_attribute_other_obj('norm_dist_y', 'FLOAT', 'POINT', this_cut_o)
        props.find_or_add_attribute_other_obj('norm_dist_z', 'FLOAT', 'POINT', this_cut_o)
        this_cut_o.select_set(True)
        lvl_cuts[lvli].append(this_cut_o)

bpy.ops.object.editmode_toggle()

bmo = bmesh.from_edit_mesh(mo.data)
bmo_vol = bmo.calc_volume()

bmox = bmesh.from_edit_mesh(axos[0].data)
bmoy = bmesh.from_edit_mesh(axos[1].data)
bmoz = bmesh.from_edit_mesh(axos[2].data)

bmvox = bmesh.from_edit_mesh(vaxos[0].data)
bmvoy = bmesh.from_edit_mesh(vaxos[1].data)
bmvoz = bmesh.from_edit_mesh(vaxos[2].data)

ol = [
    bmox,
    bmoy,
    bmoz
]

olv = [
    bmvox,
    bmvoy,
    bmvoz
]

cutbms = [[],[],[]]

for lvl_cut_i, lvl_cut_a in enumerate(lvl_cuts):
    for cut_o in lvl_cut_a:
        lvl_cut_bm = bmesh.from_edit_mesh(cut_o.data)
        cutbms[lvl_cut_i].append(lvl_cut_bm)

bpbms = []

for bp in  bake_planes:
    bake_plane_bm = bmesh.from_edit_mesh(bp.data)
    bpbms.append(bake_plane_bm)

olcovm = [
    mathutils.Vector((1,0,0)),
    mathutils.Vector((0,1,0)),
    mathutils.Vector((0,0,1))
]

current_vm = mathutils.Vector((1,0,0))

def vec_sort_x(e):
    return e[1].x
def vec_sort_y(e):
    return e[1].y
def vec_sort_z(e):
    return e[1].z

def vol_calc(this_volume, full_volume):
    if this_volume > 0:
        vol_fac = this_volume / full_volume
        return vol_fac
    else:
        return 0.0

attr_dict = [
    {"name": "volz", "type": "FLOAT", "domain": "POINT", "layer": None},
    {"name": "voly", "type": "FLOAT", "domain": "POINT", "layer": None},
    {"name": "volx", "type": "FLOAT", "domain": "POINT", "layer": None},
    {"name": "zmid", "type": "FLOAT_VECTOR", "domain": "POINT", "layer": None},
    {"name": "xmid", "type": "FLOAT_VECTOR", "domain": "POINT", "layer": None},
    {"name": "ymid", "type": "FLOAT_VECTOR", "domain": "POINT", "layer": None},
    {"name": "dist_x", "type": "FLOAT", "domain": "POINT", "layer": None},
    {"name": "dist_y", "type": "FLOAT", "domain": "POINT", "layer": None},
    {"name": "dist_z", "type": "FLOAT", "domain": "POINT", "layer": None},
    {"name": "norm_pos_x", "type": "FLOAT", "domain": "POINT", "layer": None},
    {"name": "norm_pos_y", "type": "FLOAT", "domain": "POINT", "layer": None},
    {"name": "norm_pos_z", "type": "FLOAT", "domain": "POINT", "layer": None},
    {"name": "orig_co", "type": "FLOAT_VECTOR", "domain": "POINT", "layer": None}
    ]

for nam in attr_dict:
    attr = props.find_or_add_attribute(nam["name"], nam["type"], nam["domain"])
    nam["layer"] = props.get_attribute_layer(nam["name"], nam["type"], nam["domain"], bmo)

volx = props.get_dict_layer("volx", attr_dict)
voly = props.get_dict_layer("voly", attr_dict)
volz = props.get_dict_layer("volz", attr_dict)
zmid = props.get_dict_layer("zmid", attr_dict)
xmid = props.get_dict_layer("xmid", attr_dict)
ymid = props.get_dict_layer("ymid", attr_dict)
dist_x = props.get_dict_layer("dist_x", attr_dict)
dist_y = props.get_dict_layer("dist_y", attr_dict)
dist_z = props.get_dict_layer("dist_z", attr_dict)
norm_pos_x = props.get_dict_layer("norm_pos_x", attr_dict)
norm_pos_y = props.get_dict_layer("norm_pos_y", attr_dict)
norm_pos_z = props.get_dict_layer("norm_pos_z", attr_dict)
orig_co = props.get_dict_layer("orig_co", attr_dict)

ollyr = [
    volx,
    voly,
    volz
]

olmid = [
    xmid,
    ymid,
    zmid
]

oldist = [
    dist_x,
    dist_y,
    dist_z
]

olnp = [
    norm_pos_x,
    norm_pos_y,
    norm_pos_z
]

bmo.verts.ensure_lookup_table()
bmo.edges.ensure_lookup_table()
bmo.faces.ensure_lookup_table()

bmo_flat_uv_xy = bmo.loops.layers.uv.new()
bmo_flat_uv_yz = bmo.loops.layers.uv.new()
bmo_flat_uv_zy = bmo.loops.layers.uv.new()

for bmov in bmo.verts:
    this_npx = (bmov.co.x - min_x) / dim_x
    this_npy = (bmov.co.y - min_y) / dim_y
    this_npz = (bmov.co.z - min_z) / dim_z
    bmov[norm_pos_x] = this_npx
    bmov[norm_pos_y] = this_npy
    bmov[norm_pos_z] = this_npz

    for loop in bmov.link_loops:
        loop[bmo_flat_uv_xy].uv.x = this_npx
        loop[bmo_flat_uv_xy].uv.y = this_npy
        loop[bmo_flat_uv_yz].uv.x = this_npy
        loop[bmo_flat_uv_yz].uv.y = this_npz
        loop[bmo_flat_uv_zy].uv.x = this_npz
        loop[bmo_flat_uv_zy].uv.y = this_npy

if per_vert == True:
    for obj,vec,lyr,mid,ax in zip(ol, olcovm, ollyr,olmid,olax):
        ls = []
        for v in obj.verts:
            lsi = [v, v.co, v.index]
            ls.append(lsi)
        current_vm = vec
        plane_vec = mathutils.Vector((1,1,1)) - vec
        if ax == 'x':
            ls.sort(key = vec_sort_x, reverse=True)
        if ax == 'y':
            ls.sort(key = vec_sort_y, reverse=True)
        if ax == 'z':
            ls.sort(key = vec_sort_z, reverse=True)
        last_mid = ls[0][1]

        for vert in ls:
            print('Before cut:', len(obj.verts), len([v for v in obj.verts if v.is_valid]), vert[1] * vec, 'axis: ', ax)
            vol = vol_calc(obj.calc_volume(), bmo_vol)
            bmo.verts[vert[2]][lyr] = vol
            cut, rest = bmesh.ops.bisect_plane(obj, geom=obj.verts[:] + obj.edges[:] + obj.faces[:], dist=0.0001, plane_co=vert[1] * vec, plane_no=vec, use_snap_center=False, clear_outer=True, clear_inner=False)
            new_faces = bmesh.ops.holes_fill(obj, edges=obj.edges)
            print('After cut:',len(obj.verts), len([v for v in obj.verts if v.is_valid]), 'Number of new faces: ', len(new_faces['faces']), 'New volume: ', vol)
        
        for bmov in bmo.verts:
            bmov[orig_co] = bmov.co
            inwards = -bmov.normal
            inwards_vec = inwards * plane_vec
            inwards_dir = inwards_vec.normalized()
            hit, hit_loc, hit_norm, hit_i, hit_o, hit_ab = bpy.context.scene.ray_cast(depsgraph, bmov.co + (inwards*0.01), inwards_dir)
            if hit:
                last_mid = bmov.co.lerp(hit_loc, 0.5)
            print('Hit:', hit, 'Midpoint for vert', bmov.index, 'in', ax, 'plane: ', 'x:', last_mid.x, 'y:', last_mid.y, 'z:', last_mid.z)

        bpy.context.active_object.data.update()
else:
    for obj,vobj, vec,lyr,mid,np,dst,ax,lvla,lvln,cbms,dmin,dmax,inc,bpbm,corna in zip(ol, olv, olcovm, ollyr,olmid,olnp,oldist,olax,lvls,lvlns,cutbms,mins,maxs,incs,bpbms,corners):
        current_vm = vec
        plane_vec = mathutils.Vector((1,1,1)) - vec
        
        ls = []
        for v in vobj.verts:
            lsi = [v, v.co, v.index]
            ls.append(lsi)
        if ax == 'x':
            ls.sort(key = vec_sort_x, reverse=True)
        if ax == 'y':
            ls.sort(key = vec_sort_y, reverse=True)
        if ax == 'z':
            ls.sort(key = vec_sort_z, reverse=True)

        last_mid = mathutils.Vector((0,0,0))

        for vert in ls:
            print('Before cut:', len(vobj.verts), len([v for v in vobj.verts if v.is_valid]), vert[1] * vec, 'axis: ', ax)
            cut, rest = bmesh.ops.bisect_plane(vobj, geom=vobj.verts[:] + vobj.edges[:] + vobj.faces[:], dist=0.0001, plane_co=vert[1] * vec, plane_no=vec, use_snap_center=False, clear_outer=True, clear_inner=False)
            new_faces = bmesh.ops.holes_fill(vobj, edges=vobj.edges)
            vol = vol_calc(vobj.calc_volume(), bmo_vol)
            bmo.verts[vert[2]][lyr] = vol
            print('After cut:',len(vobj.verts), len([v for v in vobj.verts if v.is_valid]), 'Number of new faces: ', len(new_faces['faces']), 'New volume: ', vol)
        
        for bmov in bmo.verts:
            bmov[orig_co] = bmov.co
        
        this_res = x_res
        if ax == 'y':
            this_res = y_res
        if ax == 'z':
            this_res = z_res
        max_dim = max(dim_x, dim_y, dim_z)

        last_mid = mathutils.Vector()

        vol_lvls = []

        lvla_rng = [i for i in range(len(lvla))]
        print(lvla_rng)
        lvla_rng.sort(reverse = True)

        lvl_vol_map = []

        for ri in lvla_rng:
            lvl = lvla[ri]
            cbm = cbms[ri]
            
            print('Before cut:', len(obj.verts), len([v for v in obj.verts if v.is_valid]), 'level: ', lvl, ' in axis: ', ax)
            
            plane_pos = vec * lvl
            cut, rest = bmesh.ops.bisect_plane(obj, geom=obj.verts[:] + obj.edges[:] + obj.faces[:], dist=0.0001, plane_co=plane_pos, plane_no=vec, use_snap_center=False, clear_outer=True, clear_inner=False)
            cut_cut, cut_rest = bmesh.ops.bisect_plane(cbm, geom=cbm.verts[:] + cbm.edges[:] + cbm.faces[:], dist=0.0001, plane_co=plane_pos, plane_no=vec, use_snap_center=False, clear_outer=True, clear_inner=True)
            if use_ray == False:
                cut_new_faces = bmesh.ops.holes_fill(cbm, edges=cbm.edges)
                new_faces = bmesh.ops.holes_fill(obj, edges=obj.edges)
            
            vol = vol_calc(obj.calc_volume(), bmo_vol)

            vol_lvls.append(vol)
            lvl_np = lvln[ri]
            lvl_vol_item = [lvl, vol, lvl_np]
            lvl_vol_map.append(lvl_vol_item)
            print('After cut:',len(obj.verts), len([v for v in obj.verts if v.is_valid]), 'Number of new faces: ', len(new_faces['faces']), 'New volume: ', vol)
        
        vol_lvls.sort(reverse=True)

        # obj.data.update()

        for cbm in cbms:
            this_mid_lyr = props.get_attribute_layer('xmid', 'FLOAT_VECTOR', 'POINT', cbm)
            this_nmid_lyr = props.get_attribute_layer('norm_xmid', 'FLOAT_VECTOR', 'POINT', cbm)
            this_dist_lyr = props.get_attribute_layer('dist_x', 'FLOAT', 'POINT', cbm)
            this_npx_lyr = props.get_attribute_layer('norm_pos_x', 'FLOAT', 'POINT', cbm)
            this_npy_lyr = props.get_attribute_layer('norm_pos_y', 'FLOAT', 'POINT', cbm)
            this_npz_lyr = props.get_attribute_layer('norm_pos_z', 'FLOAT', 'POINT', cbm)
            this_nd_lyr = props.get_attribute_layer('norm_dist_x', 'FLOAT', 'POINT', cbm)
            if ax == 'y':
                this_mid_lyr = props.get_attribute_layer('ymid', 'FLOAT_VECTOR', 'POINT', cbm)
                this_nmid_lyr = props.get_attribute_layer('norm_ymid', 'FLOAT_VECTOR', 'POINT', cbm)
                this_dist_lyr = props.get_attribute_layer('dist_y', 'FLOAT', 'POINT', cbm)
                this_nd_lyr = props.get_attribute_layer('norm_dist_y', 'FLOAT', 'POINT', cbm)
            if ax == 'z':
                this_mid_lyr = props.get_attribute_layer('zmid', 'FLOAT_VECTOR', 'POINT', cbm)
                this_nmid_lyr = props.get_attribute_layer('norm_zmid', 'FLOAT_VECTOR', 'POINT', cbm)
                this_dist_lyr = props.get_attribute_layer('dist_z', 'FLOAT', 'POINT', cbm)
                this_nd_lyr = props.get_attribute_layer('norm_dist_z', 'FLOAT', 'POINT', cbm)
            for face in cbm.faces:
                for fv in face.verts:
                    if use_ray == True:
                        inwards = -fv.normal
                        inwards_vec = inwards * plane_vec
                        inwards_dir = inwards_vec.normalized()
                        hit, hit_loc, hit_norm, hit_i, hit_o, hit_ab = bpy.context.scene.ray_cast(depsgraph, fv.co + (inwards*0.01), inwards_dir)
                        if hit:
                            last_mid = fv.co.lerp(hit_loc, 0.5)
                    else:
                        last_mid = face.calc_center_median_weighted()
                    fvec = fv.co - last_mid
                    fvdist = fvec.length
                    fv[this_mid_lyr] = last_mid
                    fv[this_dist_lyr] = fvdist
                    nmx = (last_mid.x - min_x)/dim_x
                    nmy = (last_mid.y - min_y)/dim_y
                    nmz = (last_mid.z - min_z)/dim_z
                    fv[this_nmid_lyr] = mathutils.Vector((nmx,nmy,nmz))
                    npx = (fv.co.x - min_x)/dim_x
                    npy = (fv.co.y - min_y)/dim_y
                    npz = (fv.co.z - min_z)/dim_z
                    fv[this_npx_lyr] = npx
                    fv[this_npy_lyr] = npy
                    fv[this_npz_lyr] = npz
                    fndvec = mathutils.Vector((npx,npy,npz)) - mathutils.Vector((nmx,nmy,nmz))
                    fndist = fndvec.length
                    fv[this_nd_lyr] = fndist
            if len(cbm.faces) < 1:
                for vert in cbm.verts:
                    vert[this_mid_lyr] = vert.co
                    vert[this_dist_lyr] = 0.0
                    if use_ray == True:
                        inwards = -vert.normal
                        inwards_vec = inwards * plane_vec
                        inwards_dir = inwards_vec.normalized()
                        hit, hit_loc, hit_norm, hit_i, hit_o, hit_ab = bpy.context.scene.ray_cast(depsgraph, vert.co + (inwards*0.01), inwards_dir)
                        if hit:
                            last_mid = vert.co.lerp(hit_loc, 0.5)
                            fvec = vert.co - last_mid
                            fvdist = fvec.length
                            vert[this_mid_lyr] = last_mid
                            vert[this_dist_lyr] = fvdist
        
        uv_multiplier = 1/math.ceil(math.sqrt(this_res))

        uv_layer = bpbm.loops.layers.uv.verify()

        bpbm.verts.ensure_lookup_table()
        bpbm.edges.ensure_lookup_table()
        bpbm.faces.ensure_lookup_table()

        last_face = bpbm.faces[0]

        for bpi, bp_lvl in enumerate(lvla):
            if bpi == 0:
                for bpvi, bpfv in enumerate(last_face.verts):
                    bpfv.co = corna[bpvi] + (vec * bp_lvl)
                    for bpl in bpfv.link_loops:
                        bpl[uv_layer].uv.x = uv_corners[bpvi][0] * uv_multiplier
                        bpl[uv_layer].uv.y = uv_corners[bpvi][1] * uv_multiplier
            if bpi > 0:
                nf = last_face.copy()
                uv_y_add = (math.floor(bpi * uv_multiplier)) * uv_multiplier
                uv_x_add = math.fmod(bpi * uv_multiplier, 1.0)
                for bpvi, bpfv in enumerate(nf.verts):
                    bpfv.co = corna[bpvi] + (vec * bp_lvl)
                    for bpl in bpfv.link_loops:
                        bpl[uv_layer].uv.x = (uv_corners[bpvi][0] * uv_multiplier) + uv_x_add
                        bpl[uv_layer].uv.y = (uv_corners[bpvi][1] * uv_multiplier) + uv_y_add
                last_face = nf
        
        for lvl_cut_a in lvl_cuts:
            for lvl_cut_o in lvl_cut_a:
                lvl_cut_o.data.update()

        

        # for bmov in bmo.verts:
        #     effective_co = bmov.co * vec
        #     cur_ax_val = effective_co.x + effective_co.y + effective_co.z
        #     above_or_below = False
        #     if cur_ax_val <= dmin:
        #         bmov[lyr] = 0.0
        #         bmov[np] = 0.0
        #         above_or_below = True
        #         print('Vert', bmov.index, 'is below bounds in axis ', ax)
        #     if cur_ax_val >= dmax:
        #         bmov[lyr] = 1.0
        #         bmov[np] = 1.0
        #         print('Vert', bmov.index, 'is above bounds in axis ', ax)
        #         above_or_below = True
        #     if above_or_below == False:
        #         for lvi, lv in enumerate(lvl_vol_map):
        #             # print(ax, lvli, lvla[lvli], vol_lvls[lvli])
        #             this_lvl = lvl_vol_map[lvi][0]
        #             this_vol = lvl_vol_map[lvi][1]
        #             this_np = lvl_vol_map[lvi][2]
        #             below_lvl = 0.0
        #             try:
        #                 below_lvl = lvl_vol_map[lvi + 1][0]
        #             except:
        #                 pass
        #             below_vol = 0.0
        #             try:
        #                 below_vol = lvl_vol_map[lvi + 1][1]
        #             except:
        #                 pass
        #             bwloe_np = 0.0
        #             try:
        #                 below_np = lvl_vol_map[lvi + 1][2]
        #             except:
        #                 pass
        #             if cur_ax_val < this_lvl and cur_ax_val > below_lvl:
        #                 cur_alpha = (below_lvl - cur_ax_val) / inc
        #                 new_vol = bl_math.lerp(this_vol, below_vol, cur_alpha)
        #                 new_np = bl_math.lerp(this_np, below_np, cur_alpha)
        #                 # print('Vert', bmov.index, 'is between', str(lvli), 'and', str(lvli + 1), 'at alpha', str(cur_alpha), 'with normalized position in', ax, '=', str(new_np))
        #                 bmov[lyr] = new_vol
        #                 bmov[np] = new_np
        
            
    # for bmov in bmo.verts:
    #     bmov[orig_co] = bmov.co
    
    bpy.context.active_object.data.update()

for bake_plane in bake_planes:
    bake_plane.data.update()

bpy.ops.object.editmode_toggle()
        


