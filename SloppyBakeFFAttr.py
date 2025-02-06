import bpy
import bmesh
import mathutils
import math
import bl_math
from bpy_extras import bmesh_utils

props = bpy.context.scene.sloppy_props

depsgraph = bpy.context.evaluated_depsgraph_get()

per_vert = False
# use_ray = False

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
bake_plane_bases = [None, None, None]

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

    this_bake_plane_base = bake_plane.copy()
    this_bake_plane_base.data = bake_plane.data.copy()
    this_bake_plane_base.name = mo.name + '_basePlane.' + axis
    bpy.context.collection.objects.link(this_bake_plane_base)
    this_bake_plane_base.select_set(True)
    bake_plane_bases[axi] = this_bake_plane_base

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

lvl_lerps = [[],[],[]]

for lvli, lvlsa in enumerate(lvls):
    for lvlvi, lvlv in enumerate(lvlsa):
        this_lerp_o = mo.copy()
        this_lerp_o.data = mo.data.copy()
        this_lerp_o.name = mo.name + '_lerp.' + olax[lvli] + str(lvlvi)
        bpy.context.collection.objects.link(this_lerp_o)
        props.find_or_add_attribute_other_obj('lerp_x', 'FLOAT', 'POINT', this_lerp_o)
        props.find_or_add_attribute_other_obj('lerp_y', 'FLOAT', 'POINT', this_lerp_o)
        props.find_or_add_attribute_other_obj('lerp_z', 'FLOAT', 'POINT', this_lerp_o)
        this_lerp_o.select_set(True)
        lvl_lerps[lvli].append(this_lerp_o)

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
    for cut_i, cut_o in enumerate(lvl_cut_a):
        lvl_cut_bm = bmesh.from_edit_mesh(cut_o.data)
        cutbms[lvl_cut_i].append(lvl_cut_bm)

lerpbms = [[],[],[]]

for lvl_lerp_i, lvl_lerp_a in enumerate(lvl_lerps):
    for lerp_i, lerp_o in enumerate(lvl_lerp_a):
        lvl_lerp_bm = bmesh.from_edit_mesh(lerp_o.data)
        lerpbms[lvl_lerp_i].append(lvl_lerp_bm)

bpbms = []
bpbbms = []

for bp in  bake_planes:
    bake_plane_bm = bmesh.from_edit_mesh(bp.data)
    bpbms.append(bake_plane_bm)
for bpb in  bake_plane_bases:
    bake_plane_base_bm = bmesh.from_edit_mesh(bpb.data)
    bpbbms.append(bake_plane_base_bm)

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
    {"name": "nzmid", "type": "FLOAT_VECTOR", "domain": "POINT", "layer": None},
    {"name": "nxmid", "type": "FLOAT_VECTOR", "domain": "POINT", "layer": None},
    {"name": "nymid", "type": "FLOAT_VECTOR", "domain": "POINT", "layer": None},
    {"name": "dist_x", "type": "FLOAT", "domain": "POINT", "layer": None},
    {"name": "dist_y", "type": "FLOAT", "domain": "POINT", "layer": None},
    {"name": "dist_z", "type": "FLOAT", "domain": "POINT", "layer": None},
    {"name": "norm_pos_x", "type": "FLOAT", "domain": "POINT", "layer": None},
    {"name": "norm_pos_y", "type": "FLOAT", "domain": "POINT", "layer": None},
    {"name": "norm_pos_z", "type": "FLOAT", "domain": "POINT", "layer": None},
    {"name": "norm_dist_x", "type": "FLOAT", "domain": "POINT", "layer": None},
    {"name": "norm_dist_y", "type": "FLOAT", "domain": "POINT", "layer": None},
    {"name": "norm_dist_z", "type": "FLOAT", "domain": "POINT", "layer": None},
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
nzmid = props.get_dict_layer("nzmid", attr_dict)
nxmid = props.get_dict_layer("nxmid", attr_dict)
nymid = props.get_dict_layer("nymid", attr_dict)
dist_x = props.get_dict_layer("dist_x", attr_dict)
dist_y = props.get_dict_layer("dist_y", attr_dict)
dist_z = props.get_dict_layer("dist_z", attr_dict)
norm_dist_x = props.get_dict_layer("norm_dist_x", attr_dict)
norm_dist_y = props.get_dict_layer("norm_dist_y", attr_dict)
norm_dist_z = props.get_dict_layer("norm_dist_z", attr_dict)
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

olnmid = [
    nxmid,
    nymid,
    nzmid
]

oldist = [
    dist_x,
    dist_y,
    dist_z
]

olndist = [
    norm_dist_x,
    norm_dist_y,
    norm_dist_z
]

olnp = [
    norm_pos_x,
    norm_pos_y,
    norm_pos_z
]

bmo.verts.ensure_lookup_table()
bmo.edges.ensure_lookup_table()
bmo.faces.ensure_lookup_table()

for bmov in bmo.verts:
    this_npx = (bmov.co.x - min_x) / dim_x
    this_npy = (bmov.co.y - min_y) / dim_y
    this_npz = (bmov.co.z - min_z) / dim_z
    bmov[norm_pos_x] = this_npx
    bmov[norm_pos_y] = this_npy
    bmov[norm_pos_z] = this_npz

if per_vert == False:
    for obj,vobj, vec,lyr,mid,np,dst,ax,lvla,lvln,cbms,lbms,dmin,dmax,inc,bpbm,bpbbm,corna in zip(ol, olv, olcovm, ollyr,olmid,olnp,oldist,olax,lvls,lvlns,cutbms,lerpbms,mins,maxs,incs,bpbms,bpbbms,corners):
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

        this_res = x_res
        if ax == 'y':
            this_res = y_res
        if ax == 'z':
            this_res = z_res
        max_dim = max(dim_x, dim_y, dim_z)

        # vol_lvls = []

        lvla_rng = [i for i in range(len(lvla))]
        print(lvla_rng)
        lvla_rng.sort(reverse = True)

        # lvl_vol_map = []

        ocut_geo = []

        for ri in lvla_rng:
            lvl = lvla[ri]
            lerp_lvl = lvl + inc
            cbm = cbms[ri]
            lbm = lbms[ri]

            this_lerp = props.get_attribute_layer('lerp_x', 'FLOAT', 'POINT', lbm)
            this_loco = props.get_attribute_layer('orig_co', 'FLOAT_VECTOR', 'POINT', lbm)
            this_cmid = props.get_attribute_layer('xmid', 'FLOAT_VECTOR', 'POINT', cbm)
            this_cdist = props.get_attribute_layer('dist_x', 'FLOAT', 'POINT', cbm)
            this_oco = props.get_attribute_layer('orig_co', 'FLOAT_VECTOR', 'POINT', cbm)
            if ax == 'y':
                this_lerp = props.get_attribute_layer('lerp_y', 'FLOAT', 'POINT', lbm)
                this_cmid = props.get_attribute_layer('ymid', 'FLOAT_VECTOR', 'POINT', cbm)
                this_cdist = props.get_attribute_layer('dist_y', 'FLOAT', 'POINT', cbm)
            if ax == 'z':
                this_lerp = props.get_attribute_layer('lerp_z', 'FLOAT', 'POINT', lbm)
                this_cmid = props.get_attribute_layer('zmid', 'FLOAT_VECTOR', 'POINT', cbm)
                this_cdist = props.get_attribute_layer('dist_z', 'FLOAT', 'POINT', cbm)
            
            lvl_str = '{:.3f}'.format(lvl)
            lerp_lvl_str = '{:.3f}'.format(lerp_lvl)
            print('Before cut:', len([v for v in obj.verts if v.is_valid]), 'level: ', lvl, ' in axis: ', ax)
            
            plane_pos = vec * lvl
            lerp_plane_pos = vec * lerp_lvl
            ocut = bmesh.ops.bisect_plane(obj, geom=obj.verts[:] + obj.edges[:] + obj.faces[:], dist=0.0001, plane_co=plane_pos, plane_no=vec, use_snap_center=False, clear_outer=False, clear_inner=False)
            for element in ocut['geom_cut']:
                if element not in ocut_geo:
                    ocut_geo.append(element)
            rest_verts = [i for i in ocut['geom'] if i in obj.verts and i not in ocut['geom_cut']]
            rest_edges = [i for i in ocut['geom'] if i in obj.edges and i not in ocut['geom_cut']]
            rest_faces = [i for i in ocut['geom'] if i in obj.faces and i not in ocut['geom_cut']]
            # rest_cut, rest_rest = bmesh.ops.bisect_plane(rbm, geom=rbm.verts[:] + rbm.edges[:] + rbm.faces[:], dist=0.0001, plane_co=plane_pos, plane_no=vec, use_snap_center=False, clear_outer=False, clear_inner=False)
            hicut = bmesh.ops.bisect_plane(lbm, geom=lbm.verts[:] + lbm.edges[:] + lbm.faces[:], dist=0.0001, plane_co=lerp_plane_pos, plane_no=vec, use_snap_center=False, clear_outer=True, clear_inner=False)
            locut = bmesh.ops.bisect_plane(lbm, geom=lbm.verts[:] + lbm.edges[:] + lbm.faces[:], dist=0.0001, plane_co=plane_pos, plane_no=vec, use_snap_center=False, clear_outer=False, clear_inner=True)
            
            for lv in lbm.verts:
                val_to_lerp = lv[this_loco].x
                if ax == 'y':
                    val_to_lerp = lv[this_loco].y
                if ax == 'z':
                    val_to_lerp = lv[this_loco].z
                nu_lerp = (val_to_lerp - lvl) / inc
                lv[this_lerp] = nu_lerp
                nulco = (lv.co * plane_vec) + (vec * lvl)
                lv.co = nulco

            ccut = bmesh.ops.bisect_plane(cbm, geom=cbm.verts[:] + cbm.edges[:] + cbm.faces[:], dist=0.0001, plane_co=plane_pos, plane_no=vec, use_snap_center=False, clear_outer=True, clear_inner=True)
            cut_verts = [i for i in ccut['geom_cut'] if i in cbm.verts]
            cut_edges = [i for i in ccut['geom_cut'] if i in cbm.edges]
            extruded = bmesh.ops.extrude_edge_only(cbm, edges=cbm.edges)
            new_verts = [i for i in extruded['geom'] if i in cbm.verts and i not in cut_verts]
            new_edges = [i for i in extruded['geom'] if i in cbm.edges and i not in cut_edges]
            new_faces = [i for i in extruded['geom'] if i in cbm.faces]

            for cv in cut_verts:
                cv[this_cdist] = 0.0

            for nv in new_verts:
                # nuco = nv[this_cmid] * plane_vec + (nv.co * vec)
                nv_vec = (nv[this_cmid] - nv.co) * plane_vec
                nv_dir = nv_vec.normalized()
                nuco = nv.co + (nv_dir * 0.01)
                nv.co = nuco
                nv[this_cdist] = 1.0
            
            avg_dist = 1.0
            for nv in new_verts:
                avg_dist += nv[this_cdist]
            if len(new_verts) > 0:
                avg_dist /= len(new_verts)
            else:
                avg_dist = 0.05

            # bmesh.ops.holes_fill(cbm, edges=new_edges)

            # bmesh.ops.remove_doubles(cbm, verts=new_verts, dist=0.01)
            # if use_ray == False:
            #     cut_new_faces = bmesh.ops.holes_fill(cbm, edges=cbm.edges)
            #     new_faces = bmesh.ops.holes_fill(obj, edges=obj.edges)

            # cut_cut, cut_rest = bmesh.ops.bisect_plane(cbm, geom=cbm.verts[:] + cbm.edges[:] + cbm.faces[:], dist=0.0001, plane_co=plane_pos, plane_no=vec, use_snap_center=False, clear_outer=True, clear_inner=True)
            
            # vol = vol_calc(obj.calc_volume(), bmo_vol)

            # vol_lvls.append(vol)
            # lvl_np = lvln[ri]
            # lvl_vol_item = [lvl, vol, lvl_np]
            # lvl_vol_map.append(lvl_vol_item)
            # vol_str = '{:.4f}'.format(vol)
            print('After cut:',len([v for v in obj.verts if v.is_valid]), 'Number of new faces:', len(new_faces), 'Avg. distance to midpoint:', avg_dist)

        rest_verts = [i for i in obj.verts if i not in ocut_geo]
        rest_edges = [i for i in obj.edges if i not in ocut_geo]
        rest_faces = [i for i in obj.faces if i not in ocut_geo]
        bmesh.ops.dissolve_edges(obj, edges=rest_edges, use_verts=False,use_face_split=True)

        # vol_lvls.sort(reverse=True)

        # obj.data.update()
        

        # for cbm in cbms:
        #     this_mid_lyr = props.get_attribute_layer('xmid', 'FLOAT_VECTOR', 'POINT', cbm)
        #     this_nmid_lyr = props.get_attribute_layer('norm_xmid', 'FLOAT_VECTOR', 'POINT', cbm)
        #     this_dist_lyr = props.get_attribute_layer('dist_x', 'FLOAT', 'POINT', cbm)
        #     this_npx_lyr = props.get_attribute_layer('norm_pos_x', 'FLOAT', 'POINT', cbm)
        #     this_npy_lyr = props.get_attribute_layer('norm_pos_y', 'FLOAT', 'POINT', cbm)
        #     this_npz_lyr = props.get_attribute_layer('norm_pos_z', 'FLOAT', 'POINT', cbm)
        #     this_nd_lyr = props.get_attribute_layer('norm_dist_x', 'FLOAT', 'POINT', cbm)
        #     if ax == 'y':
        #         this_mid_lyr = props.get_attribute_layer('ymid', 'FLOAT_VECTOR', 'POINT', cbm)
        #         this_nmid_lyr = props.get_attribute_layer('norm_ymid', 'FLOAT_VECTOR', 'POINT', cbm)
        #         this_dist_lyr = props.get_attribute_layer('dist_y', 'FLOAT', 'POINT', cbm)
        #         this_nd_lyr = props.get_attribute_layer('norm_dist_y', 'FLOAT', 'POINT', cbm)
        #     if ax == 'z':
        #         this_mid_lyr = props.get_attribute_layer('zmid', 'FLOAT_VECTOR', 'POINT', cbm)
        #         this_nmid_lyr = props.get_attribute_layer('norm_zmid', 'FLOAT_VECTOR', 'POINT', cbm)
        #         this_dist_lyr = props.get_attribute_layer('dist_z', 'FLOAT', 'POINT', cbm)
        #         this_nd_lyr = props.get_attribute_layer('norm_dist_z', 'FLOAT', 'POINT', cbm)
        #     if len(cbm.faces) < 1:
        #         for vert in cbm.verts:
        #             vert[this_mid_lyr] = vert.co
        #             vert[this_dist_lyr] = 0.0
        #             if use_ray == True:
        #                 inwards = -vert.normal
        #                 inwards_vec = inwards * plane_vec
        #                 inwards_dir = inwards_vec.normalized()
        #                 hit, hit_loc, hit_norm, hit_i, hit_o, hit_ab = bpy.context.scene.ray_cast(depsgraph, vert.co + (inwards*0.01), inwards_dir)
        #                 if hit:
        #                     last_mid = vert.co.lerp(hit_loc, 0.5)
        #                 fvec = vert.co - last_mid
        #                 fvdist = fvec.length
        #                 vert[this_mid_lyr] = last_mid
        #                 vert[this_dist_lyr] = fvdist
        #                 nmx = (last_mid.x - min_x)/dim_x
        #                 nmy = (last_mid.y - min_y)/dim_y
        #                 nmz = (last_mid.z - min_z)/dim_z
        #                 vert[this_nmid_lyr] = mathutils.Vector((nmx,nmy,nmz))
        #                 npx = (vert.co.x - min_x)/dim_x
        #                 npy = (vert.co.y - min_y)/dim_y
        #                 npz = (vert.co.z - min_z)/dim_z
        #                 vert[this_npx_lyr] = npx
        #                 vert[this_npy_lyr] = npy
        #                 vert[this_npz_lyr] = npz
        #                 fndvec = mathutils.Vector((npx,npy,npz)) - mathutils.Vector((nmx,nmy,nmz))
        #                 fndist = fndvec.length
        #                 vert[this_nd_lyr] = fndist

        #                 vi_str = str(vert.index).rjust(5)
        #                 col_len_a = max(len('Vert' + vi_str), len('World:'), len('Normalized:'))
        #                 col_len_bc_sub = 7
        #                 col_len_bc = (col_len_bc_sub*3) + 2
        #                 col_len_d = 9

        #                 v_str = f'Vert {vi_str}'

        #                 mid_x_str = '{:.3f}'.format(last_mid.x).rjust(col_len_bc_sub)
        #                 mid_y_str = '{:.3f}'.format(last_mid.y).rjust(col_len_bc_sub)
        #                 mid_z_str = '{:.3f}'.format(last_mid.y).rjust(col_len_bc_sub)

        #                 nmid_x_str = '{:.3f}'.format(nmx).rjust(col_len_bc_sub)
        #                 nmid_y_str = '{:.3f}'.format(nmy).rjust(col_len_bc_sub)
        #                 nmid_z_str = '{:.3f}'.format(nmz).rjust(col_len_bc_sub)
                        
        #                 pos_x_str = '{:.3f}'.format(vert.co.x).rjust(col_len_bc_sub)
        #                 pos_y_str = '{:.3f}'.format(vert.co.y).rjust(col_len_bc_sub)
        #                 pos_z_str = '{:.3f}'.format(vert.co.z).rjust(col_len_bc_sub)
                        
        #                 npos_x_str = '{:.3f}'.format(npx).rjust(col_len_bc_sub)
        #                 npos_y_str = '{:.3f}'.format(npy).rjust(col_len_bc_sub)
        #                 npos_z_str = '{:.3f}'.format(npz).rjust(col_len_bc_sub)
                        
        #                 dist_str = '{:.4f}'.format(fvdist).rjust(col_len_d)
        #                 ndist_str = '{:.4f}'.format(fndist).rjust(col_len_d)

        #                 axis_str = 'X'.center(col_len_bc_sub) + 'Y'.center(col_len_bc_sub) + 'Z'.center(col_len_bc_sub)

        #                 mid_str = f'{mid_x_str} {mid_x_str} {mid_x_str}'
        #                 nmid_str = f'{nmid_x_str} {nmid_y_str} {nmid_z_str}'
        #                 pos_str = f'{pos_x_str} {pos_y_str} {pos_z_str}'
        #                 npos_str = f'{npos_x_str} {npos_y_str} {npos_z_str}'
        #                 if hit:
        #                     print('Raycasting inwards from vertex', vert.index, 'in axis', ax, ': Hit!')
        #                     print('Hit object:', hit_o.name)
        #                     print(v_str.ljust(col_len_a) + 'Position:'.center(col_len_bc) + '   ' +  'Midpoint:'.center(col_len_bc) + 'Distance:'.rjust(col_len_d))
        #                     print(' '.ljust(col_len_a) + axis_str.center(col_len_bc) + '   | ' + axis_str.center(col_len_bc))
        #                     print('World:'.ljust(col_len_a) + pos_str.center(col_len_bc) + '   | ' + mid_str.center(col_len_bc) + dist_str.rjust(col_len_d))
        #                     print('Normalized:'.ljust(col_len_a) + npos_str.center(col_len_bc) + '   | ' + nmid_str.center(col_len_bc) + ndist_str.rjust(col_len_d))
        #                 if not hit:
        #                     print('Raycasting inwards from ', vert.index, ' in axis ', ax, ': No hit.')   
        #     else:
        #         for face in cbm.faces:
        #             for fv in face.verts:
        #                 if use_ray == True:
        #                     inwards = -fv.normal
        #                     inwards_vec = inwards * plane_vec
        #                     inwards_dir = inwards_vec.normalized()
        #                     hit, hit_loc, hit_norm, hit_i, hit_o, hit_ab = bpy.context.scene.ray_cast(depsgraph, fv.co + (inwards*0.01), inwards_dir)
        #                     if hit:
        #                         last_mid = fv.co.lerp(hit_loc, 0.5)
        #                 else:
        #                     last_mid = face.calc_center_median_weighted()
        #                 fvec = fv.co - last_mid
        #                 fvdist = fvec.length
        #                 fv[this_mid_lyr] = last_mid
        #                 fv[this_dist_lyr] = fvdist
        #                 nmx = (last_mid.x - min_x)/dim_x
        #                 nmy = (last_mid.y - min_y)/dim_y
        #                 nmz = (last_mid.z - min_z)/dim_z
        #                 fv[this_nmid_lyr] = mathutils.Vector((nmx,nmy,nmz))
        #                 npx = (fv.co.x - min_x)/dim_x
        #                 npy = (fv.co.y - min_y)/dim_y
        #                 npz = (fv.co.z - min_z)/dim_z
        #                 fv[this_npx_lyr] = npx
        #                 fv[this_npy_lyr] = npy
        #                 fv[this_npz_lyr] = npz
        #                 fndvec = mathutils.Vector((npx,npy,npz)) - mathutils.Vector((nmx,nmy,nmz))
        #                 fndist = fndvec.length
        #                 fv[this_nd_lyr] = fndist

        for lvl_cut_a in lvl_cuts:
            for lvl_cut_o in lvl_cut_a:
                lvl_cut_o.data.update()
        for lvl_lerp_a in lvl_lerps:
            for lvl_lerp_o in lvl_lerp_a:
                lvl_lerp_o.data.update()

        

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
    for obj,vobj, vec,lyr,mid,np,dst,ax,lvla,lvln,cbms,dmin,dmax,inc,bpbm,bpbbm,corna in zip(ol, olv, olcovm, ollyr,olmid,olnp,oldist,olax,lvls,lvlns,cutbms,mins,maxs,incs,bpbms,bpbbms,corners):
        this_res = x_res
        if ax == 'y':
            this_res = y_res
        if ax == 'z':
            this_res = z_res
        uv_multiplier = 1/math.ceil(math.sqrt(this_res))

        uv_layer = bpbm.loops.layers.uv.verify()
        uv_layer_b = bpbbm.loops.layers.uv.verify()

        bpbm.verts.ensure_lookup_table()
        bpbm.edges.ensure_lookup_table()
        bpbm.faces.ensure_lookup_table()
        bpbbm.verts.ensure_lookup_table()
        bpbbm.edges.ensure_lookup_table()
        bpbbm.faces.ensure_lookup_table()

        last_face = bpbm.faces[0]
        base_face = bpbbm.faces[0]

        for bpi, bp_lvl in enumerate(lvla):
            if bpi == 0:
                for bpvi, bpfv in enumerate(last_face.verts):
                    bpfv.co = corna[bpvi] + (vec * bp_lvl)
                    for bpl in bpfv.link_loops:
                        bpl[uv_layer].uv.x = uv_corners[bpvi][0] * uv_multiplier
                        bpl[uv_layer].uv.y = uv_corners[bpvi][1] * uv_multiplier
                for bpbvi, bpbfv in enumerate(base_face.verts):
                    bpbfv.co = corna[bpbvi] + (vec * bp_lvl)
                    for bpbl in bpbfv.link_loops:
                        bpbl[uv_layer_b].uv.x = uv_corners[bpbvi][0]
                        bpbl[uv_layer_b].uv.y = uv_corners[bpbvi][1]
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
    # for bmov in bmo.verts:
    #     bmov[orig_co] = bmov.co
    
    bpy.context.active_object.data.update()

for bake_plane in bake_planes:
    bake_plane.data.update()

bpy.ops.object.editmode_toggle()

bpy.ops.object.select_all(action='DESELECT')

joined_cuts = []

for cutos in lvl_cuts:
    for cuto in cutos:
        cuto.data.update()
        cuto.select_set(True)
        bpy.context.view_layer.objects.active = cuto
    bpy.ops.object.join()
    joined_cuts.append(cuto)
    bpy.ops.object.select_all(action='DESELECT')

joined_lerps = []

for lerpos in lvl_lerps:
    for lerpo in lerpos:
        lerpo.data.update()
        lerpo.select_set(True)
        bpy.context.view_layer.objects.active = lerpo
    bpy.ops.object.join()
    joined_cuts.append(lerpo)
    bpy.ops.object.select_all(action='DESELECT')

for laxo, lvaxo in zip(axos, vaxos):
    laxo.select_set(True)
    bpy.context.view_layer.objects.active = laxo
    bpy.ops.object.delete()
    lvaxo.select_set(True)
    bpy.context.view_layer.objects.active = lvaxo
    bpy.ops.object.delete()



for bapl, jcut, lax, inc in zip(bake_planes, joined_cuts, olax, incs):
    ax = lax.upper()
    bf_mat = bpy.data.materials['M_Bake']
    bt_mat = bpy.data.materials['M_BakePlane']
    bt_tex = bt_mat.node_tree.nodes['Tex_BakeTo']
    bf_bsdf = bf_mat.node_tree.nodes['Principled BSDF']
    bf_ed = bf_mat.node_tree.nodes['DistFromEdge']
    jcut.select_set(True)
    bapl.select_set(True)
    bpy.context.view_layer.objects.active = bapl
    ed_nam_str = 'T_' + mo.name + '_EdgeDist_' + ax
    img = None
    try:
        img = bpy.data.images[ed_nam_str]
    except:
        bpy.ops.image.new(name=ed_nam_str, width=4096, height=4096, color=(0,0,0,1), alpha=False, generated_type='BLANK', float=False, use_stereo_3d=False, tiled=False)
        img = bpy.data.images[ed_nam_str]
    bt_tex.select = True
    bt_tex.image = img
    bf_mat.node_tree.links.new(bf_bsdf.inputs[0], bf_ed.outputs[ax])
    bpy.context.scene.cycles.bake_type = 'DIFFUSE'
    bpy.context.scene.render.bake.margin = 1024
    bpy.context.scene.render.bake.margin_type = 'ADJACENT_FACES'
    bpy.context.scene.render.bake.use_selected_to_active = True
    bpy.context.scene.render.bake.max_ray_distance = inc/2
    bpy.context.scene.render.bake.cage_extrusion = 0.004
    bpy.context.scene.render.bake.use_pass_direct = False
    bpy.context.scene.render.bake.use_pass_indirect = False
    bpy.context.scene.render.bake.use_pass_color = True
    bpy.ops.object.bake(type='DIFFUSE')
    jcut.select_set(False)
    bapl.select_set(False)

