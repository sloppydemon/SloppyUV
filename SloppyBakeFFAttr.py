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

x_res = 24
y_res = 24
z_res = 24

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

lvl_masks = [[],[],[]]

for lvli, lvlsa in enumerate(lvls):
    for lvlvi, lvlv in enumerate(lvlsa):
        this_mask_o = mo.copy()
        this_mask_o.data = mo.data.copy()
        this_mask_o.name = mo.name + '_mask.' + olax[lvli] + str(lvlvi)
        bpy.context.collection.objects.link(this_mask_o)
        props.find_or_add_attribute_other_obj('mask_x', 'FLOAT', 'POINT', this_mask_o)
        props.find_or_add_attribute_other_obj('mask_y', 'FLOAT', 'POINT', this_mask_o)
        props.find_or_add_attribute_other_obj('mask_z', 'FLOAT', 'POINT', this_mask_o)
        this_mask_o.select_set(True)
        lvl_masks[lvli].append(this_mask_o)

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

maskbms = [[],[],[]]

for lvl_mask_i, lvl_mask_a in enumerate(lvl_masks):
    for mask_i, mask_o in enumerate(lvl_mask_a):
        lvl_mask_bm = bmesh.from_edit_mesh(mask_o.data)
        maskbms[lvl_mask_i].append(lvl_mask_bm)

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
    for obj,vobj, vec,lyr,mid,np,dst,ax,lvla,lvln,cbms,lbms,mbms,dmin,dmax,inc,bpbm,corna in zip(ol, olv, olcovm, ollyr,olmid,olnp,oldist,olax,lvls,lvlns,cutbms,lerpbms,maskbms,mins,maxs,incs,bpbms,corners):
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
            mbm = mbms[ri]

            this_lerp = props.get_attribute_layer('lerp_x', 'FLOAT', 'POINT', lbm)
            this_mask = props.get_attribute_layer('mask_x', 'FLOAT', 'POINT', mbm)
            this_loco = props.get_attribute_layer('orig_co', 'FLOAT_VECTOR', 'POINT', lbm)
            this_cmid = props.get_attribute_layer('xmid', 'FLOAT_VECTOR', 'POINT', cbm)
            this_cdist = props.get_attribute_layer('dist_x', 'FLOAT', 'POINT', cbm)
            this_oco = props.get_attribute_layer('orig_co', 'FLOAT_VECTOR', 'POINT', cbm)
            if ax == 'y':
                this_lerp = props.get_attribute_layer('lerp_y', 'FLOAT', 'POINT', lbm)
                this_mask = props.get_attribute_layer('mask_y', 'FLOAT', 'POINT', mbm)
                this_cmid = props.get_attribute_layer('ymid', 'FLOAT_VECTOR', 'POINT', cbm)
                this_cdist = props.get_attribute_layer('dist_y', 'FLOAT', 'POINT', cbm)
            if ax == 'z':
                this_lerp = props.get_attribute_layer('lerp_z', 'FLOAT', 'POINT', lbm)
                this_mask = props.get_attribute_layer('mask_z', 'FLOAT', 'POINT', mbm)
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

            mask_cut = bmesh.ops.bisect_plane(mbm, geom=mbm.verts[:] + mbm.edges[:] + mbm.faces[:], dist=0.0001, plane_co=plane_pos, plane_no=vec, use_snap_center=False, clear_outer=True, clear_inner=True)
            bmesh.ops.holes_fill(mbm, edges=mbm.edges)
            for mv in mbm.verts:
                mv[this_mask] = 1.0
            
            hicut = bmesh.ops.bisect_plane(lbm, geom=lbm.verts[:] + lbm.edges[:] + lbm.faces[:], dist=0.0001, plane_co=lerp_plane_pos, plane_no=vec, use_snap_center=False, clear_outer=True, clear_inner=False)
            locut = bmesh.ops.bisect_plane(lbm, geom=lbm.verts[:] + lbm.edges[:] + lbm.faces[:], dist=0.0001, plane_co=plane_pos, plane_no=vec, use_snap_center=False, clear_outer=False, clear_inner=True)
            
            for lv in lbm.verts:
                val_to_lerp = lv[this_loco].x
                if ax == 'y':
                    val_to_lerp = lv[this_loco].y
                if ax == 'z':
                    val_to_lerp = lv[this_loco].z
                nu_lerp = (val_to_lerp - lvl) / inc
                nulco = (lv.co * plane_vec) + (vec * lvl) - (vec * (1 - nu_lerp) * (inc/4))
                lv.co = nulco
                nu_lerp *= 1.5
                nu_lerp -= 0.25
                lv[this_lerp] = nu_lerp

            ccut = bmesh.ops.bisect_plane(cbm, geom=cbm.verts[:] + cbm.edges[:] + cbm.faces[:], dist=0.0001, plane_co=plane_pos, plane_no=vec, use_snap_center=False, clear_outer=True, clear_inner=True)
            cut_verts = [i for i in ccut['geom_cut'] if i in cbm.verts]
            cut_edges = [i for i in ccut['geom_cut'] if i in cbm.edges]
            extruded = bmesh.ops.extrude_edge_only(cbm, edges=cbm.edges)
            new_verts = [i for i in extruded['geom'] if i in cbm.verts and i not in cut_verts]
            new_edges = [i for i in extruded['geom'] if i in cbm.edges and i not in cut_edges]
            new_faces = [i for i in extruded['geom'] if i in cbm.faces]

            for cv in cut_verts:
                cv[this_cdist] = -0.25

            for nv in new_verts:
                # nuco = nv[this_cmid] * plane_vec + (nv.co * vec)
                nv_vec = (nv[this_cmid] - nv.co) * plane_vec
                nv_dir = nv_vec.normalized()
                nuco = nv.co + (nv_dir * 0.01)
                nv.co = nuco
                nv[this_cdist] = 1.25
            
            avg_dist = 1.0
            for nv in new_verts:
                avg_dist += nv[this_cdist]
            if len(new_verts) > 0:
                avg_dist /= len(new_verts)
            else:
                avg_dist = 0.05

            print('After cut:',len([v for v in obj.verts if v.is_valid]), 'Number of new faces:', len(new_faces), 'Avg. distance to midpoint:', avg_dist)

        rest_verts = [i for i in obj.verts if i not in ocut_geo]
        rest_edges = [i for i in obj.edges if i not in ocut_geo]
        rest_faces = [i for i in obj.faces if i not in ocut_geo]
        bmesh.ops.dissolve_edges(obj, edges=rest_edges, use_verts=False,use_face_split=True)

        for lvl_cut_a in lvl_cuts:
            for lvl_cut_o in lvl_cut_a:
                lvl_cut_o.data.update()
        for lvl_lerp_a in lvl_lerps:
            for lvl_lerp_o in lvl_lerp_a:
                lvl_lerp_o.data.update()

    for obj,vobj, vec,lyr,mid,np,dst,ax,lvla,lvln,cbms,dmin,dmax,inc,bpbm,corna in zip(ol, olv, olcovm, ollyr,olmid,olnp,oldist,olax,lvls,lvlns,cutbms,mins,maxs,incs,bpbms,corners):
        this_res = x_res
        if ax == 'y':
            this_res = y_res
        if ax == 'z':
            this_res = z_res
        uv_multiplier = 1/math.ceil(math.sqrt(this_res + 1))
        uv_init_add_y = (math.ceil(math.sqrt(this_res + 1)) - 1) * uv_multiplier

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
                        bpl[uv_layer].uv.y = (uv_corners[bpvi][1] * uv_multiplier) + uv_init_add_y
            if bpi > 0:
                nf = last_face.copy()
                uv_y_add = (math.floor(bpi * uv_multiplier)) * uv_multiplier
                uv_x_add = math.fmod(bpi * uv_multiplier, 1.0)
                for bpvi, bpfv in enumerate(nf.verts):
                    bpfv.co = corna[bpvi] + (vec * bp_lvl)
                    for bpl in bpfv.link_loops:
                        bpl[uv_layer].uv.x = (uv_corners[bpvi][0] * uv_multiplier) + uv_x_add
                        bpl[uv_layer].uv.y = (uv_corners[bpvi][1] * uv_multiplier) + uv_init_add_y - uv_y_add
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
    joined_lerps.append(lerpo)
    bpy.ops.object.select_all(action='DESELECT')

joined_masks = []

for maskos in lvl_masks:
    for masko in maskos:
        masko.data.update()
        masko.select_set(True)
        bpy.context.view_layer.objects.active = masko
    bpy.ops.object.join()
    joined_masks.append(masko)
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
    bt_bsdf = bt_mat.node_tree.nodes['Principled BSDF']
    bf_ed = bf_mat.node_tree.nodes['DistFromEdge']
    jcut.select_set(True)
    try:
        bpy.ops.object.shade_smooth(keep_sharp_edges=False)
    except:
        pass
    try:
        bpy.ops.mesh.customdata_custom_splitnormals_clear()
    except:
        pass
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
    bt_mat.node_tree.links.new(bt_bsdf.inputs[0], bt_tex.outputs['Color'])
    bpy.context.scene.cycles.bake_type = 'DIFFUSE'
    bpy.context.scene.render.bake.margin = 1024
    bpy.context.scene.render.bake.margin_type = 'EXTEND'
    bpy.context.scene.render.bake.use_selected_to_active = True
    bpy.context.scene.render.bake.max_ray_distance = inc/2
    bpy.context.scene.render.bake.cage_extrusion = 0.004
    bpy.context.scene.render.bake.use_pass_direct = False
    bpy.context.scene.render.bake.use_pass_indirect = False
    bpy.context.scene.render.bake.use_pass_color = True
    bpy.ops.object.bake(type='DIFFUSE')
    img.pack()
    jcut.select_set(False)
    bapl.select_set(False)

for bapl, jlerp, lax, inc in zip(bake_planes, joined_lerps, olax, incs):
    ax = lax.upper()
    bf_mat = bpy.data.materials['M_Bake']
    bt_mat = bpy.data.materials['M_BakePlane']
    bt_tex = bt_mat.node_tree.nodes['Tex_BakeTo']
    bf_bsdf = bf_mat.node_tree.nodes['Principled BSDF']
    bf_ed = bf_mat.node_tree.nodes['LevelLerp']
    jlerp.select_set(True)
    try:
        bpy.ops.object.shade_smooth(keep_sharp_edges=False)
    except:
        pass
    try:
        bpy.ops.mesh.customdata_custom_splitnormals_clear()
    except:
        pass
    bapl.select_set(True)
    bpy.context.view_layer.objects.active = bapl
    ed_nam_str = 'T_' + mo.name + '_LevelLerp_' + ax
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
    bpy.context.scene.render.bake.margin_type = 'EXTEND'
    bpy.context.scene.render.bake.use_selected_to_active = True
    bpy.context.scene.render.bake.max_ray_distance = inc/2
    bpy.context.scene.render.bake.cage_extrusion = 0.004
    bpy.context.scene.render.bake.use_pass_direct = False
    bpy.context.scene.render.bake.use_pass_indirect = False
    bpy.context.scene.render.bake.use_pass_color = True
    bpy.ops.object.bake(type='DIFFUSE')
    img.pack()
    jlerp.select_set(False)
    bapl.select_set(False)

for bapl, jmask, lax, inc in zip(bake_planes, joined_masks, olax, incs):
    ax = lax.upper()
    bf_mat = bpy.data.materials['M_Bake']
    bt_mat = bpy.data.materials['M_BakePlane']
    bt_tex = bt_mat.node_tree.nodes['Tex_BakeTo']
    bf_bsdf = bf_mat.node_tree.nodes['Principled BSDF']
    bf_ed = bf_mat.node_tree.nodes['LevelMask']
    jmask.select_set(True)
    try:
        bpy.ops.object.shade_smooth(keep_sharp_edges=False)
    except:
        pass
    try:
        bpy.ops.mesh.customdata_custom_splitnormals_clear()
    except:
        pass
    bapl.select_set(True)
    bpy.context.view_layer.objects.active = bapl
    ed_nam_str = 'T_' + mo.name + '_LevelMask_' + ax
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
    bpy.context.scene.render.bake.margin = 0
    bpy.context.scene.render.bake.margin_type = 'ADJACENT_FACES'
    bpy.context.scene.render.bake.use_selected_to_active = True
    bpy.context.scene.render.bake.max_ray_distance = inc/2
    bpy.context.scene.render.bake.cage_extrusion = 0.004
    bpy.context.scene.render.bake.use_pass_direct = False
    bpy.context.scene.render.bake.use_pass_indirect = False
    bpy.context.scene.render.bake.use_pass_color = True
    bpy.ops.object.bake(type='DIFFUSE')
    img.pack()
    jlerp.select_set(False)
    bapl.select_set(False)

