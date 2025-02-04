import bpy
import bmesh
import mathutils
import math
import bl_math
from bpy_extras import bmesh_utils

props = bpy.context.scene.sloppy_props

depsgraph = bpy.context.evaluated_depsgraph_get()

per_vert = True
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

for axi, axis in enumerate(olax):
    this_axis_o = mo.copy()
    this_axis_o.data = mo.data.copy()
    this_axis_o.name = mo.name + '_lvls.' + axis
    bpy.context.collection.objects.link(this_axis_o)
    this_axis_o.select_set(True)
    axos[axi] = this_axis_o

bpy.ops.object.editmode_toggle()

bmo = bmesh.from_edit_mesh(mo.data)
bmo_vol = bmo.calc_volume()

bmox = bmesh.from_edit_mesh(axos[0].data)
bmoy = bmesh.from_edit_mesh(axos[1].data)
bmoz = bmesh.from_edit_mesh(axos[2].data)

ol = [
    bmox,
    bmoy,
    bmoz
]

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

def dist_sort(e):
    return e[1]

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
    for obj,vec,lyr,mid,ax,np,nd,nmid,dst in zip(ol, olcovm, ollyr,olmid,olax,olnp,olndist,olnmid,oldist):
        ls = []
        for v in obj.verts:
            lsi = [v, v.co, v.index]
            ls.append(lsi)
        current_vm = vec
        plane_vec = mathutils.Vector((1,1,1)) - vec
        ray_vec_a = mathutils.Vector((0,1,0))
        ray_vec_b = mathutils.Vector((0,0,1))
        if ax == 'x':
            ls.sort(key = vec_sort_x, reverse=True)
        if ax == 'y':
            ls.sort(key = vec_sort_y, reverse=True)
            ray_vec_a = mathutils.Vector((1,0,0))
            ray_vec_b = mathutils.Vector((0,0,1))
        if ax == 'z':
            ls.sort(key = vec_sort_z, reverse=True)
            ray_vec_a = mathutils.Vector((1,0,0))
            ray_vec_b = mathutils.Vector((0,1,0))
        last_mid = ls[0][1]
        hl_a = ls[0][1]
        hl_b = ls[0][1]
        hl_c = ls[0][1]

        for vert in ls:
            print('Before cut:', len([v for v in obj.verts if v.is_valid]), vert[1] * vec, 'axis:', ax)
            vol = vol_calc(obj.calc_volume(), bmo_vol)
            bmo.verts[vert[2]][lyr] = vol
            cut, rest = bmesh.ops.bisect_plane(obj, geom=obj.verts[:] + obj.edges[:] + obj.faces[:], dist=0.0001, plane_co=vert[1] * vec, plane_no=vec, use_snap_center=False, clear_outer=True, clear_inner=False)
            new_faces = bmesh.ops.holes_fill(obj, edges=obj.edges)
            print('After cut:',len([v for v in obj.verts if v.is_valid]), 'Number of new faces:', len(new_faces['faces']), 'New volume:', vol)
        
        for bmov in bmo.verts:
            bmov[orig_co] = bmov.co
            inwards = -bmov.normal
            inwards_vec_a = inwards * ray_vec_a
            inwards_dir_a = inwards_vec_a.normalized()
            inwards_vec_b = inwards * ray_vec_b
            inwards_dir_b = inwards_vec_b.normalized()
            inwards_vec_c = inwards * plane_vec
            inwards_dir_c = inwards_vec_c.normalized()
            hit_a, hit_a_loc, hit_a_norm, hit_a_i, hit_a_o, hit_a_ab = bpy.context.scene.ray_cast(depsgraph, bmov.co + (inwards*0.002), inwards_dir_a)
            if hit_a:
                hl_a = hit_a_loc
            else:
                nuhl_a = (hl_a * plane_vec) + (bmov.co * vec)
                hl_a = nuhl_a
            hit_b, hit_b_loc, hit_b_norm, hit_b_i, hit_b_o, hit_b_ab = bpy.context.scene.ray_cast(depsgraph, bmov.co + (inwards*0.002), inwards_dir_b)
            if hit_b:
                hl_b = hit_b_loc
            else:
                nuhl_b = (hl_b * plane_vec) + (bmov.co * vec)
                hl_b = nuhl_b
            hit_c, hit_c_loc, hit_c_norm, hit_c_i, hit_c_o, hit_c_ab = bpy.context.scene.ray_cast(depsgraph, bmov.co + (inwards*0.002), inwards_dir_c)
            if hit_c:
                hl_c = hit_c_loc
            else:
                nuhl_c = (hl_c * plane_vec) + (bmov.co * vec)
                hl_c = nuhl_c

            sort_hits = [[hl_a], [hl_b], [hl_c]]

            for shi, sh in enumerate(sort_hits):
                sh_vec = sh[0] - bmov.co
                sh_len = sh_vec.length
                sort_hits[shi].append(sh_vec.length)

            sort_hits.sort(key = dist_sort)

            # avg_hl = (hl_a + hl_b + hl_c) / 3
            avg_hl = hl_a.lerp(hl_b, 0.5)
            last_mid = bmov.co.lerp(avg_hl, 0.5)
            fvec = bmov.co - last_mid
            fvdist = fvec.length
            bmov[mid] = last_mid
            bmov[dst] = fvdist
            nmx = (last_mid.x - min_x)/dim_x
            nmy = (last_mid.y - min_y)/dim_y
            nmz = (last_mid.z - min_z)/dim_z
            bmov[nmid] = mathutils.Vector((nmx,nmy,nmz))
            npx = (bmov.co.x - min_x)/dim_x
            npy = (bmov.co.y - min_y)/dim_y
            npz = (bmov.co.z - min_z)/dim_z
            bmov[norm_pos_x] = npx
            bmov[norm_pos_y] = npy
            bmov[norm_pos_z] = npz
            fndvec = mathutils.Vector((npx,npy,npz)) - mathutils.Vector((nmx,nmy,nmz))
            fndist = fndvec.length
            bmov[nd] = fndist

            vi_str = str(bmov.index).rjust(5)
            col_len_a = max(len('Vert' + vi_str), len('World:'), len('Normalized:'))
            col_len_bc_sub = 7
            col_len_bc = (col_len_bc_sub*3) + 2
            col_len_d = 12

            v_str = f'Vert {vi_str}'
            of_str = f'of {len(bmo.verts)}'

            mid_x_str = '{:.3f}'.format(last_mid.x).rjust(col_len_bc_sub)
            mid_y_str = '{:.3f}'.format(last_mid.y).rjust(col_len_bc_sub)
            mid_z_str = '{:.3f}'.format(last_mid.y).rjust(col_len_bc_sub)

            nmid_x_str = '{:.3f}'.format(nmx).rjust(col_len_bc_sub)
            nmid_y_str = '{:.3f}'.format(nmy).rjust(col_len_bc_sub)
            nmid_z_str = '{:.3f}'.format(nmz).rjust(col_len_bc_sub)
            
            pos_x_str = '{:.3f}'.format(bmov.co.x).rjust(col_len_bc_sub)
            pos_y_str = '{:.3f}'.format(bmov.co.y).rjust(col_len_bc_sub)
            pos_z_str = '{:.3f}'.format(bmov.co.z).rjust(col_len_bc_sub)
            
            npos_x_str = '{:.3f}'.format(npx).rjust(col_len_bc_sub)
            npos_y_str = '{:.3f}'.format(npy).rjust(col_len_bc_sub)
            npos_z_str = '{:.3f}'.format(npz).rjust(col_len_bc_sub)
            
            dist_str = '{:.4f}'.format(fvdist)
            ndist_str = '{:.4f}'.format(fndist)

            axis_str = 'X'.center(col_len_bc_sub) + 'Y'.center(col_len_bc_sub) + 'Z'.center(col_len_bc_sub)

            mid_str = f'{mid_x_str} {mid_x_str} {mid_x_str}'
            nmid_str = f'{nmid_x_str} {nmid_y_str} {nmid_z_str}'
            pos_str = f'{pos_x_str} {pos_y_str} {pos_z_str}'
            npos_str = f'{npos_x_str} {npos_y_str} {npos_z_str}'
            if hit_a or hit_b or hit_c:
                print('Raycasting inwards from vertex', bmov.index, 'in axis', ax, ': Hit!')
                if hit_a:
                    print('Hit object A:', hit_a_o.name)
                if hit_b:
                    print('Hit object B:', hit_b_o.name)
                if hit_c:
                    print('Hit object C:', hit_c_o.name)
                print(v_str.ljust(col_len_a) + 'Position:'.center(col_len_bc) + '     ' +  'Midpoint:'.center(col_len_bc) + 'Distance:'.rjust(col_len_d))
                print(of_str.ljust(col_len_a) + axis_str.center(col_len_bc) + '   | ' + axis_str.center(col_len_bc))
                print('World:'.ljust(col_len_a) + pos_str.center(col_len_bc) + '   | ' + mid_str.center(col_len_bc) + dist_str.rjust(col_len_d))
                print('Normalized:'.ljust(col_len_a) + npos_str.center(col_len_bc) + '   | ' + nmid_str.center(col_len_bc) + ndist_str.rjust(col_len_d))
            if not hit_a and not hit_b and not hit_c:
                print('Raycasting inwards from ', bmov.index, ' in axis ', ax, ': No hit.') 

        bpy.context.active_object.data.update()

bpy.ops.object.editmode_toggle()

bpy.ops.object.select_all(action='DESELECT')

for laxo in axos:
    laxo.select_set(True)
    bpy.context.view_layer.objects.active = laxo
    bpy.ops.object.delete()