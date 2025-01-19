import bpy
import bmesh
import mathutils
from bpy_extras import bmesh_utils

props = bpy.context.scene.sloppy_props

depsgraph = bpy.context.evaluated_depsgraph_get()

x_res = 8
y_res = 8
z_res = 8

x_div = x_res - 1
y_div = y_res - 1
z_div = z_res - 1

mo = bpy.context.object

xo = mo.copy()
xo.data = mo.data.copy()
bpy.context.collection.objects.link(xo)
xo.select_set(True)
yo = mo.copy()
yo.data = mo.data.copy()
bpy.context.collection.objects.link(yo)
yo.select_set(True)
zo = mo.copy()
zo.data = mo.data.copy()
bpy.context.collection.objects.link(zo)
zo.select_set(True)

bpy.ops.object.editmode_toggle()

bmo = bmesh.from_edit_mesh(mo.data)
bmo_vol = bmo.calc_volume()

bmox = bmesh.from_edit_mesh(xo.data)
bmoy = bmesh.from_edit_mesh(yo.data)
bmoz = bmesh.from_edit_mesh(zo.data)

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

olax = [
    'x',
    'y',
    'z'
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

bmo.verts.ensure_lookup_table()
bmo.edges.ensure_lookup_table()
bmo.faces.ensure_lookup_table()

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
        bmov[mid] = last_mid

    bpy.context.active_object.data.update()
        

bpy.ops.object.editmode_toggle()
        


