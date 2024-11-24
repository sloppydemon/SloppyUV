import bpy
import bpy_extras
from bpy_extras import bmesh_utils
import bmesh
import math
import mathutils
import bl_math
import sys

class ProceduralQuadUVUnfold(bpy.types.Operator):
    bl_idname = "operator.uv_quad_unfold"
    bl_label = "Quad Unfold UVs"
    bl_description = "Attempt to unwrap quad mesh UVs as if from flat piece of cloth, with the assumption that all face corners should be right-angled originally"
    bl_options = {'REGISTER', 'UNDO'}

    iP = bpy.props.IntProperty
    fP = bpy.props.FloatProperty
    fvP = bpy.props.FloatVectorProperty
    bP = bpy.props.BoolProperty
    eP = bpy.props.EnumProperty
    bvP = bpy.props.BoolVectorProperty
    sP = bpy.props.StringProperty

    unfold_mode : eP(
        name = "Mode",
        description = "Unfolding Mode",
        items = [
            ("A", "All Directions", "Work outwards from initial quad"),
            ("B", "X First", "Work outwards in (local) X first until no more faces are found, calculating edge lengths to use in next stage, then work outwards in Y. repeating first step when more faces are found in X"),
            ("C", "Y First", "Work outwards in (local) Y first until no more faces are found, calculating edge lengths to use in next stage, then work outwards in X. repeating first step when more faces are found in Y")
            ]
        ) # type: ignore

    initial_quad : eP(
        name = "Initial Quad(s)",
        description = "Which quad (for each UV island, if there are more than one,) to start from",
        items = [
            ("A", "Most Regular", "Start with quad that is flattest and/or that has normal most similar to island average normals"),
            ("B", "Selected", "First selected quad in each island (if an island has no faces selected, initial quad falls back to Most Regular)")
            ]
        ) # type: ignore

    reg_flat_fac : fP(
        name = "Flatness Factor",
        description = "How much flatness influences initial quad selection",
        default = 1.0,
        min = 0.0,
        max = 1.0
        ) # type: ignore

    reg_norm_fac : fP(
        name = "Normal Factor",
        description = "How much island average normal influences initial quad selection",
        default = 1.0,
        min = 0.0,
        max = 1.0
        ) # type: ignore

    quant_avg_norm : bP(
        name = "Quantize Average Normal",
        description = "Quantize island's averaged normal to up, down, right, left",
        default = False
        ) # type: ignore

    def execute(self, context):
        props = context.scene.sloppy_props
        bm = bmesh.from_edit_mesh(bpy.context.active_object.data)
        uv_layer = bm.loops.layers.uv.verify()
        islands = bmesh_utils.bmesh_linked_uv_islands(bm, uv_layer)

        current_dir_array = []
        max_avg_angle = 0.0
        min_avg_angle = 360.0
        current_avg_normal = mathutils.Vector((0,0,0))
        current_normal = mathutils.Vector((0,0,0))

        attr_dict = [
            {"name": "naboL", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "naboR", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "naboU", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "naboD", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "xi", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "yi", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "island_index", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "flatness", "type": "FLOAT", "domain": "FACE", "layer": None},
            {"name": "normal_regularity", "type": "FLOAT", "domain": "FACE", "layer": None},
            {"name": "nom", "type": "INT", "domain": "FACE", "layer": None},
            {"name": "nom_dir", "type": "FLOAT_VECTOR", "domain": "FACE", "layer": None},
            {"name": "prev_length_x", "type": "FLOAT", "domain": "FACE", "layer": None},
            {"name": "prev_length_y", "type": "FLOAT", "domain": "FACE", "layer": None},
            ]
        
        for nam in attr_dict:
            attr = props.find_or_add_attribute(nam["name"], nam["type"], nam["domain"])
            nam["layer"] = props.get_attribute_layer(nam["name"], nam["type"], nam["domain"], bm)
        
        naboL = props.get_dict_layer("naboL", attr_dict)
        naboR = props.get_dict_layer("naboR", attr_dict)
        naboU = props.get_dict_layer("naboU", attr_dict)
        flatness = props.get_dict_layer("flatness", attr_dict)
        normal_regularity = props.get_dict_layer("normal_regularity", attr_dict)
        island_index = props.get_dict_layer("island_index", attr_dict)

        quant_arr = [
            [mathutils.Vector((1,0,0)), 
                [
                    mathutils.Vector((0,0,-1)),
                    mathutils.Vector((0,-1,0)),
                    mathutils.Vector((0,0,1)),
                    mathutils.Vector((0,1,0))
                ],
                "Left"
             ],
             [mathutils.Vector((-1,0,0)),
                [
                    mathutils.Vector((0,0,-1)),
                    mathutils.Vector((0,1,0)),
                    mathutils.Vector((0,0,1)),
                    mathutils.Vector((0,-1,0))
                ],
                "Right"
             ],
             [mathutils.Vector((0,1,0)),
                [
                    mathutils.Vector((0,0,-1)),
                    mathutils.Vector((1,0,0)),
                    mathutils.Vector((0,0,1)),
                    mathutils.Vector((-1,0,0))
                ],
                "Towards"
             ],
             [mathutils.Vector((0,-1,0)),
                [
                    mathutils.Vector((0,0,-1)),
                    mathutils.Vector((-1,0,0)),
                    mathutils.Vector((0,0,1)),
                    mathutils.Vector((1,0,0))
                ],
                "Away"
             ],
             [mathutils.Vector((0,0,-1)),
                [
                    mathutils.Vector((0,-1,0)),
                    mathutils.Vector((1,0,0)),
                    mathutils.Vector((0,1,0)),
                    mathutils.Vector((-1,0,0))
                ],
                "Up"
             ],
             [mathutils.Vector((0,0,1)),
                [
                    mathutils.Vector((0,1,0)),
                    mathutils.Vector((1,0,0)),
                    mathutils.Vector((0,-1,0)),
                    mathutils.Vector((-1,0,0))
                ],
                "Down"
             ],
        ]

        def quant_sort(e):
            quant_dot = e[0].dot(current_normal)
            return quant_dot

        def regularity_sort(e):
            flat = (1.0 - e[flatness]) * -self.reg_flat_fac
            norm_reg = e[normal_regularity] * -self.reg_norm_fac
            return flat + norm_reg

        for i, island in enumerate(islands):
            faces_remain = []
            faces_selected = []
            init_face = None
            avg_norm = mathutils.Vector((0,0,0))

            for face in island:
                avg_norm += face.normal
                avg_angle = 0.0
                angle_count = 0
                for fe in face.edges:
                    angle_count += 1
                    if fe.is_boundary == True:
                        avg_angle += 90
                    else:
                        avg_angle += math.degrees(fe.calc_face_angle())
                avg_angle /= angle_count
                new_max_angle = max(avg_angle, max_avg_angle)
                new_min_angle = min(avg_angle, min_avg_angle)
                max_avg_angle = new_max_angle
                min_avg_angle = new_min_angle
                face[flatness] = avg_angle
                faces_remain.append(face)
                if face.select == True:
                    faces_selected.append(face)
            
            current_avg_normal = avg_norm.normalized()

            if self.quant_avg_norm == True:
                current_normal = current_avg_normal
                quant_arr.sort(key=quant_sort)
                current_avg_normal = quant_arr[0][0] * -1.0
                if props.verbose == True:
                    print(f"Average island normals quantized to {quant_arr[0][2]}")

            bpy.context.active_object.data.update()

            for face in island:
                new_flatness = props.remap_val(face[flatness], min_avg_angle, max_avg_angle, 0.0, 1.0)
                face[flatness] = new_flatness
                n_dot_avg = face.normal.dot(current_avg_normal)
                face[normal_regularity] = props.remap_val(n_dot_avg, -1.0, 1.0, 0.0, 1.0)
                face[island_index] = i
            
            bpy.context.active_object.data.update()

            is_selected_face = False

            if self.initial_quad == "B":
                if len(faces_selected) > 0:
                    init_face = faces_selected[0]
                    is_selected_face = True
                else:
                    faces_remain.sort(key=regularity_sort)
                    init_face = faces_remain[0]
            else:
                faces_remain.sort(key=regularity_sort)
                init_face = faces_remain[0]
            
            current_normal = init_face.normal

            quant_arr.sort(key=quant_sort)
            current_dir_array = quant_arr[0][1]

            if props.verbose == True:
                msg_start = "Most regular face for island "
                if is_selected_face == True:
                    msg_start = "Selected face for island "
                msg = "{:}{:}: {:} with a flatness of {:.3f} and a normal regularity of {:.3f}.\nFacing quantized to {:}\n".format(msg_start, i, init_face.index, init_face[flatness], init_face[normal_regularity], quant_arr[0][2])
                print(msg)
            
            init_face.select = True

        bpy.context.active_object.data.update()

        return {"FINISHED"}

# class SloppyDeTriangulate(bpy.types.Operator):
#     bl_idname = "operator.sloppy_detri"
#     bl_label = "Detriangulate From Tri Pair"
#     bl_description = "Detriangulate mesh using selected tri pair as template"
#     bl_options = {'REGISTER', 'UNDO'}