import bpy
import bmesh

bl_info = {
    "name": "Studio Nois Prep",
    "description": "Automatization of common actions",
    "author": "Kim Falck",
    "version": (0, 0, 1),
    "blender": (4, 3, 0),
    "location": "UV > Sidebar",
    "wiki_url": "...",
    "category": "Operators",
}

class StudioNoisSimplePrep(bpy.types.Operator):
    bl_idname = "operator.studionois_simpleprep"
    bl_label = "Simple FBX Prep"
    bl_description = "Travel along edges until all vertices are visited and sort vertices according to the distance traveled"
    bl_options = {'REGISTER', 'UNDO'}

    iP = bpy.props.IntProperty
    fP = bpy.props.FloatProperty
    fvP = bpy.props.FloatVectorProperty
    bP = bpy.props.BoolProperty
    eP = bpy.props.EnumProperty
    bvP = bpy.props.BoolVectorProperty
    sP = bpy.props.StringProperty

    clear_z_rot : bP(
        name = "Clear Z Rotation",
        description = "Clear rotation in Z",
        default = True
        ) # type: ignore

    apply_rot_scale : bP(
        name = "Apply Scale & Rotation",
        description = "Apply scale and rotation",
        default = True
        ) # type: ignore
    
    apply_material : bP(
        name = "Apply Material",
        description = "Apply specific material",
        default = True
        ) # type: ignore

    material_to_apply : eP(
        name = "Material:",
        description = "Material to apply"
        ) # type: ignore
    
    apply_subdiv : bP(
        name = "Apply Subdivision",
        description = "Apply 2 steps subdivision",
        default = True
        ) # type: ignore
    
    apply_seams_from_islands : bP(
        name = "Seams from Islands",
        description = "Generate seams from islands",
        default = True
        ) # type: ignore
    
    def execute(self, context):

        self.material_to_apply.items = [(str(mat_i), str(mat.name), "") for mat_i, mat in enumerate(bpy.data.materials)]
        for mat_i, mat in enumerate(bpy.data.materials):
            if mat.name == "ColorGrid" or mat.name == "colorgrid" or mat.name == "color_grid" or mat.name == "Color_Grid":
                self.material_to_apply = str(mat_i)
        mat_to_apply = bpy.data.materials[int(self.material_to_apply)]
        
        
        for obj in bpy.context.view_layer.objects:
            if self.clear_z_rot == True:
                obj.rotation_euler[2] = 0
            try:
                obj.select_set(True)
            except:
                pass
        bpy.ops.object.parent_clear(type='CLEAR_KEEP_TRANSFORM')
        for obj in bpy.context.scene.objects:
            if obj.type == 'EMPTY':
                obj.select_set(True)
            else:
                try:
                    obj.select_set(False)
                except:
                    pass
        bpy.ops.object.delete(use_global=False, confirm=False)
        for obj in bpy.context.scene.objects:
            if obj.type == 'MESH':
                if self.apply_material == True:
                    for matslot in obj.material_slots:
                        matslot.material = mat_to_apply
                try:
                    obj.select_set(True)
                    bpy.context.view_layer.objects.active = obj
                except:
                    pass
        

        if self.apply_rot_scale == True:
            bpy.ops.object.transform_apply(location=False, rotation=True, scale=True)
        if self.apply_subdiv == True:
            bpy.ops.object.subdivision_set(level=2, relative=False)
        
        # bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.object.editmode_toggle()
        # bpy.ops.uv.select_mode(type='VERTEX')
        for eobj in bpy.context.objects_in_mode:
            bm = bmesh.from_edit_mesh(eobj.data)
            uv_layer = bm.loops.layers.uv.verify()
            bm.verts.ensure_lookup_table()
            bm.edges.ensure_lookup_table()
            bm.faces.ensure_lookup_table()
            
            for face in bm.faces:
                face.select = True
                for loop in face.loops:
                    loop[uv_layer].select = True
                    loop[uv_layer].select_edge = True
            eobj.data.update()
        if self.apply_seams_from_islands == True:
            bpy.ops.uv.seams_from_islands()

        for eobj in bpy.context.objects_in_mode:
            bm = bmesh.from_edit_mesh(eobj.data)
            uv_layer = bm.loops.layers.uv.verify()
            bm.verts.ensure_lookup_table()
            bm.edges.ensure_lookup_table()
            bm.faces.ensure_lookup_table()
            
            for face in bm.faces:
                face.select = False
                for loop in face.loops:
                    loop[uv_layer].select = False
                    loop[uv_layer].select_edge = True
            
            eobj.data.update()
        return {"FINISHED"}

classes = [StudioNoisSimplePrep
           ]

def register():
    StudioNoisSimplePrep.enum_items = [(str(mat_i), str(mat.name), "") for mat_i, mat in enumerate(bpy.data.materials)]
    for cls in classes:
        try: 
            bpy.utils.register_class(cls)
        except:
            bpy.utils.unregister_class(cls)
            bpy.utils.register_class(cls)

        
def unregister():
    for cls in classes:
        bpy.utils.unregister_class(cls)
        
if __name__ == "__main__":
    register()