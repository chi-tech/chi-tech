scriptPath = chinGetScriptPath()

chiDefaultFloorObj  = chiObjectCreate("CHI_WORLD_DEFAULTFLOOR")
chiDefaultFloorSurf = chiObjectLoadSurface(scriptPath.."/obj/D01_FloorTile.obj")
chiObjectAddSurface(chiDefaultFloorObj,chiDefaultFloorSurf)

chiDefaultFloorMat  = chiMaterialCreate("CHI_WORLD_DEFAULTFLOOR")
chiDefaultFloorTex  = chiTextureLoad(scriptPath.."/textures/C01_FloorTile.png")
chiMaterialSetProperty(chiDefaultFloorMat,"DiffuseTexture",chiDefaultFloorTex)
chiMaterialSetProperty(chiDefaultFloorMat,"DiffuseTextureEnabled",true)

chiObjectSetProperty(chiDefaultFloorObj,"Material",chiDefaultFloorMat)
