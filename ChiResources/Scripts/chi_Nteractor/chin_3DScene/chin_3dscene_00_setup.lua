setupDir=chinGetScriptPath();

--========================================================= Load Floor
objname="Floor";
chiObjectCreate(objname);

--Surface
surf=chiLoadSurface(setupDir.."/Assets/Meshes/D01_FloorTile.obj");
chiObjectAddSurface(objname,surf);

--Material
chiMaterialCreate(objname .. "_Material");
newText=chiLoadTexture(setupDir.."/Assets/Textures/C01_FloorTile.png");
chiMaterialSetProperty(objname .. "_Material","DiffuseTexture",newText);

--Transform
chiTransformCreate(objname.."_Transform");

--Texture transform
chiTransformCreate(objname.."_TextureTransform");

--Assignments
chiObjectSetProperty(objname,"Material",objname .. "_Material");
chiObjectSetProperty(objname,"Transform"       ,objname.."_Transform");
chiObjectSetProperty(objname,"TextureTransform",objname.."_TextureTransform");

chiTransformSetScale(objname.."_Transform",10.0,10.0,1.0);
chiTransformSetScale(objname.."_TextureTransform",10.0,10.0,1.0);

--========================================================= Set light
chiLightCreate("DefaultLight");
ambient=0.0;
chiLightSetProperty("DefaultLight","Ambient",ambient,ambient,ambient)
chiLightSetProperty("DefaultLight","Attenuation",1.0,0.000,0.0)
chiLightSetProperty("DefaultLight","Position",0.0,15.0,10.0)





