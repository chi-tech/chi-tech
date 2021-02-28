dofile(CHI_LIBRARY)
--dofile(CHI_WORLD_DEFAULTFLOOR)

print(chiArgs)

mainCamera = chilCreateRevolverCamera("MainCamera")
chilCameraOrganizer.AddCamera(mainCamera)

mesher = chiSurfaceRemesherCreate();

if (chiArgs[3]==nil) then
    test=1;
else
    test=tonumber(chiArgs[3]);
    if (test==nil) then
        print("Invalid test argument")
        return;
    end
end

if (test==1) then
    surf=chiObjectLoadSurface("CHI_RESOURCES/TestObjects/TestSurface3_multiface.obj");
elseif (test==2) then
    surf=chiObjectLoadSurface("CHI_RESOURCES/TestObjects/TestSurface4.obj");
elseif (test==3) then
    surf=chiObjectLoadSurface("CHI_RESOURCES/TestObjects/TestSurface3_multiface_rotated.obj");
elseif (test==4) then
    surf=chiObjectLoadSurface("CHI_RESOURCES/TestObjects/TestSurface2.obj");
elseif (test==5) then
    surf=chiObjectLoadSurface("CHI_RESOURCES/TestObjects/TestSurface5_simplices.obj");
elseif (test==6) then
    surf=chiObjectLoadSurface("CHI_RESOURCES/TestObjects/TestSurface6_cylinder.obj");
elseif (test==7) then
    surf=chiObjectLoadSurface("CHI_RESOURCES/TestObjects/TestSurface7_annulus.obj");
elseif (test==8) then
    surf=chiObjectLoadSurface("CHI_RESOURCES/TestObjects/TestSurface8_sphere.obj");
end


cube=chiObjectCreate("Cube");
chiObjectAddSurface(cube,surf);
chiObjectSetProperty(cube,"Wireframe",true)
chiObjectSetProperty(cube,"Shaded",false)

newTr = chiTransformCreate("CubeTransform");

chiObjectSetProperty("Cube","Transform","CubeTransform");
if ((test==1) or (test>2)) then
    uniscale=0.8;chiTransformSetScale(newTr,uniscale,uniscale,uniscale);
else
    uniscale=3.0;chiTransformSetScale(newTr,uniscale,uniscale,uniscale);
end



if ((test==1) or (test>2)) then
    chiSurfaceRemesherSetProperty(mesher,1,0.3);
    chiSurfaceRemesherSetProperty(mesher,2,0.3);
    chiSurfaceRemesherSetProperty(mesher,3,true);
    chiSurfaceRemesherSetProperty(mesher,4,false); --Keep2D orientation
    chiSurfaceRemesherSetProperty(mesher,5,false); --Recalculate normal
    newSurf = chiSurfaceRemesherExecuteMeshing(mesher,surf);
else
    chiSurfaceRemesherSetProperty(mesher,1,0.1);
    chiSurfaceRemesherSetProperty(mesher,2,0.01);
    chiSurfaceRemesherSetProperty(mesher,3,false);
    chiSurfaceRemesherSetProperty(mesher,4,false);
    chiSurfaceRemesherSetProperty(mesher,5,true); --Recalculate normal
    newSurf = chiSurfaceRemesherExecuteMeshing(mesher,surf);
end


chiObjectAddSurface(cube,newSurf);

              