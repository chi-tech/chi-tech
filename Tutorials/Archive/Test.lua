print("############################################### LuaTest")
dofile(CHI_LIBRARY)



--############################################### Setup mesh
chiMeshHandlerCreate()

newSurfMesh = chiSurfaceMeshCreate();
chiSurfaceMeshImportFromOBJFile(newSurfMesh,"CHI_RESOURCES/TestObjects/TestSurface5_simplices_MESHED.obj")
loops,loop_count = chiSurfaceMeshGetEdgeLoops(newSurfMesh)

line_mesh = {};
line_mesh_count = 0;

for k=1,loop_count do
  split_loops,split_count = chiEdgeLoopSplitByAngle(loops,k-1);
  for m=1,split_count do
    line_mesh_count = line_mesh_count + 1;
    line_mesh[line_mesh_count] = chiLineMeshCreateFromLoop(split_loops,m-1);
  end

end

chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED);
chiVolumeMesherCreate(VOLUMEMESHER_PREDEFINED2D);


--############################################### Setup Regions

region1 = chiRegionCreate()
chiRegionAddSurfaceBoundary(region1,newSurfMesh);
for k=1,line_mesh_count do
  chiRegionAddLineBoundary(region1,line_mesh[k]);
end

--############################################### Execute meshing
chiSurfaceMesherExecute();
chiVolumeMesherExecute();


--############################################### Setup Physics

phys1 = chiDiffusionCreateSolver();
chiSolverAddRegion(phys1,region1)
chiDiffusionSetProperty(phys1,DISCRETIZATION_METHOD,CFEM2D);
chiDiffusionInitialize(phys1)

chiDiffusionExecute(phys1)


--############################################### Display interface
mainCamera = chilCreateRevolverCamera("MainCamera")
chilCameraOrganizer.AddCamera(mainCamera)
chiSetWindowProperties(400,400,3000,200)
surf = chiObjectLoadSurfaceFromSurfaceMesh(newSurfMesh);

cube=chiObjectCreate("Cube");
chiObjectAddSurface(cube,surf);
chiObjectSetProperty(cube,"Wireframe",true)
chiObjectSetProperty(cube,"Shaded",false)
chiObjectSetProperty(cube,"BackFaceCull",false)

line3d = {};
for k=1,line_mesh_count do
  name = string.format("Boundary Line %d",k);
  line3d[k] = chi3DLineCreateFromLineMesh(name,line_mesh[k]);
  chi3DLineSetStipple(line3d[k],false,1.0,1,10);
end

--================================= 2D displayscene
displayer1 = chiDisplayerCreate();
chiBindScene(0,displayer1);

mainCamera = chilCreateOrthoWindowCamera("MainCamera")
chilCameraOrganizer.AddCamera(mainCamera);

text1 = chiTextCreate("Text");
chiTextSetProperty(text1,TEXT_VALUE,"Hello World");
chiTextSetProperty(text1,TEXT_POSITION,10,10,0);

--================================= Event handler
EventHandler={}
EventHandler.callbackFunction = function(this)
  if (WM_MOUSEMOVE.occured) then
    name = WM_MOUSEMOVE.sPar0;
    chiTextSetProperty(text1,TEXT_VALUE,name);

    for k=1,line_mesh_count do
      if (WM_MOUSEMOVE.iPar5==line3d[k]) then
        chi3DLineSetStipple(line3d[k],false,1.0,4,10);
      else
        chi3DLineSetStipple(line3d[k],false,1.0,1,10);
      end
    end
  end

end

callBackObj = chilCallbacks.MakeCallback(EventHandler);
chilCallbacks.PushCallback(callBackObj);


print("############################################### Test file loaded")
