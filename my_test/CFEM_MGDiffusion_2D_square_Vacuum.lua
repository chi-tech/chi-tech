function printTable(t, f)
   local function printTableHelper(obj, cnt)
      local cnt = cnt or 0
      if type(obj) == "table" then
         io.write("\n", string.rep("\t", cnt), "{\n")
         cnt = cnt + 1
         for k,v in pairs(obj) do
            if type(k) == "string" then
               io.write(string.rep("\t",cnt), '["'..k..'"]', ' = ')
            end
            if type(k) == "number" then
               io.write(string.rep("\t",cnt), "["..k.."]", " = ")
            end
            printTableHelper(v, cnt)
            io.write(",\n")
         end
         cnt = cnt-1
         io.write(string.rep("\t", cnt), "}")
      elseif type(obj) == "string" then
         io.write(string.format("%q", obj))
      else
         io.write(tostring(obj))
      end
   end

   if f == nil then
      printTableHelper(t)
   else
      io.output(f)
      io.write("return")
      printTableHelper(t)
      io.output(io.stdout)
   end
end

--############################################### Setup mesh
chiMeshHandlerCreate()

mesh={}
N=100
L=2
xmin = -L/2
dx = L/N
for i=1,(N+1) do
    k=i-1
    mesh[i] = xmin + k*dx
end

chiMeshCreateUnpartitioned2DOrthoMesh(mesh,mesh)
chiVolumeMesherExecute();

--############################################### Set Material IDs
-- chiVolumeMesherSetMatIDToAll(0)
--   +------------------+
--   | matID=0          |
--   |       +---+      |
--   |       | 1 |      |
--   |       +---+      |
--   |                  |
--   +------------------+
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)

vol1 = chiLogicalVolumeCreate(RPP,-0.5,0.5,-0.5,0.5,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol1,1)

-- Setboundary IDs
-- xmin,xmax,ymin,ymax,zmin,zmax
e_vol = chiLogicalVolumeCreate(RPP,0.99999,1000,-1000,1000,-1000,1000)
w_vol = chiLogicalVolumeCreate(RPP,-1000,-0.9999,-1000,1000,-1000,1000)
n_vol = chiLogicalVolumeCreate(RPP,-1000,1000,0.99999,1000,-1000,1000)
s_vol = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,-0.99999,-1000,1000)

e_bndry = 0
w_bndry = 1
n_bndry = 2
s_bndry = 3

chiVolumeMesherSetProperty(BNDRYID_FROMLOGICAL,e_vol,e_bndry)
chiVolumeMesherSetProperty(BNDRYID_FROMLOGICAL,w_vol,w_bndry)
chiVolumeMesherSetProperty(BNDRYID_FROMLOGICAL,n_vol,n_bndry)
chiVolumeMesherSetProperty(BNDRYID_FROMLOGICAL,s_vol,s_bndry)

--############################################### Add materials
-- matID=0 is outer
-- matID=1 is inner
materials = {}
materials[1] = chiPhysicsAddMaterial("Mat_outer");
materials[2] = chiPhysicsAddMaterial("Mat_inner");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[2],TRANSPORT_XSECTIONS)

num_groups = 2
chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
        CHI_XSFILE,"../my_test/xs_2g_mat1.cxs")
--         CHI_XSFILE,"../my_test/xs_2g_downonly.cxs")
chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,
        CHI_XSFILE,"../my_test/xs_2g_mat2.cxs")
--         CHI_XSFILE,"../my_test/xs_2g_downonly.cxs")

chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)
chiPhysicsMaterialAddProperty(materials[2],ISOTROPIC_MG_SOURCE)

--############################################### Add external src
src={}
for g=1,num_groups do
    src[g] = 0.0
end
--  set source in fast group (1) for the inner material (=mat2)
src[1] = 1.0
printTable(src)
chiPhysicsMaterialSetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
-- reset source to  0 for outer material (=mat1)
src[1] = 0.0
printTable(src)
chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)

--############################################### Add material properties
--#### CFEM stuff
phys1 = chiCFEMMGDiffusionSolverCreate()

chiSolverSetBasicOption(phys1, "residual_tolerance", 1E-8)

chiCFEMMGDiffusionSetBCProperty(phys1,"boundary_type",e_bndry,"vacuum")
chiCFEMMGDiffusionSetBCProperty(phys1,"boundary_type",w_bndry,"vacuum")
chiCFEMMGDiffusionSetBCProperty(phys1,"boundary_type",n_bndry,"vacuum")
chiCFEMMGDiffusionSetBCProperty(phys1,"boundary_type",s_bndry,"vacuum")

chiSolverInitialize(phys1)
chiSolverExecute(phys1)

----############################################### Visualize the field function
fflist,count = chiGetFieldFunctionList(phys1)
-- export to 2 different VTK files. should be changed when new FF are in place
chiExportFieldFunctionToVTK(fflist[1],"square_dir1","Flux_Diff01")
chiExportFieldFunctionToVTK(fflist[2],"square_dir2","Flux_Diff02")