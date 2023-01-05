--############################################### Setup mesh
chiMeshHandlerCreate()

mesh={}
N=2
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
 chiVolumeMesherSetMatIDToAll(0)


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
materials = {}
materials[1] = chiPhysicsAddMaterial("Mat_outer");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)

num_groups = 2
chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
        CHI_XSFILE,"../my_test/xs_2g_mat1up.cxs")
num_groups = 168
chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
        CHI_XSFILE,"../ChiTest/xs_graphite_pure.cxs")

chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)

--############################################### Add external src
src={}
for g=1,num_groups do
    src[g] = 0.0
end
src[1] = 1.0
chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)

--############################################### Add material properties
--#### CFEM stuff
phys1 = chiCFEMMGDiffusionSolverCreate()

chiSolverSetBasicOption(phys1, "residual_tolerance", 1E-8)
chiSolverSetBasicOption(phys1, "thermal_flux_error", 1E-7)
chiSolverSetBasicOption(phys1, "max_thermal_iters", 20)
chiSolverSetBasicOption(phys1, "verbose_level", 0)

chiCFEMMGDiffusionSetBCProperty(phys1,"boundary_type",e_bndry,"reflecting")
chiCFEMMGDiffusionSetBCProperty(phys1,"boundary_type",w_bndry,"reflecting")
chiCFEMMGDiffusionSetBCProperty(phys1,"boundary_type",n_bndry,"reflecting")
chiCFEMMGDiffusionSetBCProperty(phys1,"boundary_type",s_bndry,"reflecting")

chiSolverInitialize(phys1)
chiSolverExecute(phys1)

----############################################### Visualize the field function
fflist,count = chiGetFieldFunctionList(phys1)
-- export to 2 different VTK files. should be changed when new FF are in place
--chiExportFieldFunctionToVTK(fflist[1],"square_up_flx01","Flux_Diff01")
--chiExportFieldFunctionToVTK(fflist[2],"square_up_flx02","Flux_Diff02")

--for g=1,num_groups do
--    g_string=string.format("%03d",g)
--    chiExportFieldFunctionToVTK(fflist[g],"square_up_flx"..g_string,"Flux_Diff"..g_string)
--end