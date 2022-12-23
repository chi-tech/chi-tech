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

-- vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
-- chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)
--
-- vol1 = chiLogicalVolumeCreate(RPP,-0.5,0.5,-0.5,0.5,-1000,1000)
-- chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol1,1)

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
--############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Test Material1");
-- materials[2] = chiPhysicsAddMaterial("Test Material2");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
-- chiPhysicsMaterialAddProperty(materials[2],TRANSPORT_XSECTIONS)

chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)
-- chiPhysicsMaterialAddProperty(materials[2],ISOTROPIC_MG_SOURCE)

num_groups = 2
chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
        CHI_XSFILE,"../my_test/xs_2g_downonly.cxs")
-- chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,
--         CHI_XSFILE,"../my_test/xs_2g_downonly.cxs")

src={}
for g=1,num_groups do
    src[g] = 1.0
end
chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
-- src[1] = 1.0
-- chiPhysicsMaterialSetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)

--############################################### Add material properties
--#### CFEM stuff
phys1 = chiCFEMMGDiffusionSolverCreate()

chiSolverSetBasicOption(phys1, "residual_tolerance", 1E-8)

chiCFEMMGDiffusionSetBCProperty(phys1,"boundary_type",e_bndry,"dirichlet",0.0)
chiCFEMMGDiffusionSetBCProperty(phys1,"boundary_type",w_bndry,"dirichlet",0.0)
chiCFEMMGDiffusionSetBCProperty(phys1,"boundary_type",n_bndry,"dirichlet",0.0)
chiCFEMMGDiffusionSetBCProperty(phys1,"boundary_type",s_bndry,"dirichlet",0.0)

chiSolverInitialize(phys1)
chiSolverExecute(phys1)

----############################################### Visualize the field function
fflist,count = chiGetFieldFunctionList(phys1)
chiExportFieldFunctionToVTK(fflist[1],"square_dir","Flux_Diff")