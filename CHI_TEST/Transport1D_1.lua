chiMPIBarrier()
if (chi_location_id == 0) then
    print("############################################### LuaTest")
end
--dofile(CHI_LIBRARY)



--############################################### Setup mesh
chiMeshHandlerCreate()

mesh={}
N=100
L=30.0
xmin = 0.0
dx = L/N
for i=1,(N+1) do
    k=i-1
    mesh[i] = xmin + k*dx
end
line_mesh = chiLineMeshCreateFromArray(mesh)


region1 = chiRegionCreate()
chiRegionAddLineBoundary(region1,line_mesh);


--############################################### Create meshers
chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED);
chiVolumeMesherCreate(VOLUMEMESHER_LINEMESH1D);

chiVolumeMesherSetProperty(PARTITION_Z,3)

--############################################### Execute meshing
chiSurfaceMesherExecute();
chiVolumeMesherExecute();

--############################################### Set Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)


--############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Test Material");
materials[2] = chiPhysicsAddMaterial("Test Material2");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[2],TRANSPORT_XSECTIONS)

chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)
chiPhysicsMaterialAddProperty(materials[2],ISOTROPIC_MG_SOURCE)


num_groups = 168
chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
        PDT_XSFILE,"CHI_RESOURCES/CrossSections/xs_3_170.data")
chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,
        PDT_XSFILE,"CHI_RESOURCES/CrossSections/xs_3_170.data")

--chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,SIMPLEXS0,num_groups,0.1)
--chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,SIMPLEXS0,num_groups,0.1)

src={}
for g=1,num_groups do
    src[g] = 0.0
end
--src[1] = 1.0
chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
chiPhysicsMaterialSetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)



--############################################### Setup Physics

phys1 = chiLBSCreateSolver()
chiSolverAddRegion(phys1,region1)

--========== Groups
grp = {}
for g=1,num_groups do
    grp[g] = chiLBSCreateGroup(phys1)
end

--========== ProdQuad
pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE,40)

--========== Groupset def
gs0 = chiLBSCreateGroupset(phys1)
chiLBSGroupsetAddGroups(phys1,gs0,0,62)
chiLBSGroupsetSetQuadrature(phys1,gs0,pquad)

gs1 = chiLBSCreateGroupset(phys1)
chiLBSGroupsetAddGroups(phys1,gs1,63,167)
chiLBSGroupsetSetQuadrature(phys1,gs1,pquad)

--========== Boundary conditions
bsrc={}
for g=1,num_groups do
    bsrc[g] = 0.0
end
bsrc[1] = 1.0/2
chiLBSSetProperty(phys1,BOUNDARY_CONDITION,ZMIN,INCIDENT_ISOTROPIC,bsrc);

--========== Solvers
chiLBSSetProperty(phys1,PARTITION_METHOD,FROM_SURFACE)
chiLBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD3D)
chiLBSSetProperty(phys1,SCATTERING_ORDER,5)
--chiLBSSetProperty(phys1,GROUPSET_ITERATIVEMETHOD,gs0,LBS_CLASSICRICHARDSON)
--chiLBSSetProperty(phys1,GROUPSET_ITERATIVEMETHOD,gs1,LBS_CLASSICRICHARDSON)
chiLBSSetProperty(phys1,GROUPSET_TOLERANCE,gs0,1.0e-6)
--chiLBSSetProperty(phys1,GROUPSET_MAXITERATIONS,gs0,3)
chiLBSSetProperty(phys1,GROUPSET_GMRESRESTART_INTVL,gs0,100)
chiLBSSetProperty(phys1,GROUPSET_GMRESRESTART_INTVL,gs1,100)
chiLBSSetProperty(phys1,GROUPSET_SUBSETS,gs0,8)
chiLBSSetProperty(phys1,GROUPSET_SUBSETS,gs1,8)
--chiLBSSetProperty(phys1,GROUPSET_WGDSA,gs1,true)
--chiLBSSetProperty(phys1,GROUPSET_TGDSA,gs1,true)

chiLBSInitialize(phys1)
chiLBSExecute(phys1)

fflist,count = chiLBSGetScalarFieldFunctionList(phys1)

--Testing consolidated interpolation
cline = chiFFInterpolationCreate(LINE)
chiFFInterpolationSetProperty(cline,LINE_FIRSTPOINT,0.0,0.0,0.0001+xmin)
chiFFInterpolationSetProperty(cline,LINE_SECONDPOINT,0.0,0.0, 29.999+xmin)
chiFFInterpolationSetProperty(cline,LINE_NUMBEROFPOINTS, 50)

for k=165,165 do
    chiFFInterpolationSetProperty(cline,ADD_FIELDFUNCTION,fflist[k])
end


chiFFInterpolationInitialize(cline)
chiFFInterpolationExecute(cline)


ffi1 = chiFFInterpolationCreate(VOLUME)
curffi = ffi1
chiFFInterpolationSetProperty(curffi,OPERATION,OP_MAX)
chiFFInterpolationSetProperty(curffi,LOGICAL_VOLUME,vol0)
chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fflist[1])

chiFFInterpolationInitialize(curffi)
chiFFInterpolationExecute(curffi)
maxval = chiFFInterpolationGetValue(curffi)

chiLog(LOG_0,string.format("Max-value1=%.5f", maxval))

ffi2 = chiFFInterpolationCreate(VOLUME)
curffi = ffi2
chiFFInterpolationSetProperty(curffi,OPERATION,OP_MAX)
chiFFInterpolationSetProperty(curffi,LOGICAL_VOLUME,vol0)
chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fflist[160])

chiFFInterpolationInitialize(curffi)
chiFFInterpolationExecute(curffi)
maxval = chiFFInterpolationGetValue(curffi)

chiLog(LOG_0,string.format("Max-value2=%.5e", maxval))


if (chi_location_id == 0 and master_export == nil) then
    chiFFInterpolationExportPython(cline)
    local handle = io.popen("python ZLFFI00.py")
end
