--############################################### Setup mesh
chiMeshHandlerCreate()

mesh={}
N=1000
L=4.0
xmin = 0.0
dx = L/N
for i=1,(N+1) do
    k=i-1
    mesh[i] = xmin + k*dx
end
line_mesh = chiLineMeshCreateFromArray(mesh)


region1 = chiRegionCreate()
chiRegionAddLineBoundary(region1,line_mesh);

-- Create meshers
chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED);
chiVolumeMesherCreate(VOLUMEMESHER_LINEMESH1D);

-- Execute meshing
chiSurfaceMesherExecute();
chiVolumeMesherExecute();

-- Set Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)

--############################################### Add material
material0 = chiPhysicsAddMaterial("Test Material");

chiPhysicsMaterialAddProperty(material0,TRANSPORT_XSECTIONS)

chiPhysicsMaterialSetProperty(material0,
                              TRANSPORT_XSECTIONS,
                              SIMPLEXS1,
                              1,     --Num grps
                              1.0,   --Sigma_t
                              1.0)   --Scattering ratio

--############################################### Setup Physics
phys1 = chiNPTransportCreateSolver()
chiSolverAddRegion(phys1,region1)

chiNPTCreateGroup(phys1)

pqaud = chiCreateProductQuadrature(GAUSS_LEGENDRE,32)

--========== Groupset def
gs0 = chiNPTCreateGroupset(phys1)
chiNPTGroupsetAddGroups(phys1,gs0,0,0)
chiNPTGroupsetSetQuadrature(phys1,gs0,pquad)

--========== Boundary conditions
bsrc = {0.5}
chiNPTSetProperty(phys1,BOUNDARY_CONDITION,ZMIN,INCIDENT_ISOTROPIC,bsrc);

--========== Solvers
chiNPTSetProperty(phys1,DISCRETIZATION_METHOD,PWLD1D)
chiNPTSetProperty(phys1,SCATTERING_ORDER,0)

chiNPTInitialize(phys1)
chiNPTExecute(phys1)

--############################################### Setup Output
fflist,count = chiNPTGetScalarFieldFunctionList(phys1)

cline = chiFFInterpolationCreate(LINE)
chiFFInterpolationSetProperty(cline,LINE_FIRSTPOINT,0.0,0.0,0.0+xmin)
chiFFInterpolationSetProperty(cline,LINE_SECONDPOINT,0.0,0.0, L+xmin)
chiFFInterpolationSetProperty(cline,LINE_NUMBEROFPOINTS, 50)

chiFFInterpolationSetProperty(cline,ADD_FIELDFUNCTION,fflist[1])


chiFFInterpolationInitialize(cline)
chiFFInterpolationExecute(cline)
chiFFInterpolationExportPython(cline)

if (chi_location_id == 0) then
    local handle = io.popen("python ZLFFI00.py")
end


