chiMPIBarrier()
if (chi_location_id == 0) then
    print("############################################### LuaTest")
end
--dofile(CHI_LIBRARY)

--############################################### Setup Transport mesh
    function Create1DMesh(N,L)
        tmesh_handle = chiMeshHandlerCreate()

        mesh={}
        xmin = 0.0
        dx = L/N
        for i=1,(N+1) do
            k=i-1
            mesh[i] = xmin + k*dx
        end
        line_mesh = chiLineMeshCreateFromArray(mesh)


        region0 = chiRegionCreate()
        chiRegionAddLineBoundary(region1,line_mesh);


        --############################################### Create meshers
        chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED);
        chiVolumeMesherCreate(VOLUMEMESHER_LINEMESH1D);

        if (chi_number_of_processes == 4) then
            chiVolumeMesherSetProperty(PARTITION_Z,4)
        end

        --chiVolumeMesherSetProperty(MESH_GLOBAL,true)

        --############################################### Execute meshing
        chiSurfaceMesherExecute();
        chiVolumeMesherExecute();

        --############################################### Set Material IDs
        vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
        chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)

        --############################################### Set Material IDs
        vol1 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,0.0,2.0)
        chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol1,1)

        return tmesh_handle
    end

    function Create2DMesh()
        mesh_handle = chiMeshHandlerCreate()

        newSurfMesh = chiSurfaceMeshCreate();
        chiSurfaceMeshImportFromOBJFile(newSurfMesh,
                "CHI_RESOURCES/TestObjects/SquareMesh2x2QuadsBlock.obj",true)

        region0 = chiRegionCreate()
        chiRegionAddSurfaceBoundary(region0,newSurfMesh);


        --############################################### Create meshers
        chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED);
        chiVolumeMesherCreate(VOLUMEMESHER_PREDEFINED2D);


        if (chi_number_of_processes == 4) then
            chiSurfaceMesherSetProperty(PARTITION_X,2)
            chiSurfaceMesherSetProperty(PARTITION_Y,2)
            chiSurfaceMesherSetProperty(CUT_X,0.0)
            chiSurfaceMesherSetProperty(CUT_Y,0.0)
        end

        --############################################### Execute meshing
        chiSurfaceMesherExecute();
        chiVolumeMesherExecute();

        --############################################### Set Material IDs
        vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
        chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)

        --############################################### Set Material IDs
        vol1 = chiLogicalVolumeCreate(RPP,-20,0,-20,20,-1000,1000)
        chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol1,1)

        return mesh_handle
    end

    function Create3DMesh()
        mesh_handle = chiMeshHandlerCreate()

        newSurfMesh = chiSurfaceMeshCreate();
        chiSurfaceMeshImportFromOBJFile(newSurfMesh,
                "CHI_RESOURCES/TestObjects/SquareMesh2x2QuadsBlock.obj",true)

        region0 = chiRegionCreate()
        chiRegionAddSurfaceBoundary(region0,newSurfMesh);


        --############################################### Create meshers
        chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED);
        chiVolumeMesherCreate(VOLUMEMESHER_EXTRUDER);


        if (chi_number_of_processes == 4) then
            chiSurfaceMesherSetProperty(PARTITION_X,2)
            chiSurfaceMesherSetProperty(PARTITION_Y,2)
            chiSurfaceMesherSetProperty(CUT_X,0.0)
            chiSurfaceMesherSetProperty(CUT_Y,0.0)
        end

        if (chi_number_of_processes == 8) then
            chiSurfaceMesherSetProperty(PARTITION_X,2)
            chiSurfaceMesherSetProperty(PARTITION_Y,2)
            chiSurfaceMesherSetProperty(CUT_X,0.0)
            chiSurfaceMesherSetProperty(CUT_Y,0.0)

            chiVolumeMesherSetProperty(PARTITION_Z,2)
        end

        chiVolumeMesherSetProperty(EXTRUSION_LAYER,40,8,"Main-layer height 40");--40

        --############################################### Execute meshing
        chiSurfaceMesherExecute();
        chiVolumeMesherExecute();

        --############################################### Set Material IDs
        vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
        chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)

        --############################################### Set Material IDs
        vol1 = chiLogicalVolumeCreate(RPP,-20,0,-20,20,-1000,1000)
        chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol1,1)

        return mesh_handle
    end

if (TWOD==nil and THREED==nil) then
    L=4.0
    tmesh = Create1DMesh(32,L)
    tmesh2= Create1DMesh(32,L)
elseif (TWOD==true) then
    tmesh = Create2DMesh()
    tmesh2= Create2DMesh()
elseif (THREED==true) then
    tmesh = Create3DMesh()
    tmesh2= Create3DMesh()
end



--############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Test Material");
materials[2] = chiPhysicsAddMaterial("Test Material2");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[2],TRANSPORT_XSECTIONS)

chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)
chiPhysicsMaterialAddProperty(materials[2],ISOTROPIC_MG_SOURCE)


num_groups = 1
chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,SIMPLEXS1,1,0.1,0.5)
chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,SIMPLEXS1,1,0.1,0.5)

--chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
--        PDT_XSFILE,"ChiTest/xs_3_170.data")
--chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,
--        PDT_XSFILE,"ChiTest/xs_3_170.data")

--chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,SIMPLEXS0,num_groups,0.1)


src={}
for g=1,num_groups do
    src[g] = 0.0
end
--src[1] = 1.0
chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
src[1] = 1.0
chiPhysicsMaterialSetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)






--############################################### Setup Transport Physics
chiMeshHandlerSetCurrent(tmesh)
phys1 = chiLBSCreateSolver()
chiSolverAddRegion(phys1,region0)

--========== Groups
grp = {}
for g=1,num_groups do
    grp[g] = chiLBSCreateGroup(phys1)
end

--========== ProdQuad
if (TWOD==nil and THREED==nil) then
    pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE,1)
else
    pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,1,1)
end


--========== Groupset def
gs0 = chiLBSCreateGroupset(phys1)
cur_gs = gs0
chiLBSGroupsetAddGroups(phys1,gs0,0,num_groups-1)
chiLBSGroupsetSetQuadrature(phys1,gs0,pquad)
chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,1)
chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,NPT_GMRES)
chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-6)
chiLBSGroupsetSetMaxIterations(phys1,cur_gs,300)
chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,100)
--chiLBSGroupsetSetWGDSA(phys1,cur_gs,30,1.0e-4,false," ")
--chiLBSGroupsetSetTGDSA(phys1,cur_gs,30,1.0e-4,false," ")

--========== Boundary conditions
bsrc={}
for g=1,num_groups do
    bsrc[g] = 0.0
end
bsrc[1] = 1.0
--chiLBSSetProperty(phys1,BOUNDARY_CONDITION,
--        ZMIN,LBSBoundaryTypes.INCIDENT_ISOTROPIC,bsrc);

--========== Solvers
chiLBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD3D)
chiLBSSetProperty(phys1,SCATTERING_ORDER,0)

chiLBSInitialize(phys1)
chiLBSExecute(phys1)

fflist1,count = chiLBSGetScalarFieldFunctionList(phys1)


--############################################### Setup Monte Carlo Physics
chiMeshHandlerSetCurrent(tmesh)
phys0 = chiMonteCarlonCreateSolver()
chiSolverAddRegion(phys0,region0)

--chiMonteCarlonCreateSource(phys0,MC_BNDRY_SRC,1);
--chiMonteCarlonCreateSource(phys0,MC_RESID_SRC,fflist1[1]);

chiMonteCarlonCreateSource(phys0,MCSrcTypes.RESIDUAL,99,fflist1[1],bsrc[1]);

chiMonteCarlonSetProperty(phys0,MCProperties.NUM_UNCOLLIDED_PARTICLES,5e6)
chiMonteCarlonSetProperty(phys0,MCProperties.NUM_PARTICLES,5e6)
chiMonteCarlonSetProperty(phys0,MCProperties.TFC_UPDATE_INTVL,10e3)
chiMonteCarlonSetProperty(phys0,MCProperties.TALLY_MERGE_INTVL,1e6)
chiMonteCarlonSetProperty(phys0,MCProperties.SCATTERING_ORDER,0)
chiMonteCarlonSetProperty(phys0,MCProperties.MONOENERGETIC,true)
chiMonteCarlonSetProperty(phys0,MCProperties.FORCE_ISOTROPIC,true)
chiMonteCarlonSetProperty(phys0,MCProperties.MAKE_PWLD_SOLUTION,true)
if (TWOD==nil and THREED==nil) then
    chiMonteCarlonSetProperty(phys0,MCProperties.TALLY_MULTIPLICATION_FACTOR,1.0/2)
elseif (TWOD==true) then
    chiMonteCarlonSetProperty(phys0,MCProperties.TALLY_MULTIPLICATION_FACTOR,40.0)
elseif (THREED==true) then
    chiMonteCarlonSetProperty(phys0,MCProperties.TALLY_MULTIPLICATION_FACTOR,40.0*40*2)
end

chiMonteCarlonInitialize(phys0)
chiMonteCarlonExecute(phys0)

--############################################### Setup ref Monte Carlo Physics
chiMeshHandlerSetCurrent(tmesh2)
phys2 = chiMonteCarlonCreateSolver()
chiSolverAddRegion(phys2,region1)

chiMonteCarlonCreateSource(phys2,MCSrcTypes.MATERIAL_SRC,1);
--chiMonteCarlonCreateSource(phys2,MC_RESID_SRC,fflist1[1]);

chiMonteCarlonSetProperty(phys2,MCProperties.NUM_PARTICLES,50e6)
chiMonteCarlonSetProperty(phys2,MCProperties.TFC_UPDATE_INTVL,10e3)
chiMonteCarlonSetProperty(phys2,MCProperties.TALLY_MERGE_INTVL,1e6)
chiMonteCarlonSetProperty(phys2,MCProperties.SCATTERING_ORDER,0)
chiMonteCarlonSetProperty(phys2,MCProperties.MONOENERGETIC,true)
chiMonteCarlonSetProperty(phys2,MCProperties.FORCE_ISOTROPIC,true)

chiMonteCarlonSetProperty(phys2,MCProperties.MAKE_PWLD_SOLUTION,true)

if (TWOD==nil and THREED==nil) then
    chiMonteCarlonSetProperty(phys2,MCProperties.TALLY_MULTIPLICATION_FACTOR,2.0)
elseif (TWOD==true) then
    chiMonteCarlonSetProperty(phys2,MCProperties.TALLY_MULTIPLICATION_FACTOR,20*40)
elseif (THREED==true) then
    chiMonteCarlonSetProperty(phys2,MCProperties.TALLY_MULTIPLICATION_FACTOR,40.0*40*20.0)
end


chiMonteCarlonInitialize(phys2)
chiMonteCarlonExecute(phys2)


--############################################### Post processing
fflist0,count0 = chiLBSGetScalarFieldFunctionList(phys1)
fflist1,count1 = chiGetFieldFunctionList(phys0)
fflist2,count2 = chiGetFieldFunctionList(phys2) --Fine mesh MC

if (chi_location_id == 0) then
    print(fflist0[1],count0)
    print(fflist1[1],count1)
    print(fflist2[1],count2)
end

--Testing consolidated interpolation
cline = chiFFInterpolationCreate(LINE)
if (TWOD==nil and THREED==nil) then
    chiFFInterpolationSetProperty(cline,LINE_FIRSTPOINT,0.0,0.0,0.0001)
    chiFFInterpolationSetProperty(cline,LINE_SECONDPOINT,0.0,0.0, L-0.0001)
elseif (TWOD==true) then
    chiFFInterpolationSetProperty(cline,LINE_FIRSTPOINT ,-20.0,0.0,0.0)
    chiFFInterpolationSetProperty(cline,LINE_SECONDPOINT, 20.0,0.0,0.0)
elseif (THREED==true) then
    chiFFInterpolationSetProperty(cline,LINE_FIRSTPOINT ,-20.0,0.0,20.0)
    chiFFInterpolationSetProperty(cline,LINE_SECONDPOINT, 20.0,0.0,20.0)
end

chiFFInterpolationSetProperty(cline,LINE_NUMBEROFPOINTS, 500)

--chiFFInterpolationSetProperty(cline,ADD_FIELDFUNCTION,fflist0[1])
chiFFInterpolationSetProperty(cline,ADD_FIELDFUNCTION,fflist1[1])
chiFFInterpolationSetProperty(cline,ADD_FIELDFUNCTION,fflist1[1]+count1/2)
chiFFInterpolationSetProperty(cline,ADD_FIELDFUNCTION,fflist0[1])
chiFFInterpolationSetProperty(cline,ADD_FIELDFUNCTION,fflist2[1])
chiFFInterpolationSetProperty(cline,ADD_FIELDFUNCTION,fflist2[1]+count2/2)
--chiFFInterpolationSetProperty(cline,ADD_FIELDFUNCTION,fflist1[1])
--chiFFInterpolationSetProperty(cline,ADD_FIELDFUNCTION,fflist2[1])

chiFFInterpolationInitialize(cline)
chiFFInterpolationExecute(cline)
chiFFInterpolationExportPython(cline)


--
chiExportFieldFunctionToVTKG(fflist0[1],"ZPhi","Phi")
chiExportFieldFunctionToVTKG(fflist1[1]+count1/2,"ZErr","Err")

if (chi_location_id == 0) then
    local handle = io.popen("python3 Compare3.py")
end
