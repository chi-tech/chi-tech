--############################################### Setup mesh
chiMeshHandlerCreate()
nodesx={-250.0,-220.2,-190.4,-160.6,-130.8,-101.0,-100.0,-99.0,-80.1,-61.2,-42.300000000000004,-23.400000000000006,-4.5,14.399999999999991,33.29999999999998,52.19999999999999,71.1,90.0,91.0,92.0,93.0,94.0,95.0,96.0,97.0,98.0,99.0,100.0,101.0,102.0,103.0,104.0,105.0,106.0,107.0,108.0,109.0,110.0,138.0,166.0,194.0,222.0,250.0}
nodesy={-250.0,-204.0,-158.0,-112.0,-66.0,-20.0,-19.0,-18.0,-17.0,-16.0,-15.0,-14.0,-13.0,-12.0,-11.0,-10.0,-9.0,-8.0,-7.0,-6.0,-5.0,-4.0,-3.0,-2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0,66.0,112.0,158.0,204.0,250.0}
nodesz={-250.0,-204.0,-158.0,-112.0,-66.0,-20.0,-19.0,-18.0,-17.0,-16.0,-15.0,-14.0,-13.0,-12.0,-11.0,-10.0,-9.0,-8.0,-7.0,-6.0,-5.0,-4.0,-3.0,-2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0,66.0,112.0,158.0,204.0,250.0}

chiMeshCreateUnpartitioned3DOrthoMesh(nodesx,nodesy,nodesz)
chiVolumeMesherSetProperty(PARTITION_TYPE,PARMETIS)
chiVolumeMesherExecute();

--chiRegionExportMeshToVTK(0, "../../Output/MCNP/local_angle/mesh")

--############################################### Setup materials
--############################################### Set Material IDs
vol_air = chiLogicalVolumeCreate(
    RPP,-100000,100000,-100000,100000,-100000,100000
)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol_air,0)

vol_src = chiLogicalVolumeCreate(RPP,-101,-99,-5,5,-5,5)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol_src,1)

vol_target = chiLogicalVolumeCreate(RPP,90,110,-20,20,-20,20)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol_target,2)

--############################################### Add materials
material_air = chiPhysicsAddMaterial("air");
chiPhysicsMaterialAddProperty(material_air,TRANSPORT_XSECTIONS)

material_src = chiPhysicsAddMaterial("src");
chiPhysicsMaterialAddProperty(material_src,TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(material_src,ISOTROPIC_MG_SOURCE)

material_target = chiPhysicsAddMaterial("target");
chiPhysicsMaterialAddProperty(material_target,TRANSPORT_XSECTIONS)

--########################################### Set XS values

num_groups = 1

chiPhysicsMaterialSetProperty(
    material_air,
    TRANSPORT_XSECTIONS,
    SIMPLEXS1,
    num_groups,
    3.26E-5,  -- Sigma_t
    0.99  -- Scattering ratio
)

chiPhysicsMaterialSetProperty(
    material_src,
    TRANSPORT_XSECTIONS,
    SIMPLEXS1,
    num_groups,
    3.26E-5,  -- Sigma_t
    0.99  -- Scattering ratio
)

chiPhysicsMaterialSetProperty(
    material_target,
    TRANSPORT_XSECTIONS,
    SIMPLEXS1,
    num_groups,
    4.94E-2,  -- Sigma_t
    0.992  -- Scattering ratio
)

src={}
for g=1,num_groups do
    src[g] = 1.0 -- source strength
end
chiPhysicsMaterialSetProperty(material_src,ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)

--############################################### Setup physics
phys1 = chiLBSCreateSolver()
chiSolverAddRegion(phys1,0)

--========== ProdQuad
--pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,8,16)

--========== Refined angular quadrature
pquad = chiCreateSLDFESQAngularQuadrature(2)
chiLocallyRefineSLDFESQAngularQuadrature(pquad, {1, 0, 0}, 90.0*math.pi/180, false)
chiLocallyRefineSLDFESQAngularQuadrature(pquad, {1, 0, 0}, 45.0*math.pi/180, false)
chiLocallyRefineSLDFESQAngularQuadrature(pquad, {1, 0, 0}, 30.0*math.pi/180, false)
--chiPrintToPythonSLDFESQAngularQuadrature(pquad, "../../Output/MCNP/local_angle/YQuad_")
--os.exit()

--========== Groups
grp = {}
for g=1,num_groups do
    grp[g] = chiLBSCreateGroup(phys1)
end

--========== Groupset def
gs0 = chiLBSCreateGroupset(phys1)
cur_gs = gs0
chiLBSGroupsetAddGroups(phys1,cur_gs,0,num_groups-1)
chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquad)
chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,1)
chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,NPT_GMRES_CYCLES)
chiLBSGroupsetSetAngleAggregationType(phys1,cur_gs,LBSGroupset.ANGLE_AGG_SINGLE)
chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-12)
chiLBSGroupsetSetMaxIterations(phys1,cur_gs,300)
chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,30)

--========== Boundary conditions
bsrc={}
for g=1,num_groups do
    bsrc[g] = 0.0
end
chiLBSSetProperty(
    phys1,
    BOUNDARY_CONDITION,
    XMIN,
    LBSBoundaryTypes.INCIDENT_ISOTROPIC,
    bsrc
);
chiLBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD3D)
chiLBSSetProperty(phys1,SCATTERING_ORDER,0)

chiLBSInitialize(phys1)
chiLBSExecute(phys1)

--############################################### Grab FF
fflist,count = chiLBSGetScalarFieldFunctionList(phys1)

--chiExportFieldFunctionToASCII(fflist[1], "../../Output/MCNP/local_angle/phi")

chiExportMultiFieldFunctionToVTK({fflist[1]}, "../../Output/MCNP/local_angle/phi")
