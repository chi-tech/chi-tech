--############################################### Setup mesh
chiMeshHandlerCreate()

nodes={}
N=32
ds=2.0/N
for i=0,N do
    nodes[i+1] = i*ds
end
chiMeshCreate2DOrthoMesh(nodes,nodes)
chiVolumeMesherExecute();


material = chiPhysicsAddMaterial("Test Material");

-- Set Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,material)

--############################################### Add material properties


-- Set material properties
chiPhysicsMaterialAddProperty(material,SCALAR_VALUE,"k")
chiPhysicsMaterialSetProperty(material,"k",SINGLE_VALUE,1.0)

chiPhysicsMaterialAddProperty(material,SCALAR_VALUE,"q")
chiPhysicsMaterialSetProperty(material,"q",SINGLE_VALUE,1.0)


--############################################### Setup Physics
phys1 = chiDiffusionCreateSolver();
chiSolverAddRegion(phys1,region1)
chiDiffusionSetProperty(phys1,DISCRETIZATION_METHOD,PWLC);
chiDiffusionSetProperty(phys1,RESIDUAL_TOL,1.0e-6)

--############################################### Initialize and
--                                                Execute Solver
chiDiffusionInitialize(phys1)
chiDiffusionExecute(phys1)

--############################################### Visualize the field function
fflist,count = chiGetFieldFunctionList(phys1)
chiExportFieldFunctionToVTK(fflist[1],"Tutorial1Output","Temperature")