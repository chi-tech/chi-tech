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

--chiMeshCreateUnpartitioned3DOrthoMesh(mesh,mesh,mesh)
chiMeshCreateUnpartitioned2DOrthoMesh(mesh,mesh)
--chiMeshCreateUnpartitioned1DOrthoMesh(mesh)
chiVolumeMesherExecute();

--############################################### Set Material IDs
chiVolumeMesherSetMatIDToAll(0)

chi_unit_sim_tests.chiSimTest01_FV();
chiMPIBarrier()
if (chi_location_id == 0) then
    os.execute("rm CodeTut1_FV*")
end