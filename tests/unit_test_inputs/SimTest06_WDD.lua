--############################################### Setup mesh
chiMeshHandlerCreate()

if (nmesh==nil) then nmesh = 10 end

mesh={}
N=nmesh
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

chiSimTest06_WDD();
chiMPIBarrier()
if (chi_location_id == 0) then
    os.execute("rm SimTest_06*")
end