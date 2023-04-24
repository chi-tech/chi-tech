--############################################### Setup mesh
chiMeshHandlerCreate()

mesh={}
N=25
L=2
xmin = -L/2
--xmin = 0.0
dx = L/N
for i=1,(N+1) do
    k=i-1
    mesh[i] = xmin + k*dx
end

zmesh={}
Nz=8
Lz=1.6
zmin=0.0
dz=Lz/Nz
for i=1,(Nz+1) do
    k=i-1
    zmesh[i] = zmin + k*dz
end

--chiMeshCreateUnpartitioned3DOrthoMesh(mesh,mesh,zmesh)
chiMeshCreateUnpartitioned2DOrthoMesh(mesh,mesh)
chiVolumeMesherExecute();

--############################################### Set Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)

chiSimTest91_PWLD();
chiMPIBarrier()
if (chi_location_id == 0) then
    os.execute("rm SimTest_91*")
end

--[0]  Iteration     0   1.000e+00
--[0]  Iteration     1   2.016e+02
--[0]  Iteration     2   1.941e+00
--[0]  Iteration     3   1.294e+00
--[0]  Iteration     4   3.890e-01
--[0]  Iteration     5   2.887e-02
--[0]  Iteration     6   1.239e-03
--[0]  Iteration     7   4.076e-05
--[0]  Iteration     8   1.119e-06
--[0]  Iteration     9   2.955e-08
