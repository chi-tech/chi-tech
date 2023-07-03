--############################################### Setup mesh
chiMeshHandlerCreate()

mesh={}
N=40
L=14
xmin = 0.0
dx = L/N
for i=1,(N+1) do
    k=i-1
    mesh[i] = xmin + k*dx
end

chiMeshCreateUnpartitioned2DOrthoMesh(mesh,mesh)

chiVolumeMesherExecute();

chiVolumeMesherSetMatIDToAll(0)

vol1 = chi_mesh.RPPLogicalVolume.Create
({ xmin=-1000.0,xmax=10.0,ymin=-1000.0,ymax=10.0, infz=true })
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol1,1)

