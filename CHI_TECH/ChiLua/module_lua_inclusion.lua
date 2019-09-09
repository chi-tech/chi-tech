moduleFolders={}
moduleFolders.itemCount = 0;

function Add_Folder(folderPath)
    moduleFolders.itemCount=moduleFolders.itemCount+1;
    moduleFolders[moduleFolders.itemCount]=folderPath;
end

Add_Folder("../CHI_MODULES/CHI_THERMOALPHA/lua")

Add_Folder("../ChiMath/Quadratures/lua")
Add_Folder("../ChiMath/Quadratures/LegendrePoly/lua")

Add_Folder("../ChiMesh/lua")
Add_Folder("../ChiMesh/SurfaceMesh/lua")
Add_Folder("../ChiMesh/LogicalVolume/lua")
Add_Folder("../ChiMesh/Region/lua")
Add_Folder("../ChiMesh/MeshHandler/lua")
Add_Folder("../ChiMesh/LineMesh/lua")
Add_Folder("../ChiMesh/SurfaceMesher/lua")
Add_Folder("../ChiMesh/VolumeMesher/lua")
Add_Folder("../ChiMesh/FieldFunctionInterpolation/lua")
Add_Folder("../ChiMesh/DomainDecomposition/lua")

Add_Folder("../ChiMPI/lua")
Add_Folder("../ChiLog/lua")

Add_Folder("../ChiPhysics/lua")

Add_Folder("../CHI_MODULES/CHI_MONTECARLON/lua")
Add_Folder("../CHI_MODULES/CHI_DIFFUSION/lua")
Add_Folder("../CHI_MODULES/LinearBoltzmanSolver/lua")