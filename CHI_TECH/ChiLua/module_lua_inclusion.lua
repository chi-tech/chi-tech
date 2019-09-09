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
Add_Folder("../ChiMesh/CHI_SURFACEMESH/lua")
Add_Folder("../ChiMesh/CHI_LOGICALVOLUME/lua")
Add_Folder("../ChiMesh/CHI_REGION/lua")
Add_Folder("../ChiMesh/CHI_MESHHANDLER/lua")
Add_Folder("../ChiMesh/CHI_LINEMESH/lua")
Add_Folder("../ChiMesh/CHI_SURFACEMESHER/lua")
Add_Folder("../ChiMesh/CHI_VOLUMEMESHER/lua")
Add_Folder("../ChiMesh/CHI_FFINTERPOLATION/lua")

Add_Folder("../ChiMPI/lua")
Add_Folder("../ChiLog/lua")

Add_Folder("../ChiPhysics/lua")

Add_Folder("../CHI_MODULES/CHI_MONTECARLON/lua")
Add_Folder("../CHI_MODULES/CHI_DIFFUSION/lua")
Add_Folder("../CHI_MODULES/LinearBoltzmanSolver/lua")