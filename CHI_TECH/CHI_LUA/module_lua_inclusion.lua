moduleFolders={}
moduleFolders.itemCount = 0;

function Add_Folder(folderPath)
    moduleFolders.itemCount=moduleFolders.itemCount+1;
    moduleFolders[moduleFolders.itemCount]=folderPath;
end

Add_Folder("../CHI_MODULES/CHI_PI3/lua")
Add_Folder("../CHI_MODULES/CHI_THERMOALPHA/lua")

Add_Folder("../CHI_MATH/Quadratures/lua")
Add_Folder("../CHI_MATH/Quadratures/LegendrePoly/lua")

Add_Folder("../CHI_MESH/lua")
Add_Folder("../CHI_MESH/CHI_SURFACEMESH/lua")
Add_Folder("../CHI_MESH/CHI_LOGICALVOLUME/lua")
Add_Folder("../CHI_MESH/CHI_REGION/lua")
Add_Folder("../CHI_MESH/CHI_MESHHANDLER/lua")
Add_Folder("../CHI_MESH/CHI_LINEMESH/lua")
Add_Folder("../CHI_MESH/CHI_SURFACEMESHER/lua")
Add_Folder("../CHI_MESH/CHI_VOLUMEMESHER/lua")
Add_Folder("../CHI_MESH/CHI_FFINTERPOLATION/lua")

Add_Folder("../CHI_MPI/lua")
Add_Folder("../CHI_LOG/lua")

Add_Folder("../CHI_PHYSICS/lua")

Add_Folder("../CHI_MODULES/CHI_MONTECARLON/lua")
Add_Folder("../CHI_MODULES/CHI_DIFFUSION/lua")
Add_Folder("../CHI_MODULES/LinearBoltzmanSolver/lua")