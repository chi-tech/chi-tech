-------------------------------------------------------------------------------- Shared Global Variables
--========================================================= Begin GlobVar
chinGlobal = {};


--================================= Type ____ GlobalVar
chinGlobal.xSize            = 0;
chinGlobal.ySize            = 0;
chinGlobal.dwindowxsize     = 400;
chinGlobal.dwindowysize     = 400;
chinGlobal.dwindowxposi     = 10;
chinGlobal.dwindowyposi     = 10;

CHI_DIFFUSE_COLOR    = 2;


--================================= Type ____ GlobalVar
chinGlobal.wTouch           = false;
chinGlobal.wResize          = false;


--################################## Global functions
dofile(chinBaseDir.."chin_GlobalFunctions/chin_globalfs_00_register.lua")

--################################## Modules
dofile(chinBaseDir.."chin_Modules/Table/table.lua")
dofile(chinBaseDir.."chin_Modules/Material/chinMaterial.lua")
dofile(chinBaseDir.."chin_Modules/Object/chinObject.lua")
dofile(chinBaseDir.."chin_Modules/Texture/chinTexture.lua")
dofile(chinBaseDir.."chin_Modules/PhycicsModels/physicsmodel_00_register.lua")

