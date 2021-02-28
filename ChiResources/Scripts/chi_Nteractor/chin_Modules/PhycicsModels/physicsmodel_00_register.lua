moduleFolderPath = chinGetScriptPath();
physicsRegister = {}
physicsRegisterCount = 0;
dofile(moduleFolderPath.."/ReactorDynamics1D/rdynamic1d.lua")
dofile(moduleFolderPath.."/SPH/sph.lua")