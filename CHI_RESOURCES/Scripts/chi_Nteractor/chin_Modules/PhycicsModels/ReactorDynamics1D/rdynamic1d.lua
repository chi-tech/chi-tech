PhysicsName="Reactor Point Kinetics"
ReactorDynamics1DClass = {}
ReactorDynamics1DClass.__index = ReactorDynamics1DClass
ReactorDynamics1DCount = 0;

submoduleFolderPath = chinGetScriptPath();

dofile(submoduleFolderPath.."/rdynamic1d_00_constrdestr.lua")


physicsRegisterCount = physicsRegisterCount+1;
physicsRegister[physicsRegisterCount] = PhysicsName;