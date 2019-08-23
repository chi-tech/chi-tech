newSurf = chiLoadSurface(scriptPath.."/Assets/Meshes/Manipulator.obj");

dofile(scriptPath.."/Manipulator/manipulator_00_constrdestr.lua");
dofile(scriptPath.."/Manipulator/manipulator_01_events.lua");
dofile(scriptPath.."/Manipulator/manipulator_02_hideunhide.lua");
dofile(scriptPath.."/Manipulator/manipulator_03_slavecontrol.lua");

manipulator={}
manipulator[1]=ManipulatorClass.New("MasterManip");
manipulator[1].Hide(manipulator[1]);

dRegisterFeature(manipulator[1])

function ManipulatorSelectionCallback(selectionItem)
    if (selectionItem.type==SELECTION_TYPE_OBJECT) then
        manipulator[1].AddObjectSlave(manipulator[1],selectionItem.originFeature);
    end
end

function ManipulatorDeSelectionCallback(selectionItem)
    manipulator[1].RemoveAllObjectSlaves(manipulator[1])
end

selectionStack.PushSelectionCallback(ManipulatorSelectionCallback)

selectionStack.PushDeSelectionCallback(ManipulatorDeSelectionCallback)


