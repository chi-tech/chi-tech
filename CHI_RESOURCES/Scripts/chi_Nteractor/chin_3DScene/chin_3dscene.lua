scriptPath=chinGetScriptPath();
dofile(scriptPath.."/chin_3dscene_00_setup.lua")
dofile(scriptPath.."/chin_3dscene_01_manipulator.lua");
dofile(scriptPath.."/chin_3dscene_02_camera.lua");


simulationObjectStack={};
simulationObjectCount=0;

infoLabel=chiSetLabel("Position","",5,15,0.0,0.0,0.0,1.0,3);
function main()

    CameraControl();
    chiSetLabel(infoLabel,string.format("x=%6.2f y=%6.2f z=%6.2f yaw=%6.1f pitch=%6.1f FPS=%.1f FTC=%.2f",Camera_x,Camera_y,Camera_z,Camera_alpha,Camera_azimuth,chi_frameRate,chi_physicsTimeCost),5,15,0.0,0.0,0.0,1.0,3);
end