function chinSaveScene(fileName)
    --Camera position
    folderPath = chinFilePathToFolder(fileName);
    print("Saving Scene")
    
    
    --===================================================== Write chunk to disc
    file = io.open(fileName,"w");
    io.output(file);
    
    --===================================================== Camera position
    io.write("\n--=========================== Camera\n")
    io.write("chiBindScene(1);\n")
    io.write("Camera_x          ="..string.format("%f",Camera_x      )..";\n");
    io.write("Camera_y          ="..string.format("%f",Camera_y      )..";\n");
    io.write("Camera_z          ="..string.format("%f",Camera_z      )..";\n");
    
    io.write("Camera_alpha      ="..string.format("%f",Camera_alpha  )..";\n");
    io.write("Camera_azimuth    ="..string.format("%f",Camera_azimuth)..";\n");
    io.write("Camera_tsi        ="..string.format("%f",Camera_tsi    )..";\n");
    
    io.write("chiGraphicsOrientCamera(0,Camera_alpha,Camera_azimuth,Camera_tsi);\n");
    io.write("chiGraphicsPositionCamera(Camera_x,Camera_y,Camera_z);\n"); 
    io.write("print(\"Camera loaded\");\n")
    io.write("chiBindScene(0);\n")
    --===================================================== Save timestep control
    io.write("\n--=========================== Time control\n")
    io.write("chinTimeControl.timeStep = "..string.format("%.6f",chinTimeControl.timeStep)..";\n");
    io.write("chinTimeControl.maxSteps = "..string.format("%d"  ,chinTimeControl.maxSteps)..";\n");
    io.write("chinTimeControl.endTime  = "..string.format("%.6f" ,chinTimeControl.endTime )..";\n");
    io.write("chinTimeControl.endType  = "..string.format("%d"  ,chinTimeControl.endType )..";\n");
    io.write("print(\"Time control loaded\");\n")
    
    --===================================================== Saving pkinetics objects
    SavePKinetics();
    
    --===================================================== Saving Table objects
    SaveTables();
    
    
    --===================================================== Save textures
    print("Saving Textures")
    SaveTextures();
    
    print("Scene Saved")
    
    io.close(file);
end

