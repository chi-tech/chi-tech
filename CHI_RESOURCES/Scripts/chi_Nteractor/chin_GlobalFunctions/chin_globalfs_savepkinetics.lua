function SavePKinetics()
    io.write("\n--=========================== PKinetics\n")
    numrdynamics1D = chiReactorDynamics1DQuery(0);
    for pkCount=1,numrdynamics1D do
        io.write("temp = ReactorDynamics1DClass.New(\"Temp Name\");\n");
        k=pkCount-1
        chiReactorDynamics1DQuery(4,k);
        if (chinRDynamic1D[k]~=nil) then
            io.write("chiReactorDynamics1DSetProperty(temp.moduleIndex,1,\""..chinRDynamic1D[k].name.."\");                  \n");
            io.write("\n");
            io.write("chiReactorDynamics1DSetProperty(temp.moduleIndex,2,0,"..string.format("%f",chinRDynamic1D[k].Precursor[1].betai)..");  \n");
            io.write("chiReactorDynamics1DSetProperty(temp.moduleIndex,2,1,"..string.format("%f",chinRDynamic1D[k].Precursor[2].betai)..");  \n");
            io.write("chiReactorDynamics1DSetProperty(temp.moduleIndex,2,2,"..string.format("%f",chinRDynamic1D[k].Precursor[3].betai)..");  \n");
            io.write("chiReactorDynamics1DSetProperty(temp.moduleIndex,2,3,"..string.format("%f",chinRDynamic1D[k].Precursor[4].betai)..");  \n");
            io.write("chiReactorDynamics1DSetProperty(temp.moduleIndex,2,4,"..string.format("%f",chinRDynamic1D[k].Precursor[5].betai)..");  \n");
            io.write("chiReactorDynamics1DSetProperty(temp.moduleIndex,2,5,"..string.format("%f",chinRDynamic1D[k].Precursor[6].betai)..");  \n");
            io.write("\n");
            io.write("chiReactorDynamics1DSetProperty(temp.moduleIndex,3,0,"..string.format("%f",chinRDynamic1D[k].Precursor[1].lambdai)..");  \n");
            io.write("chiReactorDynamics1DSetProperty(temp.moduleIndex,3,1,"..string.format("%f",chinRDynamic1D[k].Precursor[2].lambdai)..");  \n");
            io.write("chiReactorDynamics1DSetProperty(temp.moduleIndex,3,2,"..string.format("%f",chinRDynamic1D[k].Precursor[3].lambdai)..");  \n");
            io.write("chiReactorDynamics1DSetProperty(temp.moduleIndex,3,3,"..string.format("%f",chinRDynamic1D[k].Precursor[4].lambdai)..");  \n");
            io.write("chiReactorDynamics1DSetProperty(temp.moduleIndex,3,4,"..string.format("%f",chinRDynamic1D[k].Precursor[5].lambdai)..");  \n");
            io.write("chiReactorDynamics1DSetProperty(temp.moduleIndex,3,5,"..string.format("%f",chinRDynamic1D[k].Precursor[6].lambdai)..");  \n");
            io.write("\n");
            io.write("chiReactorDynamics1DSetProperty(temp.moduleIndex,4,"..string.format("%f",chinRDynamic1D[k].source)..");\n");
            io.write("\n");        
            io.write("chiReactorDynamics1DSetProperty(temp.moduleIndex,5,"..string.format("%d",chinRDynamic1D[k].initMode)..");\n");
            io.write("\n");       
            io.write("chiReactorDynamics1DSetProperty(temp.moduleIndex,6,"..string.format("%d",chinRDynamic1D[k].solver)..");\n");
            io.write("\n");
            io.write("chiReactorDynamics1DSetProperty(temp.moduleIndex,7,0,"..string.format("%f",chinRDynamic1D[k].Precursor[1].ci)..");  \n");
            io.write("chiReactorDynamics1DSetProperty(temp.moduleIndex,7,1,"..string.format("%f",chinRDynamic1D[k].Precursor[2].ci)..");  \n");
            io.write("chiReactorDynamics1DSetProperty(temp.moduleIndex,7,2,"..string.format("%f",chinRDynamic1D[k].Precursor[3].ci)..");  \n");
            io.write("chiReactorDynamics1DSetProperty(temp.moduleIndex,7,3,"..string.format("%f",chinRDynamic1D[k].Precursor[4].ci)..");  \n");
            io.write("chiReactorDynamics1DSetProperty(temp.moduleIndex,7,4,"..string.format("%f",chinRDynamic1D[k].Precursor[5].ci)..");  \n");
            io.write("chiReactorDynamics1DSetProperty(temp.moduleIndex,7,5,"..string.format("%f",chinRDynamic1D[k].Precursor[6].ci)..");  \n");
            io.write("\n");
            io.write("chiReactorDynamics1DSetProperty(temp.moduleIndex,8,"..string.format("%f",chinRDynamic1D[k].initPower)..");\n");
            io.write("\n");
            io.write("chiReactorDynamics1DSetProperty(temp.moduleIndex,9,"..string.format("%f",chinRDynamic1D[k].initRho)..");\n");
            
            
            
        end
    end
    io.write("print(\"Point kinetics loaded\");\n")
    
    
end