--Called from ../chiNteractor_04_treeview.lua
chinTimeControl={};
chinTimeControl.timeStep = 5.0;
chinTimeControl.maxSteps = 1000;
chinTimeControl.endTime  = 18.0;
chinTimeControl.endType  = 1; --1=EndTime,2=Steps
chinTimeControl.initialized  = false;
chinTimeControl.stepFunction = {}
chinTimeControl.stepFunctionCount = 0;
chinTimeControl.initializeFunction = {}
chinTimeControl.initializeFunctionCount = 0;
chinTimeControl.clearFunction = {}
chinTimeControl.clearFunctionCount = 0;
chinTimeControl.slaves = {}
chinTimeControl.slaveCount = 0;
chinTimeControl.stepsCompleted = 0;
chinTimeControl.timeCompleted = 0;
chinTimeControl.forceStop=false;
chinTimeControl.running = false;

chinTimeControl.time = {}
chinTimeControl.power = {}
chinTimeControl.rho = {}
chinTimeControl.powerCount = 0;

function chinTimeControl.Initialize(this)
    print("Initializing")
    chinTimeControl.stepsCompleted = 0;
    chinTimeControl.timeCompleted = 0;
    for k=1,this.slaveCount do
        if (this.slaves[k]~=nil) then
            this.slaves[k].Initialize(this.slaves[k]);
        end
    end
    
    print("Initialization complete")
end
showTime=0.0;
oldTime=0.0;
function chinTimeControl.Step(this)
    this.stepsCompleted = this.stepsCompleted+1;
    this.timeCompleted = this.timeCompleted+this.timeStep/1000.0;
    if (this.timeCompleted>oldTime+showTime) then
        print(string.format("Stepping iteration %05d, time=%10.4f",this.stepsCompleted,this.timeCompleted))
        oldTime=oldTime+showTime;
    end
    for k=1,this.slaveCount do
        if (this.slaves[k]~=nil) then
            print("Step slave ",k)
            this.slaves[k].Step(this.slaves[k],this.timeStep/1000.0);
        end
    end
end

function chinTimeControl.Run(this)
    if (this.endType==1) then
        if (this.timeCompleted<=chinTimeControl.endTime) then
            this.Step(this)
            
        end
    else
        if (this.stepsCompleted<=chinTimeControl.maxSteps) then
            this.Step(this)
            
        end
    end
end

function chinTimeControl.Dump(this)
    file = io.open("ZZZOutPutFile.txt","w");
    io.output(file);
    for k=1,this.powerCount do
        io.write(string.format("%.4e %.4e %.4e\n",this.time[k],this.power[k],this.rho[k]))
    end
    io.close(file);
    print("Dump created")
end

function chinTimeControl.AddSlave(this,feature)
    this.slaveCount=this.slaveCount+1;
    this.slaves[this.slaveCount]=feature;
end

--######################################################### UpdateTimeControl
function UpdateTimeControl()
    if (chinTimeControl.running) then
        chinTimeControl.Run(chinTimeControl);
    end
    if (not chinTimeControl.initialized) then
        chinTimeControl.initialized=true;
        
        if (TimeFolderSubFolders == nil) then
            --======================= Create selection call back for properties
            function func(this)
                GridviewHideAllItems();
                local newSel=SelectionClass.New();
                newSel.type=SELECTION_TYPE_PROPERTY;
                newSel.index=0;
                newSel.originFeature=this;
                selectionStack.PushItem(newSel);
            end
            
            chinTimeControl.prop = {}
            chinTimeControl.prop[1] = PropertyClass.New("Timestep","Scalar");
            chinTimeControl.prop[2] = PropertyClass.New("Maximum Steps","Scalar");
            chinTimeControl.prop[3] = PropertyClass.New("Maximum Time","Scalar");
            chinTimeControl.prop[4] = PropertyClass.New("Termination type","Scalar");
            
            chinTimeControl.prop[1].value = chinTimeControl.timeStep;
            chinTimeControl.prop[2].value = chinTimeControl.maxSteps;
            chinTimeControl.prop[3].value = chinTimeControl.endTime ;
            chinTimeControl.prop[4].value = chinTimeControl.endType ;
            
            chinTimeControl.prop[4].displayFormat = "%d";
            
            TimeFolderSubFolders = {};
            TimeFolderSubFolders[1] = TimeFolder.AddFolder(TimeFolder,"Timestep");
            TimeFolderSubFolders[2] = TimeFolder.AddFolder(TimeFolder,"Maximum Steps");
            TimeFolderSubFolders[3] = TimeFolder.AddFolder(TimeFolder,"Maximum Time");
            TimeFolderSubFolders[4] = TimeFolder.AddFolder(TimeFolder,"Termination type");
            
            TimeFolderSubFolders[1].iconTypeFolder = chinIconRedBall;
            TimeFolderSubFolders[2].iconTypeFolder = chinIconRedBall;
            TimeFolderSubFolders[3].iconTypeFolder = chinIconRedBall;
            TimeFolderSubFolders[4].iconTypeFolder = chinIconRedBall;
            
            TimeFolderSubFolders[1].iconTypeExpander = chinIconNothing;
            TimeFolderSubFolders[2].iconTypeExpander = chinIconNothing;
            TimeFolderSubFolders[3].iconTypeExpander = chinIconNothing;
            TimeFolderSubFolders[4].iconTypeExpander = chinIconNothing;
            
            TimeFolderSubFolders[1].label.parent = chinTimeControl.prop[1];
            TimeFolderSubFolders[2].label.parent = chinTimeControl.prop[2];
            TimeFolderSubFolders[3].label.parent = chinTimeControl.prop[3];
            TimeFolderSubFolders[4].label.parent = chinTimeControl.prop[4];
            
            TimeFolderSubFolders[1].label.CustomSelected = func;
            TimeFolderSubFolders[2].label.CustomSelected = func;
            TimeFolderSubFolders[3].label.CustomSelected = func;
            TimeFolderSubFolders[4].label.CustomSelected = func;
            
            titem=chinTimeControl.prop[1];
            function titem.CustomValueChanged(this)
                outstring = "chinTimeControl.timeStep = "..string.format("%.6f",this.value)..";\n"
                chinTimeControl.timeStep = this.value;
            end
            titem=chinTimeControl.prop[2];
            function titem.CustomValueChanged(this)
                outstring = "chinTimeControl.maxSteps = "..string.format("%d"  ,this.value)..";\n"
                chinTimeControl.maxSteps = this.value;
            end
            titem=chinTimeControl.prop[3];
            function titem.CustomValueChanged(this)
                outstring = "chinTimeControl.endTime  = "..string.format("%d"  ,this.value )..";\n"
                chinTimeControl.endTime = this.value;
            end
            titem=chinTimeControl.prop[4];
            function titem.CustomValueChanged(this)
                outstring = "chinTimeControl.endType  = "..string.format("%d"  ,this.value )..";\n"
                chinTimeControl.endType = this.value;
            end

    
                
            
        end
    end
    
    
end