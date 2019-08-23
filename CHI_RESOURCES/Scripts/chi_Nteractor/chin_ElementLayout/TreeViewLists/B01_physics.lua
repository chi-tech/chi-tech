--Called from ../chiNteractor_04_treeview.lua

numrdynamics1D    = 0;
prevnumrdynamics1D= 0;
chinRDynamic1D    = {};
--######################################################### Update modules
function UpdateReactorDynamics1D()
    
    numrdynamics1D = chiReactorDynamics1DQuery(0);
    if (numrdynamics1D>0) then
        if (ReactorDynamics1DFolder==nil) then
            ReactorDynamics1DFolder=PhysicsFolder.AddFolder(PhysicsFolder,"Reactor Dynamics 1D");
            ReactorDynamics1DFolder.iconTypeFolder=chinIconFolderBlue;
            ReactorDynamics1DFolderSubFolder={};
        end
    end
    
    if (not (numrdynamics1D == prevnumrdynamics1D)) then
        for m=1,numrdynamics1D do
            k=m-1
            if (chinRDynamic1D[k]==nil) then
                chiReactorDynamics1DQuery(2,k);
            end
        
            ReactorDynamics1DFolderSubFolder[m]=ReactorDynamics1DFolder.AddFolder(ReactorDynamics1DFolder,chinRDynamic1D[k].name);
            ReactorDynamics1DFolderSubFolder[m].iconTypeFolder=chinIconAtom;
            
            chiReactorDynamics1DQuery(4,k);
            
            --======================= Create selection call back for properties
            function func(this)
                GridviewHideAllItems();
                local newSel=SelectionClass.New();
                newSel.type=SELECTION_TYPE_PROPERTY;
                newSel.index=0;
                newSel.originFeature=this;
                selectionStack.PushItem(newSel);
            end
            
            --======================= Populating properties
            item = ReactorDynamics1DFolderSubFolder[m];
            
            --Initial conditions
            initFolder = item.AddFolder(item,"Initial Conditions");
            initFolder.iconTypeFolder=chinIconFolderBlue;
        
            chinRDynamic1D[k].initModeProp = PropertyClass.New("Initialization Mode","Droplist");
            chinRDynamic1D[k].initModeProp.index = k;
            chinRDynamic1D[k].initModeProp.displayFormat = "%d";
            chinRDynamic1D[k].initModeProp.initOption = chinRDynamic1D[k].initMode
            tempItem = chinRDynamic1D[k].initModeProp;
            tempItem.valueCount = tempItem.valueCount+4;
            tempItem.value[1] = "1 Sub-critical (need rho)";
            tempItem.value[2] = "2 Critical (need Power)";
            tempItem.value[3] = "3 SuperCritical (need both)";
            tempItem.value[4] = "4 User supplied";

            chinRDynamic1D[k].initValueProp = PropertyClass.New("Initialization Value","Vec2");
            chinRDynamic1D[k].initValueProp.index = k;
            chinRDynamic1D[k].initValueProp.displayFormat = "%.2f,%.2e";
            chinRDynamic1D[k].initValueProp.value = chinRDynamic1D[k].initValues
         
            newFolder = initFolder.AddFolder(initFolder,"Initialization mode");
            newFolder.iconTypeFolder=chinIconRedBall;
            newFolder.iconTypeExpander=chinIconNothing;
            newFolder.label.parent = chinRDynamic1D[k].initModeProp;
            newFolder.label.CustomSelected=func;
            
            newFolder = initFolder.AddFolder(initFolder,"Initialization Values");
            newFolder.iconTypeFolder=chinIconRedBall;
            newFolder.iconTypeExpander=chinIconNothing;
            newFolder.label.parent = chinRDynamic1D[k].initValueProp;
            newFolder.label.CustomSelected=func;
            
            --Callbacks
            tempitem = chinRDynamic1D[k].initModeProp;
            function tempitem.CustomValueChanged(this)
                outputS = "chiReactorDynamics1DSetProperty("..string.format("%d,%d,%d);",this.index,5,this.option)
                print(outputS)
                chunk,err=assert(load(outputS)); 
                if (chunk==nil) then print(err);
                else chunk() end
            end
            
            tempitem = chinRDynamic1D[k].initValueProp;
            function tempitem.CustomValueChanged(this)
                outputS = "chiReactorDynamics1DSetProperty("..string.format("%d,%d,%f,%f);",this.index,10,this.value[1],this.value[2])
                print(outputS)
                chunk,err=assert(load(outputS)); 
                if (chunk==nil) then print(err);
                else chunk() end
            end
            
            --Precursors
            chinRDynamic1D[k].precursorProperty={}
            for i=1,6 do
                chinRDynamic1D[k].precursorProperty[i] = PropertyClass.New("Precursor "..string.format("%d",i).."[b,l]","Vec2");
                chinRDynamic1D[k].precursorProperty[i].index = k;
                chinRDynamic1D[k].precursorProperty[i].displayFormat = "%.5f,%.4e";
                chinRDynamic1D[k].precursorProperty[i].value = {chinRDynamic1D[k].Precursor[i].betai,chinRDynamic1D[k].Precursor[i].lambdai};
                titem = chinRDynamic1D[k].precursorProperty[i];
                titem.referenceIndex2 = i;
                function titem.CustomValueChanged(this)
                    outputS="chiReactorDynamics1DSetProperty";
                    outputS=outputS.."("..string.format("%d",this.index)
                    outputS=outputS..",2,"..string.format("%d",this.referenceIndex2)
                    outputS=outputS..","..string.format("%f",this.value[1])..");"
                    
                    print(outputS)
                    chunk,err=assert(load(outputS)); 
                    if (chunk==nil) then print(err);
                    else chunk() end
                    
                    outputS="chiReactorDynamics1DSetProperty";
                    outputS=outputS.."("..string.format("%d",this.index)
                    outputS=outputS..",3,"..string.format("%d",this.referenceIndex2)
                    outputS=outputS..","..string.format("%f",this.value[2])..");"
                    
                    print(outputS)
                    chunk,err=assert(load(outputS)); 
                    if (chunk==nil) then print(err);
                    else chunk() end
                end
                
                newFolder = item.AddFolder(item,"Precursor "..string.format("%d",i));
                newFolder.iconTypeFolder=chinIconRedBall;
                newFolder.iconTypeExpander=chinIconNothing;
                newFolder.label.parent = chinRDynamic1D[k].precursorProperty[i];
                newFolder.label.CustomSelected=func;
            end
            
            
        end
        mainTree.Redraw(mainTree)
        numrdynamics1D    = 0;
        numrdynamics1D = chiReactorDynamics1DQuery(0);
        prevnumrdynamics1D = numrdynamics1D;
    end
    
end