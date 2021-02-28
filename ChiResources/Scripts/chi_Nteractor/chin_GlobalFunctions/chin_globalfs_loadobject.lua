function chinLoadObject(fileName)
    --=========================================== Create object
    numObj = chiObjectQuery(0);
    newSurf = chiLoadSurface(fileName);
    objName = "Object_"..string.format("%4d",numObj)
    chiObjectCreate(objName); 
    chiObjectAddSurface(objName,newSurf);
    
    --=========================================== Initialize lua context
    chiObjectQuery(2,numObj);
    local this=chinObject[numObj];
    this.name   =objName;
    this.obj1Num=numObj;
    this.selectable=true;
    this.selected=true;

    this.matlNum = chiMaterialCreate(objName .. "_Material");
    ambient = 1.0;
    chiMaterialSetProperty(this.matlNum,"Diffuse",ambient,ambient,ambient,1.0);
    chiObjectSetProperty(objName,"Material",objName .. "_Material");
    
    this.obj1Tra =chiTransformCreate(objName.."_Transform");
    this.obj1TTra=chiTransformCreate(objName.."_TextureTransform");
    
    chiObjectSetProperty(objName,"Transform",objName.."_Transform");
    chiObjectSetProperty(objName,"TextureTransform",objName.."_TextureTransform");
    chiObjectSetProperty(objName,"ListAble",true);
  
    --=========================================== Events function
    function this.ProcessEvents(this)
        if (WM_LBUTTONDOWN.occured) then
            
            if     ((WM_LBUTTONDOWN.iPar5==this.obj1Num) and (this.selectable)) then
                selectionStack.Clear()
                local newSel = SelectionClass.New();
                newSel.type = SELECTION_TYPE_OBJECT;
                newSel.index = this.obj1Num;
                newSel.originFeature=this;
                selectionStack.PushItem(newSel);
            end
        end        
    end
    
    --=========================================== Selection function
    function this.SelectionStackCheck(selectionItem)
        
        if (selectionItem.type==SELECTION_TYPE_OBJECT) then
            if (selectionItem.index==this.obj1Num) then
                chiObjectSetProperty(selectionItem.index,"ShowBoundingBox",true);
            end
        end
    end
    
    --=========================================== DeSelection function
    function this.DeSelectionStackCheck(selectionItem)
        
        if (selectionItem.type==SELECTION_TYPE_OBJECT) then
            if (selectionItem.index==this.obj1Num) then
                chiObjectSetProperty(selectionItem.index,"ShowBoundingBox",false);
            end
        end
    end
    
    selectionStack.PushSelectionCallback(this.SelectionStackCheck);
    selectionStack.PushDeSelectionCallback(this.DeSelectionStackCheck);
    
    --=========================================== Register feature on events
    dRegisterFeature(this);
    
    --=========================================== Add it to the stack
    simulationObjectCount=simulationObjectCount+1;
    simulationObjectStack[simulationObjectCount]=this;
end