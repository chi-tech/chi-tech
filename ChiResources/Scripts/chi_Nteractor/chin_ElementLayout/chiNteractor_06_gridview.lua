folderPath = chinGetScriptPath()
dofile(folderPath.."/GridView/materials.lua")
dofile(folderPath.."/GridView/properties.lua")
dofile(folderPath.."/GridView/objects.lua")
dofile(folderPath.."/GridView/shaders.lua")
dofile(folderPath.."/GridView/tables.lua")
dofile(folderPath.."/GridView/textures.lua")

--=============================================== Gridview
propertyGrid=GridViewClass.New("PropertyGrid");
propertyGrid.SetProperty(propertyGrid,"Master",panels[5]);
gridLtLabel={};
gridRtFeature={};

--######################################################### Gridview processing
--This function gets called everytime something gets pushed to
--the selection stack. It acts as an interface between the 
--gridview and the selected item.
function GridviewGetProperties(selectionItem)
    
    if (selectionItem.type==SELECTION_TYPE_MATERIAL) then
       GridviewGetPropertiesForMaterials(selectionItem); 
    end
    
    if (selectionItem.type==SELECTION_TYPE_PROPERTY) then
       GridviewGetPropertiesForProperties(selectionItem)
    end
    
    if (selectionItem.type==SELECTION_TYPE_OBJECT) then
       GridviewGetPropertiesForObjects(selectionItem)
    end
    if (selectionItem.type==SELECTION_TYPE_SHADER) then
        GridviewGetPropertiesForShaders(selectionItem)
    end
    if (selectionItem.type==SELECTION_TYPE_TABLE) then
        GridviewGetPropertiesForTables(selectionItem)
    end
    if (selectionItem.type==SELECTION_TYPE_TEXTURE) then
        GridviewGetPropertiesForTextures(selectionItem)
    end

end

--######################################################### GridviewHide all items
function GridviewHideAllItems()
    for k=1,100 do
        if (not (gridLtLabel[k]==nil)) then
            gridLtLabel[k].Hide(gridLtLabel[k]);
            gridLtLabel[k]=nil;
        end
        if (not (gridRtFeature[k]==nil)) then
            gridRtFeature[k].Hide(gridRtFeature[k]);
            gridRtFeature[k].CustomKeyPress=nil;
            gridRtFeature[k]=nil;
        end
    end

end

--######################################################### GridviewCreateLtLabel
--Creates a label in the left column and sets all the 
--required properties for it.
function GridviewCreateLtLabel(position)

    gridLtLabel[position]=LabelClass.New("");    
    gridLtLabel[position].SetProperty(gridLtLabel[position],"Master",propertyGrid.gridReference[position][1]);
    gridLtLabel[position].SetProperty(gridLtLabel[position],"Float",false); 
    gridLtLabel[position].SetProperty(gridLtLabel[position],"ViewportEnable",true);
    gridLtLabel[position].highlightFillsMaster=true;
    gridLtLabel[position].paddingLeft=4;
    gridLtLabel[position].paddingBot=3;
    gridLtLabel[position].selectable=false;
     
    dRegisterFeature(gridLtLabel[position]);
    panels[5].SizeChanged(panels[5]);
end

--######################################################### GridviewCreateRtFeauture
--Creates a label in the left column and sets all the 
--required properties for it.
function GridviewCreateRtFeature(position,featureType)

    --Label
    if (featureType==1) then
        gridRtFeature[position]=LabelClass.New(string.format("GridRt %d",position));    
        gridRtFeature[position].SetProperty(gridRtFeature[position],"Master",propertyGrid.gridReference[position][2]);
        gridRtFeature[position].SetProperty(gridRtFeature[position],"Float",false); 
        gridRtFeature[position].SetProperty(gridRtFeature[position],"ViewportEnable",true);
        gridRtFeature[position].highlightFillsMaster=true;
        gridRtFeature[position].paddingLeft=4;
        gridRtFeature[position].paddingBot=3;
        --gridRtFeature[position].selectable=false;
        
        dRegisterFeature(gridRtFeature[position]);
        panels[5].SizeChanged(panels[5]);
    end
    
    --TextBox
    if (featureType==2) then
        gridRtFeature[position]=TextBoxClass.New(string.format("GridRt %d",position));    
        gridRtFeature[position].SetProperty(gridRtFeature[position],"Master",propertyGrid.gridReference[position][2]);
        gridRtFeature[position].SetProperty(gridRtFeature[position],"Float",false); 
        gridRtFeature[position].SetProperty(gridRtFeature[position],"ViewportEnable",true);
        gridRtFeature[position].SetProperty(gridRtFeature[position],"ShowOutline",false);
        gridRtFeature[position].inputBox=true;
        --gridRtFeature[position].highlightFillsMaster=true;
        --gridRtFeature[position].paddingLeft=4;
        --gridRtFeature[position].paddingBot=3;
        --gridRtFeature[position].selectable=false;
        
        dRegisterFeature(gridRtFeature[position]);
        panels[5].SizeChanged(panels[5]);
    end
    
    --Checkbox
    if (featureType==3) then
        gridRtFeature[position]=CheckBoxClass.New(string.format("GridRt %d",position));    
        gridRtFeature[position].SetProperty(gridRtFeature[position],"Master",propertyGrid.gridReference[position][2]);
        gridRtFeature[position].SetProperty(gridRtFeature[position],"Float",false); 
        --gridRtFeature[position].SetProperty(gridRtFeature[position],"ViewportEnable",true);
        --gridRtFeature[position].SetProperty(gridRtFeature[position],"ShowOutline",false);
        gridRtFeature[position].iconSize=19;
        --gridRtFeature[position].highlightFillsMaster=true;
        --gridRtFeature[position].paddingLeft=4;
        --gridRtFeature[position].paddingBot=3;
        --gridRtFeature[position].selectable=false;
        
        dRegisterFeature(gridRtFeature[position]);
        panels[5].SizeChanged(panels[5]);
    end


    --Button
    if (featureType==4) then
        gridRtFeature[position]=ButtonClass.New(string.format("GridRt %d",position));
        gridRtFeature[position].SetProperty(gridRtFeature[position],"Master",propertyGrid.gridReference[position][2]);
        gridRtFeature[position].SetProperty(gridRtFeature[position],"Float",false);
        gridRtFeature[position].highlightFillsMaster=true;
        gridRtFeature[position].paddingLeft=4;
        gridRtFeature[position].paddingBot=3;
        gridRtFeature[position].selectable=true;

        dRegisterFeature(gridRtFeature[position]);
        panels[5].SizeChanged(panels[5]);

    end
    
    --Droplist
    if (featureType==5) then
        gridRtFeature[position]=DropListClass.New(string.format("GridRt %d",position));
        gridRtFeature[position].SetProperty(gridRtFeature[position],"Master",propertyGrid.gridReference[position][2]);
        gridRtFeature[position].SetProperty(gridRtFeature[position],"Float",false);
        gridRtFeature[position].zDepth=2.3
        dRegisterFeature(gridRtFeature[position]);
        panels[5].SizeChanged(panels[5]);
    end
    
end

selectionStack.PushSelectionCallback(GridviewGetProperties);

