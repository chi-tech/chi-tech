function GridviewGetPropertiesForTables(selectionItem)
    local k=selectionItem.index;
    chiTableQuery(4,k); --Uploads chinTable[k]
    if (gridLtLabel[1]==nil) then GridviewCreateLtLabel(1); end
    if (gridLtLabel[2]==nil) then GridviewCreateLtLabel(2); end
    if (gridLtLabel[3]==nil) then GridviewCreateLtLabel(3); end

    gridLtLabel[1].SetProperty(gridLtLabel[1],"Text","Name");
    gridLtLabel[2].SetProperty(gridLtLabel[2],"Text","Index");
    gridLtLabel[3].SetProperty(gridLtLabel[3],"Text","File");
    
    gridLtLabel[1].UnHide(gridLtLabel[1]);
    gridLtLabel[2].UnHide(gridLtLabel[2]);
    gridLtLabel[3].UnHide(gridLtLabel[3]);
    
    if (gridRtFeature[1]==nil) then GridviewCreateRtFeature(1,2); end
    if (gridRtFeature[2]==nil) then GridviewCreateRtFeature(2,2); end
    if (gridRtFeature[3]==nil) then GridviewCreateRtFeature(3,2); end
    if (gridRtFeature[4]==nil) then GridviewCreateRtFeature(4,4); end
    
    gridRtFeature[1].UnHide(gridRtFeature[1]);
    gridRtFeature[2].UnHide(gridRtFeature[2]);
    gridRtFeature[3].UnHide(gridRtFeature[3]);
    gridRtFeature[4].UnHide(gridRtFeature[4]);
    
    gridRtFeature[1].SetProperty(gridRtFeature[1],"Text",chinTable[k].name);
    gridRtFeature[2].SetProperty(gridRtFeature[2],"Text",string.format("%d",k));
    if (chinTable[k].fileName=="") then
        gridRtFeature[3].SetProperty(gridRtFeature[3],"Text","None")
        gridRtFeature[4].SetProperty(gridRtFeature[4],"Text","Load")
    else
        gridRtFeature[3].SetProperty(gridRtFeature[3],"Text",chinTable[k].fileName)
        gridRtFeature[4].SetProperty(gridRtFeature[4],"Text","Reload")
    end
    
    
    item=gridRtFeature[4];
    item.parentIndex = k;
    function item.ButtonDown(this)
        fileName = chiWindowFileDialogBox();
        success=chiTableLoadFromFile(item.parentIndex,fileName);
        if (success) then
            print("Table updated")
        else
            print("Table update failed")
        end
    end
    
end