function GridviewGetPropertiesForMaterials(selectionItem)

    local k=selectionItem.originFeature.parent.index;
       
    if (gridLtLabel[1]==nil) then GridviewCreateLtLabel(1); end
    if (gridLtLabel[2]==nil) then GridviewCreateLtLabel(2); end

    gridLtLabel[1].SetProperty(gridLtLabel[1],"Text","Name");
    gridLtLabel[2].SetProperty(gridLtLabel[2],"Text","Index");
    
    gridLtLabel[1].UnHide(gridLtLabel[1]);
    gridLtLabel[2].UnHide(gridLtLabel[2]);
    
    if (gridRtFeature[1]==nil) then GridviewCreateRtFeature(1,2); end
    if (gridRtFeature[2]==nil) then GridviewCreateRtFeature(2,2); end
    
    gridRtFeature[1].UnHide(gridRtFeature[1]);
    gridRtFeature[2].UnHide(gridRtFeature[2]);
    
    gridRtFeature[1].SetProperty(gridRtFeature[1],"Text",chinMaterial[k].name);
    gridRtFeature[2].SetProperty(gridRtFeature[2],"Text",string.format("%d",k));
end