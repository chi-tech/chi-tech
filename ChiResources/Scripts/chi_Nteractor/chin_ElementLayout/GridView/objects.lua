function GridviewGetPropertiesForObjects(selectionItem)
    local k=selectionItem.originFeature.parent.index;
        
    --chiObjectQuery(4,k); --Uploads chinMaterial[k]
    if (gridLtLabel[1]==nil) then GridviewCreateLtLabel(1); end
    if (gridLtLabel[2]==nil) then GridviewCreateLtLabel(2); end
    if (gridLtLabel[3]==nil) then GridviewCreateLtLabel(3); end

    gridLtLabel[1].SetProperty(gridLtLabel[1],"Text","Name");
    gridLtLabel[2].SetProperty(gridLtLabel[2],"Text","Index");
    gridLtLabel[3].SetProperty(gridLtLabel[3],"Text","Wireframe");
    
    gridLtLabel[1].UnHide(gridLtLabel[1]);
    gridLtLabel[2].UnHide(gridLtLabel[2]);
    gridLtLabel[3].UnHide(gridLtLabel[3]);
    
    if (gridRtFeature[1]==nil) then GridviewCreateRtFeature(1,2); end
    if (gridRtFeature[2]==nil) then GridviewCreateRtFeature(2,2); end
    if (gridRtFeature[3]==nil) then GridviewCreateRtFeature(3,3); end
    
    gridRtFeature[1].UnHide(gridRtFeature[1]);
    gridRtFeature[2].UnHide(gridRtFeature[2]);
    gridRtFeature[3].UnHide(gridRtFeature[3]);
    
    gridRtFeature[1].SetProperty(gridRtFeature[1],"Text",chinObject[k].name);
    gridRtFeature[2].SetProperty(gridRtFeature[2],"Text",string.format("%d",k));
    
    item=gridRtFeature[3];
    item.index=k;
    
    function item.CustomSelected(this)
        chiObjectSetProperty(this.index,"Wireframe",true);
    end
    function item.CustomDeSelected(this)
        chiObjectSetProperty(this.index,"Wireframe",false);
    end
end