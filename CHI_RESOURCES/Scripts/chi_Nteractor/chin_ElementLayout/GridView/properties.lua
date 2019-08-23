function GridviewGetPropertiesForProperties(selectionItem)
    property=selectionItem.originFeature.parent;
    if (gridLtLabel[1]==nil) then GridviewCreateLtLabel(1); end
    --if (gridRtFeature[1]==nil) then GridviewCreateRtFeature(1,2); end
    if (property.propType==7) then
        GridviewCreateRtFeature(1,5);
        gridLtLabel[1].UnHide(gridLtLabel[1]);
        gridRtFeature[1].UnHide(gridRtFeature[1]);
    
        gridLtLabel[1].SetProperty(gridLtLabel[1],"Text",property.name);
        gridRtFeature[1].parent=property;
    else
        GridviewCreateRtFeature(1,2);
        gridLtLabel[1].UnHide(gridLtLabel[1]);
        gridRtFeature[1].UnHide(gridRtFeature[1]);
    
        gridLtLabel[1].SetProperty(gridLtLabel[1],"Text",property.name);
        gridRtFeature[1].parent=property;
    end
    
    --=========================================== 1 Name
    if     (property.propType==1) then --name
        gridRtFeature[1].SetProperty(gridRtFeature[1],"Text",property.value);
    --=========================================== 2 Color
    elseif (property.propType==2) then --color   
        gridRtFeature[1].SetProperty(gridRtFeature[1],"Text",string.format("{%d,%d,%d,%d}",math.floor(property.r*255),math.floor(property.g*255),math.floor(property.b*255),math.floor(property.a*255)));
        item=gridRtFeature[1];
        function item.CustomKeyPress(this)
            if (WM_CHAR.iPar0==13) then
                chunk,err=assert(load("temp="..this.text));
                if (chunk==nil) then print(err);
                else chunk() end
                
                    gridRtFeature[1].parent.r=temp[1]/255;
                    gridRtFeature[1].parent.g=temp[2]/255;
                    gridRtFeature[1].parent.b=temp[3]/255;
                    gridRtFeature[1].parent.a=temp[4]/255;
                    gridRtFeature[1].parent.ValueChanged(gridRtFeature[1].parent);
    
                this.cursorPosition=0;
                this.cursorX=0;
                this.Redraw(this);
            end
        end
    --=========================================== 3 Scalar
    elseif (property.propType==3) then --scalar
        fmt=property.displayFormat;
        
        gridRtFeature[1].SetProperty(gridRtFeature[1],"Text",string.format(fmt,property.value));
        item=gridRtFeature[1];
        local func=function (this)
            if (WM_CHAR.iPar0==13) then
                chunk,err=assert(load("temp="..this.text));
                
                if (chunk==nil) then print(err);
                else chunk() end
                    gridRtFeature[1].parent.value=temp;
                    gridRtFeature[1].parent.ValueChanged(gridRtFeature[1].parent);
    
                this.cursorPosition=0;
                this.cursorX=0;
                this.Redraw(this);
            end
        end
        item.CustomKeyPress=func
    --=========================================== 4 Vec2
    elseif (property.propType==4) then --vec2
        fmt=property.displayFormat;
        gridRtFeature[1].SetProperty(gridRtFeature[1],"Text",string.format("{"..fmt.."}",property.value[1],property.value[2]));
        item=gridRtFeature[1];
        function func(this)
            if (WM_CHAR.iPar0==13) then
                chunk,err=assert(load("temp="..this.text));
                if (chunk==nil) then print(err);
                else chunk() end
                
                    gridRtFeature[1].parent.value=temp;
                    gridRtFeature[1].parent.ValueChanged(gridRtFeature[1].parent);
    
                this.cursorPosition=0;
                this.cursorX=0;
                this.Redraw(this);
            end
        end
        item.CustomKeyPress=func
    --=========================================== 5 Vec3
    elseif (property.propType==5) then --vec3
        fmt=property.displayFormat;
        gridRtFeature[1].SetProperty(gridRtFeature[1],"Text",string.format("{"..fmt.."}",property.value[1],property.value[2],property.value[3]));
        item=gridRtFeature[1];
        function func(this)
            if (WM_CHAR.iPar0==13) then
                chunk,err=assert(load("temp="..this.text));
                if (chunk==nil) then print(err);
                else chunk() end
                
                    gridRtFeature[1].parent.value=temp;
                    gridRtFeature[1].parent.ValueChanged(gridRtFeature[1].parent);
    
                this.cursorPosition=0;
                this.cursorX=0;
                this.Redraw(this);
            end
        end
        item.CustomKeyPress=func
    --=========================================== 6 Vec4
    elseif (property.propType==6) then --vec4
        fmt=property.displayFormat;
        gridRtFeature[1].SetProperty(gridRtFeature[1],"Text",string.format("{"..fmt.."}",property.value[1],property.value[2],property.value[3],property.value[4]));
        item=gridRtFeature[1];
        local func=function (this)
            
            if (WM_CHAR.iPar0==13) then
                chunk,err=assert(load("temp="..this.text));
                if (chunk==nil) then print(err);
                else chunk() end
                
                    gridRtFeature[1].parent.value=temp;
                    gridRtFeature[1].parent.ValueChanged(gridRtFeature[1].parent);
    
                this.cursorPosition=0;
                this.cursorX=0;
                this.Redraw(this);
            end
        end
        item.CustomKeyPress=func
    --=========================================== 7 Droplist
    elseif (property.propType==7) then
        for k=1,property.valueCount do
            gridRtFeature[1].AddListItem(gridRtFeature[1],property.value[k]);
        end
        if (property.initOption>0) then
            gridRtFeature[1].SetProperty(gridRtFeature[1],"Text",property.value[property.initOption]);
        else
            gridRtFeature[1].SetProperty(gridRtFeature[1],"Text","Choose");
        end
        item=gridRtFeature[1];
        function item.CustomSelectionMade(this)
            chunk,err=assert(load("temp="..string.format("%d",this.selectedOption)));
            if (chunk==nil) then print(err);
            else chunk() end
            
                gridRtFeature[1].parent.option=temp;
                gridRtFeature[1].parent.ValueChanged(gridRtFeature[1].parent);
        end
    end
end