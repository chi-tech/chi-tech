--#########################################################
function ContextMenuClass.ProcessEvents(this)
    if (WM_LBUTTONDOWN.occured) then
        if     ((WM_LBUTTONDOWN.iPar5==this.obj1Num)) then
            this.selected=true;
            this.Selected(this);
            
        else
            this.DeSelected(this);
            this.selected=false;
        end
        for k=1,this.listCount do
            if     ((WM_LBUTTONDOWN.iPar5==this.list[k].obj1Num)) then
                --this.text=this.list[k].text;
                --this.Redraw(this);
            end
        end
    end
    
    for k=1,this.eventCallbackCount do
        this.eventCallbacks[k](this);
    end
end

--#########################################################
function ContextMenuClass.SetProperty(this,property,value)
    if     (property=="Master") then
        this.master=value;
        
        value.slaveCount=value.slaveCount+1;
        k=value.slaveCount
        value.slaves[k]=this;
    elseif (property=="Float") then
        this.float=value;
    end
end

--#########################################################
function ContextMenuClass.SizeChanged(this)
    if (this.master~=nil) then
        
        if (not this.float) then
            this.xpos=this.master.cursorX;
            this.ypos=this.master.cursorY-this.ySize;
            
            this.Redraw(this)
        end

    end
    
    --this.Redraw(this);
end

--#########################################################
function ContextMenuClass.Redraw(this)
    chiTransformSetScale(this.transform1,this.xSize-1,this.ySize-1,1.0)
    chiTransformSetTranslation(this.transform1,this.xpos+this.paddingLeft,this.ypos+this.paddingBot,this.zDepth)

    chi3DLineChangeVertex(this.lineNum,0,this.xpos+1+this.paddingLeft,this.ypos+this.paddingBot,this.zDepth+0.1);
    chi3DLineChangeVertex(this.lineNum,1,this.xpos+1+this.paddingLeft,this.ypos+this.ySize-1+this.paddingBot,this.zDepth+0.1);
    chi3DLineChangeVertex(this.lineNum,2,this.xpos+1+this.xSize+this.paddingLeft,this.ypos+this.ySize-1+this.paddingBot,this.zDepth+0.1);
    chi3DLineChangeVertex(this.lineNum,3,this.xpos+1+this.xSize+this.paddingLeft,this.ypos+this.paddingBot,this.zDepth+0.1);
    chi3DLineChangeVertex(this.lineNum,4,this.xpos+this.paddingLeft,this.ypos+this.paddingBot,this.zDepth+0.1);

    if (this.showOutline) then
        chi3DLineChangeColor(this.lineNum,0.0,0.0,0.0,1.0);
    else
        chi3DLineChangeColor(this.lineNum,0.0,0.0,0.0,0.0);
    end
    
    r=this.fontColor[1];
    g=this.fontColor[2];
    b=this.fontColor[3];
    a=this.fontColor[4];

    chiSetLabel3D(this.textNum,this.text,this.xpos+this.paddingLeft+this.marginLeft,this.ypos+this.paddingBot+this.marginBot,r,g,b,0,this.fontType)
    chiSetLabelProperty(this.textNum,"ViewportEnable",true);
    chiSetLabelProperty(this.textNum,"Viewport",this.xpos+this.paddingLeft,this.ypos+this.paddingBot,this.xpos+this.xSize-1+this.paddingLeft-this.iconSize,this.ypos+this.ySize-1+this.paddingBot);
    chiSetLabelProperty3D(this.textNum,"Depth",this.zDepth+0.01);
    
    --chiTransformSetTranslation(this.obj2Tra , this.xpos+this.xSize-this.iconSize,  this.ypos-2, this.zDepth+0.01)
        --print(this.panel1.zDepth)
    this.panel1.zDepth=this.zDepth;
    this.panel1.xmin=this.xpos;
    this.panel1.ymax=this.ypos+1;
    this.panel1.xmax=this.panel1.xmin+this.panelxSize+1;
    this.panel1.ymin=this.panel1.ymax-this.ySize*this.listCount-1;
    
    for k=1,this.listCount do
        this.list[k].xSize=this.panelxSize;
    end
    
    this.panel1.SizeChanged(this.panel1);
end

--######################################################### Add list item
function ContextMenuClass.AddListItem(this,text)
    this.listCount=this.listCount+1;
    k=this.listCount;
    this.list[k]=LabelClass.New(text);
    local newLabel=this.list[k];
    this.list[k].SetProperty(this.list[k],"Float",false);
    this.list[k].selectable=true;
    this.list[k].zDepth=this.zDepth+0.1;
    this.list[k].SetProperty(this.list[k],"Master",this.panel1);
    this.panel1.SizeChanged(this.panel1);
    this.list[k].Hide(this.list[k]);
    dRegisterFeature(newLabel)
    return newLabel;
end

--######################################################### Selected event
function ContextMenuClass.Selected(this)
    this.panel1.UnHide(this.panel1);
    
    for k=1,this.listCount do
        this.list[k].UnHide(this.list[k]);
    end
    
    this.Redraw(this);
end

--######################################################### Selected event
function ContextMenuClass.DeSelected(this)
    this.panel1.Hide(this.panel1);
    
    for k=1,this.listCount do
        this.list[k].Hide(this.list[k]);
    end
    --currentScene=chiGetScene();
    --chiBindScene(1);
    --chiSetWindowProperties("RESTORE_REGION");
    --chiBindScene(currentScene);
    this.Redraw(this);
end












