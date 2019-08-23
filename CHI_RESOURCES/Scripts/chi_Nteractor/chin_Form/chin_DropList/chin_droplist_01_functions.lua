--#########################################################
function DropListClass.ProcessEvents(this)
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
                this.text=this.list[k].text;
                this.selectedOption = k;
                this.SelectionMade(this);
                this.Redraw(this);
                print("Redrawn")
            end
        end
    end
    
    for k=1,this.eventCallbackCount do
        this.eventCallbacks[k](this);
    end
end

--#########################################################
function DropListClass.SetProperty(this,property,value)
    if     (property=="Master") then
        this.master=value;
        
        value.slaveCount=value.slaveCount+1;
        k=value.slaveCount
        value.slaves[k]=this;
    elseif (property=="Float") then
        this.float=value;
    elseif (property=="Text") then
        this.text=value;
        r=this.fontColor[1];
        g=this.fontColor[2];
        b=this.fontColor[3];
        a=this.fontColor[4];
        chiSetLabel3D(this.textNum,this.text,this.xpos+this.marginLeft+this.paddingLeft,this.ypos+this.marginBot+this.paddingBot,r,g,b,a,this.fontType)
    elseif (property=="Color") then
        this.fontColor[1]=value[1];
        this.fontColor[2]=value[2];
        this.fontColor[3]=value[3];
        this.fontColor[4]=value[4];
        r=value[1];
        g=value[2];
        b=value[3];
        a=value[4];
        chiSetLabel3D(this.textNum,this.text,this.xpos+this.marginLeft+this.paddingLeft,this.ypos+this.marginBot+this.paddingBot,r,g,b,a,this.fontType)
    end
end

--#########################################################
function DropListClass.SizeChanged(this)
    if (this.master~=nil) then
        
        if (not this.float) then
            this.xpos=this.master.cursorX;
            this.ypos=this.master.cursorY-this.ySize;
            
            if (this.xSize>(this.master.xmax-this.master.xmin)) then
                this.xSize = this.master.xmax-this.master.xmin;
            end
            
            this.Redraw(this)
        end

    end
    
    this.Redraw(this);
end

--#########################################################
function DropListClass.Redraw(this)
    chiTransformSetScale(this.transform1,this.xSize-1,this.ySize-1,1.0)
    chiTransformSetTranslation(this.transform1,this.xpos+this.paddingLeft,this.ypos+this.paddingBot,this.zDepth)

    chi3DLineChangeVertex(this.lineNum,0,this.xpos+this.paddingLeft,this.ypos+this.paddingBot,this.zDepth+0.1);
    chi3DLineChangeVertex(this.lineNum,1,this.xpos+this.paddingLeft,this.ypos+this.ySize-1+this.paddingBot,this.zDepth+0.1);
    chi3DLineChangeVertex(this.lineNum,2,this.xpos+this.xSize-1+this.paddingLeft,this.ypos+this.ySize-1+this.paddingBot,this.zDepth+0.1);
    chi3DLineChangeVertex(this.lineNum,3,this.xpos+this.xSize-1+this.paddingLeft,this.ypos+this.paddingBot,this.zDepth+0.1);
    chi3DLineChangeVertex(this.lineNum,4,this.xpos-1+this.paddingLeft,this.ypos+this.paddingBot,this.zDepth+0.1);
    
    r=this.fontColor[1];
    g=this.fontColor[2];
    b=this.fontColor[3];
    a=this.fontColor[4];
    chiSetLabel3D(this.textNum,this.text,this.xpos+this.paddingLeft+this.marginLeft,this.ypos+this.paddingBot+this.marginBot,r,g,b,a,this.fontType)
    chiSetLabelProperty(this.textNum,"ViewportEnable",true);
    chiSetLabelProperty(this.textNum,"Viewport",this.xpos+this.paddingLeft,this.ypos+this.paddingBot,this.xpos+this.xSize-1+this.paddingLeft-this.iconSize,this.ypos+this.ySize-1+this.paddingBot);
    chiSetLabelProperty3D(this.textNum,"Depth",this.zDepth+0.01);
    
    chiTransformSetTranslation(this.obj2Tra , this.xpos+this.xSize-this.iconSize,  this.ypos-2, this.zDepth+0.01)
        
    this.panel1.zDepth=this.zDepth;
    this.panel1.xmin=this.xpos;
    this.panel1.ymax=this.ypos-2;
    this.panel1.xmax=this.panel1.xmin+this.panelxSize+1;
    this.panel1.ymin=this.panel1.ymax-this.ySize*this.listCount-1;
    
    for k=1,this.listCount do
        this.list[k].xSize=this.panelxSize+2;
    end
    
    this.panel1.SizeChanged(this.panel1);
end

--######################################################### Add list item
function DropListClass.AddListItem(this,text)
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
    
    return newLabel;
end

--######################################################### Selected event
function DropListClass.Selected(this)
    this.panel1.UnHide(this.panel1);
    
    
    for k=1,this.listCount do
        this.list[k].UnHide(this.list[k]);
    end
    
    this.Redraw(this);
end

--######################################################### Selected event
function DropListClass.DeSelected(this)
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

function DropListClass.Hide(this)
    chi3DLineChangeColor(this.lineNum,0.0,0.0,0.0,0.0);
    this.SetProperty(this,"Color",{0,0,0,0});
    chiObjectSetProperty(this.obj1Num,"Hidden",true);
    chiObjectSetProperty(this.obj2Num,"Hidden",true);
end

function DropListClass.UnHide(this)
    
    chi3DLineChangeColor(this.lineNum,0.0,0.0,0.0,1.0);
    r=0;
    g=0;
    b=0;
    a=1;
    this.SetProperty(this,"Color",{r,g,b,a});
    chiObjectSetProperty(this.obj1Num,"Hidden",false);
    chiObjectSetProperty(this.obj2Num,"Hidden",false);
end

function DropListClass.SelectionMade(this)
    if (this.CustomSelectionMade~=nil) then
        this.CustomSelectionMade(this);
    end
end











