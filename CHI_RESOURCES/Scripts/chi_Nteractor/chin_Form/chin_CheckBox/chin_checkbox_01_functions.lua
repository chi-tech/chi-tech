--#########################################################
function CheckBoxClass.ProcessEvents(this)
    if (WM_MOUSEMOVE.occured) then
    --
        --if     (WM_MOUSEMOVE.iPar5==this.obj1Num) then
        --   
        --else
        --    
        --end
    end
    
    if (WM_LBUTTONDOWN.occured) then
        if     (WM_LBUTTONDOWN.iPar5==this.obj1Num) then
            if (this.selected) then
                
                this.DeSelected(this);
            else
                
                this.Selected(this);
            end
            this.Redraw(this)
        end
    end
    
    for k=1,this.eventCallbackCount do
        this.eventCallbacks[k](this);
    end
end

--#########################################################
function CheckBoxClass.SetProperty(this,property,value)
    if (property=="Master") then
        this.master=value;
        
        value.slaveCount=value.slaveCount+1;
        k=value.slaveCount
        value.slaves[k]=this;
     elseif (property=="Float") then
        this.float=value;
    end
end


--#########################################################
function CheckBoxClass.SizeChanged(this)
    
    if (this.master~=nil) then
        if (not this.float) then
            
            this.xpos=this.master.cursorX+1;
            this.ypos=this.master.cursorY-this.ySize;
            chiTransformSetTranslation(this.obj1Tra,this.xpos, this.ypos, this.zpos);
        end
    end
    this.Redraw(this);
end

--#########################################################
function CheckBoxClass.Redraw(this)
    chiTransformSetScale(this.obj1Tra ,this.iconSize,this.iconSize,1.0);
    chiTransformSetTranslation(this.obj1Tra,this.xpos, this.ypos, this.zpos);
    dx = this.iconOffset*this.iconType[1]/this.iconTextureSize
    dy = this.iconOffset*this.iconType[2]/this.iconTextureSize
    chiTransformSetTranslation(this.obj1TTra,this.iconPadding/this.iconTextureSize+dx,this.iconPadding/this.iconTextureSize+dy,1.0);
end

function CheckBoxClass.Hide(this)
    chiObjectSetProperty(this.obj1Num,"Hidden",true);
end

function CheckBoxClass.UnHide(this)
    chiObjectSetProperty(this.obj1Num,"Hidden",false);
end

function CheckBoxClass.Selected(this)
    this.selected=true;
    this.iconType=chinIconCheckboxC;
    if (not (this.CustomSelected==nil)) then
        this.CustomSelected(this);
    end
end
function CheckBoxClass.DeSelected(this)
    this.selected=false;  
    this.iconType=chinIconCheckboxU;
    if (not (this.CustomDeSelected==nil)) then
        this.CustomDeSelected(this);
    end
end