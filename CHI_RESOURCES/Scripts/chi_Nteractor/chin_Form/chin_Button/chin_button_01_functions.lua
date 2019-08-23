--######################################################### Process Events
function ButtonClass.ProcessEvents(this)
     if (WM_MOUSEMOVE.occured) then
    --
        if     (WM_MOUSEMOVE.iPar5==this.obj1Num) then
            chiMaterialSetProperty(this.matlNum,"Diffuse",0.8,0.8,0.9,1.0); 
            chi3DLineChangeColor(this.lineNum,0.5,0.5,0.8,this.fontColor[4]);    
            chi3DLineSetStipple(this.lineNum,false,1.0,2.0);          
        else
            if (not this.selected) then
                ambient=0.8;
                chiMaterialSetProperty(this.matlNum,"Diffuse",ambient,ambient,ambient,1.0);
                chi3DLineChangeColor(this.lineNum,0.5,0.5,0.5,this.fontColor[4]); 
                chi3DLineSetStipple(this.lineNum,false,1.0,1.0);  
            end
        end
    end
    
    if (WM_LBUTTONDOWN.occured) then
        if     ((WM_LBUTTONDOWN.iPar5==this.obj1Num) and (this.selectable)) then
            --this.selected=true;
            chiMaterialSetProperty(this.matlNum,"Diffuse",0.6,0.6,0.9,1.0); 
            chi3DLineChangeColor(this.lineNum,0.4,0.4,0.8,this.fontColor[4]);    
            chi3DLineSetStipple(this.lineNum,false,1.0,2.0);
            this.ButtonDown(this);
        else
            --this.selected=false;
        end
    end
    
    if (WM_LBUTTONUP.occured) then
        ambient=0.8;
        chiMaterialSetProperty(this.matlNum,"Diffuse",ambient,ambient,ambient,1.0);
        chi3DLineChangeColor(this.lineNum,0.5,0.5,0.5,this.fontColor[4]); 
        chi3DLineSetStipple(this.lineNum,false,1.0,1.0); 
        this.ButtonUp(this);
    end
    
    for k=1,this.eventCallbackCount do
        this.eventCallbacks[k](this);
    end
end

--#########################################################
function ButtonClass.SetProperty(this,property,value)
    if     (property=="Text") then
        this.text=value;
    elseif (property=="Position") then
        if (this.float) then
            this.xpos=value[1];
            this.ypos=value[2];
        end
    elseif (property=="Master") then
        this.master=value;
        
        value.slaveCount=value.slaveCount+1;
        k=value.slaveCount
        value.slaves[k]=this;
     elseif (property=="Float") then
        this.float=value;
    end
    this.Redraw(this);
end

--######################################################### SizeChanged
function ButtonClass.SizeChanged(this)
    
    if (this.master~=nil) then
        if (not this.float) then
            this.xpos=this.master.cursorX;
            this.ypos=this.master.cursorY-this.ySize+this.yoffset;
            
            this.Redraw(this);
        end
        
    end
end


--######################################################### Redraw
function ButtonClass.Redraw(this)
    stringWidth=chinGetStringPixelWidth(this.text);
    textxpos=this.xpos+0.5*this.xSize-0.5*stringWidth-2*this.paddingLeft;
    r=this.fontColor[1];
    g=this.fontColor[2];
    b=this.fontColor[3];
    a=this.fontColor[4];

    chiSetLabel(this.textNum,this.text,textxpos,this.ypos+this.paddingBot+this.textPadding,r,g,b,a,this.fontType)
    
    chiTransformSetTranslation(this.obj1Tra ,this.xpos+this.paddingLeft,this.ypos+this.paddingBot,this.zDepth);
    chiTransformSetScale(this.obj1Tra ,this.xSize-2*this.paddingLeft-1,this.ySize-2*this.paddingBot+1,1.0);
    
    chi3DLineChangeVertex(this.lineNum,0,this.xpos+this.paddingLeft,this.ypos+this.paddingBot,this.zDepth);
    chi3DLineChangeVertex(this.lineNum,1,this.xpos+this.paddingLeft,this.ypos+this.ySize-2*this.paddingBot+2,this.zDepth);
    chi3DLineChangeVertex(this.lineNum,2,this.xpos+this.xSize-2*this.paddingLeft+2,this.ypos+this.ySize-2*this.paddingBot+2,this.zDepth);
    chi3DLineChangeVertex(this.lineNum,3,this.xpos+this.xSize-2*this.paddingLeft+2,this.ypos+this.paddingBot,this.zDepth);
    chi3DLineChangeVertex(this.lineNum,4,this.xpos-1+this.paddingLeft,this.ypos+this.paddingBot,this.zDepth);
end

function ButtonClass.Hide(this)
    this.hidden =  true;
    this.fontColor[4]=0.0;
    chi3DLineChangeColor(this.lineNum,0.5,0.5,0.5,0.0); 
    chiObjectSetProperty(this.obj1Num,"Hidden",true);
    this.Redraw(this);
end

function ButtonClass.UnHide(this)
    this.hidden = false;
    this.fontColor[4]=1.0;
    chi3DLineChangeColor(this.lineNum,0.5,0.5,0.5,1.0); 
    chiObjectSetProperty(this.obj1Num,"Hidden",false);
    this.Redraw(this);
end

--######################################################### ButtonDown
function ButtonClass.ButtonDown(this)
    if (not (this.CustomButtonDown==nil)) then
        this.CustomButtonDown(this);
    end
end

--######################################################### ButtonUp
function ButtonClass.ButtonUp(this)
    if (not (this.CustomButtonUp==nil)) then
        this.CustomButtonUp(this);
    end
end