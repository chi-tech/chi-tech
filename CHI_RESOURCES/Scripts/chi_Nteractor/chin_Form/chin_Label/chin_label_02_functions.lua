--#########################################################
function LabelClass.ProcessEvents(this)
    if (WM_MOUSEMOVE.occured) then
    --
        if     (WM_MOUSEMOVE.iPar5==this.obj1Num) then
            
            if ((not this.selected) and (this.selectable)) then
                ambient=0.95;
                chiMaterialSetProperty(this.matlNum,"Diffuse",ambient,ambient,ambient,1.0);
                chiObjectSetProperty(this.obj1Num,"Renderable",true);
            end
                
        else
            if (not this.selected) then
                chiObjectSetProperty(this.obj1Num,"Renderable",false);
            end
        end
    end
    
    if (WM_LBUTTONDOWN.occured) then
        
        if     ((WM_LBUTTONDOWN.iPar5==this.obj1Num) and (this.selectable)) then
        
            this.selected=true;
            ambient=0.9;
            chiMaterialSetProperty(this.matlNum,"Diffuse",ambient,ambient,ambient,1.0);
            chiObjectSetProperty(this.obj1Num,"Renderable",true);
            this.Selected(this);
            
        else
        
            if (this.selectable) then
                if (this.selected) then
                    this.DeSelected(this);
                end
                this.selected=false;
            end
        end
    end
    
    for k=1,this.eventCallbackCount do
        this.eventCallbacks[k](this);
    end
end

--#########################################################
function LabelClass.SizeChanged(this)
    
    if ((this.master~=nil) ) then
        
        if (this.highlightFillsMaster) then
            xmin=this.master.xmin;
            xmax=this.master.xmax;
            ymin=this.master.ymin;
            ymax=this.master.ymax;
        else
            xmin=this.xpos;
            xmax=this.xpos+this.xSize;
            ymin=this.ypos;
            ymax=this.ypos+this.ySize;
        end
       
        chiSetLabelProperty3D(this.textNum,"Viewport",xmin,ymin,xmax,ymax);
        chiViewportSetProperty(this.obj1View,xmin,ymin,xmax,ymax);
        chiTransformSetScale(this.obj1Tra ,xmax-xmin-1,ymax-ymin-1,1.0);
        chiTransformSetTranslation(this.obj1Tra,xmin, ymin, this.zDepth-0.1);
        
        r=0.0; --this.fontColor[1];
        g=0.0; --this.fontColor[2];
        b=0.0; --this.fontColor[3];
        a=this.fontColor[4];
        chiSetLabel3D(this.textNum,this.text,this.xpos+this.paddingLeft,this.ypos+this.paddingBot,r,g,b,a,this.fontType)
        chiSetLabelProperty3D(this.textNum,"Depth",this.zDepth+0.01);

        if (not this.float) then
            this.xpos=this.master.cursorX;
            this.ypos=this.master.cursorY-this.ySize;
            r=this.fontColor[1];
            g=this.fontColor[2];
            b=this.fontColor[3];
            a=this.fontColor[4];
            chiSetLabel3D(this.textNum,this.text,this.xpos+this.paddingLeft,this.ypos+this.paddingBot,r,g,b,a,this.fontType)
            chiSetLabelProperty3D(this.textNum,"Depth",this.zDepth+0.01);
        end
        
    end
    

end

--######################################################### Redraw
function LabelClass.Redraw(this)
    this.SizeChanged(this);
end

--#########################################################
function LabelClass.Selected(this)
    if (not (this.CustomSelected==nil)) then
        this.CustomSelected(this);
    end
end

--#########################################################
function LabelClass.DeSelected(this)
    if (not (this.CustomDeSelected==nil)) then
        this.CustomDeSelected(this);
    end
end

--#########################################################
function LabelClass.Hide(this)
    this.SetProperty(this,"Color",{0,0,0,0});
    chiObjectSetProperty(this.obj1Num,"Hidden",true);
    this.hidden=true;
    --this.Redraw(this)
end

--#########################################################
function LabelClass.UnHide(this)
    r=0;
    g=0;
    b=0;
    a=1;
    this.SetProperty(this,"Color",{r,g,b,a});
    chiObjectSetProperty(this.obj1Num,"Hidden",false);
    this.hidden=false;
    --this.Redraw(this)
end