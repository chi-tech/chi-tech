--######################################################### Events
function PanelClass.ProcessEvents(this)
    --======================= Size changed
    if (WM_SIZE.occured) then
        this.UpdateBoundaryLocks(this);
        this.WindowSizeChanged(this);

    end
    --this.UpdateBoundaryLocks(this);
    --======================= Mouse Move
    if (WM_MOUSEMOVE.occured) then
        --print("WINDOW MOUSE MOVE:",WM_MOUSEMOVE.iPar4)

        -- ==================== Main Panel
        if     (WM_MOUSEMOVE.iPar5==this.obj1Num) then
            chiWindowSetCursor(0);
        end

        --Scrollbar cursor change
        if     (WM_MOUSEMOVE.iPar5==this.objSNum) then
                chiWindowSetCursor(2);
        end

        --Scroll icon movement
        if      (this.scrollActive) then
            this.scrollyf=WM_MOUSEMOVE.iPar1-this.scrollHeight*0.5;
            this.scrollOffset=this.scrollyf-(chinGlobal.dwindowysize-this.ymax);
            --print(this.scrollOffset,this.contentYSize,this.ymax-this.ymin,this.scrollScale);
            if (this.scrollOffset<0.0) then this.scrollOffset=0.0; end
            if (this.scrollOffset>(this.ymax-this.ymin-this.scrollHeight-1)) then this.scrollOffset=this.ymax-this.ymin-this.scrollHeight-1; end
            this.SizeChanged(this);
        end

        -- ==================== Window Resize
        if ((not this.sizeChanging) and (this.ownResize)) then
            if     ((WM_MOUSEMOVE.iPar5==this.obj2Num) or (WM_MOUSEMOVE.iPar5==this.obj3Num)) then
                chiWindowSetCursor(1);
            elseif ((WM_MOUSEMOVE.iPar5==this.obj4Num) or (WM_MOUSEMOVE.iPar5==this.obj5Num)) then
                chiWindowSetCursor(2);
            else
                chiWindowSetCursor(0);
            end
        elseif (this.ownResize) then
            if (WM_MOUSEMOVE.iPar4==0) then
                if (this.currentGripper==2) then
                    this.xmax=WM_MOUSEMOVE.iPar0 --- WM_MOUSEMOVE.iPar2;
                    this.SizeChanged(this);--
                    this.Redraw(this);   
                elseif (this.currentGripper==3) then       --
                    this.xmin=WM_MOUSEMOVE.iPar0 --- WM_MOUSEMOVE.iPar2;
                    this.SizeChanged(this);--
                    this.Redraw(this);                     --
                elseif (this.currentGripper==4) then       --
                    this.ymax=chinGlobal.dwindowysize-WM_MOUSEMOVE.iPar1 --+ WM_MOUSEMOVE.iPar3;
                    this.SizeChanged(this);--
                    this.Redraw(this);                     --
                elseif (this.currentGripper==5) then       --
                    this.ymin=chinGlobal.dwindowysize-WM_MOUSEMOVE.iPar1 --+ WM_MOUSEMOVE.iPar3;
                    this.SizeChanged(this);--
                    this.Redraw(this);
                end
               
            else
                if (this.currentGripper==2) then
                    this.xmax=WM_MOUSEMOVE.iPar0 + this.xmin--- WM_MOUSEMOVE.iPar2;
                    this.SizeChanged(this);--
                    this.Redraw(this);   
                elseif (this.currentGripper==3) then       --
                    this.xmin=WM_MOUSEMOVE.iPar0 + this.xmin--- WM_MOUSEMOVE.iPar2;
                    this.SizeChanged(this);--
                    this.Redraw(this);                     --
                elseif (this.currentGripper==4) then       --
                    this.ymax=this.ymax - WM_MOUSEMOVE.iPar1 --+ WM_MOUSEMOVE.iPar3;
                    this.SizeChanged(this);--
                    this.Redraw(this);                     --
                elseif (this.currentGripper==5) then       --
                    this.ymin=this.ymax - WM_MOUSEMOVE.iPar1 --+ WM_MOUSEMOVE.iPar3;
                    this.SizeChanged(this);--
                    this.Redraw(this);
                end
            end
          
        end
    end
    
    --======================= Left Mouse Button
    if (WM_LBUTTONDOWN.occured) then
        if     (WM_LBUTTONDOWN.iPar5==this.obj2Num) then
            this.currentGripper=2;
            this.sizeChanging=true;
        elseif (WM_LBUTTONDOWN.iPar5==this.obj3Num) then
            this.currentGripper=3;
            this.sizeChanging=true;
        elseif (WM_LBUTTONDOWN.iPar5==this.obj4Num) then
            this.currentGripper=4;
            this.sizeChanging=true;
        elseif (WM_LBUTTONDOWN.iPar5==this.obj5Num) then
            this.currentGripper=5;
            this.sizeChanging=true;
        end
        if     (WM_LBUTTONDOWN.iPar5==this.objSNum) then
            this.SizeChanged(this);
            this.scrollActive       = true;
            this.scrollyi=WM_MOUSEMOVE.iPar1-this.scrollHeight;
            --this.scrollScale=(this.ymax-this.cursorY)/(this.ymax-this.ymin);
            this.scrollScale=(this.contentYSize)/(this.ymax-this.ymin);
            --print(this.contentYSize,,this.ymax-this.ymin,this.cursorY)
            if (this.scrollScale<1.0) then this.scrollScale=0; end
        end
    end

    if (WM_LBUTTONUP.occured) then
        this.currentGripper=0;
        this.sizeChanging=false;
        chiWindowSetCursor(0);
        this.scrollActive       = false;
    end
    
    for k=1,this.eventCallbackCount do
        this.eventCallbacks[k](this);
    end
end

--######################################################### Attach Window
function PanelClass.AttachWindow(this,windowNumber)
    this.windowAttached=true;
    this.window=windowNumber;
    
    curScene=chiGetScene();
    
    parentWindowx,parentWindowy=chiGetWindowProperties();
    --print(parentWindowx,parentWindowy)
    
    chiBindScene(windowNumber);
    
    chiSetWindowProperties(this.xmax-this.xmin-1,this.ymax-this.ymin-1,this.xmin,parentWindowy-this.ymax+1);
    
    chiBindScene(curScene);
end


--######################################################### Redraw 
function PanelClass.Redraw(this)
    chiTransformSetScale(this.transform1,this.xmax-this.xmin-1-this.scrollWidth,this.ymax-this.ymin-1,1.0)
    chiTransformSetTranslation(this.transform1,this.xmin,this.ymin,this.zDepth)
    
    girth =this.xmax-this.xmin-1;
    length=this.ymax-this.ymin-1;
    chiTransformSetScale(this.transform2,this.gripWidth,length,1); --E
    chiTransformSetScale(this.transform3,this.gripWidth,length,1); --W
    chiTransformSetScale(this.transform4,girth,this.gripWidth,1); --N
    chiTransformSetScale(this.transform5,girth,this.gripWidth,1); --S

    chiTransformSetTranslation(this.transform2,this.xmax                 -0.5*this.gripWidth,this.ymin,this.zDepth+0.1); --E
    chiTransformSetTranslation(this.transform3,this.xmin-this.gripWidth-1+0.5*this.gripWidth,this.ymin,this.zDepth+0.1); --W
    chiTransformSetTranslation(this.transform4,this.xmin                 ,this.ymax-0.5*this.gripWidth,this.zDepth+0.1); --N
    chiTransformSetTranslation(this.transform5,this.xmin,this.ymin-this.gripWidth-1+0.5*this.gripWidth,this.zDepth+0.1); --S

    chi3DLineChangeVertex(this.lineNum,0,this.xmin  ,this.ymin  ,this.zDepth);
    chi3DLineChangeVertex(this.lineNum,1,this.xmin  ,this.ymax  ,this.zDepth);
    chi3DLineChangeVertex(this.lineNum,2,this.xmax-this.scrollWidth,this.ymax  ,this.zDepth);
    chi3DLineChangeVertex(this.lineNum,3,this.xmax-this.scrollWidth,this.ymin  ,this.zDepth);
    chi3DLineChangeVertex(this.lineNum,4,this.xmin-1,this.ymin,this.zDepth);

    chi3DLineChangeVertex(this.scrollBorder, 0, this.xmax-this.scrollWidth,this.ymin, this.zDepth+this.scrollDepthTest);
    chi3DLineChangeVertex(this.scrollBorder, 1, this.xmax-this.scrollWidth,this.ymax, this.zDepth+this.scrollDepthTest);
    chi3DLineChangeVertex(this.scrollBorder, 2, this.xmax                 ,this.ymax, this.zDepth+this.scrollDepthTest);
    chi3DLineChangeVertex(this.scrollBorder, 3, this.xmax                 ,this.ymin, this.zDepth+this.scrollDepthTest);
    chi3DLineChangeVertex(this.scrollBorder, 4, this.xmax-this.scrollWidth,this.ymin, this.zDepth+this.scrollDepthTest);
    
    chiTransformSetScale(this.transformS,this.scrollWidth-1,this.scrollHeight,1.0)
    chiTransformSetScale(this.transformS2,this.scrollWidth-1,this.ymax-this.ymin-1,1.0)
    chiTransformSetTranslation(this.transformS,this.xmax-this.scrollWidth,this.ymax-this.scrollHeight-1-this.scrollOffset,this.zDepth+1.0)
    chiTransformSetTranslation(this.transformS2,this.xmax-this.scrollWidth,this.ymin,this.zDepth+0.9)
    --print(this.name,this.cutWindow,this.cutPanel,this.hidden)
    if (not this.hidden) then
        
        if ((this.cutWindow>0) and (not (this.cutPanel==nil))) then
            currentScene=chiGetScene();
            chiBindScene(this.cutWindow);
            xmin=this.xmin-this.cutPanel.xmin-1; if (xmin<=0) then xmin=0; end
            xmax=this.xmax-this.cutPanel.xmin; if (xmax<=0) then xmax=0; end
            ymin=this.cutPanel.ymax-this.ymax; if (ymin<=0) then ymax=0; end
            ymax=this.cutPanel.ymax-this.ymin;
            --print(xmin,ymin,xmax,ymax)
            chiSetWindowProperties("CUT_REGION",xmin,ymin,xmax,ymax);
            chiBindScene(currentScene);
        end
    else
        if ((this.cutWindow>0) and (not (this.cutPanel==nil))) then
            currentScene=chiGetScene();
            chiBindScene(this.cutWindow);
            chiSetWindowProperties("RESTORE_REGION");
            chiBindScene(currentScene);
        end
    end
end





--######################################################### Dimension request
function PanelClass.RequestDimension(this,side,value)
    if     (side=="EAST") then
        if ((value-this.xmin)<this.widthMin) then
            return false;
        else
            this.xmax=value;
        end
    elseif (side=="WEST") then
        if ((this.xmax-value)<this.widthMin) then
            return false;
        else
            this.xmin=value;
        end
    elseif (side=="NORTH") then
        if ((value-this.ymin)<this.lengthMin) then
            return false;
        else
            this.ymax=value;
        end
    elseif (side=="SOUTH") then
        if ((this.ymax-value)<this.lengthMin) then
            return false;
        else
            this.ymin=value;
        end
    end
    this.SizeChanged(this);
    this.Redraw(this);
    return true;
end

--######################################################### Dimension query
function PanelClass.QueryDimension(this,side,value)
    if     (side=="EAST") then
        if ((value-this.xmin)<this.widthMin) then
            return false;
        end
    elseif (side=="WEST") then
        if ((this.xmax-value)<this.widthMin) then
            return false;
        end
    elseif (side=="NORTH") then
        if ((value-this.ymin)<this.lengthMin) then
            return false;
        end
    elseif (side=="SOUTH") then
        if ((this.ymax-value)<this.lengthMin) then
            return false;
        end
    end

    return true;
end


--######################################################### Boundary lock
function PanelClass.AddBoundaryLock(this,side,offset,feature,lockopposite)
    this.boundaryLockCount=this.boundaryLockCount+1;
    
    this.boundaryLocks[this.boundaryLockCount]={}
    this.boundaryLocks[this.boundaryLockCount].side=side;
    if (lockopposite~=nil) then
        this.boundaryLocks[this.boundaryLockCount].lockOpposite=lockopposite;
    else
        this.boundaryLocks[this.boundaryLockCount].lockOpposite=false;
    end
    
    this.boundaryLocks[this.boundaryLockCount].offset=0;
    if (offset~=nil) then
        this.boundaryLocks[this.boundaryLockCount].offset=offset;
    end
    
    
    if (feature==nil) then
        this.boundaryLocks[this.boundaryLockCount].window=true;
    else
        this.boundaryLocks[this.boundaryLockCount].window=false;
        this.boundaryLocks[this.boundaryLockCount].feature=feature;
        feature.bslaveCount = feature.bslaveCount + 1;
        feature.bslaves[feature.bslaveCount]=this;
    end
    this.UpdateBoundaryLocks(this);
    this.Redraw(this);
end

--######################################################### Update boundary locks
function PanelClass.UpdateBoundaryLocks(this)
    height=this.ymax-this.ymin;
    width =this.xmax-this.xmin;
    
    for k=1,this.boundaryLockCount do
        
        if (this.boundaryLocks[k].window) then
            parentWindowx,parentWindowy=chiGetWindowProperties();
            --print(parentWindowx,parentWindowy);
            if     (this.boundaryLocks[k].side=="NORTH") then
                this.ymax=parentWindowy-this.boundaryLocks[k].offset;
                if (this.fixedHeight) then this.ymin=this.ymax-height; end
            elseif (this.boundaryLocks[k].side=="SOUTH") then
                this.ymin=0+this.boundaryLocks[k].offset;
                if (this.fixedHeight) then this.ymax=this.ymin+height; end
            elseif (this.boundaryLocks[k].side=="EAST") then
                this.xmax=parentWindowx-this.boundaryLocks[k].offset;
                if (this.fixedWidth) then this.xmin=this.xmax-width; end
            elseif (this.boundaryLocks[k].side=="WEST") then
                this.xmin=0+this.boundaryLocks[k].offset;
                if (this.fixedWidth) then this.xmax=this.xmin+width; end
            end
            
            this.SizeChanged(this);
        else
            
            theFeature=this.boundaryLocks[k].feature;
            if     (this.boundaryLocks[k].side=="NORTH") then
                if (this.boundaryLocks[k].lockOpposite) then
                    this.ymax=theFeature.ymax-this.boundaryLocks[k].offset;
                else
                    this.ymax=theFeature.ymin-this.boundaryLocks[k].offset;
                end
                if (this.fixedHeight) then this.ymin=this.ymax-height; end
            elseif (this.boundaryLocks[k].side=="SOUTH") then
                if (this.boundaryLocks[k].lockOpposite) then
                    this.ymin=theFeature.ymin+this.boundaryLocks[k].offset;
                else
                    this.ymin=theFeature.ymax+this.boundaryLocks[k].offset;
                end
                if (this.fixedHeight) then this.ymax=this.ymin+height; end
            elseif (this.boundaryLocks[k].side=="EAST") then
                if (this.boundaryLocks[k].lockOpposite) then
                    this.xmax=theFeature.xmax-this.boundaryLocks[k].offset;
                else
                    this.xmax=theFeature.xmin-this.boundaryLocks[k].offset;
                end
                if (this.fixedWidth) then this.xmin=this.xmax-width; end
            elseif (this.boundaryLocks[k].side=="WEST") then
                if (this.boundaryLocks[k].lockOpposite) then
                    this.xmin=theFeature.xmin+this.boundaryLocks[k].offset;
                else
                    this.xmin=theFeature.xmax+this.boundaryLocks[k].offset;
                end
                if (this.fixedWidth) then this.xmax=this.xmin+width; end
            end
 
            this.SizeChanged(this);
        end
    end
    
end

--#########################################################
function PanelClass.Hide(this)
    this.hidden=true;
    chiObjectSetProperty(this.obj1Num,"Hidden",true);
    chi3DLineChangeColor(this.lineNum,0.0,0.0,0.0,0.0);
    chi3DLineChangeColor(this.scrollBorder,0.0,0.0,0.0,0.0);
    --this.Redraw(this);
end

--#########################################################
function PanelClass.UnHide(this)
    this.hidden=false;
    chiObjectSetProperty(this.obj1Num,"Hidden",false);

    if (this.showOutline) then
        chi3DLineChangeColor(this.lineNum,0.0,0.0,0.0,1.0);
    else
        chi3DLineChangeColor(this.lineNum,0.0,0.0,0.0,0.0);
    end
    chi3DLineChangeColor(this.lineNum,0.0,0.0,0.0,1.0);
    if (this.scrollActive) then
        chi3DLineChangeColor(this.scrollBorder,0.0,0.0,0.0,1.0);
    end
    --this.Redraw(this);
end

--######################################################### Add Scrollbar
function PanelClass.AddScrollbar(this,width,offset,active)
    if (width == nil) then
        this.scrollWidth = 4.0;
    else
        this.scrollWidth = width;
    end

    if (offset == nil) then
        this.scrollOffset = 0.0;
    else
        this.scrollOffset = offset;
    end



    if (active == false) then
        chi3DLineChangeColor(this.scrollBorder,0.0,0.0,0.0,0.0);
        chiObjectSetProperty(this.objSNum,"Hidden",true);
        chiObjectSetProperty(this.objS2Num,"Hidden",true);
    else
        chi3DLineChangeColor(this.scrollBorder,0.0,0.0,0.7,1.0);
        chiObjectSetProperty(this.objSNum,"Hidden",false);
        chiObjectSetProperty(this.objS2Num,"Hidden",false);
    end
end

--######################################################### Disable Scrollbar
function PanelClass.DisableScrollbar(this,width,offset,active)
    if (width == nil) then
        this.scrollWidth = 4.0;
    else
        this.scrollWidth = width;
    end

    if (offset == nil) then
        this.scrollOffset = 0.0;
    else
        this.scrollOffset = offset;
    end

    if (active == false) then
        chi3DLineChangeColor(this.scrollBorder,0.0,0.0,0.0,0.0);
    else
        chi3DLineChangeColor(this.scrollBorder,0.0,0.0,0.7,1.0);

    end
end

function PanelClass.Click(this)
end

function PanelClass.MouseMove(this)
end

function PanelClass.MouseEnter(this)
end

function PanelClass.MouseLeave(this)
end

function PanelClass.MouseHover(this)
end

function PanelClass.KeyPress(this)
end

function PanelClass.KeyDn(this)
end     

function PanelClass.KeyUp(this)
end


function PanelClass.WindowSizeChanged(this)
    this.SizeChanged(this)
end

function PanelClass.SizeChanged(this)
    curScene=chiGetScene();
    
    if (this.window>=0) then
        parentWindowx,parentWindowy=chiGetWindowProperties();
    
        chiBindScene(this.window);
        
        chiSetWindowProperties(this.xmax-this.xmin-1,this.ymax-this.ymin-1,this.xmin,parentWindowy-this.ymax+1);
    end
    
    chiBindScene(curScene);
    
    this.cursorX=this.xmin;
    
    this.cursorY=this.ymax+this.scrollOffset*this.scrollScale;
    
    newxmax=0;
    this.cursorHeight=0;
    this.contentYSize=0;
    for k=1,this.slaveCount do
        --Check overflow
        if (this.slaves[k].xSize~=nil) then
            newxmax=this.cursorX+this.slaves[k].xSize;
            if (newxmax>this.xmax) then
                this.cursorX=this.xmin;
                this.cursorY=this.cursorY-this.cursorHeight;
                this.cursorHeight=0;
            end
        end
        
        --Call slave's event
        this.slaves[k].SizeChanged(this.slaves[k]);
        
        --Calculate new cursor
        if ((this.slaves[k].xSize~=nil) and (this.slaves[k].ySize~=nil)) then
            this.cursorX=this.cursorX+this.slaves[k].xSize+this.slaves[k].paddingLeft+this.slaves[k].paddingRight;
            if (this.slaves[k].ySize>this.cursorHeight) then
                this.cursorHeight=this.slaves[k].ySize+this.slaves[k].paddingBot+this.slaves[k].paddingTop;
                this.contentYSize=this.contentYSize+this.cursorHeight;
            end
        end
        
    end
    this.cursorY=this.cursorY-this.cursorHeight;
    
    --================================= Update slave boundary locks
    for lk=1,this.bslaveCount do
        this.bslaves[lk].UpdateBoundaryLocks(this.bslaves[lk]);
    end

    this.Redraw(this)
end
