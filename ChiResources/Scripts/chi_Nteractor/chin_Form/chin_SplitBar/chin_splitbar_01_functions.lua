

--######################################################### Events
function SplitBarClass.ProcessEvents(this)
    --======================= Size changed
    if (WM_SIZE.occured) then
        this.WindowSizeChanged(this);
         --print("MUFF")
    end    
    
    --======================= Left click
    if (WM_LBUTTONDOWN.occured) then
        if     (WM_LBUTTONDOWN.iPar5==this.obj1Num) then
            this.sizeChanging=true;
        end
    end
    if (WM_LBUTTONUP.occured) then
        this.sizeChanging=false;
    end
    
    --======================= Mouse move
    if (WM_MOUSEMOVE.occured) then
        if (not this.sizeChanging) then
            if     (WM_MOUSEMOVE.iPar5==this.obj1Num) then
                if (this.vertical) then
                    chiWindowSetCursor(1);
                else
                    chiWindowSetCursor(2);
                end
            end
        else
            --Default window
            if (WM_MOUSEMOVE.iPar4==0) then
                if (this.vertical) then
                    oldxmin=this.xmin;
                    oldxmax=this.xmax;
                    this.xmin=WM_MOUSEMOVE.iPar0-0.5*this.gripWidth;
                    this.xmax=WM_MOUSEMOVE.iPar0+0.5*this.gripWidth;
                    
                    if (not this.QuerySlaves(this)) then
                        this.xmin=oldxmin;
                        this.xmax=oldxmax;
                    else
                        this.SizeChanged(this);--
                        this.Redraw(this); 
                    end  
                else
                    oldymin=this.ymin;
                    oldymax=this.ymax;
                    this.ymin=chinGlobal.dwindowysize-WM_MOUSEMOVE.iPar1-0.5*this.gripWidth;
                    this.ymax=chinGlobal.dwindowysize-WM_MOUSEMOVE.iPar1+0.5*this.gripWidth;
                    
                    if (not this.QuerySlaves(this)) then
                        this.ymin=oldymin;
                        this.ymax=oldymax;
                    else
                        this.SizeChanged(this);--
                        this.Redraw(this); 
                    end 
                end
            --Other windows
            else
                if (this.vertical) then
                    oldxmin=this.xmin;
                    oldxmax=this.xmax;
                    this.xmin=this.xmin+WM_MOUSEMOVE.iPar0-0.5*this.gripWidth;
                    this.xmax=this.xmin+this.gripWidth;
                    
                    if (not this.QuerySlaves(this)) then
                        this.xmin=oldxmin
                        this.xmax=oldxmax
                    else
                        this.SizeChanged(this);--
                        this.Redraw(this);
                    end
                else
                    curScene=chiGetScene();
                        chiBindScene(WM_MOUSEMOVE.iPar4);
                        parentWindowx,parentWindowy=chiGetWindowProperties();
                        --print(parentWindowx,parentWindowy,WM_MOUSEMOVE.iPar1)
                    chiBindScene(curScene);
                
                    oldymin=this.ymin;
                    oldymax=this.ymax;
                    this.ymin=this.ymin+parentWindowy-WM_MOUSEMOVE.iPar1-0.5*this.gripWidth;
                    this.ymax=this.ymin+this.gripWidth;
                    
                    if (not this.QuerySlaves(this)) then
                        this.ymin=oldymin
                        this.ymax=oldymax
                    else
                        this.SizeChanged(this);--
                        this.Redraw(this);
                    end
                end
            end
        end
    end
   
    for k=1,this.eventCallbackCount do
        this.eventCallbacks[k](this);
    end
end

--######################################################### Redraw 
function SplitBarClass.Redraw(this)

    this.QuerySlaves(this);
    if (    this.vertical) then this.xmax=this.xmin+this.gripWidth; end
    if (not this.vertical) then this.ymax=this.ymin+this.gripWidth; end
    girth =this.xmax-this.xmin-1;
    length=this.ymax-this.ymin-1;
    
    
    chiTransformSetScale(this.transform1,girth,length,1); 
    chiTransformSetTranslation(this.transform1,this.xmin,this.ymin,1.2); --E
    
    this.UpdateSlaves(this);
    
    
end

--######################################################### AddSlave
function SplitBarClass.AddSlave(this,newSlave,typeOfSlave,offset)
    this.slaveCount=this.slaveCount+1;
    this.slaves[this.slaveCount]={}
    this.slaves[this.slaveCount].feature=newSlave
    this.slaves[this.slaveCount].typeOfSlave=typeOfSlave;
    this.slaves[this.slaveCount].offset=offset;
    
    this.QuerySlaves(this);
    this.Redraw(this);
end

--######################################################### Update Slaves
function SplitBarClass.UpdateSlaves(this)
    success=true;
    for k=1,this.slaveCount do
        theFeature=this.slaves[k].feature;
        if (this.vertical) then
            success=theFeature.RequestDimension(theFeature,this.slaves[k].typeOfSlave,this.xmin+this.slaves[k].offset)
        else
            success=theFeature.RequestDimension(theFeature,this.slaves[k].typeOfSlave,this.ymin+this.slaves[k].offset)
        end
        if (not success) then
            return false;
        end
    end
    
    
    return true;
end

--######################################################### Query slaves
function SplitBarClass.QuerySlaves(this)
    success=true;
    minx=10000;
    maxx=0;
    miny=10000;
    maxy=0;
    
    
    
    for k=1,this.slaveCount do
        theFeature=this.slaves[k].feature;
        if (this.vertical) then
            success=theFeature.QueryDimension(theFeature,this.slaves[k].typeOfSlave,this.xmin+this.slaves[k].offset)
            
            if (theFeature.ymin < miny) then miny=theFeature.ymin-1; end
            if (theFeature.ymax > maxy) then maxy=theFeature.ymax+1; end

        else
            success=theFeature.QueryDimension(theFeature,this.slaves[k].typeOfSlave,this.ymin+this.slaves[k].offset)
            
            if (theFeature.xmin < minx) then minx=theFeature.xmin-1; end
            if (theFeature.xmax > maxx) then maxx=theFeature.xmax+1; end
        end
        if (not success) then
            return false;
        end
    end
    
    if (this.vertical) then
        this.ymin=miny;
        this.ymax=maxy;
    else    
        this.xmin=minx;
        this.xmax=maxx;
    end
    
    return true;
end



function SplitBarClass.Click(this)
end

function SplitBarClass.MouseMove(this)
end

function SplitBarClass.MouseEnter(this)
end

function SplitBarClass.MouseLeave(this)
end

function SplitBarClass.MouseHover(this)
end

function SplitBarClass.KeyPress(this)
print("PRESS")
end

function SplitBarClass.KeyDn(this)
end     

function SplitBarClass.KeyUp(this)

end

function SplitBarClass.WindowSizeChanged(this)
    this.Redraw(this)
end

function SplitBarClass.SizeChanged(this)
    
    for k=1,this.slaveCount do
        theFeature=this.slaves[k].feature;
        theFeature.SizeChanged(theFeature);
       
    end
end