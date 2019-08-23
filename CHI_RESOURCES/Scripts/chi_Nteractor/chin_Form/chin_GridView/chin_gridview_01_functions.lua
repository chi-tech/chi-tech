--#########################################################
function GridViewClass.ProcessEvents(this)
    for k=1,this.eventCallbackCount do
        this.eventCallbacks[k](this);
    end
end


--#########################################################
function GridViewClass.SetProperty(this,property,value)
    if     (property=="Master") then
           this.master=value;
           
           value.slaveCount=value.slaveCount+1;
           k=value.slaveCount
           value.slaves[k]=this;
           this.SizeChanged(this)
    end
end

--#########################################################
function GridViewClass.Redraw(this)
    chi3DLineChangeVertex(this.outline,0,this.xpos           ,this.ypos           ,this.zDepth);
    chi3DLineChangeVertex(this.outline,1,this.xpos+this.xSize,this.ypos           ,this.zDepth);
    chi3DLineChangeVertex(this.outline,2,this.xpos+this.xSize,this.ypos-this.ySize,this.zDepth);
    chi3DLineChangeVertex(this.outline,3,this.xpos           ,this.ypos-this.ySize,this.zDepth);
    chi3DLineChangeVertex(this.outline,4,this.xpos           ,this.ypos           ,this.zDepth);
    
    for k=1,(this.numberOfColumns-1) do
        deltaX=this.xSize/this.numberOfColumns;
        chi3DLineChangeVertex(this.hdivider[k].line,0,this.xpos+k*deltaX,this.ypos           ,this.zDepth);
        chi3DLineChangeVertex(this.hdivider[k].line,1,this.xpos+k*deltaX,this.ypos-this.ySize,this.zDepth);
    end
    
    for k=1,(this.numberOfRows-1) do        
        deltaY=this.ySize/this.numberOfRows;
        chi3DLineChangeVertex(this.vdivider[k].line,0,this.xpos           ,this.ypos-k*deltaY,this.zDepth);
        chi3DLineChangeVertex(this.vdivider[k].line,1,this.xpos+this.xSize,this.ypos-k*deltaY,this.zDepth);
    end
end

--#########################################################
function GridViewClass.SizeChanged(this)
    if (this.master~=nil) then
        this.xpos=this.master.cursorX;
        this.ypos=this.master.cursorY;
        this.xSize=this.master.xmax-this.master.xmin;
        if (this.master.scrollWidth~=nil) then
            this.xSize=this.xSize-this.master.scrollWidth-1;
        end
        
        chiViewportSetProperty(this.viewportNum,this.master.xmin,this.master.ymin,this.master.xmax,this.master.ymax);

        deltaX=this.xSize/this.numberOfColumns;
        deltaY=this.ySize/this.numberOfRows;
        for y=1,this.numberOfRows do
            for x=1,this.numberOfColumns do
                for k=1,this.gridReference[y][x].slaveCount do
                    this.gridReference[y][x].cursorX=this.xpos+deltaX*(x-1);
                    this.gridReference[y][x].cursorY=this.ypos-deltaY*(y-1);
                    this.gridReference[y][x].xmin=this.xpos+deltaX*(x-1);
                    this.gridReference[y][x].ymax=this.ypos-deltaY*(y-1);
                    this.gridReference[y][x].xmax=this.xpos+deltaX*(x);
                    this.gridReference[y][x].ymin=this.ypos-deltaY*(y);
                    --Overflow
                    if (this.gridReference[y][x].xmax>this.master.xmax) then
                        this.gridReference[y][x].xmax=this.master.xmax;
                    end
                    if (this.gridReference[y][x].xmin<this.master.xmin) then
                        this.gridReference[y][x].xmin=this.master.xmin;
                    end
                    if (this.gridReference[y][x].ymax>this.master.ymax) then
                        this.gridReference[y][x].ymax=this.master.ymax;
                    end
                    if (this.gridReference[y][x].ymin<this.master.ymin) then
                        this.gridReference[y][x].ymin=this.master.ymin;
                    end
                    this.gridReference[y][x].slaves[k].SizeChanged(this.gridReference[y][x].slaves[k]);
                end
            end
        end
        

    end
    
    this.Redraw(this);
end
