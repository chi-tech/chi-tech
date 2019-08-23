--############################################################################# Constructor
function GridViewClass.New(name)
    local this=setmetatable({},GridViewClass)
    GridViewCount=GridViewCount+1;
    
    objName="GridView"..string.format("%03d",GridViewCount);
    this.name=name;
    this.numberOfColumns=2;
    this.numberOfRows=50;
    this.xpos=100;
    this.ypos=100;
    this.xSize=200;
    this.ySize=this.numberOfRows*20;
    this.zDepth=1.1;
    
    this.eventCallbacks     = {}
    this.eventCallbackCount = 0;
    
    this.paddingTop=0;
    this.paddingBot=0;
    this.paddingLeft=0;
    this.paddingRight=0;
    
    --=========================================== Viewport
    this.viewportNum=chiViewportCreate(objName.."Viewport");
    chiViewportSetProperty(this.viewportNum,0,0,4000,4000);
    
    --=========================================== Outline
    this.outline=chi3DLineCreate(objName.."Outline");
    chi3DLineChangeColor(this.outline,0.6,0.6,0.6,1.0);
    chi3DLineSetviewport(this.outline,this.viewportNum);
    
    chi3DLineAddVertex(this.outline,this.xpos           ,this.ypos           ,this.zDepth);
    chi3DLineAddVertex(this.outline,this.xpos+this.xSize,this.ypos           ,this.zDepth);
    chi3DLineAddVertex(this.outline,this.xpos+this.xSize,this.ypos-this.ySize,this.zDepth);
    chi3DLineAddVertex(this.outline,this.xpos           ,this.ypos-this.ySize,this.zDepth);
    chi3DLineAddVertex(this.outline,this.xpos           ,this.ypos           ,this.zDepth);
    
    --=========================================== Horizontal Dividers
    this.hdivider={};
    for k=1,(this.numberOfColumns-1) do
        this.hdivider[k]={}
        this.hdivider[k].line=chi3DLineCreate(objName.."hdivider"..string.format("%02d",k).."Line");
        chi3DLineChangeColor(this.hdivider[k].line,0.6,0.6,0.6,1.0);
        chi3DLineSetviewport(this.hdivider[k].line,this.viewportNum);
        
        deltaX=this.xSize/this.numberOfColumns;
        chi3DLineAddVertex(this.hdivider[k].line,this.xpos+k*deltaX,this.ypos           ,this.zDepth);
        chi3DLineAddVertex(this.hdivider[k].line,this.xpos+k*deltaX,this.ypos-this.ySize,this.zDepth);
    end
    
    --=========================================== Vertical Dividers
    this.vdivider={};
    for k=1,(this.numberOfRows-1) do
        this.vdivider[k]={};
        this.vdivider[k].line=chi3DLineCreate(objName.."vdivider"..string.format("%02d",k).."Line");
        chi3DLineChangeColor(this.vdivider[k].line,0.6,0.6,0.6,1.0);
        chi3DLineSetviewport(this.vdivider[k].line,this.viewportNum);
        
        deltaY=this.ySize/this.numberOfRows;
        chi3DLineAddVertex(this.vdivider[k].line,this.xpos           ,this.ypos-k*deltaY,this.zDepth);
        chi3DLineAddVertex(this.vdivider[k].line,this.xpos+this.xSize,this.ypos-k*deltaY,this.zDepth);
    end
    
    --=========================================== GridPositions
    this.gridReference={}
    deltaX=this.xSize/this.numberOfColumns;
    deltaY=this.ySize/this.numberOfRows;
    for y=1,this.numberOfRows do
        this.gridReference[y]={};
        for x=1,this.numberOfColumns do
            this.gridReference[y][x]={};
            this.gridReference[y][x].cursorX=this.xpos+deltaX*(x-1);
            this.gridReference[y][x].cursorY=this.ypos-deltaY*(y-1);
            this.gridReference[y][x].xmin=this.xpos+deltaX*(x-1);
            this.gridReference[y][x].ymax=this.ypos-deltaY*(y-1);
            this.gridReference[y][x].xmax=this.xpos+deltaX*(x);
            this.gridReference[y][x].ymin=this.ypos-deltaY*(y);
            this.gridReference[y][x].slaves={}
            this.gridReference[y][x].slaveCount=0;
        end
    end
    
    return this;
end