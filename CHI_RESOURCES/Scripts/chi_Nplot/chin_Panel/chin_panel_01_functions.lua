--What are panels? Panels are essentially viewports. Stuff drawn within panels 
--are clipped on the edges. Panels can be free floating or bound by edges.
--Panels should have different edge styles and should always have a resize 
--handler.
--Panels can also be constrained to a minimum size.
--For connected panels, code should be implemented that dynamically moves all
--panels that are attached together.

PanelClass={}
PanelClass.__index = PanelClass
PanelCount=0;

PanelSurface=chiLoadSurface(dBaseDir .. "Panel.obj");
PanelTexture=chiLoadTexture("CHI_RESOURCES/Textures/LogoTest.tga");

--######################################################### Constructor
function PanelClass.New(name)
    local this=setmetatable({},PanelClass)
    PanelCount=PanelCount+1;
    
    this.name=name;
    this.xmin=100;
    this.ymin=100;
    this.xmax=200;
    this.ymax=200;
    this.gripWidth=10;
    this.sizeChanging=false;
    this.currentGripper=0;
    
    
    objName="Panel"..string.format("%02d",PanelCount);
    

    
    --=========================================== Creating Handles
    this.obj2Num=chiObjectCreate(objName.."_Handle1"); --E
    this.obj3Num=chiObjectCreate(objName.."_Handle2"); --W
    this.obj4Num=chiObjectCreate(objName.."_Handle3"); --N
    this.obj5Num=chiObjectCreate(objName.."_Handle4"); --S
    
    chiObjectAddSurface(objName.."_Handle1",PanelSurface); --E
    chiObjectAddSurface(objName.."_Handle2",PanelSurface); --W
    chiObjectAddSurface(objName.."_Handle3",PanelSurface); --N
    chiObjectAddSurface(objName.."_Handle4",PanelSurface); --S
    
    this.transform2=chiTransformCreate(objName.."_Handle1".."_Transform"); --E
    this.transform3=chiTransformCreate(objName.."_Handle2".."_Transform"); --W
    this.transform4=chiTransformCreate(objName.."_Handle3".."_Transform"); --N
    this.transform5=chiTransformCreate(objName.."_Handle4".."_Transform"); --S
    
    chiObjectSetProperty(this.obj2Num,"Transform",objName.."_Handle1".."_Transform"); --E
    chiObjectSetProperty(this.obj3Num,"Transform",objName.."_Handle2".."_Transform"); --W
    chiObjectSetProperty(this.obj4Num,"Transform",objName.."_Handle3".."_Transform"); --N
    chiObjectSetProperty(this.obj5Num,"Transform",objName.."_Handle4".."_Transform"); --S
    
    chiObjectSetProperty(this.obj2Num,"Renderable",false); --E
    chiObjectSetProperty(this.obj3Num,"Renderable",false); --W
    chiObjectSetProperty(this.obj4Num,"Renderable",false); --N
    chiObjectSetProperty(this.obj5Num,"Renderable",false); --S
    
    girth =this.xmax-this.xmin-1;
    length=this.ymax-this.ymin-1;
    chiTransformSetScale(this.transform2,this.gripWidth,length,1); --E
    chiTransformSetScale(this.transform3,this.gripWidth,length,1); --W
    chiTransformSetScale(this.transform4,girth,this.gripWidth,1); --N
    chiTransformSetScale(this.transform5,girth,this.gripWidth,1); --S
    
    chiTransformSetTranslation(this.transform2,this.xmax                 -0.5*this.gripWidth,this.ymin,1.1); --E
    chiTransformSetTranslation(this.transform3,this.xmin-this.gripWidth-1+0.5*this.gripWidth,this.ymin,1.1); --W
    chiTransformSetTranslation(this.transform4,this.xmin                 ,this.ymax-0.5*this.gripWidth,1.1); --N
    chiTransformSetTranslation(this.transform5,this.xmin,this.ymin-this.gripWidth-1+0.5*this.gripWidth,1.1); --S
    
    --=========================================== Creating base object
    this.obj1Num=chiObjectCreate(objName);
    chiObjectAddSurface(this.obj1Num,PanelSurface);
    
    this.transform1=chiTransformCreate(objName.."_Transform");
    chiObjectSetProperty(this.obj1Num,"Transform",objName.."_Transform");
    chiTransformSetScale(this.transform1,this.xmax-this.xmin-1,this.ymax-this.ymin-1,1.5)
    chiTransformSetTranslation(this.transform1,this.xmin,this.ymin,1.0)
    
    this.matlNum=chiMaterialCreate(objName .. "_Material");
    --chiMaterialSetProperty(this.matlNum,"AmbientTexture","CHI_RESOURCES/Textures/LogoTest.tga");
    --chiMaterialSetProperty(this.matlNum,"AmbientTextureEnabled",true);
    ambient=0.8;
    chiMaterialSetProperty(this.matlNum,"Ambient",ambient,ambient,ambient,1.0);
    chiObjectSetProperty(this.obj1Num,"Material",objName .. "_Material");
    
    --=========================================== Creating outline    
    this.lineNum=chiLineCreate("Yes");
    print("LineCreated",this.lineNum)
    chiLineAddVertex(this.lineNum,this.xmin,this.ymin,0.0);
    chiLineAddVertex(this.lineNum,this.xmin,this.ymax,0.0);
    chiLineAddVertex(this.lineNum,this.xmax,this.ymax,0.0);
    chiLineAddVertex(this.lineNum,this.xmax,this.ymin,0.0);
    chiLineAddVertex(this.lineNum,this.xmin-1,this.ymin,0.0);
    chiLineChangeColor(this.lineNum,0.0,0.0,0.0,1.0);
    
    return this;
end

--######################################################### Set bounds
function PanelClass.SetBounds(xmin,ymin,xmax,ymax)
    --Update bounding line
    --If theres a window attached update it
end



function PanelClass.Redraw(this)
    chiTransformSetScale(this.transform1,this.xmax-this.xmin-1,this.ymax-this.ymin-1,1.5)
    chiTransformSetTranslation(this.transform1,this.xmin,this.ymin,1.0)
    
    girth =this.xmax-this.xmin-1;
    length=this.ymax-this.ymin-1;
    chiTransformSetScale(this.transform2,this.gripWidth,length,1); --E
    chiTransformSetScale(this.transform3,this.gripWidth,length,1); --W
    chiTransformSetScale(this.transform4,girth,this.gripWidth,1); --N
    chiTransformSetScale(this.transform5,girth,this.gripWidth,1); --S

    chiTransformSetTranslation(this.transform2,this.xmax                 -0.5*this.gripWidth,this.ymin,1.1); --E
    chiTransformSetTranslation(this.transform3,this.xmin-this.gripWidth-1+0.5*this.gripWidth,this.ymin,1.1); --W
    chiTransformSetTranslation(this.transform4,this.xmin                 ,this.ymax-0.5*this.gripWidth,1.1); --N
    chiTransformSetTranslation(this.transform5,this.xmin,this.ymin-this.gripWidth-1+0.5*this.gripWidth,1.1); --S

    chiLineChangeVertex(this.lineNum,0,this.xmin,this.ymin,0.0);
    chiLineChangeVertex(this.lineNum,1,this.xmin,this.ymax,0.0);
    chiLineChangeVertex(this.lineNum,2,this.xmax,this.ymax,0.0);
    chiLineChangeVertex(this.lineNum,3,this.xmax,this.ymin,0.0);
    chiLineChangeVertex(this.lineNum,4,this.xmin-1,this.ymin,0.0);
end



--######################################################### Events
function PanelClass.ProcessEvents(this)
    --======================= Size changed
    if (WM_SIZE.occured) then
        --print("WindowSizeChanged()");
        this.WindowSizeChanged();
    end
    
    --======================= Mouse Move
    if (WM_MOUSEMOVE.occured) then
        if (not this.sizeChanging) then
            if     ((WM_MOUSEMOVE.iPar5==this.obj2Num) or (WM_MOUSEMOVE.iPar5==this.obj3Num)) then
                chiWindowSetCursor(1);
            elseif ((WM_MOUSEMOVE.iPar5==this.obj4Num) or (WM_MOUSEMOVE.iPar5==this.obj5Num)) then
                chiWindowSetCursor(2);
            else
                chiWindowSetCursor(0);
            end
        else
            if     (this.currentGripper==2) then
                this.xmax=WM_MOUSEMOVE.iPar0 --- WM_MOUSEMOVE.iPar2;
                this.Redraw(this);                     --
            elseif (this.currentGripper==3) then       --
                this.xmin=WM_MOUSEMOVE.iPar0 --- WM_MOUSEMOVE.iPar2;
                this.Redraw(this);                     --
            elseif (this.currentGripper==4) then       --
                this.ymax=chinGlobal.dwindowysize-WM_MOUSEMOVE.iPar1 --+ WM_MOUSEMOVE.iPar3;
                this.Redraw(this);                     --
            elseif (this.currentGripper==5) then       --
                this.ymin=chinGlobal.dwindowysize-WM_MOUSEMOVE.iPar1 --+ WM_MOUSEMOVE.iPar3;
                this.Redraw(this);
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
    end
    if (WM_LBUTTONUP.occured) then
        this.currentGripper=0;
        this.sizeChanging=false;
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
print("PRESS")
end

function PanelClass.KeyDn(this)
end     

function PanelClass.KeyUp(this)
print("UP")
end
print("Panel code loaded")

function PanelClass.WindowSizeChanged(this)
print("SizeChanged")
end
