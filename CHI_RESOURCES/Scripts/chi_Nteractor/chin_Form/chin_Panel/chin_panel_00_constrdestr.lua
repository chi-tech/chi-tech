--######################################################### Constructor
function PanelClass.New(name)
    local this = setmetatable({},PanelClass)
    PanelCount = PanelCount+1;

    -- ===== Panel Variables =====
    this.name               = name;
    this.xmin               = 100;
    this.ymin               = 100;
    this.xmax               = 200;
    this.ymax               = 200;
    this.cursorX            = 0;
    this.cursorY            = 0;
    this.cursorHeight       = 0;
    this.widthMin           = 100;
    this.lengthMin          = 100;
    this.gripWidth          = 10;
    this.sizeChanging       = false;
    this.currentGripper     = 0;
    this.showOutline        = true;
    this.contentYSize       = 0;

    this.boundaryLocks      = {}
    this.boundaryLockCount  = 0;
    this.fixedHeight        = false;
    this.fixedWidth         = false;
    this.ownResize          = false;
    this.zDepth             = 1.0;
    this.hidden             = false;

    this.windowAttached     = false;
    this.window             = -1;


    this.cutWindow          = -1;
    this.cutPanel           = nil;

    -- =====  Slave Variables =====
    this.slaves             = {}
    this.slaveCount         = 0;
    
    this.bslaves            = {}
    this.bslaveCount        = 0;
    
    this.eventCallbacks     = {}
    this.eventCallbackCount = 0;


    this.scrollActive       = false;
    this.scrollWidth        = 0.0;
    this.scrollOffset       = 0.0;
    this.scrollDepthTest    = 1.1;
    this.scrollHeight       = 50;
    this.scrollyi           = 0;
    this.scrollyf           = 0;
    this.scrollScale        = 1.0;
    this.scrollReference    = 0;


   
    
    
    objName = "Panel"..string.format("%02d",PanelCount);
    

    --=========================================== Creating Handles
    this.obj2Num = chiObjectCreate(objName.."_Handle1"); --E
    this.obj3Num = chiObjectCreate(objName.."_Handle2"); --W
    this.obj4Num = chiObjectCreate(objName.."_Handle3"); --N
    this.obj5Num = chiObjectCreate(objName.."_Handle4"); --S
    
    chiObjectAddSurface(objName.."_Handle1",PanelSurface); --E
    chiObjectAddSurface(objName.."_Handle2",PanelSurface); --W
    chiObjectAddSurface(objName.."_Handle3",PanelSurface); --N
    chiObjectAddSurface(objName.."_Handle4",PanelSurface); --S
    
    this.transform2 = chiTransformCreate(objName.."_Handle1".."_Transform"); --E
    this.transform3 = chiTransformCreate(objName.."_Handle2".."_Transform"); --W
    this.transform4 = chiTransformCreate(objName.."_Handle3".."_Transform"); --N
    this.transform5 = chiTransformCreate(objName.."_Handle4".."_Transform"); --S
    
    chiObjectSetProperty(this.obj2Num,"Transform",objName.."_Handle1".."_Transform"); --E
    chiObjectSetProperty(this.obj3Num,"Transform",objName.."_Handle2".."_Transform"); --W
    chiObjectSetProperty(this.obj4Num,"Transform",objName.."_Handle3".."_Transform"); --N
    chiObjectSetProperty(this.obj5Num,"Transform",objName.."_Handle4".."_Transform"); --S
    
    chiObjectSetProperty(this.obj2Num,"Renderable",false); --E
    chiObjectSetProperty(this.obj3Num,"Renderable",false); --W
    chiObjectSetProperty(this.obj4Num,"Renderable",false); --N
    chiObjectSetProperty(this.obj5Num,"Renderable",false); --S

    chiObjectSetProperty(this.obj2Num,"Hidden",true); --E
    chiObjectSetProperty(this.obj3Num,"Hidden",true); --W
    chiObjectSetProperty(this.obj4Num,"Hidden",true); --N
    chiObjectSetProperty(this.obj5Num,"Hidden",true); --S
    
    chiObjectSetProperty(this.obj2Num,"donotList",true);
    chiObjectSetProperty(this.obj3Num,"donotList",true);
    chiObjectSetProperty(this.obj4Num,"donotList",true);
    chiObjectSetProperty(this.obj5Num,"donotList",true);

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
    
    --=========================================== Creating base object
    this.obj1Num = chiObjectCreate(objName);
    chiObjectAddSurface(this.obj1Num,PanelSurface);
    
    this.transform1 = chiTransformCreate(objName.."_Transform");
    chiObjectSetProperty(this.obj1Num,"Transform",objName.."_Transform");
    chiTransformSetScale(this.transform1,this.xmax-this.xmin-1,this.ymax-this.ymin-1,1.0)
    chiTransformSetTranslation(this.transform1,this.xmin,this.ymin,this.zDepth)
    
    this.matlNum = chiMaterialCreate(objName .. "_Material");
    chiMaterialSetProperty(this.matlNum,"DisableShading",true);
    ambient = 1.0;
    chiMaterialSetProperty(this.matlNum,"Diffuse",ambient,ambient,ambient,1.0);
    chiObjectSetProperty(this.obj1Num,"Material",objName .. "_Material");
    chiObjectSetProperty(this.obj1Num,"donotList",true);
    
    --=========================================== Creating outline    
    this.lineNum = chi3DLineCreate("Yes");
    chi3DLineAddVertex(this.lineNum,this.xmin,this.ymin,this.zDepth);
    chi3DLineAddVertex(this.lineNum,this.xmin,this.ymax,this.zDepth);
    chi3DLineAddVertex(this.lineNum,this.xmax,this.ymax,this.zDepth);
    chi3DLineAddVertex(this.lineNum,this.xmax,this.ymin,this.zDepth);
    chi3DLineAddVertex(this.lineNum,this.xmin-1,this.ymin,this.zDepth);

    if (this.showOutline) then
        chi3DLineChangeColor(this.lineNum,0.0,0.0,0.0,1.0);
    else
        chi3DLineChangeColor(this.lineNum,0.0,0.0,0.0,0.0);
    end

    --=========================================== Creating Scrollbar
    -- Scrollbar Border
    this.scrollBorder = chi3DLineCreate(1);
    chi3DLineAddVertex(this.scrollBorder, this.xmax-this.scrollWidth, this.ymin, this.zDepth+this.scrollDepthTest); -- SW
    chi3DLineAddVertex(this.scrollBorder, this.xmax-this.scrollWidth, this.ymax, this.zDepth+this.scrollDepthTest); -- NW
    chi3DLineAddVertex(this.scrollBorder, this.xmax                 , this.ymax, this.zDepth+this.scrollDepthTest); -- NE
    chi3DLineAddVertex(this.scrollBorder, this.xmax                 , this.ymin, this.zDepth+this.scrollDepthTest); -- SE
    chi3DLineAddVertex(this.scrollBorder, this.xmax-this.scrollWidth, this.ymin, this.zDepth+this.scrollDepthTest); -- ??
    chi3DLineChangeColor(this.scrollBorder,0.0,0.0,0.7,0.0);                                                        -- Temp Blue

    --=========================================== Creating scroll Icon
    this.objSNum = chiObjectCreate(objName.."Scroll Icon");
    chiObjectAddSurface(this.objSNum,PanelSurface);

    this.transformS = chiTransformCreate(objName.."Scroll Icon".."_Transform");
    chiObjectSetProperty(this.objSNum,"Transform",objName.."Scroll Icon".."_Transform");
    chiTransformSetScale(this.transformS,this.scrollWidth,this.scrollHeight,1.0)
    chiTransformSetTranslation(this.transformS,this.xmax-this.scrollWidth,this.ymax-this.scrollHeight,this.zDepth+3.0)

    this.matSNum = chiMaterialCreate(objName.."Scroll Icon".. "_Material");
    chiMaterialSetProperty(this.matSNum,"DisableShading",true);

    ambient = 0.8;
    chiMaterialSetProperty(this.matSNum,"Diffuse",ambient,ambient,ambient,1.0);
    chiObjectSetProperty(this.objSNum,"Material",objName.."Scroll Icon".. "_Material");
    chiObjectSetProperty(this.objSNum,"Hidden",true);

    --=========================================== Creating scroll Base
    this.objS2Num = chiObjectCreate(objName.."Scroll Base");
    chiObjectAddSurface(this.objS2Num,PanelSurface);

    this.transformS2 = chiTransformCreate(objName.."Scroll Base".."_Transform");
    chiObjectSetProperty(this.objS2Num,"Transform",objName.."Scroll Base".."_Transform");
    chiTransformSetScale(this.transformS2,this.scrollWidth-1,this.ymax-this.ymin-1,1.0)
    chiTransformSetTranslation(this.transformS2,this.xmax-this.scrollWidth,this.ymin,this.zDepth+2.9)

    this.matS2Num = chiMaterialCreate(objName.."Scroll Base".. "_Material");
    chiMaterialSetProperty(this.matS2Num,"DisableShading",true);

    ambient = 0.97;
    chiMaterialSetProperty(this.matS2Num,"Diffuse",ambient,ambient,ambient,1.0);
    chiObjectSetProperty(this.objS2Num,"Material",objName.."Scroll Base".. "_Material");
    chiObjectSetProperty(this.objS2Num,"Hidden",true);

    return this;
end


