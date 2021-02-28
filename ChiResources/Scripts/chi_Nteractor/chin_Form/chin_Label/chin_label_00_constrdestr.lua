--######################################################### Constructor
function LabelClass.New(name)
    local this=setmetatable({},LabelClass)
    LabelCount=LabelCount+1;
    
    local objName="Label"..string.format("%03d",LabelCount);
    this.name=name;
    this.xpos=100;
    this.ypos=100;
    this.xSize=100;
    this.ySize=20;
    this.zDepth=1.1;
    this.paddingTop=0;
    this.paddingBot=0;
    this.paddingLeft=0;
    this.paddingRight=0;
    this.float=true;
    this.selectable=true;
    this.selected=false;
    this.highlightFillsMaster=false;
    this.hidden=false;
    
    this.eventCallbacks     = {}
    this.eventCallbackCount = 0;
    
    --Label
    this.fontType=2;
    this.fontColor={0,0,0,1};
    this.text=name;
    r=this.fontColor[1];
    g=this.fontColor[2];
    b=this.fontColor[3];
    a=this.fontColor[4];
    this.textNum=chiSetLabel3D(objName,this.text,this.xpos+this.paddingLeft,this.ypos+this.paddingBot,r,g,b,a,this.fontType)
    chiSetLabelProperty3D(this.textNum,"Depth",this.zDepth+0.01);
    this.master=nil;
    this.parent=nil
    
    --Selection box
    subObjName=objName.."_selectionBox";
    this.obj1Num=chiObjectCreate(subObjName);
    chiObjectAddSurface(subObjName,PanelSurface);
    
    this.matlNum=chiMaterialCreate(objName .. "_Material");
    chiMaterialSetProperty(this.matlNum,"DisableShading",true);
    ambient=0.8;
    chiMaterialSetProperty(this.matlNum,"Diffuse",ambient,ambient,ambient,1.0);
    --chiMaterialSetProperty(this.matlNum,"Diffuse",0.0,0.0,0.0,0.0);
    
    this.obj1Tra =chiTransformCreate(subObjName.."_Transform");
    this.obj1TTra=chiTransformCreate(subObjName.."_TextureTransform");
    this.obj1View=chiViewportCreate(subObjName.."_Viewport");
    
    chiObjectSetProperty(this.obj1Num,"Transform",subObjName.."_Transform");
    chiObjectSetProperty(this.obj1Num,"TextureTransform",subObjName.."_TextureTransform");
    chiObjectSetProperty(this.obj1Num,"Material",objName .. "_Material");
    chiObjectSetProperty(this.obj1Num,"Renderable",false);
    chiObjectSetProperty(this.obj1Num, "Viewport", this.obj1View);
    chiObjectSetProperty(this.obj1Num, "ViewportEnabled", true);
    chiObjectSetProperty(this.obj1Num,"donotList",true);
    
    chiTransformSetScale(this.obj1Tra ,this.xSize-1,this.ySize-1,1.0);
    chiTransformSetTranslation(this.obj1Tra,this.xpos, this.ypos, this.zDepth-0.1);
    
    return this;
end