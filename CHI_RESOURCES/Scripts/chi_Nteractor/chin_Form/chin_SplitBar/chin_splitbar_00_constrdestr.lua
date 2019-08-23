
--######################################################### Constructor
function SplitBarClass.New(name)
    local this=setmetatable({},SplitBarClass)
    SplitBarCount=SplitBarCount+1;
    
    this.name=name;
    this.xmin=2;
    this.ymin=0;
    this.xmax=10;
    this.ymax=500;
    this.gripWidth=10;
    this.sizeChanging=false;
    this.vertical=true;
    this.slaves={}
    this.slaveCount=0;
    
    this.eventCallbacks     = {}
    this.eventCallbackCount = 0;
    
    objName="SplitBar"..string.format("%02d",SplitBarCount);
    
    this.obj1Num=chiObjectCreate(objName.."_Handle1");
    
    chiObjectAddSurface(objName.."_Handle1",SplitBarSurface);
    
    this.transform1=chiTransformCreate(objName.."_Handle1".."_Transform");
    chiObjectSetProperty(this.obj1Num,"Transform",objName.."_Handle1".."_Transform");
    chiObjectSetProperty(this.obj1Num,"Renderable",false);

    this.matlNum=chiMaterialCreate(objName .. "_Material");
    chiMaterialSetProperty(this.matlNum,"DisableShading",true);
    ambient=0.82;
    chiMaterialSetProperty(this.matlNum,"Diffuse",ambient,ambient,ambient,1.0);
    chiObjectSetProperty(this.obj1Num,"Material",objName .. "_Material");
    
    
    
    girth =this.xmax-this.xmin-1;
    length=this.ymax-this.ymin-1;
    chiTransformSetScale(this.transform1,girth,length,1); 
    chiTransformSetTranslation(this.transform1,this.xmin,this.ymin,1.2); 
    return this;
end