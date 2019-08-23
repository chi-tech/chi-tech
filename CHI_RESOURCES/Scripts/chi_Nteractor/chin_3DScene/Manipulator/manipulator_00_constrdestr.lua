ManipulatorClass = {}
ManipulatorClass.__index = ManipulatorClass
ManipulatorCount = 0;



--######################################################### Constructor
function ManipulatorClass.New(name)
    local this=setmetatable({},ManipulatorClass)
    ManipulatorCount=ManipulatorCount+1;
    
    local objName="Manipulator"..string.format("%03d",ManipulatorCount);
    
    this.name   = name;
    this.objectSlaves = {};
    this.objectSlaveCount = 0;
    
    this.lbuttondn = false;
    this.lbuttonscope = 0;
    this.rbuttondn = false;
    this.hidden = false;
    
    this.sceneScope = chiGetScene();
    
    --=============================================== Master transform
    this.manipTransform = chiTransformCreate(objName.."_Master_Transform");
    
    --=============================================== Load Z-translate
    this.zobjName ="Z-translate";
    this.zobjNum  = chiObjectCreate(this.zobjName);
    chiObjectAddSurface(this.zobjNum,newSurf);
    
    this.zmatNum = chiMaterialCreate(this.zobjName .. "_Material");
    chiMaterialSetProperty(this.zmatNum,"Diffuse",0,0,1,1.0);
    chiMaterialSetProperty(this.zmatNum,"DisableShading",true);
    chiObjectSetProperty(this.zobjNum,"Material",this.zobjName .. "_Material");
    chiObjectSetProperty(this.zobjNum,"Layer",1);
    
    this.zobjTra =chiTransformCreate(this.zobjName.."_Transform");
    chiObjectSetProperty(this.zobjNum,"Transform",this.zobjName.."_Transform");
    chiTransformSetParent(this.zobjTra,this.manipTransform);
    
    --=============================================== Load X-translate
    this.xobjName ="X-translate";
    this.xobjNum  = chiObjectCreate(this.xobjName);
    chiObjectAddSurface(this.xobjNum,newSurf);
    
    this.xmatNum = chiMaterialCreate(this.xobjName .. "_Material");
    chiMaterialSetProperty(this.xmatNum,"Diffuse",1,0,0,1.0);
    chiMaterialSetProperty(this.xmatNum,"DisableShading",true);
    chiObjectSetProperty(this.xobjNum,"Material",this.xobjName .. "_Material");
    chiObjectSetProperty(this.xobjNum,"Layer",1);
    
    this.xobjTra =chiTransformCreate(this.xobjName.."_Transform");
    chiObjectSetProperty(this.xobjNum,"Transform",this.xobjName.."_Transform");
    chiTransformSetParent(this.xobjTra,this.manipTransform);
    chiTransformSetRotation(this.xobjTra,0.0,90.0,0.0);
    
    --=============================================== Load Y-translate
    this.yobjName ="Y-translate";
    this.yobjNum  = chiObjectCreate(this.yobjName);
    chiObjectAddSurface(this.yobjNum,newSurf);
    
    this.ymatNum = chiMaterialCreate(this.yobjName .. "_Material");
    chiMaterialSetProperty(this.ymatNum,"Diffuse",0,1,0,1.0);
    chiMaterialSetProperty(this.ymatNum,"DisableShading",true);
    chiObjectSetProperty(this.yobjNum,"Material",this.yobjName .. "_Material");
    chiObjectSetProperty(this.yobjNum,"Layer",1);
    
    this.yobjTra =chiTransformCreate(this.yobjName.."_Transform");
    chiObjectSetProperty(this.yobjNum,"Transform",this.yobjName.."_Transform");
    chiTransformSetParent(this.yobjTra,this.manipTransform);
    chiTransformSetRotation(this.yobjTra,-90.0,0.0,0.0);
    
    return this;
end