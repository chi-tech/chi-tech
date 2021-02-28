--############################################### Constructor
function CheckBoxClass.New(name)
    local this=setmetatable({},CheckBoxClass)
    CheckBoxCount = CheckBoxCount + 1;
    
    local objName="CheckBox"..string.format("%03d",CheckBoxCount);
    this.name=name;
    this.xpos=100;
    this.ypos=100;
    this.zpos=1.2;
    this.xSize=200;
    this.ySize=20;
    this.paddingTop=0;
    this.paddingBot=0;
    this.paddingLeft=0;
    this.paddingRight=0;
    this.iconSize=20;
    this.float=true;
    this.selectable=false;
    this.selected=false;
    this.iconType=chinIconCheckboxU;
    
    this.eventCallbacks     = {}
    this.eventCallbackCount = 0;
    
    this.iconTextureSize    = chinGlobal.Icons.iconTextureSize;
    this.iconCutOutSize     = chinGlobal.Icons.iconsCutSize;
    this.iconPadding        = chinGlobal.Icons.iconPadding;
    this.iconOffset         = this.iconPadding*1+ this.iconCutOutSize;
    this.iconScale          = this.iconCutOutSize/this.iconTextureSize;
    
    --Material 
    this.matlNum=chiMaterialCreate(objName .. "_Material");
    chiMaterialSetProperty(this.matlNum,"DisableShading",true);
    chiMaterialSetProperty(objName .. "_Material","DiffuseTexture",PanelTexture);
    
    --General Icon
    subObjName=objName.."_icon";
    this.obj1Num=chiObjectCreate(subObjName);
    chiObjectAddSurface(subObjName,PanelSurface);
    
    this.obj1Tra =chiTransformCreate(subObjName.."_Transform");
    this.obj1TTra=chiTransformCreate(subObjName.."_TextureTransform");
    
    chiObjectSetProperty(this.obj1Num,"Transform",subObjName.."_Transform");
    chiObjectSetProperty(this.obj1Num,"TextureTransform",subObjName.."_TextureTransform");
    chiObjectSetProperty(this.obj1Num,"Material",objName .. "_Material");
    chiObjectSetProperty(this.obj1Num,"donotList",true);
    
    chiTransformSetScale(this.obj1Tra ,this.iconSize,this.iconSize,1.0);
    chiTransformSetScale(this.obj1TTra,this.iconScale,this.iconScale,1.0);

    dx = this.iconOffset*this.iconType[1]/this.iconTextureSize
    dy = this.iconOffset*this.iconType[2]/this.iconTextureSize
    chiTransformSetTranslation(this.obj1TTra,this.iconPadding/this.iconTextureSize+dx,this.iconPadding/this.iconTextureSize+dy,1.0);
    
    chiTransformSetTranslation(this.obj1Tra ,this.xpos,this.ypos,this.zpos)
    
    return this;
end