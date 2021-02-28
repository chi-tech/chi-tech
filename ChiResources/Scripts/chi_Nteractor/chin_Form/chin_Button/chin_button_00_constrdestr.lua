--######################################################### Constructor
function ButtonClass.New(name)
    local this=setmetatable({},ButtonClass);
    ButtonCount=ButtonCount+1;
    
    local objName="Button"..string.format("%03d",ButtonCount);
    this.name=name;
    this.xpos=100;
    this.ypos=100;
    this.xSize=50;
    this.ySize=20;
    this.zDepth=1.2;
    this.paddingTop=0;
    this.paddingBot=2;
    this.paddingLeft=2;
    this.paddingRight=2;
    this.textPadding=3;
    this.float=true;
    this.selectable=true;
    this.selected=false;
    this.highlightFillsMaster=false;
    this.yoffset=0;
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
    stringWidth=chinGetStringPixelWidth(this.text);
    textxpos=this.xpos+0.5*this.xSize-0.5*stringWidth;
    this.textNum=chiSetLabel(objName,this.text,textxpos,this.ypos+this.paddingBot+this.textPadding,r,g,b,a,this.fontType)
    this.master=nil;
    
    --Object
    subObjName=objName.."_Button";
    this.obj1Num=chiObjectCreate(subObjName);
    --print(subObjName);
    chiObjectAddSurface(subObjName,PanelSurface);
    
    this.matlNum=chiMaterialCreate(objName .. "_Material");
    chiMaterialSetProperty(this.matlNum,"DisableShading",true);
    ambient=0.8;
    chiMaterialSetProperty(this.matlNum,"Diffuse",ambient,ambient,ambient,1.0);
    
    this.obj1Tra =chiTransformCreate(subObjName.."_Transform");
    this.obj1TTra=chiTransformCreate(subObjName.."_TextureTransform");
    
    chiObjectSetProperty(this.obj1Num,"Transform",subObjName.."_Transform");
    chiObjectSetProperty(this.obj1Num,"TextureTransform",subObjName.."_TextureTransform");
    chiObjectSetProperty(this.obj1Num,"Material",objName .. "_Material");
    chiObjectSetProperty(this.obj1Num,"donotList",true);
    
    chiTransformSetScale(this.obj1Tra ,this.xSize-2*this.paddingLeft-1,this.ySize-2*this.paddingBot-1,1.0);
    chiTransformSetTranslation(this.obj1Tra ,this.xpos+this.paddingLeft,this.ypos+this.paddingBot,this.zDepth);
    
    --Line
    this.lineNum=chi3DLineCreate(objName.."outline");
    chi3DLineAddVertex(this.lineNum,this.xpos+this.paddingLeft,this.ypos+this.paddingBot,this.zDepth);
    chi3DLineAddVertex(this.lineNum,this.xpos+this.paddingLeft,this.ypos+this.ySize-2*this.paddingBot+2,this.zDepth);
    chi3DLineAddVertex(this.lineNum,this.xpos+this.xSize-2*this.paddingLeft+2,this.ypos+this.ySize-2*this.paddingBot+2,this.zDepth);
    chi3DLineAddVertex(this.lineNum,this.xpos+this.xSize-2*this.paddingLeft+2,this.ypos+this.paddingBot,this.zDepth);
    chi3DLineAddVertex(this.lineNum,this.xpos-1+this.paddingLeft,this.ypos+this.paddingBot,this.zDepth);
    chi3DLineChangeColor(this.lineNum,0.5,0.5,0.5,1.0);
    
    
    return this;
end


