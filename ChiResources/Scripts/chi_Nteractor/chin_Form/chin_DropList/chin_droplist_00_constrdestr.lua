--######################################################### Constructor
function DropListClass.New(name)
    local this=setmetatable({},DropListClass);
    DropListCount=DropListCount+1;
    
    local objName="DropList"..string.format("%03d",DropListCount);
    this.name=name;
    this.xpos=100;
    this.ypos=100;
    this.xSize=200;
    this.ySize=20;
    this.zDepth=1.3;
    
    this.paddingTop=0;
    this.paddingBot=0;
    this.paddingLeft=0;
    this.paddingRight=0;
    
    this.marginTop=0;
    this.marginBot=3;
    this.marginLeft=3;
    this.marginRight=0;
    
    this.eventCallbacks     = {}
    this.eventCallbackCount = 0;
    
    this.cursorPosition=0;
    this.cursorX=0;
    this.lastBlinkUpdated=0.0;
    this.cursorShown=false;
    
    this.float=true;
    this.inputBox=false;
    this.fontType=2;
    this.fontColor={0,0,0,1};
    this.inputfontColor={0.5,0.5,0.5,1}
    this.selectfontColor={0,0,0.8,1}
    this.text=name;
    this.inputText="Input"..string.format("%02d",DropListCount);
    r=this.fontColor[1];
    g=this.fontColor[2];
    b=this.fontColor[3];
    a=this.fontColor[4];
    this.master=nil;
    
    this.selected=false;
    this.previousSelected=false;
    this.shiftDown=false;
    this.selectedOption = 0;
    
    this.iconSize           = 21;
    this.iconTextureSize    = chinGlobal.Icons.iconTextureSize;
    this.iconCutOutSize     = chinGlobal.Icons.iconsCutSize;
    this.iconPadding        = chinGlobal.Icons.iconPadding;
    this.iconOffset         = chinGlobal.Icons.iconOffset;
    this.iconScale          = chinGlobal.Icons.iconScale;
    this.iconTypeFolder     = chinIconExpander;
    
    this.panelxSize         = 250;
    
    --=========================================== Creating base object
    this.obj1Num=chiObjectCreate(objName);
    chiObjectAddSurface(this.obj1Num,PanelSurface);
    
    this.transform1=chiTransformCreate(objName.."_Transform");
    chiObjectSetProperty(this.obj1Num,"Transform",objName.."_Transform");
    chiTransformSetScale(this.transform1,this.xSize-1,this.ySize-1,1.0)
    chiTransformSetTranslation(this.transform1,this.xpos,this.ypos,this.zDepth)
    
    this.matlNum=chiMaterialCreate(objName .. "_Material");
    chiMaterialSetProperty(this.matlNum,"DisableShading",true);
    ambient=1.0;
    chiMaterialSetProperty(this.matlNum,"Diffuse",ambient,ambient,ambient,1.0);
    chiObjectSetProperty(this.obj1Num,"Material",objName .. "_Material");
    chiObjectSetProperty(this.obj1Num,"donotList",true);
    
    --=========================================== Creating outline    
    this.lineNum=chi3DLineCreate(objName.."outline");
    chi3DLineAddVertex(this.lineNum,this.xpos+this.paddingLeft,this.ypos+this.paddingBot,this.zDepth+0.01);
    chi3DLineAddVertex(this.lineNum,this.xpos+this.paddingLeft,this.ypos+this.ySize-1+this.paddingBot,this.zDepth+0.01);
    chi3DLineAddVertex(this.lineNum,this.xpos+this.xSize-1+this.paddingLeft,this.ypos+this.ySize-1+this.paddingBot,this.zDepth+0.01);
    chi3DLineAddVertex(this.lineNum,this.xpos+this.xSize-1+this.paddingLeft,this.ypos+this.paddingBot,this.zDepth+0.01);
    chi3DLineAddVertex(this.lineNum,this.xpos-1+this.paddingLeft,this.ypos+this.paddingBot,this.zDepth+0.01);
    chi3DLineChangeColor(this.lineNum,0.0,0.0,0.0,1.0);
    
    --=========================================== Creating text
    this.textNum=chiSetLabel3D(objName,this.text,this.xpos+this.paddingLeft+this.marginLeft,this.ypos+this.paddingBot+this.marginBot,r,g,b,a,this.fontType);
    chiSetLabelProperty(this.textNum,"ViewportEnable",true);
    chiSetLabelProperty(this.textNum,"Viewport",this.xpos+this.paddingLeft,this.ypos+this.paddingBot,this.xpos+this.xSize-1+this.paddingLeft-this.iconSize,this.ypos+this.ySize-1+this.paddingBot);
    chiSetLabelProperty3D(this.textNum,"Depth",this.zDepth+0.01);
    
    --=========================================== Dropdown icon
    subObjName      = objName.."_icon";
    this.obj2Num    = chiObjectCreate(subObjName);

    chiObjectAddSurface(subObjName, PanelSurface);
    
    this.obj2Tra    = chiTransformCreate(subObjName.."_Transform");
    this.obj2TTra   = chiTransformCreate(subObjName.."_TextureTransform");
    
    this.matlNum=chiMaterialCreate(subObjName .. "_Material");
    chiMaterialSetProperty(this.matlNum,"DisableShading",true);
    chiMaterialSetProperty(subObjName .. "_Material","DiffuseTexture",PanelTexture);
    
    chiObjectSetProperty(this.obj2Num, "Transform", subObjName.."_Transform");
    chiObjectSetProperty(this.obj2Num, "TextureTransform", subObjName.."_TextureTransform");
    chiObjectSetProperty(this.obj2Num, "Material", subObjName .. "_Material");
    chiObjectSetProperty(this.obj2Num,"donotList",true);
    
    chiTransformSetScale(this.obj2Tra , this.iconSize,  this.iconSize , 1.0)
    chiTransformSetScale(this.obj2TTra, this.iconScale, this.iconScale, 1.0)

    dx = chinGlobal.Icons.iconShift*this.iconTypeFolder[1]
    dy = chinGlobal.Icons.iconShift*this.iconTypeFolder[2]

    chiTransformSetTranslation(this.obj2Tra , this.xpos+this.xSize-this.iconSize,  this.ypos-2, this.zDepth+0.01)
    chiTransformSetTranslation(this.obj2TTra, chinGlobal.Icons.iconTranslation +dx, chinGlobal.Icons.iconTranslation +dy, 0.0);
    
    --=========================================== Panel for list
    this.panel1=PanelClass.New(objName.."Panel");
    this.panel1.zDepth=this.zDepth;
    this.panel1.xmin=this.xpos;
    this.panel1.ymax=this.ypos;
    this.panel1.xmax=this.panel1.xmin+this.panelxSize+1;
    this.panel1.ymin=this.panel1.ymax-this.ySize;
    this.panel1.SizeChanged(this.panel1);
    
    this.list={};
    this.listCount=0;
    
    chiObjectSetProperty(this.panel1.obj1Num,"Hidden",true);
    chi3DLineChangeColor(this.panel1.lineNum,0.0,0.0,0.0,0.0);
    
    for k=1,this.listCount do
        this.list[k].Hide(this.list[k]);
    end
    
    
    return this;
end