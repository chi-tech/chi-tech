

--#########################################################
function TextBoxClass.New(name)
    local this=setmetatable({},TextBoxClass)
    TextBoxCount=TextBoxCount+1
    
    local objName="TextBox"..string.format("%02d",TextBoxCount);
    
    this.name=name;
    this.xpos=100;
    this.ypos=100;
    this.xSize=200;
    this.ySize=20;
    this.zDepth=1.1;
    
    this.paddingTop=0;
    this.paddingBot=0;
    this.paddingLeft=0;
    this.paddingRight=0;
    
    this.marginTop=0;
    this.marginBot=3;
    this.marginLeft=4;
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
    this.outlineColor={0,0,0,1};
    this.text=name;
    this.inputText="Input"..string.format("%02d",TextBoxCount);
    r=this.fontColor[1];
    g=this.fontColor[2];
    b=this.fontColor[3];
    a=this.fontColor[4];
    this.master=nil;
    this.parent=nil;
    
    this.selected=false;
    this.previousSelected=false;
    this.shiftDown=false;
    
    
    
    --=========================================== Creating base object
    this.obj1Num=chiObjectCreate(objName);
    chiObjectAddSurface(this.obj1Num,PanelSurface);
    
    this.transform1=chiTransformCreate(objName.."_Transform");
    chiObjectSetProperty(this.obj1Num,"Transform",objName.."_Transform");
    chiTransformSetScale(this.transform1,this.xSize-1,this.ySize-1,1.0)
    chiTransformSetTranslation(this.transform1,this.xpos,this.ypos,this.zDepth);
    this.obj1View=chiViewportCreate(subObjName.."_Viewport");
    chiObjectSetProperty(this.obj1Num, "Viewport", this.obj1View);
    chiObjectSetProperty(this.obj1Num, "ViewportEnabled", true);
    
    this.matlNum=chiMaterialCreate(objName .. "_Material");
    chiMaterialSetProperty(this.matlNum,"DisableShading",true);
    ambient=1.0;
    chiMaterialSetProperty(this.matlNum,"Diffuse",ambient,ambient,ambient,1.0);
    chiObjectSetProperty(this.obj1Num,"Material",objName .. "_Material");
    
    --=========================================== Creating outline    
    this.lineNum=chi3DLineCreate(objName.."outline");
    chi3DLineAddVertex(this.lineNum,this.xpos+this.paddingLeft,this.ypos+this.paddingBot,this.zDepth);
    chi3DLineAddVertex(this.lineNum,this.xpos+this.paddingLeft,this.ypos+this.ySize+this.paddingBot,this.zDepth);
    chi3DLineAddVertex(this.lineNum,this.xpos+this.xSize+this.paddingLeft,this.ypos+this.ySize+this.paddingBot,this.zDepth);
    chi3DLineAddVertex(this.lineNum,this.xpos+this.xSize+this.paddingLeft,this.ypos+this.paddingBot,this.zDepth);
    chi3DLineAddVertex(this.lineNum,this.xpos-1+this.paddingLeft,this.ypos+this.paddingBot,this.zDepth);
    col=this.outlineColor;
    chi3DLineChangeColor(this.lineNum,col[1],col[2],col[3],col[4]);
    chi3DLineSetviewport(this.lineNum,this.obj1View);
    
    --=========================================== Add cursor line
    this.cursorlineNum=chiLineCreate(objName.."cursor");
    chiLineAddVertex(this.cursorlineNum,this.xpos+this.paddingLeft+this.marginLeft,this.ypos+this.paddingBot+2,0.0);
    chiLineAddVertex(this.cursorlineNum,this.xpos+this.paddingLeft+this.marginLeft,this.ypos+this.paddingBot+this.ySize,0.0);
    chiLineChangeColor(this.cursorlineNum,0.0,0.0,0.0,0.0);
    
    --=========================================== Creating text
    this.textNum=chiSetLabel3D(objName,this.text,this.xpos+this.paddingLeft+this.marginLeft,this.ypos+this.paddingBot+this.marginBot,r,g,b,a,this.fontType)
    chiSetLabelProperty3D(this.textNum,"ViewportEnable",true);
    chiSetLabelProperty3D(this.textNum,"Depth",this.zDepth+0.01);
    
    --=========================================== Creating input label
    r=this.inputfontColor[1];
    g=this.inputfontColor[2];
    b=this.inputfontColor[3];
    a=this.inputfontColor[4];
    this.inputNum=chiSetLabel3D(objName.."input",this.inputText,this.xpos+this.paddingLeft+this.marginLeft,this.ypos+this.paddingBot+this.marginBot,r,g,b,a,this.fontType)
    chiSetLabelProperty3D(this.inputNum,"Depth",this.zDepth+0.01);
    --=========================================== Creating selection label
    
    r=this.selectfontColor[1];
    g=this.selectfontColor[2];
    b=this.selectfontColor[3];
    a=this.selectfontColor[4];
    this.selectNum=chiSetLabel3D(objName.."select",this.text,this.xpos+this.paddingLeft+this.marginLeft,this.ypos+this.paddingBot+this.marginBot,r,g,b,a,this.fontType)
    chiSetLabelProperty3D(this.selectNum,"Depth",this.zDepth+0.01);
    return this;
end