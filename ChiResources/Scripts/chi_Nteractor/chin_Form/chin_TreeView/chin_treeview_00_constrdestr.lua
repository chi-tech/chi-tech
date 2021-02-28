

--#########################################################
function TreeviewClass.New(name)
    local this              = setmetatable({},TreeviewClass)
    TreeviewCount           = TreeviewCount+1
    
    this.name               = name;

    this.subfolders         = {}
    this.subfolderCount     = 0;

    this.xpos               = 100;
    this.ypos               = 100;
    this.xmin               = 0;
    this.ymin               = 0;
    this.xmax               = 4000;
    this.ymax               = 4000;
    this.xSize              = 200;
    this.ySize              = 200;
    this.paddingTp          = 20;
    this.paddingLt          = 10;
    this.paddingLeft        = 0;
    this.paddingRight       = 0;
    this.paddingBot         = 0;
    this.paddingTop         = 0;
    this.lineSpacing        = 18;
    this.lineOffset         = 20;

    return this;
end


--#########################################################
function TreeviewFolderClass.New(name)
    local this              = setmetatable({},TreeviewFolderClass)
    TreeviewFolderCount     = TreeviewFolderCount+1;
    objName                 = "TreeFolder"..string.format("%04d",TreeviewFolderCount);
    
    this.name               = name;
    this.xpos               = 0;
    this.ypos               = 0;
    this.zpos               = 1.2;
    this.xmin               = 0;
    this.ymin               = 0;
    this.xmax               = 4000;
    this.ymax               = 4000;
    this.xSize              = 200;
    this.ySize              = 200;
    this.iconSize           = chinGlobal.Icons.iconSizeSmall;           --Overall size of icons
    this.iconSpacing1       = 2;                                        --Space between 1st icon and second
    this.iconSpacing2       = 0;                                        --Space between 2nd icon and text
    this.iconSpacing3       = 2;                                        --Y offset of icon
    this.iconSpacing4       = 3;                                        --Y offset of label
    this.iconTextureSize    = chinGlobal.Icons.iconTextureSize;
    this.iconCutOutSize     = chinGlobal.Icons.iconsCutSize;
    this.iconPadding        = chinGlobal.Icons.iconPadding;
    this.iconOffset         = chinGlobal.Icons.iconOffset;
    this.iconScale          = chinGlobal.Icons.iconScale;
    this.iconTypeExpander   = chinIconExpander;
    this.iconTypeReducer    = chinIconReducer;
    this.iconTypeFolder     = chinIconFolder;
    this.expanded           = false;
    this.nextLine           = 0;
    this.lineSpacing        = 15;
    this.lineOffset         = 35;

    this.parent             = nil;
    this.subfolders         = {}
    this.subfolderCount     = 0;

    --Label
    this.label      = LabelClass.New(name);
    this.label.SetProperty(this.label,"ViewportEnable",true)
    this.label.selectable = true;
    
    this.eventCallbacks     = {}
    this.eventCallbackCount = 0;
    
    --
    this.matlNum = chiMaterialCreate(objName .. "_Material");
    chiMaterialSetProperty(this.matlNum,"DisableShading",true);
    chiMaterialSetProperty(objName .. "_Material","DiffuseTexture",PanelTexture);
    --chiMaterialSetProperty(objName .. "_Material","Diffuse",0.0,0.0,0.0,0.0);


    --================================================ Expansion Icon
    subObjName      = objName.."_expander";
    this.obj1Num    = chiObjectCreate(subObjName);

    chiObjectAddSurface(subObjName,PanelSurface);
    
    this.obj1Tra    = chiTransformCreate(subObjName.."_Transform");
    this.obj1TTra   = chiTransformCreate(subObjName.."_TextureTransform");
    this.obj1View   = chiViewportCreate(subObjName.."_Viewport");
    
    chiObjectSetProperty(this.obj1Num, "Transform", subObjName.."_Transform");
    chiObjectSetProperty(this.obj1Num, "TextureTransform", subObjName.."_TextureTransform");
    chiObjectSetProperty(this.obj1Num, "Material", objName .. "_Material");
    chiObjectSetProperty(this.obj1Num, "Viewport", this.obj1View);
    chiObjectSetProperty(this.obj1Num, "ViewportEnabled", true);
    
    chiTransformSetScale(this.obj1Tra , this.iconSize,  this.iconSize,1.0);
    chiTransformSetScale(this.obj1TTra, this.iconScale, this.iconScale,1.0);

    dx = chinGlobal.Icons.iconShift*this.iconTypeExpander[1]
    dy = chinGlobal.Icons.iconShift*this.iconTypeExpander[2]

    chiTransformSetTranslation(this.obj1TTra, chinGlobal.Icons.iconTranslation +dx, chinGlobal.Icons.iconTranslation +dy, this.zpos);
    --================================================


    --================================================ Folder Icon
    subObjName      = objName.."_item";
    this.obj2Num    = chiObjectCreate(subObjName);

    chiObjectAddSurface(subObjName, PanelSurface);
    
    this.obj2Tra    = chiTransformCreate(subObjName.."_Transform");
    this.obj2TTra   = chiTransformCreate(subObjName.."_TextureTransform");
    this.obj2View   = chiViewportCreate(subObjName.."_Viewport");
    
    chiObjectSetProperty(this.obj2Num, "Transform", subObjName.."_Transform");
    chiObjectSetProperty(this.obj2Num, "TextureTransform", subObjName.."_TextureTransform");
    chiObjectSetProperty(this.obj2Num, "Material", objName .. "_Material");
    chiObjectSetProperty(this.obj2Num, "Viewport", this.obj2View);
    chiObjectSetProperty(this.obj2Num, "ViewportEnabled", true);
    
    chiTransformSetScale(this.obj2Tra , this.iconSize,  this.iconSize , 1.0)
    chiTransformSetScale(this.obj2TTra, this.iconScale, this.iconScale, 1.0)

    dx = chinGlobal.Icons.iconShift*this.iconTypeFolder[1]
    dy = chinGlobal.Icons.iconShift*this.iconTypeFolder[2]

    chiTransformSetTranslation(this.obj2TTra, chinGlobal.Icons.iconTranslation +dx, chinGlobal.Icons.iconTranslation +dy, this.zpos);
    --================================================


    --================================================ Highlight and selection
    subObjName      = objName.."_highlight";
    this.obj3Num    = chiObjectCreate(subObjName);

    chiObjectAddSurface(subObjName,PanelSurface);

    this.obj3Tra    = chiTransformCreate(subObjName.."_Transform");
    this.obj3TTra   = chiTransformCreate(subObjName.."_TextureTransform");
    this.obj3View   = chiViewportCreate(subObjName.."_Viewport");

    chiObjectSetProperty(this.obj3Num, "Transform", subObjName.."_Transform");
    chiObjectSetProperty(this.obj3Num, "TextureTransform", subObjName.."_TextureTransform");
    chiObjectSetProperty(this.obj3Num, "Material", objName .. "_Material");
    chiObjectSetProperty(this.obj3Num, "Viewport", this.obj3View);
    chiObjectSetProperty(this.obj3Num, "ViewportEnabled", true);

    chiTransformSetScale(this.obj3Tra , this.iconSize,  this.iconSize,1.0);
    chiTransformSetScale(this.obj3TTra, this.iconScale, this.iconScale,1.0);

    dx = chinGlobal.Icons.iconShift*this.iconTypeReducer[1]
    dy = chinGlobal.Icons.iconShift*this.iconTypeReducer[2]

    chiTransformSetTranslation(this.obj1TTra, chinGlobal.Icons.iconTranslation +dx, chinGlobal.Icons.iconTranslation +dy, this.zpos);
    --================================================

   
    
    --Default transforms
    chiTransformSetTranslation(this.obj1Tra ,this.xpos,this.ypos+this.iconSpacing3,this.zpos)
    chiTransformSetTranslation(this.obj2Tra ,this.xpos+this.iconSize+this.iconSpacing1,this.ypos+this.iconSpacing3,this.zpos)
    this.label.SetProperty(this.label,"Position",{this.xpos+2*this.iconSize+this.iconSpacing1+this.iconSpacing2+4,this.ypos});
    this.label.paddingBot = this.iconSpacing4;
    return this;
end
