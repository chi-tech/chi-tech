--#########################################################
function TreeviewClass.ProcessEvents(this)

    TreeviewUpdateFlag = false;

    for k = 1,this.subfolderCount do
        this.subfolders[k].ProcessEvents(this.subfolders[k]);
    end

    if (TreeviewUpdateFlag) then
        TreeviewUpdateFlag = false;
        this.Redraw(this);
    end
    
    
end


--#########################################################
function TreeviewFolderClass.ProcessEvents(this)

    --======================= Left Mouse Button
    if (WM_LBUTTONDOWN.occured) then
        if (WM_LBUTTONDOWN.iPar5 == this.obj1Num) then
            if (this.expanded) then
                this.expanded       = false;
                TreeviewUpdateFlag  = true;
                this.iconTypeExpander=chinIconExpander;
                this.HideSubFolders(this);
            else
                this.expanded       = true;
                this.iconTypeExpander=chinIconReducer;
                this.UnHideSubFolders(this);
                TreeviewUpdateFlag  = true;
            end
        end
    end

    if (WM_LBUTTONUP.occured) then
        
    end
    
    for k = 1,this.subfolderCount do
        this.subfolders[k].ProcessEvents(this.subfolders[k]);
    end
    
    this.label.ProcessEvents(this.label);
    
    for k=1,this.eventCallbackCount do
        this.eventCallbacks[k](this);
    end
end


--#########################################################
function TreeviewClass.AddFolder(this,name)

    newFolder = TreeviewFolderClass.New(name);
    this.subfolderCount = this.subfolderCount+1;
    
    local k = this.subfolderCount
    this.subfolders[k] = newFolder;
    
    this.Redraw(this);
    
    return newFolder;
end


--#########################################################
function TreeviewClass.SetProperty(this,property,value)

    if (property == "Master") then
        this.master = value;
        
        value.slaveCount = value.slaveCount+1;
        k = value.slaveCount
        value.slaves[k] = this;

        --print(value.slaveCount)
        --v-alue.slaves[value.slaveCount]=this;
        --print("Master set as",value.name,value.slaves[value.slaveCount].name)
    end
end


--#########################################################
function TreeviewClass.SizeChanged(this)
    
    if (this.master ~= nil) then
        
        this.xmin = this.master.xmin;
        this.ymin = this.master.ymin;
        this.xmax = this.master.xmax-1;
        this.ymax = this.master.ymax-1;
        this.xSize = this.xmax - this.xmin;

        this.xpos = this.xmin+this.paddingLt;
        this.ypos = this.master.cursorY-this.paddingTp;
        --this.ypos = this.ymax-this.paddingTp;
        
        this.Redraw(this);
    end

end


-- ============================================================ Treeview Folder Classes
--#########################################################
function TreeviewFolderClass.AddFolder(this,name)

    local newFolder = TreeviewFolderClass.New(name);

    this.subfolderCount = this.subfolderCount+1;
    
    local k = this.subfolderCount
    this.subfolders[k] = newFolder;
    
    
    if (not this.expanded) then
        this.HideSubFolders(this);
    end

    this.Redraw(this);
    return newFolder;
end


--#########################################################
function TreeviewClass.Redraw(this)

    ypos = this.ypos;

    for k = 1,this.subfolderCount do
        this.subfolders[k].lineSpacing = this.lineSpacing;
        this.subfolders[k].lineOffset  = this.lineOffset;
        this.subfolders[k].xmin = this.xmin;
        this.subfolders[k].ymin = this.ymin;
        this.subfolders[k].xmax = this.xmax;
        this.subfolders[k].ymax = this.ymax;
        this.subfolders[k].xpos = this.xpos;
        this.subfolders[k].ypos = ypos;
        this.subfolders[k].Redraw(this.subfolders[k]);
        ypos = this.subfolders[k].nextLine;
    end
    this.ySize = this.ymax - ypos;
    if (this.master~=nil) then
        this.ySize = this.master.cursorY - ypos;
    end
end


--#########################################################
function TreeviewFolderClass.Redraw(this)

    --chiObjectSetProperty(this.obj1Num,"Hidden",false);
    --chiObjectSetProperty(this.obj2Num,"Hidden",false);
    chiTransformSetTranslation(this.obj1Tra ,this.xpos,this.ypos+this.iconSpacing3,this.zpos)
    chiTransformSetTranslation(this.obj2Tra ,this.xpos+this.iconSize+this.iconSpacing1,this.ypos+this.iconSpacing3,this.zpos)

    r = this.label.fontColor[1];
    g = this.label.fontColor[2];
    b = this.label.fontColor[3];
    a = this.label.fontColor[4];

    chiSetLabel(this.label.textNum,this.label.text,this.label.xpos,this.label.ypos,r,g,b,a,this.label.fontType)
    this.label.SetProperty(this.label,"Position",{this.xpos+2*this.iconSize+this.iconSpacing1+this.iconSpacing2+4,this.ypos});

    dx = this.iconOffset*this.iconTypeExpander[1]/this.iconTextureSize
    dy = this.iconOffset*this.iconTypeExpander[2]/this.iconTextureSize
    chiTransformSetTranslation(this.obj1TTra,this.iconPadding/this.iconTextureSize+dx,this.iconPadding/this.iconTextureSize+dy,this.zpos);

    dx = this.iconOffset*this.iconTypeFolder[1]/this.iconTextureSize
    dy = this.iconOffset*this.iconTypeFolder[2]/this.iconTextureSize
    chiTransformSetTranslation(this.obj2TTra,this.iconPadding/this.iconTextureSize+dx,this.iconPadding/this.iconTextureSize+dy,this.zpos);
    --
    xmin = this.xmin;
    --ymin = this.label.ypos;
    ymin = this.ymin;
    xmax = this.xmax;
    --ymax = this.label.ypos+this.label.ySize;
    ymax = this.ymax;

    chiViewportSetProperty(this.obj1View,xmin,ymin,xmax,ymax);
    chiViewportSetProperty(this.obj2View,xmin,ymin,xmax,ymax);
    chiViewportSetProperty(this.obj3View,xmin,ymin,xmax,ymax);
    this.label.SetProperty(this.label,"Viewport",{xmin,ymin,xmax,ymax});
    --print(xmin,ymin,xmax,ymax);

    ymin = this.label.ypos;
    ymax = this.label.ypos+this.label.ySize;

    chiTransformSetScale(this.label.obj1Tra ,xmax-xmin,ymax-ymin,1.0);
    chiTransformSetTranslation(this.label.obj1Tra,xmin, ymin, this.label.zDepth);


    ypos = this.ypos;
    this.nextLine = ypos-this.lineSpacing;
    if (this.expanded) then
        for k = 1,this.subfolderCount do
            this.subfolders[k].lineSpacing = this.lineSpacing;
            this.subfolders[k].lineOffset  = this.lineOffset;
            this.subfolders[k].xmin = this.xmin;
            this.subfolders[k].ymin = this.ymin;
            this.subfolders[k].xmax = this.xmax;
            this.subfolders[k].ymax = this.ymax;
            this.subfolders[k].xpos = this.xpos+this.lineOffset;
            this.subfolders[k].ypos = ypos-this.lineSpacing;
            this.subfolders[k].Redraw(this.subfolders[k]);
            this.nextLine = this.subfolders[k].nextLine;
        end
    end
end





--#########################################################
function TreeviewFolderClass.HideSubFolders(this)

    for k=1,this.subfolderCount do
        chiObjectSetProperty(this.subfolders[k].obj1Num,"Hidden",true);
        chiObjectSetProperty(this.subfolders[k].obj2Num,"Hidden",true);

        this.subfolders[k].label.Hide(this.subfolders[k].label);
        
        this.subfolders[k].HideSubFolders(this.subfolders[k]);
    end
end


--#########################################################
function TreeviewFolderClass.UnHideSubFolders(this)

    for k = 1,this.subfolderCount do
        
        chiObjectSetProperty(this.subfolders[k].obj1Num,"Hidden",false);
        chiObjectSetProperty(this.subfolders[k].obj2Num,"Hidden",false);

        this.subfolders[k].label.UnHide(this.subfolders[k].label);
        
        if (this.subfolders[k].expanded) then
            this.subfolders[k].UnHideSubFolders(this.subfolders[k]);
        end
    end
end