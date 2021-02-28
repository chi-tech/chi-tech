folderPath = chinGetScriptPath()
dofile(folderPath.."/TreeViewLists/A01_materials.lua")
dofile(folderPath.."/TreeViewLists/A02_textures.lua")
dofile(folderPath.."/TreeViewLists/A03_objects.lua")
dofile(folderPath.."/TreeViewLists/A04_transforms.lua")
dofile(folderPath.."/TreeViewLists/A05_shaders.lua")
dofile(folderPath.."/TreeViewLists/A99_dxdevices.lua")
dofile(folderPath.."/TreeViewLists/B01_physics.lua")
dofile(folderPath.."/TreeViewLists/B02_timecontrol.lua")
dofile(folderPath.."/TreeViewLists/D01_tables.lua")

--=============================================== Treeview update function
--called from chiNteractor_00_layout.lua:CustomProcessing()
function TreeViewUpdate()

    UpdateObjects();
    UpdateTextures();
    UpdateMaterials();
    
    UpdateTransforms();
    
    UpdateShaders();
    UpdateDXDevices();
    UpdateReactorDynamics1D();
    UpdateTimeControl();
    
    UpdateTables()
end




--=============================================== TreeView bare skeleton
mainTree=TreeviewClass.New("MainTree");
mainTree.SetProperty(mainTree,"Master",panels[3])

dRegisterFeature(mainTree);
--mainTree.xpos=10;
--mainTree.ypos=600;

X_techElementFolder=mainTree.AddFolder(mainTree,"X-Tech Elements")
    X_techObjectsFolder=X_techElementFolder.AddFolder(X_techElementFolder,"Objects")
    X_techObjectsSubFolders={}

    X_techMaterialsFolder=X_techElementFolder.AddFolder(X_techElementFolder,"Materials")
    X_techMaterialsSubFolders={}

    X_techTexturesFolder=X_techElementFolder.AddFolder(X_techElementFolder,"Textures")
    X_techTexturesSubFolders={}

    X_techTransformationsFolder=X_techElementFolder.AddFolder(X_techElementFolder,"Transformations")
    X_techTransformationsSubFolders={}

    X_techShadersFolder=X_techElementFolder.AddFolder(X_techElementFolder,"Shaders")
    X_techShadersSubFolders={}

PhysicsFolder = mainTree.AddFolder(mainTree,"Physics");
PhysicsFolder.iconTypeFolder = chinIconFolderBlue;

TimeFolder = mainTree.AddFolder(mainTree,"Time step control");
TimeFolder.iconTypeFolder = chinIconFolderBlue;

ToolsFolder = mainTree.AddFolder(mainTree,"Tools");
ToolsFolder.iconTypeFolder = chinIconFolderBlue;
    
    TablesFolder = ToolsFolder.AddFolder(ToolsFolder,"Tables");
    TablesFolder.iconTypeFolder = chinIconFolderBlue;

--=============================================== Context menus for bare skeleton
function PhysicsFolderContextMenuControl(this)
    if (WM_RBUTTONDOWN.occured) then
        if (WM_RBUTTONDOWN.iPar5==this.label.obj1Num) then
            parentWindowx,parentWindowy=chiGetWindowProperties();
            physicsFolderCMenu.xpos = WM_RBUTTONDOWN.iPar0;
            physicsFolderCMenu.ypos = parentWindowy-WM_RBUTTONDOWN.iPar1;
            physicsFolderCMenu.Redraw(physicsFolderCMenu);
            physicsFolderCMenu.Selected(physicsFolderCMenu);
            
        else
            physicsFolderCMenu.Redraw(physicsFolderCMenu);
            physicsFolderCMenu.DeSelected(physicsFolderCMenu);
            
        end
    end

end
PhysicsFolder.eventCallbackCount=PhysicsFolder.eventCallbackCount+1;
PhysicsFolder.eventCallbacks[PhysicsFolder.eventCallbackCount]=PhysicsFolderContextMenuControl;

function TablesFolderContextMenuControl(this)
    if (WM_RBUTTONDOWN.occured) then
        if (WM_RBUTTONDOWN.iPar5==this.label.obj1Num) then
            parentWindowx,parentWindowy=chiGetWindowProperties();
            TablesFolderCMenu.xpos = WM_RBUTTONDOWN.iPar0;
            TablesFolderCMenu.ypos = parentWindowy-WM_RBUTTONDOWN.iPar1;
            TablesFolderCMenu.Redraw(TablesFolderCMenu);
            TablesFolderCMenu.Selected(TablesFolderCMenu);
        else
            TablesFolderCMenu.Redraw(TablesFolderCMenu);
            TablesFolderCMenu.DeSelected(TablesFolderCMenu);
        end
    end

end
TablesFolder.eventCallbackCount=TablesFolder.eventCallbackCount+1;
TablesFolder.eventCallbacks[TablesFolder.eventCallbackCount]=TablesFolderContextMenuControl;
