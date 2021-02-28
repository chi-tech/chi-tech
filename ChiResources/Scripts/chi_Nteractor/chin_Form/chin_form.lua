folderPath = chinGetScriptPath()

chinCameraDir           = folderPath.."/chin_Camera/";
chinEventDir            = folderPath.."/chin_Event/";
chinPanelDir            = folderPath.."/chin_Panel/";
chinSplitBarDir         = folderPath.."/chin_SplitBar/";
chinTreeViewDir         = folderPath.."/chin_TreeView/";
chinWindowSetupDir      = folderPath.."/chin_WindowSetup/";
chinLabelDir            = folderPath.."/chin_Label/";
chinConsoleDir          = folderPath.."/chin_Console/";
chinTextBoxDir          = folderPath.."/chin_TextBox/";
chinGridViewDir         = folderPath.."/chin_GridView/";
chinSelectionStackDir   = folderPath.."/chin_SelectionStack/";
chinCheckBoxDir         = folderPath.."/chin_CheckBox/";
chinButtonDir           = folderPath.."/chin_Button/";
chinDropListDir         = folderPath.."/chin_DropList/";
chinMenuDir             = folderPath.."/chin_Menu/";
chinPropertyDir         = folderPath.."/chin_Property/";
chinMatPrevDir          = folderPath.."/chin_MaterialPreview/";
chinContextMenuDir      = folderPath.."/chin_ContextMenu/";



dofile(chinEventDir.."chin_event.lua");
dofile(chinCameraDir.."chin_camera.lua");
dofile(chinPanelDir.."chin_panel.lua");
dofile(chinSplitBarDir.."chin_splitbar.lua");
dofile(chinLabelDir.."chin_label.lua");
dofile(chinTreeViewDir.."chin_treeview.lua");
dofile(chinConsoleDir.."chin_console.lua");
dofile(chinTextBoxDir.."chin_textbox.lua");
dofile(chinGridViewDir.."chin_gridview.lua");
dofile(chinSelectionStackDir.."chin_selectionstack.lua");
dofile(chinCheckBoxDir.."chin_checkbox.lua");
dofile(chinButtonDir.."chin_button.lua");
dofile(chinDropListDir.."chin_droplist.lua");
dofile(chinMenuDir.."chin_menu.lua");
dofile(chinPropertyDir.."chin_property.lua");
dofile(chinMatPrevDir.."chin_materialpreview.lua");
dofile(chinContextMenuDir.."chin_contextmenu.lua");

FormClass = {}
FormClass.__index = FormClass
FormCount = 0;


--######################################################### Constructor
function FormClass.New(name,newWindow)
    local this=setmetatable({},FormClass);
    
    
    this.dwindowxsize=1280;
    this.dwindowysize=1024;
    this.dwindowxposi=10;
    this.dwindowyposi=10;
    
    this.feature={}
    this.feature.count=0;
    
    this.cycleCount=0;
    currentScene=chiGetScene();
    if (newWindow==nil) then
        this.sceneNumber=chiGetScene();
    else
        this.sceneNumber=chiWindowCreate(name);
    end
    
    chiBindScene(this.sceneNumber);
    
    parentWindowx,parentWindowy=chiGetWindowProperties();
    chiGraphicsCameraType(2);
    chiGraphicsCameraOrthoWidth(parentWindowx);
    chiGraphicsPositionCamera(parentWindowx/2,parentWindowy/2,100.0);
    chiBindScene(currentScene); 
    return this;
end


--######################################################### ProcessEvents
function FormClass.ProcessEvents(this)
    this.cycleCount=this.cycleCount+1;
    currentScene=chiGetScene();
    chiBindScene(this.sceneNumber);
    
    if (this.cycleCount==1) then
        parentWindowx,parentWindowy=chiGetWindowProperties();
        chiSetWindowProperties(parentWindowx,parentWindowy); --Flush window
    end
    
    if (WM_SIZE.occured) then
        parentWindowx,parentWindowy=chiGetWindowProperties();
        chiGraphicsCameraType(2);
        chiGraphicsCameraOrthoWidth(parentWindowx);
        chiGraphicsPositionCamera(parentWindowx/2,parentWindowy/2,100.0);
    end
    if (WM_CLOSE.occured) then
        if (WM_CLOSE.iPar4==this.sceneNumber) then
            this.Hide(this);
        end
    end
    
    for k=1,this.feature.count do
        this.feature[k].ProcessEvents(this.feature[k]);
    end  
    
    chiBindScene(currentScene);    
end

--######################################################### RegisterFeature
function FormClass.RegisterFeature(this,thefeature)
    this.feature.count=this.feature.count+1;
    k=this.feature.count;
    this.feature[k]=thefeature;
end

--########################################################## Show
function FormClass.Show(this)
    currentScene=chiGetScene();
    chiBindScene(this.sceneNumber);
    chiSetWindowProperties("RESTORED");
    chiBindScene(currentScene);
end

--########################################################## Hide
function FormClass.Hide(this)
    currentScene=chiGetScene();
    chiBindScene(this.sceneNumber);
    chiSetWindowProperties("HIDE");
    chiBindScene(currentScene);
end