physicsFolderCMenu = ContextMenuClass.New("RDynamics")
physicsFolderCMenu.xpos = 200;
physicsFolderCMenu.ypos = 500;
physicsFolderCMenu.panel1.cutWindow=1;
physicsFolderCMenu.panel1.cutPanel=panels[4];
physicsFolderCMenu.Redraw(physicsFolderCMenu)

item  =physicsFolderCMenu.AddListItem(physicsFolderCMenu,"Select Models")

function item.CustomSelected(this)
    selectPhysicsModelsForm.Show(selectPhysicsModelsForm)
end

dRegisterFeature(physicsFolderCMenu)