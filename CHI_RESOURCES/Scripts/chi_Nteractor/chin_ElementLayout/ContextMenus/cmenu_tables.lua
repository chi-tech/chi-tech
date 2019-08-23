TablesFolderCMenu = ContextMenuClass.New("Tables")
TablesFolderCMenu.xpos = 200;
TablesFolderCMenu.ypos = 500;
TablesFolderCMenu.panel1.cutWindow=1;
TablesFolderCMenu.panel1.cutPanel=panels[4];
TablesFolderCMenu.Redraw(TablesFolderCMenu)

item  =TablesFolderCMenu.AddListItem(TablesFolderCMenu,"Add Table")

function item.CustomSelected(this)
    print("Creating table")
    temp = TableClass.New("Reactivity vs Time")
    
end

dRegisterFeature(TablesFolderCMenu)