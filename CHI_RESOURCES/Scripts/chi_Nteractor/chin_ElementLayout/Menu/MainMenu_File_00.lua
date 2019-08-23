--=============================================== Open
item = menuItem[1].subItems[2];
function item.CustomSelected(this)
    fileName = chiWindowFileDialogBox();


    chinOpenScene(fileName)
    
    
end

--=============================================== Save As
item = menuItem[1].subItems[4];
function item.CustomSelected(this)
    fileName = chiWindowFileDialogBox();
    chiBindScene(1);
    chinSaveScene(fileName)
    chiBindScene(0);
end

--=============================================== Import object
item = menuItem[1].subItems[5];
function item.CustomSelected(this)
    fileName = chiWindowFileDialogBox();
    chiBindScene(1);
    chinLoadObject(fileName)
    chiBindScene(0);
end