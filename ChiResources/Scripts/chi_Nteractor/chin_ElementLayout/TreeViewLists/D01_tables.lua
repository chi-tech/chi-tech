--Called from ../chiNteractor_04_treeview.lua

numTables       = 0;
prevnumTables   = 0;
chinTable       = {};
--######################################################### Update Tables
function UpdateTables()
    numTables = chiTableQuery(0);
    if (not (numTables == prevnumTables)) then
        for k = 1,numTables do
            chiTableQuery(4,k-1);
            
            TabelFolderSubFolder = TablesFolder.AddFolder(TablesFolder,chinTable[k-1].name);
            item = TabelFolderSubFolder;
            
            --======================= Create selection call back for properties
            function func(this)
                GridviewHideAllItems();
                local newSel=SelectionClass.New();
                newSel.type=SELECTION_TYPE_TABLE;
                newSel.index=k-1;
                newSel.originFeature=this;
                selectionStack.PushItem(newSel);
            end
            item.label.CustomSelected = func;
            
        end
        numTables = 0;
        numTables = chiTableQuery(0);
        prevnumTables = numTables;
    end
end