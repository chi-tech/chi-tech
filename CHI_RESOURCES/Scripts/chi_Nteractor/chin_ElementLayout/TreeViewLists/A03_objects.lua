--Called from ../chiNteractor_04_treeview.lua

numObj          = 0;
prevnumObj      = 0;
--######################################################### Update objects
function UpdateObjects()
    numObj = chiObjectQuery(0);
    if (not (numObj == prevnumObj)) then

        sk = 0;
        for k = 1,numObj do
            chiObjectQuery(4,k-1);
            --print("k ",sceneNum)
            --print(chinObject[k-1].name,chinObject[k-1].listable)
            if ( (chinObject[k-1].listable)) then
                sk = sk + 1;
                if (X_techObjectsSubFolders[sk] == nil) then
                
                    if (chinObject[k-1]==nil) then 
                        chiObjectQuery(2,k-1); 
                    end
                    
                    X_techObjectsSubFolders[sk] = X_techObjectsFolder.AddFolder(X_techObjectsFolder,chinObject[k-1].name)
                    
                    item = X_techObjectsSubFolders[sk]
                    item.iconTypeFolder = chinIconObject;
                    item.label.parent = chinObject[k-1];
                    
                    --================= Adding callback function to object main folder
                    function item.label.CustomSelected(this)
                        local newSel = SelectionClass.New();
                        newSel.type = SELECTION_TYPE_OBJECT;
                        newSel.originFeature=this;
                        selectionStack.PushItem(newSel);
                    end
                else

                end
            end
        end
        mainTree.Redraw(mainTree)
        numObj = 0;
        numObj = chiObjectQuery(0);
        prevnumObj = numObj;
    end
end