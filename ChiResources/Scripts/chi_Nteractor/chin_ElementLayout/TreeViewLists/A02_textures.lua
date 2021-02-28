--Called from ../chiNteractor_04_treeview.lua

numTextures=0;
prevnumTextures=0;
--######################################################### Update textures
function UpdateTextures()
    numTextures=chiTextureQuery(0);
    if (not (numTextures==prevnumTextures)) then
        sk=0;
        for k=1,numTextures do
            sceneNum = chiTextureQuery(3,k-1);
            if (sceneNum > 0) then
                sk=sk+1;
                if (X_techTexturesSubFolders[sk]==nil) then

                    name=chinTrimPath(chinTexture[k-1].name);
                    X_techTexturesSubFolders[sk]=X_techTexturesFolder.AddFolder(X_techTexturesFolder,name);
                    item=X_techTexturesSubFolders[sk];
                    X_techTexturesSubFolders[sk].iconTypeFolder=chinIconTexture;
                    X_techTexturesSubFolders[sk].label.parent = chinTexture[k-1];
                    
                    --===================== Assigning callback to texture folder
                    function item.label.CustomSelected(this)
                        local newSel=SelectionClass.New();
                        
                        newSel.type=SELECTION_TYPE_TEXTURE;
                        newSel.index=this.parent.index;
                        newSel.originFeature=this;
                        selectionStack.PushItem(newSel);
                    end
                else
                end
            end
        end
    end
end



--######################################################### Callback for texture preview
function TextureSetPreview(selectionItem)

    if (selectionItem.type==SELECTION_TYPE_TEXTURE) then
        chiMaterialSetProperty(mfdDisplayPanel.matlNum,8,selectionItem.index);
        
    end
end

selectionStack.PushCallBackCount=selectionStack.PushCallBackCount+1;
k=selectionStack.PushCallBackCount;
selectionStack.PushCallBacks[k]=TextureSetPreview;