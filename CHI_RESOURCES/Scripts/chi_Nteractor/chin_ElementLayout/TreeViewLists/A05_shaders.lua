numShaders=0;
prevnumShaders=0;
--######################################################### Update Transforms
function UpdateShaders()
    numShaders = chiShaderQuery(0);
    if (not (numShaders == prevnumShaders)) then
        sk = 0;
        for k = 1,numShaders do
            sk = sk + 1;

            if (X_techShadersSubFolders[sk] == nil) then
                chiShaderQuery(2,k-1);
                X_techShadersSubFolders[sk] = X_techShadersFolder.AddFolder(X_techShadersFolder,chinShader[k-1].name)
                item = X_techShadersSubFolders[sk]
                item.iconTypeFolder = chinIconTransform;

                function item.label.CustomSelected(this)
                    local newSel = SelectionClass.New();
                    newSel.type = SELECTION_TYPE_SHADER;
                    newSel.index = k-1;
                    selectionStack.PushItem(newSel);
                end
            else
                chiShaderQuery(2,k-1);
                X_techShadersSubFolders[sk].name = chinShader[k-1].name;
            end
        end
        mainTree.Redraw(mainTree)
        numShaders = 0;
        numShaders = chiShaderQuery(0);
        prevnumShaders = numShaders;
    end
end