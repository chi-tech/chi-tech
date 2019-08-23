numTrans=0;
prevnumTrans=0;
chinTransform={}
--######################################################### Update Transforms
function UpdateTransforms()
    numTrans = chiTransformQuery(0);
    if (not (numTrans == prevnumTrans)) then

        sk = 0;
        for k = 1,numTrans do

            if (chinTransform[k-1]==nil) then
                chinTransform[k-1]={}
            end

            chiTransformQuery(5,k-1,chinTransform[k-1]);

            if (chinTransform[k-1].listable) then
                sk = sk + 1;
                if (X_techTransformationsSubFolders[sk] == nil) then
                    chiTransformQuery(2,k-1);
                    X_techTransformationsSubFolders[sk] = X_techTransformationsFolder.AddFolder(X_techTransformationsFolder,chinTransform[k-1].name)
                    item = X_techTransformationsSubFolders[sk]
                    item.iconTypeFolder = chinIconTransform;
                    function item.label.CustomSelected(this)
                        local newSel = SelectionClass.New();
                        newSel.type = SELECTION_TYPE_TRANSFORMATION;
                        newSel.index = k-1;
                        selectionStack.PushItem(newSel);
                    end
                else
                    chiTransformQuery(2,k-1);
                    X_techTransformationsSubFolders[sk].name = chinTransform[k-1].name;
                end
            end
        end
        mainTree.Redraw(mainTree)
        numTrans = 0;
        numTrans = chiTransformQuery(0);
        prevnumTrans = numTrans;
    end
end