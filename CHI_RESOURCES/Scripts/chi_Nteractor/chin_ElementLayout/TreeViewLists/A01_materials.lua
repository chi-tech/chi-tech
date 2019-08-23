--Called from ../chiNteractor_04_treeview.lua


numMat          = 0;
prevnumMat      = 0;
--######################################################### Update materials
function UpdateMaterials()
    --================================= Update materials
    
    numMat=chiMaterialQuery(0);
    if (not (numMat==prevnumMat)) then
        
        sk=0;
        for k=1,numMat do
            
            
            if (chinMaterial[k-1].listable) then
            --print(chinMaterial[k-1].listable)
                sk=sk+1;
                if (X_techMaterialsSubFolders[sk]==nil) then
                   
                    X_techMaterialsSubFolders[sk]=X_techMaterialsFolder.AddFolder(X_techMaterialsFolder,chinMaterial[k-1].name)
                    
                    item=X_techMaterialsSubFolders[sk]
                    item.iconTypeFolder = chinIconMaterial;
                    item.label.parent = chinMaterial[k-1];
                    
                    --================= Adding callback function to material main folder
                    function item.label.CustomSelected(this)
                        GridviewHideAllItems();
                        local newSel=SelectionClass.New();
                        
                        newSel.type=SELECTION_TYPE_MATERIAL;
                        newSel.originFeature=this;
                        selectionStack.PushItem(newSel);
                    end
                    
                    --================= Setting property sub-folders
                    af=item.AddFolder(item,"Diffuse Color");  
                    bf=item.AddFolder(item,"Ambient Color");  
                    cf=item.AddFolder(item,"Specular Color"); 
                    df=item.AddFolder(item,"Shininess Color");
                    
                    af.iconTypeExpander={0,14};
                    bf.iconTypeExpander={0,14};
                    cf.iconTypeExpander={0,14};
                    df.iconTypeExpander={0,14};
                     
                    af.iconTypeFolder=chinIconRedBall;
                    bf.iconTypeFolder=chinIconRedBall;
                    cf.iconTypeFolder=chinIconRedBall;
                    df.iconTypeFolder=chinIconRedBall;
                    
                    --================= Connecting folder label to property
                    af.label.parent=chinMaterial[k-1].property[0];
                    bf.label.parent=chinMaterial[k-1].property[1];
                    cf.label.parent=chinMaterial[k-1].property[2];
                    df.label.parent=chinMaterial[k-1].property[3];
                    
                    --================== Adding callback function to subfolders
                    function func(this)
                        GridviewHideAllItems();
                        local newSel=SelectionClass.New();
                        newSel.type=SELECTION_TYPE_PROPERTY;
                        newSel.originFeature=this;
                        selectionStack.PushItem(newSel);
                    end
                    
                    af.label.CustomSelected=func;
                    bf.label.CustomSelected=func;
                    cf.label.CustomSelected=func;
                    df.label.CustomSelected=func;
                    
                    --================= Adding textures folder
                    ef=item.AddFolder(item,"Textures");  
                    tef={}
                    tef[1] = ef.AddFolder(ef,"Ambient");  
                    tef[2] = ef.AddFolder(ef,"Diffuse");
                    tef[3] = ef.AddFolder(ef,"Specular");
                    tef[4] = ef.AddFolder(ef,"Normal");
                    tef[5] = ef.AddFolder(ef,"Environment");
                    tef[6] = ef.AddFolder(ef,"Emissive");
                    
                    tindex={}
                    tindex[1]=chinMaterial[k-1].ambientTexture.index;
                    tindex[2]=chinMaterial[k-1].diffuseTexture.index;
                    tindex[3]=chinMaterial[k-1].specularTexture.index;
                    tindex[4]=chinMaterial[k-1].normalTexture.index;
                    tindex[5]=chinMaterial[k-1].environmentTexture.index;
                    tindex[6]=chinMaterial[k-1].emissiveTexture.index;
                    
                    tenabled={}
                    tenabled[1]=chinMaterial[k-1].ambientTexture.enabled;
                    tenabled[2]=chinMaterial[k-1].diffuseTexture.enabled;
                    tenabled[3]=chinMaterial[k-1].specularTexture.enabled;
                    tenabled[4]=chinMaterial[k-1].normalTexture.enabled;
                    tenabled[5]=chinMaterial[k-1].environmentTexture.enabled;
                    tenabled[6]=chinMaterial[k-1].emissiveTexture.enabled;

                    
                    --================== Creating subfolders for textures
                    for tcount=1,6 do
                        if (tindex[tcount]>-1) then
                            
                            name=chinTrimPath(chinTexture[tindex[tcount]].name);
                            aaef=tef[tcount].AddFolder(tef[tcount],name)
                            aaef.iconTypeFolder=chinIconTexture;
                            aaef.label.parent = chinTexture[tindex[tcount]];
                            aaef.label.parente= tenabled[tcount]
                            aaef.label.mparent = chinMaterial[k-1];
                            aaef.label.tindex = tcount;
                        
                            --===================== Assigning callback to texture folder
                            function aaef.label.CustomSelected(this)
                                GridviewHideAllItems();
                                local newSel=SelectionClass.New();
                                
                                newSel.type=SELECTION_TYPE_TEXTURE;
                                newSel.index=this.parent.index;
                                newSel.mindex=this.mparent.index;
                                newSel.tindex = this.tindex;
                                newSel.tenabled = this.parente;
                                newSel.originFeature=this;
                                selectionStack.PushItem(newSel);
                            end
                        end
                    end
                    

                    
                    
                else
                    --DO SOMETHING HERE
                end
            end
        end
        mainTree.Redraw(mainTree)
        numMat=0;
        numMat=chiMaterialQuery(0);
        prevnumMat=numMat;
    end
end

