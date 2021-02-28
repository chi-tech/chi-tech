savedTextures={}

function SaveTextures()

    io.write("\n--=========================== Textures\n")
    chiDeleteDirectoryContents(folderPath.."/Textures");
    chiCreateDirectoryS(folderPath.."/Textures");
    numTextures=chiTextureQuery(0);
    for k=1,numTextures do
        sceneNum = chiTextureQuery(3,k-1);
        sk=0;
        if (sceneNum > 0) then
            chiTextureQuery(4,k-1);
            sk=sk+1;
            savedTextures[k-1]={}
            savedTextures[k-1].savedIndex=sk;
            savedTextures[k-1].name = chinFilePathToFileName(chinTexture[k-1].name);
            
            
            if (sk>1) then
                chiTextureSave(k-1,"Textures/"..savedTextures[k-1].name);
                io.write("chiLoadTexture(\"Textures/"..savedTextures[k-1].name.."\");\n")
            end
        end
    end
end