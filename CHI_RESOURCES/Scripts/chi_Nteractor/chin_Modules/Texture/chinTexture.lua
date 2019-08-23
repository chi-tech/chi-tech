TextureClass = {}
TextureClass.__index = TextureClass
TextureCount = 0;

chinTexture = {}
textureStackCount = 0;
prevtextureStackCount = 0;

folderPath = chinGetScriptPath()
dofile(folderPath.."/chinTexture_00_constrdestr.lua")
dofile(folderPath.."/chinTexture_01_functions.lua")

--######################################################### Element manager
function UpdateTextureStack()

    textureStackCount = chiTextureQuery(0);
    if (not (textureStackCount==prevtextureStackCount)) then
    
        for k=1,textureStackCount do
        
            if (chinTexture[k-1]==nil) then
                chinTexture[k-1]=TextureClass.New(k-1);
            end
            
        end
        prevtextureStackCount=textureStackCount;
    end
    
end