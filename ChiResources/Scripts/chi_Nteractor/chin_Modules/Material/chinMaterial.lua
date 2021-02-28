MaterialClass = {}
MaterialClass.__index = MaterialClass
MaterialCount = 0;

chinMaterial = {}
materialStackCount = 0;
prevmaterialStackCount = 0;

folderPath = chinGetScriptPath()
dofile(folderPath.."/chinMaterial_00_constrdestr.lua")
dofile(folderPath.."/chinMaterial_01_functions.lua")


--######################################################### Element manager
function UpdateMaterialStack()

    materialStackCount = chiMaterialQuery(0);
    if (not (materialStackCount==prevmaterialStackCount)) then
    
        for k=1,materialStackCount do
        
            if (chinMaterial[k-1]==nil) then
                chinMaterial[k-1]=MaterialClass.New(k-1);
            end
            
        end
        prevmaterialStackCount=materialStackCount;
    end
    
end