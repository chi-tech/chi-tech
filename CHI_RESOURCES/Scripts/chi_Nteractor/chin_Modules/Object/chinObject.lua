ObjectClass = {}
ObjectClass.__index = ObjectClass
ObjectCount = 0;

chinObject = {}
objectStackCount = 0;
prevobjectStackCount = 0;

folderPath = chinGetScriptPath()
dofile(folderPath.."/chinObject_00_constrdestr.lua")
dofile(folderPath.."/chinObject_01_functions.lua")


--######################################################### Element manager
function UpdateObjectStack()

    objectStackCount = chiObjectQuery(0);
    if (not (objectStackCount==prevobjectStackCount)) then
    
        for k=1,objectStackCount do
        
            if (chinObject[k-1]==nil) then
                chinObject[k-1]=ObjectClass.New(k-1);
            end
            
        end
        prevobjectStackCount=objectStackCount;
    end
    
end
