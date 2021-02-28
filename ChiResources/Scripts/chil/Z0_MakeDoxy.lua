




os.execute("find . -name \"*.lua\" > Z1_FileList.txt");

fFile = io.open('Z1_FileList.txt',"r");
cFile = io.open('Z2_CPPequivalent.cpp',"w");

cFile:write("/** \\defgroup Lua_chil Y0 CHIL scripts \n \\ingroup LuaGeneralUtilities */\n");

for filename in fFile:lines() do

    if (string.find(filename,"Z0_MakeDoxy.lua") == nil) then
        extractGate = false;
        iFile = io.open(filename,"r");
        for fileLine in iFile:lines() do


            if     (string.find(fileLine,"--#{") ~= nil) then
                extractGate = true;
            elseif (string.find(fileLine,"--#}") ~= nil) then
                extractGate = false;
            elseif (extractGate==true) then
                --print(string.sub(fileLine,3))
                cFile:write(string.sub(fileLine,3).."\n");
            end

        end
        iFile:close();
    end

end

fFile:close();

cFile:write("};".."\n");
cFile:close();