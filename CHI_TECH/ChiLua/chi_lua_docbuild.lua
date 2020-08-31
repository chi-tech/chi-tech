
dofile("module_lua_inclusion.lua")

--[[
Function that reads in the file into memory and loads it into a table for parsing.
]]--
function readCXXFile()
	local fTable, fName, fContent;
	local index, value;

	fTable = scanDirectory();

	for index,value in pairs(fTable) do
		fName = io.open(value, 'r+');
		fContent = fName:read("*all");
		print("Searching: "..value)
		searchCXXFile(fContent);
		fName:close();
	end

	print("=============== Finished Reading lua files ===============")
end


--[[
Funtion that goes link by line looking for file contents.
]]--
function searchCXXFile(fileContent)
	local line, wIndex, pCount, count 		= 0, 0, 0, 1;
	local inComment, inFunction, inParam  	= false, false, false;
	local wTable, pTable, wTable2			= {}, {}, {};
	local fLine, fParam						= "","";

	for line in string.gmatch(fileContent, "[^\r\n]+") do
		fParam = "";


		--============================ Searching for doxygen style comments
		if (string.find(line, "/%*%*") ~= nil) then
			wIndex 			= wIndex + 1;
			inComment 		= true;
			wTable[wIndex] 	= "";
			wTable2[wIndex] = "";
		end

		if (inComment == true) then
			wTable[wIndex] 	= wTable[wIndex]..line.."\n";
		end

		if (string.find(line, "%*/") ~= nil) then
			inComment 		= false;
			if (wTable[wIndex] ~= nil) then
				--print(wTable[wIndex]);
			end
		end


		--============================ Searching for doxygen param
		if (string.find(line, "\\param") ~= nil) then
			wIndex 			= wIndex + 1;
			pCount 			= pCount + 1;
			inParam 		= true;
			wTable[wIndex] 	= "";
			wTable2[wIndex] = "";
            
		end

		if (inParam == true) then
			pTable[pCount] 	= extractParameter(line);
            inParam 		= false;
		end

		if (string.find(line, "|") ~= nil) then
			inParam 		= false;
			if (wTable[wIndex] ~= nil) then
				--print(wTable[wIndex]);
			end
		end


		--============================ Searching for function name
		if ((string.find(line,"lua_State") ~= nil) and (string.find(line,"chi") ~= nil)) then
			wIndex 			= wIndex + 1;
			inFunction 		= true;
			wTable[wIndex] 	= "";
			wTable2[wIndex] = "";
		end

		if (inFunction == true) then

			while (count <= pCount) do
				local item = pTable[count];
				if (count == pCount) then
					fParam = fParam..item[1].." "..item[2];
				else
					fParam = fParam..item[1].." "..item[2]..", ";
				end
				count = count + 1;
			end

			count 			= 1;
			pCount 			= 0;

			local lBegin, lEnd = string.find(line,"chi");

			--fLine = string.sub(line,1,lBegin-1);
			--fLine = fLine..string.sub(line,lEnd+1);
               --fLine = string.sub(line,lEnd+1);
               fLine = string.sub(line,lBegin);
            --fLine = line;
            --print(fLine)

			local cBegin, cEnd = string.find(fLine,"%(");
			local fLine2 = string.sub(fLine,1,cBegin-1);

            
			wTable2[wIndex] 	= wTable2[wIndex].."int "..fLine2.."("..fParam..");\n";
            


			local cBegin, cEnd = string.find(fLine,"%(");
			local fLine3 = string.sub(fLine,1,cBegin-1);

            
			wTable[wIndex] 	= wTable[wIndex].."int chi_lua::"..fLine3.."("..fParam..")\n".."{return;} \n";
            
			--print(fParam)
			--print(wTable[wIndex])
			inFunction 		= false;

			lBegin, lEnd, cBegin = 0,0,0;
		end
	end
	writeDocoFile(wTable);
	writeDocoFile2(wTable2);
end


--[[
Function that fractures strings into words, and extracts the parameter followed by the value
]]--
function extractParameter(fString)
	local sIndex = 0;
	local sTable = {};
	local sParam = {};

	for sValue in string.gmatch(fString, "%S+") do
		sIndex = sIndex + 1;
		sTable[sIndex] = sValue;
	end

	sParam[1] = sTable[3];
	sParam[2] = sTable[2];

	return sParam;
end


--[[
Function that writes the extracted strings into a seperate file.
]]--
function writeDocoFile(wordTable)
	local index, value;

	for index,value in pairs(wordTable) do
		newFile:write(value);
	end
end


--[[
Function that writes the extracted strings into a seperate file.
]]--
function writeDocoFile2(wordTable)
	local index, value;

	for index,value in pairs(wordTable) do
		newFile2:write(value);
	end
end


--[[
Windows specific funtion that goes through directories
recursively and creates a table with the various filenames.

** Note:
The function already contains the absolute path.
]]--
function scanDirectory()
	local fIndex, fTable, popen = 0, {}, io.popen;
	local fFile;

	--os.execute("dir /B/S *.cpp > LUA_DOCUMENTATION/Z_GeneratedFileList.txt");
	print("Generating filelist");
	os.execute("mkdir LUA_DOCUMENTATION")
	os.execute(">LUA_DOCUMENTATION/Z_GeneratedFileList.txt")
	os.execute("find . -name \"*.cc\" >> LUA_DOCUMENTATION/Z_GeneratedFileList.txt");
	os.execute("find . -name \"*.cpp\">> LUA_DOCUMENTATION/Z_GeneratedFileList.txt");


    
    if (moduleFolders.itemCount>0) then
        for folderDex=1,moduleFolders.itemCount do 
            --os.execute("dir \""..moduleFolders[folderDex].."\" /B/S *.cpp >> LUA_DOCUMENTATION/Z_GeneratedFileList.txt");
						os.execute("find \""..moduleFolders[folderDex].."\" -name \"*.cpp\" >> LUA_DOCUMENTATION/Z_GeneratedFileList.txt");
						os.execute("find \""..moduleFolders[folderDex].."\" -name \"*.cc\" >> LUA_DOCUMENTATION/Z_GeneratedFileList.txt");
        end
        
    end

	fFile = io.open('LUA_DOCUMENTATION/Z_GeneratedFileList.txt',"r");

	for filename in fFile:lines() do
		fIndex = fIndex + 1;
		fTable[fIndex] = filename;
	end

	fFile:close();

	return fTable
end


--[[
To keep the lua code a little bit cleaner a simple function was
written which takes a table and outputs the size of the table.

** Note:
Used for testing purposes.
]]--
function tableSize(tTable)
	local counter = 0;
	for _ in pairs(tTable) do counter = counter + 1; end
	return counter;
end

--[[
This functions writes a table in HTML format
--]]
function WriteTableHTML(newFile, table_name, table_ref, table_vals)
	if (rawlen(table_vals) > 0) then
		newFile:write("### "..table_name.."\n")
		newFile:write(table_ref.."\n")
		newFile:write("<table>\n")
		colcnt = 0
		for k=1,(rawlen(table_vals)) do
			if (colcnt == 0) then
				newFile:write("<tr>")
			end
			colcnt = colcnt + 1
			newFile:write("<td width=\"33%\">"..table_vals[k].."()".."</td>")

			if ((colcnt == 3) or (k == rawlen(table_vals))) then
				colcnt = 0
				newFile:write("</tr>\n")
			end
		end
		newFile:write("</table>\n")
	end
end

--[[
This function reads the lua register and creates a formatted table
containing all the function names.
]]--
function generateFunctionList()
	print("Generating function list for main page")
	newFile = io.open("../../CHI_DOC/PAGES/MainPage.h","w");
	headFile = io.open("../../CHI_DOC/PAGES/mp_header.txt","r")
	reg = io.open("../../CHI_DOC/LUA_DOCUMENTATION/lua_register.txt","r");

	--=============================== Write header
	line = headFile:read("L")
	while (line ~= nil) do
		newFile:write(line)
		line = headFile:read("L")
	end

	block = false

	linecount = 0
	line = reg:read()
	linecount = linecount + 1
	table_vals = {}
	table_name = " "
	table_ref = " "
	strword = nil
	strline = nil
	while (line ~= nil) do

		ifdef = string.find(line,"#ifdef",0,true)
		if (ifdef ~= nil) then
			block = true
		end

		endif = string.find(line,"#endif",0,true)
		if (endif ~= nil) then
			block = false
		end

		--======================================== "module:" keyword
		modulestr = string.find(line,"module:", 0,true)
		if (modulestr ~= nil) then
			if (rawlen(table_vals) > 0) then
				WriteTableHTML(newFile, table_name, table_ref, table_vals)
			end
			table_name = string.sub(line,modulestr+7,-1)
			table_vals = {}
			table_ref = " "
		end

		if (block ~= true) then
			--=================================== "RegisterFunction" keyword
			open_paren = string.find(line,"RegisterFunction(",0,true)

			if (open_paren ~= nil) then
				clos_paren = string.find(line,")")

				if ((clos_paren ~= nil) and (linecount>5)) then
					k = rawlen(table_vals)+1
					table_vals[k] = string.sub(line,open_paren+17,clos_paren-1)
				end
			end

			--=================================== "\ref " keyword
			ref_start = string.find(line,"//\\ref",0,true)

			if (ref_start ~= nil) then
				table_ref = string.sub(line,ref_start+2)
			end

			--=================================== "string:" keyword
			strword = string.find(line,"string:",0,true)
			if (strword ~= nil) then
				strline = line
			end
		end





		line = reg:read()
		linecount = linecount + 1
	end

	--======================================== Last table in buffer
	if (rawlen(table_vals) > 0) then
		WriteTableHTML(newFile, table_name, table_ref, table_vals)
	end
	--======================================== END OF WRITES A TABLE

	if (strword ~= nil) then
		newFile:write(string.sub(strline,10).."\n")
	end

	newFile:write("\n\n*/\n")

	headFile:close()
	reg:close();
	newFile:close();

	print("Done generating function list for main page")
end

-- ============================================== Recursive includes
--[[
This function can be called recursively to find #include statements
and put them into the master file.
]]--
function IncludeIntoMaster(master_file,slave_file_name)
	print("Opening Slave file: "..slave_file_name)
	local in_file = io.open(slave_file_name)

	if (in_file == nil) then
		print("ERROR (non-fatal): Could not open file \""..slave_file_name.."\"")
		return 0;
	end

	line = in_file:read("L")
	while (line ~= nil) do

		incl_statement = string.find(line,"#include",0,true)
		if (incl_statement ~= nil) then

			name_start  = string.find(line,"\"",0,true)
			name_end  = string.find(line,"\"",name_start+1,true)

			new_slave_name = string.sub(line,name_start+1,name_end-1)
			master_file:write("\n\n// Consolidating "..new_slave_name)
			master_file:write("\n")

			IncludeIntoMaster(master_file,new_slave_name)
		else
			master_file:write(line)
		end

		line = in_file:read("L")
	end
	in_file:close()
end

-- ============================================== Consolidate registers
--[[
This function consolidates all registers into a single one. This is
mostly because registers are c macro-fied files that link together with
#include statements
]]--
function ConsolidateRegisters()
	print("Consolidating lua register")
	consolidated_file_name = "../../CHI_DOC/LUA_DOCUMENTATION/lua_register.txt"
	consolidated_file = io.open(consolidated_file_name,"w+");

	IncludeIntoMaster(consolidated_file,"chi_lua_register.h")

	consolidated_file:close()
	print("Done consolidating lua register")
end

newFile = io.open("../../CHI_DOC/LUA_DOCUMENTATION/lua_functions.c","w");
newFile2 = io.open("../../CHI_DOC/LUA_DOCUMENTATION/lua_namespace.hpp","w");
newFile2:write("namespace chi_lua \n {\n");
readCXXFile()
newFile2:write("}\n");
newFile2:close();
newFile:close();

ConsolidateRegisters()
generateFunctionList()

