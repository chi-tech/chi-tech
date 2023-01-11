--[[
To keep the lua code a little bit cleaner a simple function was
written which takes a table and outputs the size of the table.

** Note:
Used for testing purposes.
]]--
function TableSize(tTable)
	local counter = 0;
	for _ in pairs(tTable) do counter = counter + 1; end
	return counter;
end

--[[
Searches directories for .cc and .cpp files and adds these
file names to a list.

** Note:
The function already contains the absolute path.
]]--
function ScanDirectories(dirs_to_scan)
	print("ScanDirectories: Generating list of cc and cpp files.");

	local fIndex, fTable = 0, {};
	local fFile;

	--================================= Make a clean list
	os.execute("mkdir LUA_DOCUMENTATION")
	os.execute(">LUA_DOCUMENTATION/Z_GeneratedFileList.txt")

	--================================= Scan each directory
    if (dirs_to_scan.itemCount>0) then
        for folderDex=1,dirs_to_scan.itemCount do
			os.execute("find \""..dirs_to_scan[folderDex].."\" -name \"*.cpp\""..
					   " >> LUA_DOCUMENTATION/Z_GeneratedFileList.txt");
			os.execute("find \""..dirs_to_scan[folderDex].."\" -name \"*.cc\""..
					   " >> LUA_DOCUMENTATION/Z_GeneratedFileList.txt");
        end

    end

	--================================= Convert file-list to table
	fFile = io.open("LUA_DOCUMENTATION/Z_GeneratedFileList.txt","r");

	for filename in fFile:lines() do
		fIndex = fIndex + 1;
		fTable[fIndex] = filename;
	end

	fFile:close();

	return fTable
end

--[[
Function that fractures strings into words, and extracts the parameter
followed by the value.

This expects a parameter line to look like this:
   \param ParamName ParamType Description

The returned table is [ParamType, ParamName]
]]--
function ExtractParameter(fString)
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
Function that scans through a cxx file looking for certain
keywords/identifiers that will identify a lua wrapper function.
]]--
function ProcessCXXFileLuaFunctions(fileContent)
	local inComment, inFunction, inParam  	           = false, false, false;
	local function_line, comment_string, param_strings = "", "", {};

	for line in string.gmatch(fileContent, "[^\r\n]+") do

		--============================ Searching for doxygen style comments
		--Looking for opening comment block /**
		if (string.find(line, "/%*%*") ~= nil) then
			inComment 		= true;
			comment_string = ""
			param_strings = {};
		end

		if (inComment == true) then
			comment_string 	= comment_string..line.."\n";
		end

		--Looking for closing comment block */
		if (string.find(line, "%*/") ~= nil) then
			inComment 		= false;
		end


		--============================ Searching for doxygen param
		--Looking for the \param keyword
		if (string.find(line, "\\param") ~= nil) then
			inParam 		= true;
		end

		if (inParam == true) then
			params = ExtractParameter(line);
			table.insert(param_strings, params[1].." "..params[2])
            inParam 		= false;
		end

		if (string.find(line, "|") ~= nil) then
			inParam 		= false;
		end

		--============================ Searching for function name
		-- Looking for keywords chi and lua_State
		if ((string.find(line,"lua_State") ~= nil) and
			(string.find(line,"chi") ~= nil)) then
			inFunction 		= true;
		end

		if (inFunction == true) then
			local fParam = ""
			tlen = TableSize(param_strings)
			if (tlen > 0) then
				for k=1,tlen-1 do
					fParam = fParam..param_strings[k]..", "
				end
				fParam = fParam..param_strings[tlen]
			end

			local lBegin = string.find(line,"chi");

			function_line = string.sub(line,lBegin);

			local cBegin = string.find(function_line,"%(");
			local fName = string.sub(function_line,1,cBegin-1);

			inFunction 		= false;

			function_definition=""
			function_declaration=""
			if (fName ~= "") then
				function_definition = comment_string..
					"int chi_lua::"..fName.."("..fParam..")\n".."{return;} \n";
				function_declaration = "int ".. fName .."("..fParam..");\n";
			end


			definition_file:write(function_definition);
			declaration_file:write(function_declaration);
		end
	end
end

--[[
Function that reads in the file into memory and loads it into a table for parsing.
]]--
function ProcessCXXFiles(dirs_to_scan)
	local fTable, fName, fContent;

	fTable = ScanDirectories(dirs_to_scan);

	for _,value in pairs(fTable) do
		fName = io.open(value, 'r+');
		fContent = fName:read("*all");
		print("Processing: "..value)
		ProcessCXXFileLuaFunctions(fContent);
		fName:close();
	end

	print("=============== Finished Reading lua files ===============")
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
function GenerateFunctionList()
	print("Generating function list for main page")
	main_page_file = io.open("../../doc/PAGES/MainPage.h","w");
	headFile = io.open("../../doc/PAGES/mp_header.txt","r")
	reg = io.open("../../doc/LUA_DOCUMENTATION/lua_register.txt","r");

	--=============================== Write header
	line = headFile:read("L")
	while (line ~= nil) do
		main_page_file:write(line)
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
				WriteTableHTML(main_page_file, table_name, table_ref, table_vals)
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
		WriteTableHTML(main_page_file, table_name, table_ref, table_vals)
	end
	--======================================== END OF WRITES A TABLE

	if (strword ~= nil) then
		main_page_file:write(string.sub(strline,10).."\n")
	end

	main_page_file:write("\n\n*/\n")

	headFile:close()
	reg:close();
	main_page_file:close();

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
	consolidated_file_name = "../../doc/LUA_DOCUMENTATION/lua_register.txt"
	consolidated_file = io.open(consolidated_file_name,"w+");

	IncludeIntoMaster(consolidated_file,"chi_lua_register.h")

	consolidated_file:close()
	print("Done consolidating lua register")
end

dofile("module_lua_inclusion.lua") --Defines moduleFolders

definition_file = io.open("../../doc/LUA_DOCUMENTATION/lua_functions.c","w");
declaration_file = io.open("../../doc/LUA_DOCUMENTATION/lua_namespace.hpp","w");
declaration_file:write("namespace chi_lua \n {\n");
ProcessCXXFiles(moduleFolders)
declaration_file:write("}\n");
declaration_file:close();
definition_file:close();

ConsolidateRegisters()
GenerateFunctionList()
