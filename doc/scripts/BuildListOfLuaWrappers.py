"""
This script walks the modules- and framework-folders and
looks for directories named "lua" or "LuaTest". When found
it adds all the .cc files in those folders to a list.
"""
import os
import re

folders_to_walk = ["framework", "modules"]

wrapper_sources = []
for search_folder in folders_to_walk:
    for dir_path, sub_dirs, files in os.walk(search_folder):
        dir_name = os.path.basename(dir_path)
        if dir_name == "lua" or dir_name == "LuaTest":
            for file_name in files:
                base_name, extension = os.path.splitext(file_name)
                if extension == ".cc":
                    print(dir_path + "/" + file_name)
                    wrapper_sources.append(dir_path + "/" + file_name)

"""
Next each file is scanned. Doxygen comments are extracted as well
as functions with the signature `int chiFunction(lua_State* L)`. The 
latter is identified with 3 items, i.e., the `int`, the prefix `chi` and
the presence of `lua_State`.

The functions only (without the doxy comments) are wrapped in a
namespace (chi_lua) and written to 
`doc/generated_files/lua_namespace.h` (the declaration file).

The doxy comments AND the functions are written to 
`doc/generated_files/lua_functions.c` (the definition file). 

For the definitions the doxy comments of a function is scanned for the presence 
of `param` and `return`. If `param` is encountered, as word0, then word1 is 
taken as the parameter name and word2 is taken as the parameter type. 
If `return` is encountered, as word0, then word1 is taken as the return type.

The function (declaration and definition) is then reformatted when written, such
that the return type replaces the classical `int` return and the parameters
replace the `lua_State* L` parameter.
"""
definition_file = open("doc/generated_files/lua_functions.c", "w")
declaration_file = open("doc/generated_files/lua_namespace.h", "w")
declaration_file.write("namespace chi_lua\n{\n")
for src_path in wrapper_sources:
    file = open(src_path, "r")

    definition_file.write("//" + src_path + "\n")
    lines = file.readlines()

    in_comment = False
    params = []
    returns = []
    for line in lines:

        if line.find("/**") >= 0:
            in_comment = True

        if in_comment:
            definition_file.write(line)

        if in_comment:
            param_start = line.find(r"\param")
            if param_start >= 0:
                words = line[param_start:].split()
                if len(words) >= 3:
                    params.append(words[2] + " " + words[1])
            return_start = line.find(r"\return")
            if return_start >= 0:
                words = line[return_start:].split()
                if len(words) >= 2:
                    returns.append(words[1])

        if line.find("*/") >= 0:
            in_comment = False

        if line.find("int") >= 0 and line.find("chi") >= 0 and \
                line.find("lua_State") >= 0:
            words = re.split(r"\(|\)|\s", line.strip())
            print(src_path, words)

            function_name = words[1]

            # Defintion file
            if len(returns) > 0:
                definition_file.write(returns[0] + " ")
            else:
                definition_file.write("void ")

            definition_file.write("chi_lua::" + function_name + "(")
            for k in range(0, len(params)):
                definition_file.write(params[k])
                if k < (len(params)-1):
                    definition_file.write(", ")
            definition_file.write("){}\n")

            # Declaration file
            if len(returns) > 0:
                declaration_file.write(returns[0] + " ")
            else:
                declaration_file.write("void ")

            declaration_file.write(function_name + "(")
            for k in range(0, len(params)):
                declaration_file.write(params[k])
                if k < (len(params)-1):
                    declaration_file.write(", ")

            declaration_file.write(");\n")

            params = []
            returns = []

    file.close()

definition_file.close()

declaration_file.write("}//namespace chi_lua\n")
declaration_file.close()
