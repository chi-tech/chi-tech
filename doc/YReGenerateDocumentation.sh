cd ../

./bin/ChiTech --suppress_color --dump-object-registry -v 2 > doc/generated_files/test.txt

python3 "doc/scripts/MakeInputParamsDocs.py"

#============================== Transform lua wrappers for documentation
# Lua wrapper functions are normally int chiFunction(lua_State* L)
# The scipt below uses their doc-strings to transform to doxy style
# functions
python3 "doc/scripts/BuildListOfLuaWrappers.py"

#============================== Build mainpage quick reference
python3 "doc/scripts/MakeMainPage.py"

#============================== Making main documentation
echo "Running DoxyGen"
doxygen "doc/DoxyfileLua"
