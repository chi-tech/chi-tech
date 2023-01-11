cd ../

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
