cd ../

LUA="lua"

python3 "doc/scripts/BuildListOfLuaWrappers.py"
python3 "doc/scripts/MakeMainPage.py"

##============================== Making lua documentation
#cd framework/ChiLua
#${LUA} chi_lua_docbuild.lua
#cd ../..
#
#============================== Making main documentation
echo "Running DoxyGen"
doxygen "doc/DoxyfileLua"
