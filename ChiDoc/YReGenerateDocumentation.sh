cd ../

LUA="lua"

#============================== Making lua documentation
cd ChiTech/ChiLua
${LUA} chi_lua_docbuild.lua
cd ../..

#============================== Making chil documentation
cd ChiResources/Scripts/chil
${LUA} Z0_MakeDoxy.lua
cd ../../..

#============================== Making main documentation
doxygen "ChiDoc/DoxyfileLua"
cd CHI_BUILD
