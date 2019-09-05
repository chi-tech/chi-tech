cd ../

LUA="$PWD/CHI_RESOURCES/Dependencies/lua-5.3.5/install/bin/lua"

#============================== Making lua documentation
cd CHI_TECH/CHI_LUA
${LUA} chi_lua_docbuild.lua
cd ../..

#============================== Making chil documentation
cd CHI_RESOURCES/Scripts/chil
${LUA} Z0_MakeDoxy.lua
cd ../../..

#============================== Making main documentation
doxygen "CHI_DOC/DoxyfileLua"
cd CHI_BUILD
