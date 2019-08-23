cd ../

#============================== Making lua documentation
cd CHI_TECH/CHI_LUA
lua chi_lua_docbuild.lua
cd ../..

#============================== Making chil documentation
cd CHI_RESOURCES/Scripts/chil
lua Z0_MakeDoxy.lua
cd ../../..

#============================== Making main documentation
doxygen "CHI_DOC/DoxyfileLua"
cd CHI_BUILD
