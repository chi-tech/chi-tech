#define RegisterFunction(x) \
        int x(lua_State *L); \
        lua_register(L, #x, x);

#define RegisterConstant(x,y) \
        lua_pushnumber(L,y); \
        lua_setglobal(L, #x);

#define RegisterNamespace(x) \
        lua_newtable(L); \
        lua_setglobal(L, #x);

#define AddNamedConstantToNamespace(const_name,const_value,namespace_name) \
        lua_getglobal(L,#namespace_name); \
        lua_pushstring(L,#const_name); \
        lua_pushnumber(L,const_value); \
        lua_settable(L,-3);
