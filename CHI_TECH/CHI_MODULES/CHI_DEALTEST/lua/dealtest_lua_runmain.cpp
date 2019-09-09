#include "../../../ChiLua/chi_lua.h"

int dealtest();

int chiRunDealTest(lua_State *L)
{
    dealtest();

    return 0;
}
