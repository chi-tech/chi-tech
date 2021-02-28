#include "material_property_transportxsections.h"

//###################################################################
/**Pushes all of the relevant items of the transport xs to a lua table.*/
void chi_physics::TransportCrossSections::PushLuaTable(lua_State *L)
{

  lua_newtable(L);
  lua_pushstring(L,"is_empty");
  lua_pushboolean(L,false);
  lua_settable(L,-3);

  lua_pushstring(L,"G");
  lua_pushnumber(L,G);
  lua_settable(L,-3);

  lua_pushstring(L,"L");
  lua_pushnumber(L,this->L);
  lua_settable(L,-3);

  int g=0;

  lua_pushstring(L,"sigma_tg");
  lua_newtable(L);
  g=0;
  for (auto val : sigma_tg)
  {
    ++g;
    lua_pushnumber(L,g);
    lua_pushnumber(L,val);
    lua_settable(L,-3);
  }
  lua_settable(L,-3);

  lua_pushstring(L,"sigma_fg");
  lua_newtable(L);
  g=0;
  for (auto val : sigma_fg)
  {
    ++g;
    lua_pushnumber(L,g);
    lua_pushnumber(L,val);
    lua_settable(L,-3);
  }
  lua_settable(L,-3);

  lua_pushstring(L,"sigma_captg");
  lua_newtable(L);
  g=0;
  for (auto val : sigma_captg)
  {
    ++g;
    lua_pushnumber(L,g);
    lua_pushnumber(L,val);
    lua_settable(L,-3);
  }
  lua_settable(L,-3);

  lua_pushstring(L,"chi_g");
  lua_newtable(L);
  g=0;
  for (auto val : chi_g)
  {
    ++g;
    lua_pushnumber(L,g);
    lua_pushnumber(L,val);
    lua_settable(L,-3);
  }
  lua_settable(L,-3);

  lua_pushstring(L,"nu_sigma_fg");
  lua_newtable(L);
  g=0;
  for (auto val : nu_sigma_fg)
  {
    ++g;
    lua_pushnumber(L,g);
    lua_pushnumber(L,val);
    lua_settable(L,-3);
  }
  lua_settable(L,-3);

  lua_pushstring(L,"nu_p_sigma_fg");
  lua_newtable(L);
  g=0;
  for (auto val : nu_p_sigma_fg)
  {
    ++g;
    lua_pushnumber(L,g);
    lua_pushnumber(L,val);
    lua_settable(L,-3);
  }
  lua_settable(L,-3);

  lua_pushstring(L,"nu_d_sigma_fg");
  lua_newtable(L);
  g=0;
  for (auto val : nu_d_sigma_fg)
  {
    ++g;
    lua_pushnumber(L,g);
    lua_pushnumber(L,val);
    lua_settable(L,-3);
  }
  lua_settable(L,-3);

  lua_pushstring(L,"diffg");
  lua_newtable(L);
  g=0;
  for (auto val : diffg)
  {
    ++g;
    lua_pushnumber(L,g);
    lua_pushnumber(L,val);
    lua_settable(L,-3);
  }
  lua_settable(L,-3);

  lua_pushstring(L,"sigma_rg");
  lua_newtable(L);
  g=0;
  for (auto val : sigma_rg)
  {
    ++g;
    lua_pushnumber(L,g);
    lua_pushnumber(L,val);
    lua_settable(L,-3);
  }
  lua_settable(L,-3);

  lua_pushstring(L,"sigma_ag");
  lua_newtable(L);
  g=0;
  for (auto val : sigma_ag)
  {
    ++g;
    lua_pushnumber(L,g);
    lua_pushnumber(L,val);
    lua_settable(L,-3);
  }
  lua_settable(L,-3);

  lua_pushstring(L,"sigma_s_gtog");
  lua_newtable(L);
  g=0;
  for (auto val : sigma_s_gtog)
  {
    ++g;
    lua_pushnumber(L,g);
    lua_pushnumber(L,val);
    lua_settable(L,-3);
  }
  lua_settable(L,-3);

  lua_pushstring(L,"xi_Jfull_g");
  lua_newtable(L);
  g=0;
  for (auto val : xi_Jfull_g)
  {
    ++g;
    lua_pushnumber(L,g);
    lua_pushnumber(L,val);
    lua_settable(L,-3);
  }
  lua_settable(L,-3);

  lua_pushstring(L,"xi_Jpart_g");
  lua_newtable(L);
  g=0;
  for (auto val : xi_Jpart_g)
  {
    ++g;
    lua_pushnumber(L,g);
    lua_pushnumber(L,val);
    lua_settable(L,-3);
  }
  lua_settable(L,-3);

  lua_pushstring(L,"D_jfull");
  lua_pushnumber(L,D_jfull);
  lua_settable(L,-3);

  lua_pushstring(L,"D_jpart");
  lua_pushnumber(L,D_jpart);
  lua_settable(L,-3);

  lua_pushstring(L,"sigma_a_jfull");
  lua_pushnumber(L,sigma_a_jfull);
  lua_settable(L,-3);

  lua_pushstring(L,"sigma_a_jpart");
  lua_pushnumber(L,sigma_a_jpart);
  lua_settable(L,-3);
}