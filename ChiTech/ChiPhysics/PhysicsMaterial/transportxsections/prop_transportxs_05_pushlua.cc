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
  lua_pushinteger(L,static_cast<lua_Integer>(num_groups));
  lua_settable(L,-3);

  lua_pushstring(L,"L");
  lua_pushinteger(L,static_cast<lua_Integer>(scattering_order));
  lua_settable(L,-3);

  lua_pushstring(L,"J");
  lua_pushinteger(L,static_cast<lua_Integer>(num_precursors));
  lua_settable(L,-3);

  lua_pushstring(L,"is_fissile");
  lua_pushboolean(L,is_fissile);
  lua_settable(L,-3);

  auto Push1DXS = [L](const std::vector<double>& xs, const std::string& name)
  {
    lua_pushstring(L,name.c_str());
    lua_newtable(L);
    {
      int g=0;
      for (auto val : xs)
      {
        ++g;
        lua_pushinteger(L,g);
        lua_pushnumber(L,val);
        lua_settable(L,-3);
      }
    }
    lua_settable(L,-3);
  };

  Push1DXS(sigma_tg,"sigma_tg");
  Push1DXS(sigma_fg,"sigma_fg");
  Push1DXS(sigma_ag,"sigma_ag");
  Push1DXS(chi_g,"chi_g");
  Push1DXS(nu,"nu");
  Push1DXS(nu_prompt,"nu_prompt");
  Push1DXS(nu_delayed,"nu_delayed");
  Push1DXS(nu_sigma_fg,"nu_sigma_fg");
  Push1DXS(nu_p_sigma_fg,"nu_p_sigma_fg");
  Push1DXS(nu_d_sigma_fg,"nu_d_sigma_fg");
  Push1DXS(ddt_coeff,"ddt_coeff");

  lua_pushstring(L,"chi_d");
  lua_newtable(L);
  {
    int g = 0;
    for (auto& row : chi_d)
    {
      ++g;
      lua_pushinteger(L, g);
      lua_newtable(L);
        int j=0;
        for (auto val : row)
        {
          ++j;
          lua_pushinteger(L,j);
          lua_pushnumber(L,val);
          lua_settable(L,-3);
        }
      lua_settable(L, -3);
    }
  }
  lua_settable(L,-3);

  lua_pushstring(L,"lambda");
  lua_newtable(L);
  {
    int j = 0;
    for (auto val : lambda)
    {
      ++j;
      lua_pushinteger(L, j);
      lua_pushnumber(L, val);
      lua_settable(L, -3);
    }
  }
  lua_settable(L,-3);

  lua_pushstring(L,"gamma");
  lua_newtable(L);
  {
    int j = 0;
    for (auto val : gamma)
    {
      ++j;
      lua_pushinteger(L, j);
      lua_pushnumber(L, val);
      lua_settable(L, -3);
    }
  }
  lua_settable(L,-3);

  //============================================= Transfer matrices
  lua_pushstring(L,"transfer_matrix");
  lua_newtable(L);
  {
    int ell = 0;
    for (const auto& matrix : transfer_matrix)
    {
      ++ell;
      lua_pushinteger(L, ell);
      lua_newtable(L);
      {
        for (int g=0; g<matrix.NumRows(); ++g)
        {
          const auto& col_indices = matrix.rowI_indices[g];
          const auto& col_values  = matrix.rowI_values[g];

          size_t num_vals = col_values.size();
          lua_pushinteger(L,g+1);
          lua_newtable(L);
          for (int gg=0; gg<num_vals; ++gg)
          {
            lua_pushinteger(L, static_cast<long long>(col_indices[gg])+1);
            lua_pushnumber(L, col_values[gg]);
            lua_settable(L, -3);
          }
          lua_settable(L, -3);
        }//for g
      }
      lua_settable(L, -3);
    }
  }
  lua_settable(L,-3);


  //============================================= Diffusion quantities
  Push1DXS(diffg,"diffg");
  Push1DXS(sigma_rg,"sigma_rg");
  Push1DXS(sigma_s_gtog,"sigma_s_gtog");
  Push1DXS(xi_Jfull_g,"xi_Jfull_g");
  Push1DXS(xi_Jpart_g,"xi_Jpart_g");

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