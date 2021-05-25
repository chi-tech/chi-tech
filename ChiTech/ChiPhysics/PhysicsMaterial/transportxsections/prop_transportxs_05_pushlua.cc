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

  lua_pushstring(L,"sigma_tg");
  lua_newtable(L);
  {
    int g=0;
    for (auto val : sigma_tg)
    {
      ++g;
      lua_pushinteger(L,g);
      lua_pushnumber(L,val);
      lua_settable(L,-3);
    }
  }
  lua_settable(L,-3);

  lua_pushstring(L,"sigma_fg");
  lua_newtable(L);
  {
    int g=0;
    for (auto val : sigma_fg)
    {
      ++g;
      lua_pushinteger(L,g);
      lua_pushnumber(L,val);
      lua_settable(L,-3);
    }
  }
  lua_settable(L,-3);

  lua_pushstring(L,"sigma_captg");
  lua_newtable(L);
  {
    int g=0;
    for (auto val : sigma_captg)
    {
      ++g;
      lua_pushinteger(L,g);
      lua_pushnumber(L,val);
      lua_settable(L,-3);
    }
  }
  lua_settable(L,-3);

  lua_pushstring(L,"chi_g");
  lua_newtable(L);
  {
    int g = 0;
    for (auto val : chi_g)
    {
      ++g;
      lua_pushinteger(L, g);
      lua_pushnumber(L, val);
      lua_settable(L, -3);
    }
  }
  lua_settable(L,-3);

  lua_pushstring(L,"nu_sigma_fg");
  lua_newtable(L);
  {
    int g = 0;
    for (auto val : nu_sigma_fg)
    {
      ++g;
      lua_pushinteger(L, g);
      lua_pushnumber(L, val);
      lua_settable(L, -3);
    }
  }
  lua_settable(L,-3);

  lua_pushstring(L,"nu_p_sigma_fg");
  lua_newtable(L);
  {
    int g = 0;
    for (auto val : nu_p_sigma_fg)
    {
      ++g;
      lua_pushinteger(L, g);
      lua_pushnumber(L, val);
      lua_settable(L, -3);
    }
  }
  lua_settable(L,-3);

  lua_pushstring(L,"nu_d_sigma_fg");
  lua_newtable(L);
  {
    int g = 0;
    for (auto val : nu_d_sigma_fg)
    {
      ++g;
      lua_pushinteger(L, g);
      lua_pushnumber(L, val);
      lua_settable(L, -3);
    }
  }
  lua_settable(L,-3);

  lua_pushstring(L,"ddt_coeff");
  lua_newtable(L);
  {
    int g = 0;
    for (auto val : ddt_coeff)
    {
      ++g;
      lua_pushinteger(L, g);
      lua_pushnumber(L, val);
      lua_settable(L, -3);
    }
  }
  lua_settable(L,-3);

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
  lua_pushstring(L,"diffg");
  lua_newtable(L);
  {
    int g = 0;
    for (auto val : diffg)
    {
      ++g;
      lua_pushinteger(L, g);
      lua_pushnumber(L, val);
      lua_settable(L, -3);
    }
  }
  lua_settable(L,-3);

  lua_pushstring(L,"sigma_rg");
  lua_newtable(L);
  {
    int g = 0;
    for (auto val : sigma_rg)
    {
      ++g;
      lua_pushinteger(L, g);
      lua_pushnumber(L, val);
      lua_settable(L, -3);
    }
  }
  lua_settable(L,-3);

  lua_pushstring(L,"sigma_ag");
  lua_newtable(L);
  {
    int g = 0;
    for (auto val : sigma_ag)
    {
      ++g;
      lua_pushinteger(L, g);
      lua_pushnumber(L, val);
      lua_settable(L, -3);
    }
  }
  lua_settable(L,-3);

  lua_pushstring(L,"sigma_s_gtog");
  lua_newtable(L);
  {
    int g = 0;
    for (auto val : sigma_s_gtog)
    {
      ++g;
      lua_pushinteger(L, g);
      lua_pushnumber(L, val);
      lua_settable(L, -3);
    }
  }
  lua_settable(L,-3);

  lua_pushstring(L,"xi_Jfull_g");
  lua_newtable(L);
  {
    int g = 0;
    for (auto val : xi_Jfull_g)
    {
      ++g;
      lua_pushinteger(L, g);
      lua_pushnumber(L, val);
      lua_settable(L, -3);
    }
  }
  lua_settable(L,-3);

  lua_pushstring(L,"xi_Jpart_g");
  lua_newtable(L);
  {
    int g = 0;
    for (auto val : xi_Jpart_g)
    {
      ++g;
      lua_pushinteger(L, g);
      lua_pushnumber(L, val);
      lua_settable(L, -3);
    }
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