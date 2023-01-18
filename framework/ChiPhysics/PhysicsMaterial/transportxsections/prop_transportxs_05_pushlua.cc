#include "material_property_transportxsections.h"

//###################################################################
/**Pushes all of the relevant items of the transport xs to a lua table.*/
void chi_physics::TransportCrossSections::PushLuaTable(lua_State *L)
{
  lua_newtable(L);
  lua_pushstring(L, "is_empty");
  lua_pushboolean(L, false);
  lua_settable(L, -3);

  lua_pushstring(L, "num_groups");
  lua_pushinteger(L, static_cast<lua_Integer>(num_groups));
  lua_settable(L, -3);

  lua_pushstring(L, "scattering_order");
  lua_pushinteger(L, static_cast<lua_Integer>(scattering_order));
  lua_settable(L, -3);

  lua_pushstring(L, "num_precursors");
  lua_pushinteger(L, static_cast<lua_Integer>(num_precursors));
  lua_settable(L, -3);

  lua_pushstring(L, "is_fissionable");
  lua_pushboolean(L, is_fissionable);
  lua_settable(L, -3);

  auto Push1DXS =
      [L](const std::vector<double>& xs,
          const std::string& name)
  {
    lua_pushstring(L, name.c_str());
    lua_newtable(L);
    {
      unsigned int g = 0;
      for (const auto& val : xs)
      {
        lua_pushinteger(L, g++);
        lua_pushnumber(L ,val);
        lua_settable(L, -3);
      }
    }
    lua_settable(L, -3);
  };

  Push1DXS(sigma_t, "sigma_t");
  Push1DXS(sigma_a, "sigma_a");
  Push1DXS(sigma_f, "sigma_f");
  Push1DXS(nu_sigma_f, "nu_sigma_f");
  Push1DXS(nu_prompt_sigma_f, "nu_prompt_sigma_f");
  Push1DXS(nu_delayed_sigma_f, "nu_delayed_sigma_f");
  Push1DXS(inv_velocity, "inv_velocity");

  lua_pushstring(L, "chi_delayed");
  lua_newtable(L);
  {
    for (unsigned int g = 0; g < num_groups; ++g)
    {
      lua_pushinteger(L, g + 1);
      lua_newtable(L);
      {
        for (unsigned int j = 0; j < num_precursors; ++j)
        {
          lua_pushinteger(L, j + 1);
          lua_pushnumber(L, precursors[j].emission_spectrum[g]);
          lua_settable(L, -3);
        }
      }
      lua_settable(L, -3);
    }//for g
  }
  lua_settable(L,-3);

  lua_pushstring(L,"precursor_decay_constants");
  lua_newtable(L);
  {
    unsigned int j = 0;
    for (const auto& precursor : precursors)
    {
      lua_pushinteger(L, j++);
      lua_pushnumber(L, precursor.decay_constant);
      lua_settable(L, -3);
    }
  }
  lua_settable(L, -3);

  lua_pushstring(L, "precursor_fractional_yields");
  lua_newtable(L);
  {
    unsigned int j = 0;
    for (const auto& precursor : precursors)
    {
      lua_pushinteger(L, j++);
      lua_pushnumber(L, precursor.fractional_yield);
      lua_settable(L, -3);
    }
  }
  lua_settable(L, -3);

  //============================================= Transfer matrices
  lua_pushstring(L, "transfer_matrix");
  lua_newtable(L);
  {
    unsigned int ell = 0;
    for (const auto& matrix : transfer_matrices)
    {
      ++ell;
      lua_pushinteger(L, ell);
      lua_newtable(L);
      {
        for (unsigned int g=0; g < matrix.NumRows(); ++g)
        {
          const auto& col_indices = matrix.rowI_indices[g];
          const auto& col_values  = matrix.rowI_values[g];
          size_t num_vals = col_values.size();

          lua_pushinteger(L, g + 1);
          lua_newtable(L);
          {
            for (unsigned int gg = 0; gg < num_vals; ++gg)
            {
              lua_pushinteger(L, static_cast<long long>(col_indices[gg]) + 1);
              lua_pushnumber(L, col_values[gg]);
              lua_settable(L, -3);
            }
            lua_settable(L, -3);
          }
        }//for g
      }
      lua_settable(L, -3);
    }
  }
  lua_settable(L, -3);


  lua_pushstring(L, "production_matrix");
  lua_newtable(L);
  {
    for (unsigned int g = 0; g < num_groups; ++g)
    {
      const auto& prod = production_matrix[g];

      lua_pushinteger(L, g + 1);
      lua_newtable(L);
      {
        for (unsigned int gp = 0; gp < num_groups; ++gp)
        {
          lua_pushinteger(L, gp + 1);
          lua_pushnumber(L, prod[gp]);
          lua_settable(L, -3);
        }
        lua_settable(L, -3);
      }
    }//for g
  }
  lua_settable(L, -3);

  //============================================= Diffusion quantities
  Push1DXS(diffusion_coeff, "diffusion_coeff");
  Push1DXS(sigma_removal, "sigma_removal");
  Push1DXS(sigma_s_gtog, "sigma_s_gtog");
}
