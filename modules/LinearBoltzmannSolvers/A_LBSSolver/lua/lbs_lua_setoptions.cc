#include "lbs_lua_utils.h"
#include "chi_lua.h"

#include "LinearBoltzmannSolvers/A_LBSSolver/lbs_solver.h"

#include "ChiConsole/chi_console.h"
#include "chi_runtime.h"
#include "chi_log.h"

namespace lbs::common_lua_utils
{

RegisterLuaFunctionAsIs(chiLBSSetOptions);

/**Sets an LBS option.
 *
 * \param handle int Handle to the lbs-based object.
 * \param specs Table Various options in a table.

##_

### Example usage
Example:
\code
chiLBSSetOptions(phys1,
{
  spatial_discretization = "pwld",
  scattering_order = 2
})
\endcode

### Table keys
`spatial_discretization`\n
<I>type=</I><span style="color: blue;"><TT>STRING</TT></span>
What spatial discretization to use. Currently only `"pwld"` is supported.
Default: `"pwld"`.\n\n

`scattering_order`\n
<I>type=</I><span style="color: blue;"><TT>INTEGER</TT></span>
Defines the level of harmonic expansion for the scattering source. Default: `1`.
\n\n

`sweep_eager_limit`\n
<I>type=</I><span style="color: blue;"><TT>INTEGER</TT></span>
The eager limit to be used in message size during sweep initialization.
 This expects to be followed by a size in bytes (Max 64,0000).Default: `32,000`.
 See note below.\n\n

 ###Note on the Eager limit
The eager limit is the message size limit before which non-blocking MPI send
calls will execute without waiting for a matching receive call. The limit is
platform dependent but in general 64 kb. Some systems have 32 kb as a limit
and therefore we use that as a default limit in ChiTech. There is a fine
interplay between message size and the shear amount of messages that will be
sent. In general smaller messages tend to be more efficient, however, when
there are too many small messages being sent around the communication system
on the given platform will start to suffer. One can gain a small amount of
parallel efficiency by lowering this limit, however, there is a point where
the parallel efficiency will actually get worse so use with caution.
\n\n

`read_restart_data`\n
<I>type=</I><span style="color: blue;"><TT>BOOLEAN</TT></span>
Flag indicating whether restart data is to be read. Default: `false`.\n\n

`read_restart_folder_name`\n
<I>type=</I><span style="color: blue;"><TT>STRING</TT></span>
Folder name to use when reading restart data. Default: `"YRestart"`.\n\n

`read_restart_file_base`\n
<I>type=</I><span style="color: blue;"><TT>STRING</TT></span>
File base name to use when reading restart data. Default: `"restart"`.\n\n

`write_restart_data`\n
<I>type=</I><span style="color: blue;"><TT>BOOLEAN</TT></span>
Flag indicating whether restart data is to be written. Default: `false`.\n\n

`write_restart_folder_name`\n
<I>type=</I><span style="color: blue;"><TT>STRING</TT></span>
Folder name to use when writing restart data. Default: `"YRestart"`.\n\n

`write_restart_file_base`\n
<I>type=</I><span style="color: blue;"><TT>STRING</TT></span>
File base name to use when writing restart data. Default: `"restart"`.\n\n

`write_restart_interval`\n
<I>type=</I><span style="color: blue;"><TT>FLOAT</TT></span>
Interval at which restart data is to be writtin. Default: `30.0`. Currently
not implemented.\n\n

`use_precursors`\n
<I>type=</I><span style="color: blue;"><TT>BOOLEAN</TT></span>
Flag for using delayed neutron precursors. Default: `false`.\n\n

`use_source_moments`\n
<I>type=</I><span style="color: blue;"><TT>BOOLEAN</TT></span>
Flag for ignoring fixed sources and selectively using source moments
obtained elsewhere. Default: `false`.\n\n

`save_angular_flux`\n
<I>type=</I><span style="color: blue;"><TT>BOOLEAN</TT></span>
Flag indicating whether angular fluxes are to be stored or not. Default:
`false`. \n\n

`verbose_inner_iterations`\n
<I>type=</I><span style="color: blue;"><TT>BOOLEAN</TT></span>
Flag to control verbosity of inner iterations. Default: `true`.\n\n

`verbose_ags_iterations`\n
<I>type=</I><span style="color: blue;"><TT>BOOLEAN</TT></span>
Flag to control verbosity of across-groupset iterations. Default: `false`.\n\n

`verbose_outer_iterations`\n
<I>type=</I><span style="color: blue;"><TT>BOOLEAN</TT></span>
Flag to control verbosity of outer iterations. Default: `true`.\n\n

`power_field_function_on`\n
<I>type=</I><span style="color: blue;"><TT>BOOLEAN</TT></span>
Flag to control the creation of the power generation field function. If set to
`true` then a field function will be created with the general name
`<solver_name>_power_generation`.
Default: `false`.\n\n

`power_default_kappa`\n
<I>type=</I><span style="color: blue;"><TT>FLOAT</TT></span>
Default `kappa` value (Energy released per fission) to use for power generation
when cross sections do not have `kappa` values.
Default: `3.20435e-11` (Joules in 200 MeV).\n\n

`power_normalization`\n
<I>type=</I><span style="color: blue;"><TT>FLOAT</TT></span>
Power normalization factor to use. Supply a negative or zero number to turn
this off. Default: `-1.0` (Disabled).\n\n

`field_function_prefix_option`\n
<I>type=</I><span style="color: blue;"><TT>STRING</TT></span>
Prefix option on field function names. Default: `"prefix"`. Can be `"prefix"` or
`"solver_name"`.
By default this option is `"prefix"` which means it uses the designated
"prefix" (another option), however, that is defaulted to nothing. Therefore,
default behavior is export flux moment fields functions as `phi_gXXX_mYYY`
where `XXX` is the zero padded 3 digit group number and similarly for `YYY`.\n\n

`field_function_prefix`\n
<I>type=</I><span style="color: blue;"><TT>STRING</TT></span>
Prefix to use on all field functions. Default: `""`.
By default this option is empty but if specified then flux moments will exported
as `prefix_phi_gXXX_mYYY` where `XXX` is the zero padded 3 digit group number
and similarly for `YYY`. The underscore after "prefix" is added automatically.
\n\n

`boundary_conditions`\n
<I>type=</I><span style="color: blue;"><TT>BLOCK</TT></span>
A table contain sub-tables for each boundary specification.\n

Each sub-block:
- `name` <I>type=</I><span style="color: blue;"><TT>STRING</TT></span>. Name
 of the boundary as it exists in the mesh.
- `type` <I>type=</I><span style="color: blue;"><TT>STRING</TT></span>. Type
 name of the boundary. Can be any of the following:
  - `"vacuum"`
  - `"incident_isotropic"`
  - `"reflecting"`
  - `"incident_anisotropic_heterogeneous"`
- `group_strength` <I>type=</I><span style="color: blue;"><TT>ARRAY</TT></span>.
  Required only if `type` is `"incident_isotropic"`. An array of isotropic
  strength per group.
- `function_name` <I>type=</I><span style="color: blue;"><TT>STRING</TT></span>.
  Text name of the lua function to be called for this boundary condition. For
  more on this boundary condition type see \ref LBSBCs

Example:
\code
bsrc={}
for g=1,num_groups do
    bsrc[g] = 0.0
end
bsrc[1] = 1.0

chiLBSSetOptions(phys1,
{
  spatial_discretization = "pwld",
  scattering_order = 2,
  boundary_conditions =
  {
    { name = "xmin", type = "vacuum" },
    { name = "xmax", type = "incident_isotropic", group_strength = bsrc},
    { name = "ymin", type = "reflecting" },
    { name = "zmax", type = "incident_anisotropic_heterogeneous",
      function_name = "TopBoundaryFunction"},
  }
})
\endcode

\ingroup LuaLBS*/
int chiLBSSetOptions(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 2)
    LuaPostArgAmountError(fname, /*expected=*/2, /*given=*/num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckTableValue(fname, L, 2);

  const size_t handle = lua_tointeger(L, 1);

  auto& lbs_solver =
    chi::GetStackItem<lbs::LBSSolver>(chi::object_stack, handle, fname);

  auto specs = chi_lua::TableParserAsParameterBlock::ParseTable(L, 2);

  specs.SetErrorOriginScope(fname);
  for (const auto& spec : specs.Parameters())
  {
    if (spec.Name() == "spatial_discretization")
    {
      auto sdm_name = spec.GetValue<std::string>();
      if (sdm_name == "pwld")
        lbs_solver.Options().sd_type =
          chi_math::SpatialDiscretizationType::PIECEWISE_LINEAR_DISCONTINUOUS;
      else
        ChiInvalidArgument("spatial_discretization only supports \"pwld\"");
    }

    else if (spec.Name() == "scattering_order")
      lbs_solver.Options().scattering_order = spec.GetValue<int>();

    else if (spec.Name() == "sweep_eager_limit")
      lbs_solver.Options().sweep_eager_limit = spec.GetValue<int>();

    else if (spec.Name() == "read_restart_data")
      lbs_solver.Options().read_restart_data = spec.GetValue<bool>();

    else if (spec.Name() == "read_restart_folder_name")
      lbs_solver.Options().read_restart_folder_name =
        spec.GetValue<std::string>();

    else if (spec.Name() == "read_restart_file_base")
      lbs_solver.Options().read_restart_file_base =
        spec.GetValue<std::string>();

    else if (spec.Name() == "write_restart_data")
      lbs_solver.Options().write_restart_data = spec.GetValue<bool>();

    else if (spec.Name() == "write_restart_folder_name")
      lbs_solver.Options().write_restart_folder_name =
        spec.GetValue<std::string>();

    else if (spec.Name() == "write_restart_file_base")
      lbs_solver.Options().write_restart_file_base =
        spec.GetValue<std::string>();

    else if (spec.Name() == "write_restart_interval")
      lbs_solver.Options().write_restart_interval = spec.GetValue<double>();

    else if (spec.Name() == "use_precursors")
      lbs_solver.Options().use_precursors = spec.GetValue<bool>();

    else if (spec.Name() == "use_source_moments")
      lbs_solver.Options().use_src_moments = spec.GetValue<bool>();

    else if (spec.Name() == "save_angular_flux")
      lbs_solver.Options().save_angular_flux = spec.GetValue<bool>();

    else if (spec.Name() == "verbose_inner_iterations")
      lbs_solver.Options().verbose_inner_iterations = spec.GetValue<bool>();

    else if (spec.Name() == "verbose_ags_iterations")
      lbs_solver.Options().verbose_ags_iterations = spec.GetValue<bool>();

    else if (spec.Name() == "verbose_outer_iterations")
      lbs_solver.Options().verbose_outer_iterations = spec.GetValue<bool>();

    else if (spec.Name() == "power_field_function_on")
      lbs_solver.Options().power_field_function_on = spec.GetValue<bool>();

    else if (spec.Name() == "power_default_kappa")
      lbs_solver.Options().power_default_kappa = spec.GetValue<double>();

    else if (spec.Name() == "power_normalization")
      lbs_solver.Options().power_normalization = spec.GetValue<double>();

    else if (spec.Name() == "field_function_prefix_option")
    {
      const std::string value = spec.GetValue<std::string>();

      ChiInvalidArgumentIf(not (value == "solver_name" or value == "prefix"),
                           std::string("Option field_function_prefix_option "
                                       "can only be \"solver_name\" or "
                                       "\"prefix\". \"" +
                                       value + " is not supported."));

      lbs_solver.Options().field_function_prefix_option =
        spec.GetValue<std::string>();
    }

    else if (spec.Name() == "field_function_prefix")
      lbs_solver.Options().field_function_prefix =
        spec.GetValue<std::string>();

    else if (spec.Name() == "boundary_conditions")
    {
      spec.RequireBlockTypeIs(chi_objects::ParameterBlockType::ARRAY);

      const size_t num_spec_params = spec.NumParameters();
      for (size_t p = 0; p < num_spec_params; ++p)
      {
        const auto& bndry_block = spec.GetParam(p);

        bndry_block.RequireParameter("name");
        bndry_block.RequireParameter("type");

        const auto boundary_name =
          bndry_block.GetParamValue<std::string>("name");
        const std::map<std::string, uint64_t> supported_bndry_names = {
          {"xmin", 1},
          {"xmax", 0},
          {"ymin", 3},
          {"ymax", 2},
          {"zmin", 5},
          {"zmax", 4}};

        if (supported_bndry_names.count(boundary_name) == 0)
        {
          std::string message = fname;
          message += ":boundary_conditions:"
                     "name can only be \"xmin\", \"xmax\","
                     " \"ymin\", \"ymax\", \"zmin\", and \"zmax\". \"";
          message += boundary_name + "\" is not supported.";
          throw std::invalid_argument(message);
        }

        const auto bndry_type = bndry_block.GetParamValue<std::string>("type");
        const std::map<std::string, lbs::BoundaryType> type_list = {
          {"vacuum", BoundaryType::VACUUM},
          {"incident_isotropic", BoundaryType::INCIDENT_ISOTROPIC},
          {"reflecting", BoundaryType::REFLECTING},
          {"incident_anisotropic_heterogeneous",
           BoundaryType::INCIDENT_ANISTROPIC_HETEROGENEOUS}};

        if (type_list.count(bndry_type) == 0)
        {
          std::string message = fname;
          message += ":boundary_conditions:"
                     "type can only be \"vacuum\", \"incident_isotropic\","
                     " \"reflecting\", and "
                     "\"incident_anisotropic_heterogeneous\". \"";
          message += bndry_type + "\" is not supported.";
          throw std::invalid_argument(message);
        }

        const auto bid = supported_bndry_names.at(boundary_name);
        const auto type = type_list.at(bndry_type);
        switch (type)
        {
          case BoundaryType::VACUUM:
          case BoundaryType::REFLECTING:
          {
            lbs_solver.BoundaryPreferences()[bid] = {type};
            break;
          }
          case BoundaryType::INCIDENT_ISOTROPIC:
          {
            if (not bndry_block.Has("group_strength"))
            {
              std::string message = fname;
              message += ":boundary_conditions:"
                         "type=\"incident_isotropic\" requires parameter "
                         "\"group_strength\".";

              throw std::invalid_argument(message);
            }
            bndry_block.RequireParameterBlockTypeIs(
              "group_strength", chi_objects::ParameterBlockType::ARRAY);

            const auto group_strength =
              bndry_block.GetParamVectorValue<double>("group_strength");
            lbs_solver.BoundaryPreferences()[bid] = {type, group_strength};
            break;
          }
          case BoundaryType::INCIDENT_ANISTROPIC_HETEROGENEOUS:
          {
            if (not bndry_block.Has("function_name"))
            {
              std::string message = fname;
              message += ":boundary_conditions:"
                         "type=\"incident_anisotropic_heterogeneous\" requires "
                         "parameter \"function_name\".";

              throw std::invalid_argument(message);
            }
            const auto bndry_function_name =
              bndry_block.GetParamValue<std::string>("function_name");

            lbs_solver.BoundaryPreferences()[bid] = {
              type, {}, bndry_function_name};
            break;
          }
        }
      }
    } // "boundary_conditions"

    else
      ChiInvalidArgument(std::string("Unsupported option ") + spec.Name());
  }

  return 0;
}

} // namespace lbs::common_lua_utils