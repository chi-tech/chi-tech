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
 * \param specs Table Various options in a table. Detailed below.

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

\ingroup LBSLuaFunctions*/
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

  auto options_params = LBSSolver::OptionsBlock();
  options_params.AssignParameters(specs);

  lbs_solver.SetOptions(options_params);

  return 0;
}

} // namespace lbs::common_lua_utils