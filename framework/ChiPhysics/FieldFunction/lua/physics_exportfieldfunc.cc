#include "ChiLua/chi_lua.h"

#include "ChiPhysics/FieldFunction/fieldfunction_gridbased.h"

#include "chi_runtime.h"
#include "chi_log.h"

// #############################################################################
/** Exports a field function to VTK format.
 *
\param FFHandle int Global handle to the field function.
\param BaseName char Base name for the exported file.

\ingroup LuaFieldFunc
\author Jan*/
int chiExportFieldFunctionToVTK(lua_State* L)
{
  const std::string fname = "chiExportFieldFunctionToVTK";
  const int num_args = lua_gettop(L);
  if (num_args != 2) LuaPostArgAmountError(fname, 2, num_args);

  int ff_handle = lua_tonumber(L, 1);
  const char* base_name = lua_tostring(L, 2);

  typedef chi_physics::FieldFunctionGridBased FFGridBased;

  auto ff_base =
    chi::GetStackItemPtr(chi::field_function_stack, ff_handle, fname);
  auto ff = std::dynamic_pointer_cast<FFGridBased>(ff_base);

  ChiLogicalErrorIf(not ff, "Only grid-based field functions can be exported");

  //  ff->ExportToVTK(base_name);
  chi_physics::FieldFunctionGridBased::ExportMultipleToVTK(base_name, {ff});

  return 0;
}

// #############################################################################
/** Exports all the field functions in a list to VTK format.
 *  *
\param listFFHandles table Global handles to the field functions
\param BaseName char Base name for the exported file.

\ingroup LuaFieldFunc
\author Jan*/
int chiExportMultiFieldFunctionToVTK(lua_State* L)
{
  const std::string fname = "chiExportMultiFieldFunctionToVTK";
  const int num_args = lua_gettop(L);
  if (num_args != 2) LuaPostArgAmountError(fname, 2, num_args);

  const char* base_name = lua_tostring(L, 2);

  LuaCheckTableValue(fname, L, 1);

  const size_t table_size = lua_rawlen(L, 1);
  std::vector<std::shared_ptr<const chi_physics::FieldFunctionGridBased>> ffs;
  ffs.reserve(table_size);
  for (int i = 0; i < table_size; ++i)
  {
    lua_pushnumber(L, i + 1);
    lua_gettable(L, 1);

    int ff_handle = lua_tonumber(L, -1);
    lua_pop(L, 1);

    auto ff_base =
      chi::GetStackItemPtr(chi::field_function_stack, ff_handle, fname);

    typedef chi_physics::FieldFunctionGridBased FFGridBased;
    auto ff = std::dynamic_pointer_cast<FFGridBased>(ff_base);

    ChiLogicalErrorIf(not ff,
                      "Only grid-based field functions can be exported");

    ffs.push_back(ff);
  }

  chi_physics::FieldFunctionGridBased::ExportMultipleToVTK(base_name, ffs);

  return 0;
}
