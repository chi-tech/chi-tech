#include "chi_volumemesher.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiMesh/LogicalVolume/chi_mesh_logicalvolume.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "ChiConsole/chi_console.h"

#include "ChiTimer/chi_timer.h"
#include "ChiLua/chi_lua.h"


//###################################################################
/**Sets material id's using a logical volume.*/
void chi_mesh::VolumeMesher::
  SetMatIDFromLogical(const chi_mesh::LogicalVolume& log_vol,
                      bool sense, int mat_id)
{
  chi::log.Log0Verbose1()
    << chi::program_timer.GetTimeString()
    << " Setting material id from logical volume.";
  //============================================= Get current mesh handler
  auto& handler = chi_mesh::GetCurrentHandler();

  //============================================= Get back mesh
  chi_mesh::MeshContinuumPtr vol_cont = handler.GetGrid();

  int num_cells_modified = 0;
  for (auto& cell : vol_cont->local_cells)
  {
    if (log_vol.Inside(cell.centroid_) && sense){
      cell.material_id_ = mat_id;
      ++num_cells_modified;
    }
  }

  const auto& ghost_ids = vol_cont->cells.GetGhostGlobalIDs();
  for (uint64_t ghost_id : ghost_ids)
  {
    auto& cell = vol_cont->cells[ghost_id];
    if (log_vol.Inside(cell.centroid_) && sense)
      cell.material_id_ = mat_id;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  chi::log.Log0Verbose1()
    << chi::program_timer.GetTimeString()
    << " Done setting material id from logical volume. "
    << "Number of cells modified = " << num_cells_modified << ".";
}

//###################################################################
/**Sets material id's using a logical volume.*/
void chi_mesh::VolumeMesher::
  SetBndryIDFromLogical(const chi_mesh::LogicalVolume& log_vol,
                        bool sense, int bndry_id)
{
  chi::log.Log()
    << chi::program_timer.GetTimeString()
    << " Setting boundary id from logical volume.";
  //============================================= Get current mesh handler
  auto& handler = chi_mesh::GetCurrentHandler();

  //============================================= Get back mesh
  chi_mesh::MeshContinuumPtr vol_cont = handler.GetGrid();

  int num_faces_modified = 0;
  for (auto& cell : vol_cont->local_cells)
  {
    for (auto& face : cell.faces_)
    {
      if (face.has_neighbor_) continue;
      if (log_vol.Inside(face.centroid_) && sense){
        face.neighbor_id_ = abs(bndry_id);
        ++num_faces_modified;
      }
    }
  }

  chi::log.Log()
    << chi::program_timer.GetTimeString()
    << " Done setting boundary id from logical volume. "
    << "Number of faces modified = " << num_faces_modified << ".";
}

//###################################################################
/**Sets material id's for all cells to the specified material id.*/
void chi_mesh::VolumeMesher::SetMatIDToAll(int mat_id)
{
  chi::log.Log()
    << chi::program_timer.GetTimeString()
    << " Setting material id " << mat_id << " to all cells.";

  //============================================= Get current mesh handler
  auto& handler = chi_mesh::GetCurrentHandler();

  //============================================= Get back mesh
  auto vol_cont = handler.GetGrid();

  for (auto& cell : vol_cont->local_cells)
    cell.material_id_ = mat_id;

  const auto& ghost_ids = vol_cont->cells.GetGhostGlobalIDs();
  for (uint64_t ghost_id : ghost_ids)
    vol_cont->cells[ghost_id].material_id_ = mat_id;

  MPI_Barrier(MPI_COMM_WORLD);
  chi::log.Log()
    << chi::program_timer.GetTimeString()
    << " Done setting material id " << mat_id << " to all cells";
}


//###################################################################
/**Sets material id's using a lua function. The lua function is called
with for each cell with 4 arguments, the cell's centroid x,y,z values
and the cell's current material id.

The lua function's prototype should be:
\code
function LuaFuncName(x,y,z,id)
 --stuff
end
\endcode
*/
void chi_mesh::VolumeMesher::
  SetMatIDFromLuaFunction(const std::string& lua_fname)
{
  const std::string fname = "chi_mesh::VolumeMesher::SetMatIDFromLuaFunction";

  chi::log.Log0Verbose1()
    << chi::program_timer.GetTimeString()
    << " Setting material id from lua function.";

  //============================================= Define console call
  auto L = chi::console.GetConsoleState();
  auto CallLuaXYZFunction = [&L,&lua_fname,&fname](const chi_mesh::Cell& cell)
  {
    //============= Load lua function
    lua_getglobal(L, lua_fname.c_str());

    //============= Error check lua function
    if (not lua_isfunction(L, -1))
      throw std::logic_error(fname + " attempted to access lua-function, " +
                               lua_fname + ", but it seems the function"
                                             " could not be retrieved.");

    const auto& xyz = cell.centroid_;

    //============= Push arguments
    lua_pushnumber(L, xyz.x);
    lua_pushnumber(L, xyz.y);
    lua_pushnumber(L, xyz.z);
    lua_pushinteger(L, cell.material_id_);

    //============= Call lua function
    //4 arguments, 1 result (double), 0=original error object
    int lua_return;
    if (lua_pcall(L,4,1,0) == 0)
    {
      LuaCheckNumberValue(fname, L, -1);
      lua_return = lua_tointeger(L,-1);
    }
    else
      throw std::logic_error(fname + " attempted to call lua-function, " +
                             lua_fname + ", but the call failed.");

    lua_pop(L,1); //pop the int, or error code

    return lua_return;
  };

  //============================================= Get current mesh handler
  auto& handler = chi_mesh::GetCurrentHandler();

  //============================================= Get back mesh
  chi_mesh::MeshContinuum& grid = *handler.GetGrid();

  int local_num_cells_modified = 0;
  for (auto& cell : grid.local_cells)
  {
    int new_matid = CallLuaXYZFunction(cell);

    if (cell.material_id_ != new_matid)
    {
      cell.material_id_ = new_matid;
      ++local_num_cells_modified;
    }
  }//for local cell

  const auto& ghost_ids = grid.cells.GetGhostGlobalIDs();
  for (uint64_t ghost_id : ghost_ids)
  {
    auto& cell = grid.cells[ghost_id];
    int new_matid = CallLuaXYZFunction(cell);

    if (cell.material_id_ != new_matid)
    {
      cell.material_id_ = new_matid;
      ++local_num_cells_modified;
    }
  }//for ghost cell id

  int globl_num_cells_modified;
  MPI_Allreduce(&local_num_cells_modified, //sendbuf
                &globl_num_cells_modified, //recvbuf
                1, MPI_INT,                //count+datatype
                MPI_SUM,                   //operation
                MPI_COMM_WORLD);           //comm

  chi::log.Log0Verbose1()
    << chi::program_timer.GetTimeString()
    << " Done setting material id from lua function. "
    << "Number of cells modified = " << globl_num_cells_modified << ".";
}


//###################################################################
/**Sets boundary id's using a lua function. The lua function is called
for each boundary face with 7 arguments, the face's centroid x,y,z values,
the face's normal x,y,z values and the face's current boundary id.

The lua function's prototype should be:
\code
function LuaFuncName(x,y,z,nx,ny,nz,id)
 --stuff
end
\endcode
*/
void chi_mesh::VolumeMesher::
  SetBndryIDFromLuaFunction(const std::string& lua_fname)
{
  const std::string fname = "chi_mesh::VolumeMesher::SetBndryIDFromLuaFunction";

  chi::log.Log0Verbose1()
    << chi::program_timer.GetTimeString()
    << " Setting boundary id from lua function.";

  //============================================= Define console call
  auto L = chi::console.GetConsoleState();
  auto CallLuaXYZFunction = [&L,&lua_fname,&fname]
    (const chi_mesh::CellFace& face)
  {
    //============= Load lua function
    lua_getglobal(L, lua_fname.c_str());

    //============= Error check lua function
    if (not lua_isfunction(L, -1))
      throw std::logic_error(fname + " attempted to access lua-function, " +
                             lua_fname + ", but it seems the function"
                                         " could not be retrieved.");

    const auto& xyz = face.centroid_;
    const auto& n   = face.normal_;

    //============= Push arguments
    lua_pushnumber(L, xyz.x);
    lua_pushnumber(L, xyz.y);
    lua_pushnumber(L, xyz.z);
    lua_pushnumber(L, n.x);
    lua_pushnumber(L, n.y);
    lua_pushnumber(L, n.z);
    lua_pushinteger(L, static_cast<lua_Integer>(face.neighbor_id_));

    //============= Call lua function
    //7 arguments, 1 result (double), 0=original error object
    int lua_return;
    if (lua_pcall(L,7,1,0) == 0)
    {
      LuaCheckNumberValue(fname, L, -1);
      lua_return = lua_tointeger(L,-1);
    }
    else
      throw std::logic_error(fname + " attempted to call lua-function, " +
                             lua_fname + ", but the call failed.");

    lua_pop(L,1); //pop the int, or error code

    return lua_return;
  };

  //============================================= Get current mesh handler
  auto& handler = chi_mesh::GetCurrentHandler();

  //============================================= Get back mesh
  chi_mesh::MeshContinuum& grid = *handler.GetGrid();

  int local_num_faces_modified = 0;
  for (auto& cell : grid.local_cells)
    for (auto& face : cell.faces_)
      if (not face.has_neighbor_)
      {
        int new_bndryid = CallLuaXYZFunction(face);

        if (face.neighbor_id_ != new_bndryid)
        {
          face.neighbor_id_ = new_bndryid;
          ++local_num_faces_modified;
        }
      }//for bndry face

  const auto& ghost_ids = grid.cells.GetGhostGlobalIDs();
  for (uint64_t ghost_id : ghost_ids)
  {
    auto& cell = grid.cells[ghost_id];
    for (auto& face : cell.faces_)
      if (not face.has_neighbor_)
      {
        int new_bndryid = CallLuaXYZFunction(face);

        if (face.neighbor_id_ != new_bndryid)
        {
          face.neighbor_id_ = new_bndryid;
          ++local_num_faces_modified;
        }
      }//for bndry face
  }//for ghost cell id

  int globl_num_faces_modified;
  MPI_Allreduce(&local_num_faces_modified, //sendbuf
                &globl_num_faces_modified, //recvbuf
                1, MPI_INT,                //count+datatype
                MPI_SUM,                   //operation
                MPI_COMM_WORLD);           //comm

  chi::log.Log0Verbose1()
    << chi::program_timer.GetTimeString()
    << " Done setting boundary id from lua function. "
    << "Number of cells modified = " << globl_num_faces_modified << ".";
}