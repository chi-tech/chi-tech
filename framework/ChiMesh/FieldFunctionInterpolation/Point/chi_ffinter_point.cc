#include "chi_ffinter_point.h"

#include "ChiPhysics/FieldFunction/fieldfunction.h"
#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_mpi.h"

//###################################################################
/**Initializes the point interpolator.*/
void chi_mesh::FieldFunctionInterpolationPoint::Initialize()
{
  const std::string fname = "FieldFunctionInterpolationPoint::Initialize";
  const auto& grid = *field_functions.front()->SDM().ref_grid;

  std::vector<uint64_t> cells_potentially_owning_point;
  for (const auto& cell : grid.local_cells)
  {
    const auto& vcc = cell.centroid;
    const auto& poi = m_point_of_interest;
    const auto nudged_point = poi + 1.0e-6*(vcc-poi);
    if (grid.CheckPointInsideCell(cell, nudged_point))
      cells_potentially_owning_point.push_back(cell.global_id);
  }

  const int local_count = static_cast<int>(cells_potentially_owning_point.size());
  std::vector<int> locI_count(chi::mpi.process_count,0);
  MPI_Allgather(&local_count,      //sendbuf
                1, MPI_INT,        //sendcount + sendtype
                locI_count.data(), //recvbuf + recvtype
                1, MPI_INT,        //recvcount + recvtype
                MPI_COMM_WORLD);   //communicator

  std::vector<int> recvdispls(chi::mpi.process_count,0);

  int running_count = 0;
  for (int locI=0; locI<chi::mpi.process_count; ++locI)
  {
    recvdispls[locI] = running_count;
    running_count += locI_count[locI];
  }

  const auto& sendbuf = cells_potentially_owning_point;
  std::vector<uint64_t> recvbuf(running_count, 0);
  MPI_Allgatherv(sendbuf.data(),            //sendbuf
                 local_count, MPI_UINT64_T, //sendcount + sendtype
                 recvbuf.data(),            //recvbuf
                 locI_count.data(),         //recvcount
                 recvdispls.data(),         //recvdispl
                 MPI_UINT64_T,              //recvtype
                 MPI_COMM_WORLD);           //communicator

 if (recvbuf.empty())
   throw std::logic_error(fname + ": No cell identified containing the point.");

 uint64_t owning_cell_gid = recvbuf.front();
 for (const uint64_t gid : recvbuf)
   owning_cell_gid = std::min(owning_cell_gid, gid);

 m_locally_owned = false;
 for (const uint64_t gid : cells_potentially_owning_point)
   if (gid == owning_cell_gid)
   {
     m_locally_owned = true;
     m_owning_cell_gid = owning_cell_gid;
     break;
   }
}

//###################################################################
/**Executes the point interpolator.*/
void chi_mesh::FieldFunctionInterpolationPoint::Execute()
{
  if (not m_locally_owned) return;

  const auto& ref_ff = *field_functions.front();
  const auto& sdm    = ref_ff.SDM();
  const auto& grid   = *sdm.ref_grid;

  const auto& uk_man = ref_ff.UnkManager();
  const auto uid = 0;
  const auto cid = m_ref_component;

  using namespace chi_mesh::ff_interpolation;
  const auto field_data = ref_ff.GetGhostedFieldVector();

  const auto& cell = grid.cells[m_owning_cell_gid];
  const auto& cell_mapping = sdm.GetCellMapping(cell);
  const size_t num_nodes = cell_mapping.NumNodes();

  std::vector<double> node_dof_values(num_nodes, 0.0);
  for (size_t i=0; i<num_nodes; ++i)
  {
    const int64_t imap = sdm.MapDOFLocal(cell,i,uk_man,uid,cid);
    node_dof_values[i] = field_data[imap];
  }//for i

  std::vector<double> shape_values(num_nodes, 0.0);
  cell_mapping.ShapeValues(m_point_of_interest, shape_values);

  m_point_value = 0.0;
  for (size_t i=0; i<num_nodes; ++i)
    m_point_value += node_dof_values[i]*shape_values[i];
}

//###################################################################
/**Gets the value of the field function evaluation at the point.*/
double chi_mesh::FieldFunctionInterpolationPoint::GetPointValue() const
{
  double global_point_value;
  MPI_Allreduce(&m_point_value,      //sendbuf
                &global_point_value, //recvbuf
                1, MPI_DOUBLE,       //count + datatype
                MPI_SUM,             //operation
                MPI_COMM_WORLD);     //communicator

  return global_point_value;
}
