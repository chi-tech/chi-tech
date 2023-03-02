#include "lbsadj_solver.h"

#include "ChiMesh/LogicalVolume/chi_mesh_logicalvolume.h"

#include "chi_runtime.h"
#include "chi_log.h"


void lbs::DiscOrdSteadyStateAdjointSolver::Initialize()
{
  lbs::DiscOrdSteadyStateSolver::Initialize();

  //============================================= Create adjoint cross sections
  using AdjXSPtr = std::shared_ptr<chi_physics::AdjointMultiGroupXS>;

  std::map<int, XSPtr> matid_to_adj_xs_map;
  for (const auto& matid_xs_pair : matid_to_xs_map_)
  {
    const auto matid = matid_xs_pair.first;
    const auto fwd_xs = std::dynamic_pointer_cast<
        chi_physics::AdjointMultiGroupXS>(matid_xs_pair.second);

    AdjXSPtr adj_xs_ptr(new chi_physics::AdjointMultiGroupXS(*fwd_xs));
    matid_to_adj_xs_map[matid] = adj_xs_ptr;
  }//for each mat
  matid_to_xs_map_ = std::move(matid_to_adj_xs_map);

  //============================================= Initialize QOIs
  for (auto& qoi_pair : response_functions_)
  {
    const auto& qoi_designation = qoi_pair.first;
    auto& qoi_cell_subscription = qoi_pair.second;

    for (const auto& cell : grid_ptr_->local_cells)
      if (qoi_designation.logical_volume->Inside(cell.centroid_))
        qoi_cell_subscription.push_back(cell.local_id_);

    size_t num_local_subs = qoi_cell_subscription.size();
    size_t num_globl_subs = 0;

    MPI_Allreduce(&num_local_subs,           //sendbuf
                  &num_globl_subs,           //recvbuf
                  1, MPI_UNSIGNED_LONG_LONG, //count + datatype
                  MPI_SUM,                   //operation
                  MPI_COMM_WORLD );          //communicator

    chi::log.Log() << "LBAdjointSolver: Number of cells subscribed to "
                  << qoi_designation.name << " = "
                  << num_globl_subs;
  }

  //================================================== Initialize source func
  using namespace std::placeholders;
  active_set_source_function_ =
    std::bind(&DiscOrdSteadyStateAdjointSolver::SetAdjointSource, this, _1, _2, _3, _4);

}