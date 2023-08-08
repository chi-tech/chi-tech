#include "lbsadj_solver.h"

#include "mesh/LogicalVolume/LogicalVolume.h"
#include "physics/PhysicsMaterial/MultiGroupXS/adjoint_mgxs.h"
#include "A_LBSSolver/SourceFunctions/adjoint_src_function.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "chi_log.h"


void lbs::DiscreteOrdinatesAdjointSolver::Initialize()
{
  lbs::DiscreteOrdinatesSolver::Initialize();

  //============================================= Create adjoint cross sections
  using AdjXS = chi_physics::AdjointMGXS;

  // define the actual cross sections
  std::map<int, XSPtr> matid_to_adj_xs_map;
  for (const auto& matid_xs_pair : matid_to_xs_map_)
  {
    const auto matid = matid_xs_pair.first;
    const auto fwd_xs = std::dynamic_pointer_cast<
        chi_physics::MultiGroupXS>(matid_xs_pair.second);
    matid_to_adj_xs_map[matid] = std::make_shared<AdjXS>(*fwd_xs);
  }//for each mat
  matid_to_xs_map_ = std::move(matid_to_adj_xs_map);

  // reassign transport view to adjoint cross sections
  if (grid_ptr_->local_cells.size() == cell_transport_views_.size())
    for (const auto& cell : grid_ptr_->local_cells)
    {
      const auto& xs_ptr = matid_to_xs_map_[cell.material_id_];
      auto& transport_view = cell_transport_views_[cell.local_id_];

      transport_view.ReassingXS(*xs_ptr);
    }

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
                  Chi::mpi.comm );          //communicator

    Chi::log.Log() << "LBAdjointSolver: Number of cells subscribed to "
                  << qoi_designation.name << " = "
                  << num_globl_subs;
  }

  //================================================== Initialize source func
  auto src_function = std::make_shared<AdjointSourceFunction>(*this);

  using namespace std::placeholders;
  active_set_source_function_ =
    std::bind(&SourceFunction::operator(), src_function, _1, _2, _3, _4);
}