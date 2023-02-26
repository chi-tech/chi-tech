#include "lbsadj_solver.h"

#include "ChiMesh/LogicalVolume/chi_mesh_logicalvolume.h"

#include "chi_runtime.h"
#include "chi_log.h"

void lbs::DiscOrdSteadyStateAdjointSolver::Initialize()
{
  lbs::DiscOrdSteadyStateSolver::Initialize();

  //============================================= Transpose Scattering operators
  for (const auto& matid_xs_pair : matid_to_xs_map_)
  {
    const auto  matid = matid_xs_pair.first;
    const auto& S = matid_xs_pair.second->transfer_matrices_;

    std::vector<chi_math::SparseMatrix> S_transpose;
    for (const auto& S_ell : S)
    {
      chi_math::SparseMatrix S_ell_transpose(S_ell.NumRows(), S_ell.NumCols());

      for (size_t g=0; g<S_ell.NumRows(); ++g)
      {
        size_t num_col_indices = S_ell.rowI_indices[g].size();
        for (size_t j=0; j<num_col_indices; ++j)
        {
          size_t gprime = S_ell.rowI_indices[g][j];
          double value  = S_ell.rowI_values[g][j];

          S_ell_transpose.Insert(gprime, g, value);
        }
      }

      S_transpose.push_back(std::move(S_ell_transpose));
    }//for each S

    matid_to_S_transpose_[matid] = std::move(S_transpose);
  }//for each mat

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