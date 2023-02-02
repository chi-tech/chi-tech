#include "lbs_linear_boltzmann_solver.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "ChiConsole/chi_console.h"



#include <iomanip>

void lbs::SteadyStateSolver::InitializeSpatialDiscretization()
{
  using namespace chi_math::finite_element;
  chi::log.Log() << "Initializing spatial discretization.\n";
  discretization_ =
    chi_math::SpatialDiscretization_PWLD::New(grid_ptr_, COMPUTE_CELL_MAPPINGS |
                                                         COMPUTE_UNIT_INTEGRALS);

  ComputeUnitIntegrals();
}

void lbs::SteadyStateSolver::ComputeUnitIntegrals()
{
  chi::log.Log() << "Computing unit integrals.\n";
  const auto& sdm = *discretization_;

  const size_t num_local_cells = grid_ptr_->local_cells.size();
  unit_cell_matrices_.resize(num_local_cells);

  for (const auto& cell : grid_ptr_->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t cell_num_faces = cell.faces.size();
    const size_t cell_num_nodes = cell_mapping.NumNodes();
    const auto vol_qp_data = cell_mapping.MakeVolumeQuadraturePointData();

    MatDbl  IntV_gradshapeI_gradshapeJ(cell_num_nodes, VecDbl(cell_num_nodes));
    MatVec3 IntV_shapeI_gradshapeJ(cell_num_nodes, VecVec3(cell_num_nodes));
    MatDbl  IntV_shapeI_shapeJ(cell_num_nodes, VecDbl(cell_num_nodes));
    VecDbl  IntV_shapeI(cell_num_nodes);

    std::vector<MatDbl>  IntS_shapeI_shapeJ(cell_num_faces);
    std::vector<MatVec3> IntS_shapeI_gradshapeJ(cell_num_faces);
    std::vector<VecDbl>  IntS_shapeI(cell_num_faces);

    //Volume integrals
    for (unsigned int i = 0; i < cell_num_nodes; ++i)
    {
      for (unsigned int j = 0; j < cell_num_nodes; ++j)
      {
        for (const auto& qp : vol_qp_data.QuadraturePointIndices())
        {
          IntV_gradshapeI_gradshapeJ[i][j]
            += vol_qp_data.ShapeGrad(i, qp).Dot(vol_qp_data.ShapeGrad(j, qp)) *
               vol_qp_data.JxW(qp);  //K-matrix

          IntV_shapeI_gradshapeJ[i][j]
            += vol_qp_data.ShapeValue(i, qp) *
               vol_qp_data.ShapeGrad(j, qp) *
               vol_qp_data.JxW(qp);  //G-matrix

          IntV_shapeI_shapeJ[i][j]
            += vol_qp_data.ShapeValue(i, qp) *
               vol_qp_data.ShapeValue(j, qp) *
               vol_qp_data.JxW(qp);  //M-matrix
        }// for qp
      }// for j

      for (const auto& qp : vol_qp_data.QuadraturePointIndices())
      {
        IntV_shapeI[i]
          += vol_qp_data.ShapeValue(i, qp) * vol_qp_data.JxW(qp);
      }// for qp
    }//for i


    //  surface integrals
    for (size_t f = 0; f < cell_num_faces; ++f)
    {
      const auto faces_qp_data = cell_mapping.MakeFaceQuadraturePointData(f);
      IntS_shapeI_shapeJ[f].resize(cell_num_nodes, VecDbl(cell_num_nodes));
      IntS_shapeI[f].resize(cell_num_nodes);
      IntS_shapeI_gradshapeJ[f].resize(cell_num_nodes, VecVec3(cell_num_nodes));

      for (unsigned int i = 0; i < cell_num_nodes; ++i)
      {
        for (unsigned int j = 0; j < cell_num_nodes; ++j)
        {
          for (const auto& qp : faces_qp_data.QuadraturePointIndices())
          {
            IntS_shapeI_shapeJ[f][i][j]
              += faces_qp_data.ShapeValue(i, qp) *
                 faces_qp_data.ShapeValue(j, qp) *
                 faces_qp_data.JxW(qp);
            IntS_shapeI_gradshapeJ[f][i][j]
              += faces_qp_data.ShapeValue(i, qp) *
                 faces_qp_data.ShapeGrad(j, qp) *
                 faces_qp_data.JxW(qp);
          }// for qp
        }//for j

        for (const auto& qp : faces_qp_data.QuadraturePointIndices())
        {
          IntS_shapeI[f][i]
            += faces_qp_data.ShapeValue(i, qp) * faces_qp_data.JxW(qp);
        }// for qp
      }//for i
    }//for f

    unit_cell_matrices_[cell.local_id] =
      UnitCellMatrices{IntV_gradshapeI_gradshapeJ, //K-matrix
                       IntV_shapeI_gradshapeJ,     //G-matrix
                       IntV_shapeI_shapeJ,         //M-matrix
                       IntV_shapeI,                //Vi-vectors

                       IntS_shapeI_shapeJ,         //face M-matrices
                       IntS_shapeI_gradshapeJ,     //face G-matrices
                       IntS_shapeI};               //face Si-vectors
  }//for cell


  MPI_Barrier(MPI_COMM_WORLD);
  chi::log.Log()
    << "Cell matrices computed.                   Process memory = "
    << std::setprecision(3)
    << chi_objects::ChiConsole::GetMemoryUsageInMB() << " MB";
}