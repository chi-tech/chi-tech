#include "chi_lua.h"

#include "DFEMDiffusion/dfem_diffusion_solver.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/finite_element.h"
#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace chi_unit_sim_tests
{

  int chiSimTest_IP_MMS_L2error(lua_State* L)
  {
    chi::log.Log() << "chiSimTest_IP_MMS_L2error";

    //============================================= Make solver
    dfem_diffusion::Solver solver("SimTest_IP_MMS_L2err");

    solver.Initialize();

    const auto& sdm = *solver.sdm_ptr;
    const auto& OneDofPerNode = sdm.UNITARY_UNKNOWN_MANAGER;
    const size_t num_local_dofs = sdm.GetNumLocalDOFs(OneDofPerNode);
    const size_t num_globl_dofs = sdm.GetNumGlobalDOFs(OneDofPerNode);
    chi::log.Log() << "Num local DOFs: " << num_local_dofs;
    chi::log.Log() << "Num globl DOFs: " << num_globl_dofs;
    chi::log.Log() << "Done constructing solver" << std::endl;

    solver.Execute();

    //============================================= Compute error
    //First get ghosted values
    const auto& field = solver.field;

    double local_error = 0.0;
    for (const auto& cell : sdm.ref_grid->local_cells)
    {
      const int mat_id = cell.material_id;
      const auto& cell_mapping = sdm.GetCellMapping(cell);
      const size_t num_nodes = cell_mapping.NumNodes();
      const auto qp_data = cell_mapping.MakeVolumeQuadraturePointData();

      //======================= Grab nodal phi values
      std::vector<double> nodal_phi(num_nodes,0.0);
      for (size_t j=0; j < num_nodes; ++j)
      {
        const size_t jmap = sdm.MapDOFLocal(cell, j);
        nodal_phi[j] = field[jmap];
      }//for j

      //======================= Quadrature loop
      for (size_t qp : qp_data.QuadraturePointIndices())
      {
        double phi_fem = 0.0;
        for (size_t j=0; j < num_nodes; ++j)
          phi_fem += nodal_phi[j] * qp_data.ShapeValue(j, qp);

        double phi_true = dfem_diffusion::Solver::
          CallLua_iXYZFunction(L,"MMS_phi",mat_id, qp_data.QPointXYZ(qp));

        local_error += std::pow(phi_true - phi_fem,2.0) * qp_data.JxW(qp);
      }
    }//for cell

    double global_error = 0.0;
    MPI_Allreduce(&local_error,     //sendbuf
                  &global_error,    //recvbuf
                  1, MPI_DOUBLE,    //count+datatype
                  MPI_SUM,          //operation
                  MPI_COMM_WORLD);  //communicator

    global_error = std::sqrt(global_error);

    chi::log.Log() << "Error: " << std::scientific << global_error
                   << " Num-cells: " << sdm.ref_grid->GetGlobalNumberOfCells();

    auto stl_vector = new std::vector<double>();
    sdm.LocalizePETScVector(solver.x, *stl_vector, OneDofPerNode);

    //Create ff but use stl_vector instead of &solver.x

    return 0;
  }

}//namespace chi_unit_sim_tests