#include "LBSCurvilinear/lbs_curvilinear_solver.h"

#include <iomanip>

#include "chi_log.h"
#include "chi_mpi.h"
#include "ChiConsole/chi_console.h"
#include "ChiMath/Quadratures/cylindrical_angular_quadrature.h"
#include "ChiMath/Quadratures/spherical_angular_quadrature.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"
#include "LinearBoltzmannSolver/lbs_structs.h"
#include "LBSCurvilinear/lbs_curvilinear_sweepchunk_pwl.h"

typedef chi_mesh::sweep_management::SweepChunk SweepChunk;

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;
extern ChiConsole& chi_console;


void
LBSCurvilinear::Solver::PerformInputChecks()
{
  chi_log.Log(LOG_0) << "LBSCurvilinear::Solver::PerformInputChecks : enter";

  //  --------------------------------------------------------------------------
  //  perform all verifications of Cartesian LBS
  //  --------------------------------------------------------------------------

  LinearBoltzmann::Solver::PerformInputChecks();

  //  --------------------------------------------------------------------------
  //  perform additional verifications for curvilinear LBS
  //  --------------------------------------------------------------------------

  //  coordinate system must be curvilinear
  if (coord_system_type != chi_math::CoordinateSystemType::CYLINDRICAL &&
      coord_system_type != chi_math::CoordinateSystemType::SPHERICAL)
  {
    chi_log.Log(LOG_ALLERROR)
      << "LBSCurvilinear::Solver::PerformInputChecks : "
      << "invalid coordinate system, static_cast<int>(type) = "
      << static_cast<int>(coord_system_type);
    std::exit(EXIT_FAILURE);
  }

  //  re-interpret geometry type to curvilinear
  switch (options.geometry_type)
  {
    case LinearBoltzmann::GeometryType::ONED_SLAB:
    {
      switch (coord_system_type)
      {
      //case chi_math::CoordinateSystemType::CYLINDRICAL:
      //{
      //  options.geometry_type = LinearBoltzmann::GeometryType::ONED_CYLINDRICAL;
      //  break;
      //}
      //case chi_math::CoordinateSystemType::SPHERICAL:
      //{
      //  options.geometry_type = LinearBoltzmann::GeometryType::ONED_SPHERICAL;
      //  break;
      //}
        default:
        {
          chi_log.Log(LOG_ALLERROR)
            << "LBSCurvilinear::Solver::PerformInputChecks : "
            << "invalid geometry, static_cast<int>(type) = "
            << static_cast<int>(options.geometry_type) << " "
            << "for curvilinear coordinate system, static_cast<int>(type) = "
            << static_cast<int>(coord_system_type);
          std::exit(EXIT_FAILURE);
        }
      }
      break;
    }
    case LinearBoltzmann::GeometryType::TWOD_CARTESIAN:
    {
      switch (coord_system_type)
      {
        case chi_math::CoordinateSystemType::CYLINDRICAL:
        {
          options.geometry_type = LinearBoltzmann::GeometryType::TWOD_CYLINDRICAL;
          break;
        }
        default:
        {
          chi_log.Log(LOG_ALLERROR)
            << "LBSCurvilinear::Solver::PerformInputChecks : "
            << "invalid geometry, static_cast<int>(type) = "
            << static_cast<int>(options.geometry_type) << " "
            << "for curvilinear coordinate system, static_cast<int>(type) = "
            << static_cast<int>(coord_system_type);
          std::exit(EXIT_FAILURE);
        }
      }
      break;
    }
    default:
    {
      chi_log.Log(LOG_ALLERROR)
        << "LBSCurvilinear::Solver::PerformInputChecks : "
        << "invalid geometry, static_cast<int>(type) = "
        << static_cast<int>(options.geometry_type) << " "
        << "for curvilinear coordinate system";
      std::exit(EXIT_FAILURE);
    }
  }

  for (size_t gs = 0; gs < groupsets.size(); ++gs)
  {
    //  angular quadrature type must be compatible with coordinate system
    const auto angular_quad_ptr = groupsets[gs].quadrature;
    switch (coord_system_type)
    {
      case chi_math::CoordinateSystemType::CYLINDRICAL:
      {
        const auto curvilinear_angular_quad_ptr =
          std::dynamic_pointer_cast<chi_math::CylindricalAngularQuadrature>(angular_quad_ptr);
        if (!curvilinear_angular_quad_ptr)
        {
          chi_log.Log(LOG_ALLERROR)
            << "LBSCurvilinear::Solver::PerformInputChecks : "
            << "invalid angular quadrature, static_cast<int>(type) = "
            << static_cast<int>(angular_quad_ptr->type)
            << ", for groupset = " << gs;
          std::exit(EXIT_FAILURE);
        }
        break;
      }
      case chi_math::CoordinateSystemType::SPHERICAL:
      {
        const auto curvilinear_angular_quad_ptr =
          std::dynamic_pointer_cast<chi_math::SphericalAngularQuadrature>(angular_quad_ptr);
        if (!curvilinear_angular_quad_ptr)
        {
          chi_log.Log(LOG_ALLERROR)
            << "LBSCurvilinear::Solver::PerformInputChecks : "
            << "invalid angular quadrature, static_cast<int>(type) = "
            << static_cast<int>(angular_quad_ptr->type)
            << ", for groupset = " << gs;
          std::exit(EXIT_FAILURE);
        }
        break;
      }
      default:
      {
        chi_log.Log(LOG_ALLERROR)
          << "LBSCurvilinear::Solver::PerformInputChecks : "
          << "invalid curvilinear coordinate system, static_cast<int>(type) = "
          << static_cast<int>(coord_system_type);
        std::exit(EXIT_FAILURE);
      }
    }

    //  angle aggregation type must be compatible with coordinate system
    const auto angleagg_method = groupsets[gs].angleagg_method;
    switch (coord_system_type)
    {
      case chi_math::CoordinateSystemType::CYLINDRICAL:
      {
        if (angleagg_method != LinearBoltzmann::AngleAggregationType::AZIMUTHAL)
        {
          chi_log.Log(LOG_ALLERROR)
            << "LBSCurvilinear::Solver::PerformInputChecks : "
            << "invalid angle aggregation type, static_cast<int>(type) = "
            << static_cast<int>(angleagg_method)
            << ", for groupset = " << gs;
          std::exit(EXIT_FAILURE);
        }
        break;
      }
      case chi_math::CoordinateSystemType::SPHERICAL:
      {
        if (angleagg_method != LinearBoltzmann::AngleAggregationType::POLAR)
        {
          chi_log.Log(LOG_ALLERROR)
            << "LBSCurvilinear::Solver::PerformInputChecks : "
            << "invalid angle aggregation type, static_cast<int>(type) = "
            << static_cast<int>(angleagg_method)
            << ", for groupset = " << gs;
          std::exit(EXIT_FAILURE);
        }
        break;
      }
      default:
      {
        chi_log.Log(LOG_ALLERROR)
          << "LBSCurvilinear::Solver::PerformInputChecks : "
          << "invalid curvilinear coordinate system, static_cast<int>(type) = "
          << static_cast<int>(coord_system_type);
        std::exit(EXIT_FAILURE);
      }
    }
  }

  //  boundary of mesh must be rectangular with origin at (0, 0, 0)
  const std::vector<chi_mesh::Vector3> unit_normal_vectors =
    { chi_mesh::Vector3(1.0, 0.0, 0.0),
      chi_mesh::Vector3(0.0, 1.0, 0.0),
      chi_mesh::Vector3(0.0, 0.0, 1.0), };
  for (const auto& cell : grid->local_cells)
  {
    for (const auto& face : cell.faces)
    {
      if (!face.has_neighbor)
      {
        bool face_orthogonal = false;
        for (size_t d = 0; d < unit_normal_vectors.size(); ++d)
        {
          const auto n_dot_e = face.normal.Dot(unit_normal_vectors[d]);
          if      (n_dot_e >  0.999999)
          {
            face_orthogonal = true;
            break;
          }
          else if (n_dot_e < -0.999999)
          {
            for (const auto& v_id : face.vertex_ids)
            {
              const auto& vertex = grid->vertices[v_id];
              if (std::abs(vertex[d]) > 1.0e-12)
              {
                chi_log.Log(LOG_ALLERROR)
                  << "LBSCurvilinear::Solver::PerformInputChecks : "
                  << "mesh contains boundary faces with outward-oriented unit "
                  << "normal vector " << (-1*unit_normal_vectors[d]).PrintS()
                  << "with vertices characterised by v(" << d << ") != 0.";
                std::exit(EXIT_FAILURE);
              }
            }
            face_orthogonal = true;
            break;
          }
        }
        if (!face_orthogonal)
        {
          chi_log.Log(LOG_ALLERROR)
            << "LBSCurvilinear::Solver::PerformInputChecks : "
            << "mesh contains boundary faces not orthogonal with respect to "
            << "Cartesian reference frame.";
          std::exit(EXIT_FAILURE);
        }
      }
    }
  }

  chi_log.Log(LOG_0) << "LBSCurvilinear::Solver::PerformInputChecks : exit";
}


void
LBSCurvilinear::Solver::InitializeSpatialDiscretization()
{
  chi_log.Log(LOG_0) << "Initializing spatial discretization.\n";

  const auto setup_flags =
    chi_math::finite_element::COMPUTE_CELL_MAPPINGS |
    chi_math::finite_element::COMPUTE_UNIT_INTEGRALS;
  auto qorder = chi_math::QuadratureOrder::INVALID_ORDER;
  auto system = chi_math::CoordinateSystemType::UNDEFINED;

  //  primary discretisation
  switch (options.geometry_type)
  {
    case LinearBoltzmann::GeometryType::ONED_SPHERICAL:
    {
      qorder = chi_math::QuadratureOrder::FOURTH;
      system = chi_math::CoordinateSystemType::SPHERICAL;
      break;
    }
    case LinearBoltzmann::GeometryType::ONED_CYLINDRICAL:
    {
      qorder = chi_math::QuadratureOrder::THIRD;
      system = chi_math::CoordinateSystemType::CYLINDRICAL;
      break;
    }
    case LinearBoltzmann::GeometryType::TWOD_CYLINDRICAL:
    {
      qorder = chi_math::QuadratureOrder::THIRD;
      system = chi_math::CoordinateSystemType::CYLINDRICAL;
      break;
    }
    default:
    {
      chi_log.Log(LOG_ALLERROR)
        << "LBSCurvilinear::Solver::InitializeSpatialDiscretization : "
        << "invalid geometry, static_cast<int>(type) = "
        << static_cast<int>(options.geometry_type);
      std::exit(EXIT_FAILURE);
    }
  }

  discretization =
    SpatialDiscretization_PWLD::New(grid, setup_flags, qorder, system);

  //  secondary discretisation
  //  system - manipulated such that the spatial discretisation returns
  //  a cell view of the same type but with weighting of degree one less
  //  than the primary discretisation
  switch (options.geometry_type)
  {
    case LinearBoltzmann::GeometryType::ONED_SPHERICAL:
    {
      qorder = chi_math::QuadratureOrder::THIRD;
      system = chi_math::CoordinateSystemType::CYLINDRICAL;
      break;
    }
    case LinearBoltzmann::GeometryType::ONED_CYLINDRICAL:
    {
      qorder = chi_math::QuadratureOrder::SECOND;
      system = chi_math::CoordinateSystemType::CARTESIAN;
      break;
    }
    case LinearBoltzmann::GeometryType::TWOD_CYLINDRICAL:
    {
      qorder = chi_math::QuadratureOrder::SECOND;
      system = chi_math::CoordinateSystemType::CARTESIAN;
      break;
    }
    default:
    {
      chi_log.Log(LOG_ALLERROR)
        << "LBSCurvilinear::Solver::InitializeSpatialDiscretization : "
        << "invalid geometry, static_cast<int>(type) = "
        << static_cast<int>(options.geometry_type);
      std::exit(EXIT_FAILURE);
    }
  }

  discretization_secondary =
    SpatialDiscretization_PWLD::New(grid, setup_flags, qorder, system);


  MPI_Barrier(MPI_COMM_WORLD);
  chi_log.Log(LOG_0)
    << "Cell matrices computed.                   Process memory = "
    << std::setprecision(3)
    << chi_console.GetMemoryUsageInMB() << " MB";
}


std::shared_ptr<SweepChunk>
LBSCurvilinear::Solver::SetSweepChunk(LBSGroupset& groupset)
{
  auto pwld_sdm_primary =
    std::dynamic_pointer_cast<SpatialDiscretization_PWLD>(discretization);
  auto pwld_sdm_secondary =
    std::dynamic_pointer_cast<SpatialDiscretization_PWLD>(discretization_secondary);

   auto sweep_chunk =
     std::make_shared<SweepChunkPWL>(grid,
                                     *pwld_sdm_primary,
                                     *pwld_sdm_secondary,
                                     cell_transport_views,
                                     phi_new_local,
                                     psi_new_local[groupset.id],
                                     q_moments_local,
                                     groupset,
                                     material_xs,
                                     num_moments,
                                     max_cell_dof_count);

  return sweep_chunk;
}
