#include "LBSCurvilinear/lbs_curvilinear_solver.h"

#include <iomanip>

#include "chi_runtime.h"
#include "chi_log.h"
#include "chi_mpi.h"
#include "ChiConsole/chi_console.h"
#include "ChiMath/Quadratures/cylindrical_angular_quadrature.h"
#include "ChiMath/Quadratures/spherical_angular_quadrature.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"
#include "LinearBoltzmannSolver/lbs_structs.h"
#include "LBSCurvilinear/lbs_curvilinear_sweepchunk_pwl.h"
#include "LinearBoltzmannSolver/Groupset/lbs_groupset.h"

typedef chi_mesh::sweep_management::SweepChunk SweepChunk;

void
lbs_curvilinear::Solver::PerformInputChecks()
{
  chi::log.Log() << "LBSCurvilinear::Solver::PerformInputChecks : enter";

  //  --------------------------------------------------------------------------
  //  perform all verifications of Cartesian LBS
  //  --------------------------------------------------------------------------

  lbs::SteadySolver::PerformInputChecks();

  //  --------------------------------------------------------------------------
  //  perform additional verifications for curvilinear LBS
  //  --------------------------------------------------------------------------

  //  coordinate system must be curvilinear
  if (coord_system_type != chi_math::CoordinateSystemType::CYLINDRICAL &&
      coord_system_type != chi_math::CoordinateSystemType::SPHERICAL)
  {
    chi::log.LogAllError()
      << "LBSCurvilinear::Solver::PerformInputChecks : "
      << "invalid coordinate system, static_cast<int>(type) = "
      << static_cast<int>(coord_system_type);
    chi::Exit(EXIT_FAILURE);
  }

  //  re-interpret geometry type to curvilinear
  switch (options.geometry_type)
  {
    case lbs::GeometryType::ONED_SLAB:
    {
      switch (coord_system_type)
      {
        default:
        {
          chi::log.LogAllError()
            << "LBSCurvilinear::Solver::PerformInputChecks : "
            << "invalid geometry, static_cast<int>(type) = "
            << static_cast<int>(options.geometry_type) << " "
            << "for curvilinear coordinate system, static_cast<int>(type) = "
            << static_cast<int>(coord_system_type);
          chi::Exit(EXIT_FAILURE);
        }
      }
      break;
    }
    case lbs::GeometryType::TWOD_CARTESIAN:
    {
      switch (coord_system_type)
      {
        case chi_math::CoordinateSystemType::CYLINDRICAL:
        {
          options.geometry_type = lbs::GeometryType::TWOD_CYLINDRICAL;
          break;
        }
        default:
        {
          chi::log.LogAllError()
            << "LBSCurvilinear::Solver::PerformInputChecks : "
            << "invalid geometry, static_cast<int>(type) = "
            << static_cast<int>(options.geometry_type) << " "
            << "for curvilinear coordinate system, static_cast<int>(type) = "
            << static_cast<int>(coord_system_type);
          chi::Exit(EXIT_FAILURE);
        }
      }
      break;
    }
    default:
    {
      chi::log.LogAllError()
        << "LBSCurvilinear::Solver::PerformInputChecks : "
        << "invalid geometry, static_cast<int>(type) = "
        << static_cast<int>(options.geometry_type) << " "
        << "for curvilinear coordinate system";
      chi::Exit(EXIT_FAILURE);
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
        typedef chi_math::CylindricalAngularQuadrature CylAngQuad;
        const auto curvilinear_angular_quad_ptr =
          std::dynamic_pointer_cast<CylAngQuad>(angular_quad_ptr);
        if (!curvilinear_angular_quad_ptr)
        {
          chi::log.LogAllError()
            << "LBSCurvilinear::Solver::PerformInputChecks : "
            << "invalid angular quadrature, static_cast<int>(type) = "
            << static_cast<int>(angular_quad_ptr->type)
            << ", for groupset = " << gs;
          chi::Exit(EXIT_FAILURE);
        }
        break;
      }
      case chi_math::CoordinateSystemType::SPHERICAL:
      {
        typedef chi_math::SphericalAngularQuadrature SphAngQuad;
        const auto curvilinear_angular_quad_ptr =
          std::dynamic_pointer_cast<SphAngQuad>(angular_quad_ptr);
        if (!curvilinear_angular_quad_ptr)
        {
          chi::log.LogAllError()
            << "LBSCurvilinear::Solver::PerformInputChecks : "
            << "invalid angular quadrature, static_cast<int>(type) = "
            << static_cast<int>(angular_quad_ptr->type)
            << ", for groupset = " << gs;
          chi::Exit(EXIT_FAILURE);
        }
        break;
      }
      default:
      {
        chi::log.LogAllError()
          << "LBSCurvilinear::Solver::PerformInputChecks : "
          << "invalid curvilinear coordinate system, static_cast<int>(type) = "
          << static_cast<int>(coord_system_type);
        chi::Exit(EXIT_FAILURE);
      }
    }

    //  angle aggregation type must be compatible with coordinate system
    const auto angleagg_method = groupsets[gs].angleagg_method;
    switch (coord_system_type)
    {
      case chi_math::CoordinateSystemType::CYLINDRICAL:
      {
        if (angleagg_method != lbs::AngleAggregationType::AZIMUTHAL)
        {
          chi::log.LogAllError()
            << "LBSCurvilinear::Solver::PerformInputChecks : "
            << "invalid angle aggregation type, static_cast<int>(type) = "
            << static_cast<int>(angleagg_method)
            << ", for groupset = " << gs;
          chi::Exit(EXIT_FAILURE);
        }
        break;
      }
      case chi_math::CoordinateSystemType::SPHERICAL:
      {
        if (angleagg_method != lbs::AngleAggregationType::POLAR)
        {
          chi::log.LogAllError()
            << "LBSCurvilinear::Solver::PerformInputChecks : "
            << "invalid angle aggregation type, static_cast<int>(type) = "
            << static_cast<int>(angleagg_method)
            << ", for groupset = " << gs;
          chi::Exit(EXIT_FAILURE);
        }
        break;
      }
      default:
      {
        chi::log.LogAllError()
          << "LBSCurvilinear::Solver::PerformInputChecks : "
          << "invalid curvilinear coordinate system, static_cast<int>(type) = "
          << static_cast<int>(coord_system_type);
        chi::Exit(EXIT_FAILURE);
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
                chi::log.LogAllError()
                  << "LBSCurvilinear::Solver::PerformInputChecks : "
                  << "mesh contains boundary faces with outward-oriented unit "
                  << "normal vector " << (-1*unit_normal_vectors[d]).PrintS()
                  << "with vertices characterised by v(" << d << ") != 0.";
                chi::Exit(EXIT_FAILURE);
              }
            }
            face_orthogonal = true;
            break;
          }
        }
        if (!face_orthogonal)
        {
          chi::log.LogAllError()
            << "LBSCurvilinear::Solver::PerformInputChecks : "
            << "mesh contains boundary faces not orthogonal with respect to "
            << "Cartesian reference frame.";
          chi::Exit(EXIT_FAILURE);
        }
      }
    }
  }

  chi::log.Log() << "LBSCurvilinear::Solver::PerformInputChecks : exit";
}


void
lbs_curvilinear::Solver::InitializeSpatialDiscretization()
{
  chi::log.Log() << "Initializing spatial discretization.\n";

  const auto setup_flags =
    chi_math::finite_element::COMPUTE_CELL_MAPPINGS |
    chi_math::finite_element::COMPUTE_UNIT_INTEGRALS;
  auto qorder = chi_math::QuadratureOrder::INVALID_ORDER;
  auto system = chi_math::CoordinateSystemType::UNDEFINED;

  //  primary discretisation
  switch (options.geometry_type)
  {
    case lbs::GeometryType::ONED_SPHERICAL:
    {
      qorder = chi_math::QuadratureOrder::FOURTH;
      system = chi_math::CoordinateSystemType::SPHERICAL;
      break;
    }
    case lbs::GeometryType::ONED_CYLINDRICAL:
    case lbs::GeometryType::TWOD_CYLINDRICAL:
    {
      qorder = chi_math::QuadratureOrder::THIRD;
      system = chi_math::CoordinateSystemType::CYLINDRICAL;
      break;
    }
    default:
    {
      chi::log.LogAllError()
        << "LBSCurvilinear::Solver::InitializeSpatialDiscretization : "
        << "invalid geometry, static_cast<int>(type) = "
        << static_cast<int>(options.geometry_type);
      chi::Exit(EXIT_FAILURE);
    }
  }

  typedef chi_math::SpatialDiscretization_PWLD SDM_PWLD;
  discretization = SDM_PWLD::New(grid, setup_flags, qorder, system);

  ComputeUnitIntegrals();

  //  secondary discretisation
  //  system - manipulated such that the spatial discretisation returns
  //  a cell view of the same type but with weighting of degree one less
  //  than the primary discretisation
  switch (options.geometry_type)
  {
    case lbs::GeometryType::ONED_SPHERICAL:
    {
      qorder = chi_math::QuadratureOrder::THIRD;
      system = chi_math::CoordinateSystemType::CYLINDRICAL;
      break;
    }
    case lbs::GeometryType::ONED_CYLINDRICAL:
    case lbs::GeometryType::TWOD_CYLINDRICAL:
    {
      qorder = chi_math::QuadratureOrder::SECOND;
      system = chi_math::CoordinateSystemType::CARTESIAN;
      break;
    }
    default:
    {
      chi::log.LogAllError()
        << "LBSCurvilinear::Solver::InitializeSpatialDiscretization : "
        << "invalid geometry, static_cast<int>(type) = "
        << static_cast<int>(options.geometry_type);
      chi::Exit(EXIT_FAILURE);
    }
  }

  discretization_secondary = SDM_PWLD::New(grid, setup_flags, qorder, system);


  MPI_Barrier(MPI_COMM_WORLD);
  chi::log.Log()
    << "Cell matrices computed.                   Process memory = "
    << std::setprecision(3)
    << chi_objects::ChiConsole::GetMemoryUsageInMB() << " MB";
}


std::shared_ptr<SweepChunk>
lbs_curvilinear::Solver::SetSweepChunk(lbs::LBSGroupset& groupset)
{
  auto pwld_sdm_primary =
    std::dynamic_pointer_cast<chi_math::SpatialDiscretization_PWLD>(discretization);
  auto pwld_sdm_secondary =
    std::dynamic_pointer_cast<chi_math::SpatialDiscretization_PWLD>(discretization_secondary);

   auto sweep_chunk =
     std::make_shared<SweepChunkPWL>(grid,
                                     *pwld_sdm_primary,
                                     *pwld_sdm_secondary,
                                     cell_transport_views,
                                     phi_new_local,
                                     psi_new_local[groupset.id],
                                     q_moments_local,
                                     groupset,
                                     matid_to_xs_map,
                                     num_moments,
                                     max_cell_dof_count);

  return sweep_chunk;
}
