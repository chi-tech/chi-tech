#include "C_DO_RZ_SteadyState/lbs_curvilinear_solver.h"

#include "ChiMath/Quadratures/cylindrical_angular_quadrature.h"
#include "ChiMath/Quadratures/spherical_angular_quadrature.h"

#include "A_LBSSolver/lbs_structs.h"
#include "C_DO_RZ_SteadyState/lbs_curvilinear_sweepchunk_pwl.h"
#include "LinearBoltzmannSolvers/A_LBSSolver/Groupset/lbs_groupset.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "ChiConsole/chi_console.h"

#include <iomanip>

typedef chi_mesh::sweep_management::SweepChunk SweepChunk;

void
lbs_curvilinear::DiscOrdSteadyStateSolver::PerformInputChecks()
{
  chi::log.Log() << "C_DO_RZ_SteadyState::SteadyStateSolver::PerformInputChecks : enter";

  //  --------------------------------------------------------------------------
  //  perform all verifications of Cartesian LBS
  //  --------------------------------------------------------------------------

  lbs::DiscOrdSteadyStateSolver::PerformInputChecks();

  //  --------------------------------------------------------------------------
  //  perform additional verifications for curvilinear LBS
  //  --------------------------------------------------------------------------

  //  coordinate system must be curvilinear
  if (coord_system_type_ != chi_math::CoordinateSystemType::CYLINDRICAL &&
      coord_system_type_ != chi_math::CoordinateSystemType::SPHERICAL)
  {
    chi::log.LogAllError()
      << "C_DO_RZ_SteadyState::SteadyStateSolver::PerformInputChecks : "
      << "invalid coordinate system, static_cast<int>(type) = "
      << static_cast<int>(coord_system_type_);
    chi::Exit(EXIT_FAILURE);
  }

  //  re-interpret geometry type to curvilinear
  switch (options_.geometry_type)
  {
    case lbs::GeometryType::ONED_SLAB:
    {
      switch (coord_system_type_)
      {
        default:
        {
          chi::log.LogAllError()
            << "C_DO_RZ_SteadyState::SteadyStateSolver::PerformInputChecks : "
            << "invalid geometry, static_cast<int>(type) = "
            << static_cast<int>(options_.geometry_type) << " "
            << "for curvilinear coordinate system, static_cast<int>(type) = "
            << static_cast<int>(coord_system_type_);
          chi::Exit(EXIT_FAILURE);
        }
      }
      break;
    }
    case lbs::GeometryType::TWOD_CARTESIAN:
    {
      switch (coord_system_type_)
      {
        case chi_math::CoordinateSystemType::CYLINDRICAL:
        {
          options_.geometry_type = lbs::GeometryType::TWOD_CYLINDRICAL;
          break;
        }
        default:
        {
          chi::log.LogAllError()
            << "C_DO_RZ_SteadyState::SteadyStateSolver::PerformInputChecks : "
            << "invalid geometry, static_cast<int>(type) = "
            << static_cast<int>(options_.geometry_type) << " "
            << "for curvilinear coordinate system, static_cast<int>(type) = "
            << static_cast<int>(coord_system_type_);
          chi::Exit(EXIT_FAILURE);
        }
      }
      break;
    }
    default:
    {
      chi::log.LogAllError()
        << "C_DO_RZ_SteadyState::SteadyStateSolver::PerformInputChecks : "
        << "invalid geometry, static_cast<int>(type) = "
        << static_cast<int>(options_.geometry_type) << " "
        << "for curvilinear coordinate system";
      chi::Exit(EXIT_FAILURE);
    }
  }

  for (size_t gs = 0; gs < groupsets_.size(); ++gs)
  {
    //  angular quadrature type must be compatible with coordinate system
    const auto angular_quad_ptr = groupsets_[gs].quadrature;
    switch (coord_system_type_)
    {
      case chi_math::CoordinateSystemType::CYLINDRICAL:
      {
        typedef chi_math::CylindricalAngularQuadrature CylAngQuad;
        const auto curvilinear_angular_quad_ptr =
          std::dynamic_pointer_cast<CylAngQuad>(angular_quad_ptr);
        if (!curvilinear_angular_quad_ptr)
        {
          chi::log.LogAllError()
            << "C_DO_RZ_SteadyState::SteadyStateSolver::PerformInputChecks : "
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
            << "C_DO_RZ_SteadyState::SteadyStateSolver::PerformInputChecks : "
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
          << "C_DO_RZ_SteadyState::SteadyStateSolver::PerformInputChecks : "
          << "invalid curvilinear coordinate system, static_cast<int>(type) = "
          << static_cast<int>(coord_system_type_);
        chi::Exit(EXIT_FAILURE);
      }
    }

    //  angle aggregation type must be compatible with coordinate system
    const auto angleagg_method = groupsets_[gs].angleagg_method;
    switch (coord_system_type_)
    {
      case chi_math::CoordinateSystemType::CYLINDRICAL:
      {
        if (angleagg_method != lbs::AngleAggregationType::AZIMUTHAL)
        {
          chi::log.LogAllError()
            << "C_DO_RZ_SteadyState::SteadyStateSolver::PerformInputChecks : "
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
            << "C_DO_RZ_SteadyState::SteadyStateSolver::PerformInputChecks : "
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
          << "C_DO_RZ_SteadyState::SteadyStateSolver::PerformInputChecks : "
          << "invalid curvilinear coordinate system, static_cast<int>(type) = "
          << static_cast<int>(coord_system_type_);
        chi::Exit(EXIT_FAILURE);
      }
    }
  }

  //  boundary of mesh must be rectangular with origin at (0, 0, 0)
  const std::vector<chi_mesh::Vector3> unit_normal_vectors =
    { chi_mesh::Vector3(1.0, 0.0, 0.0),
      chi_mesh::Vector3(0.0, 1.0, 0.0),
      chi_mesh::Vector3(0.0, 0.0, 1.0), };
  for (const auto& cell : grid_ptr_->local_cells)
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
              const auto& vertex = grid_ptr_->vertices[v_id];
              if (std::abs(vertex[d]) > 1.0e-12)
              {
                chi::log.LogAllError()
                  << "C_DO_RZ_SteadyState::SteadyStateSolver::PerformInputChecks : "
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
            << "C_DO_RZ_SteadyState::SteadyStateSolver::PerformInputChecks : "
            << "mesh contains boundary faces not orthogonal with respect to "
            << "Cartesian reference frame.";
          chi::Exit(EXIT_FAILURE);
        }
      }
    }
  }

  chi::log.Log() << "C_DO_RZ_SteadyState::SteadyStateSolver::PerformInputChecks : exit";
}


std::shared_ptr<SweepChunk>
lbs_curvilinear::DiscOrdSteadyStateSolver::SetSweepChunk(lbs::LBSGroupset& groupset)
{
   auto sweep_chunk =
     std::make_shared<SweepChunkPWL>(grid_ptr_,
                                     *discretization_,
                                     unit_cell_matrices_,
                                     secondary_unit_cell_matrices_,
                                     cell_transport_views_,
                                     phi_new_local_,
                                     psi_new_local_[groupset.id],
                                     q_moments_local_,
                                     groupset,
                                     matid_to_xs_map_,
                                     num_moments_,
                                     max_cell_dof_count_);

  return sweep_chunk;
}
