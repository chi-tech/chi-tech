#include "pwlc.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "chi_mpi.h"

#include "ChiMath/PETScUtils/petsc_utils.h"
#define sc_int64 static_cast<int64_t>

//###################################################################
/**Get the number of local degrees-of-freedom.*/
size_t chi_math::SpatialDiscretization_PWLC::
  GetNumLocalDOFs(const chi_math::UnknownManager& unknown_manager) const
{
  unsigned int N = unknown_manager.GetTotalUnknownStructureSize();

  return local_base_block_size*N;
}

//###################################################################
/**Get the number of global degrees-of-freedom.*/
size_t chi_math::SpatialDiscretization_PWLC::
  GetNumGlobalDOFs(const chi_math::UnknownManager& unknown_manager) const
{
  unsigned int N = unknown_manager.GetTotalUnknownStructureSize();

  return globl_base_block_size*N;
}

//###################################################################
/**Get the number of ghost degrees-of-freedom.*/
size_t chi_math::SpatialDiscretization_PWLC::
  GetNumGhostDOFs(const chi_math::UnknownManager& unknown_manager) const
{
  unsigned int N = unknown_manager.GetTotalUnknownStructureSize();

  return m_ghost_node_mapping.size()*N;
}

//###################################################################
/**Returns the ghost DOF indices.*/
std::vector<int64_t> chi_math::SpatialDiscretization_PWLC::
  GetGhostDOFIndices(const chi_math::UnknownManager& unknown_manager) const
{
  std::vector<int64_t> dof_ids;
  dof_ids.reserve(GetNumGhostDOFs(unknown_manager));

  size_t num_unknowns = unknown_manager.GetTotalUnknownStructureSize();
  auto   storage      = unknown_manager.dof_storage_type;

  for (const auto& vid_gnid : m_ghost_node_mapping)
  {
    const int64_t global_id = vid_gnid.second;

    for (size_t u=0; u<num_unknowns; ++u)
    {
      const auto& unkn = unknown_manager.unknowns[u];
      const size_t num_comps = unkn.num_components;
      for (size_t c=0; c<num_comps; ++c)
      {
        size_t block_id     = unknown_manager.MapUnknown(u, c);
        int64_t address=-1;
        if (storage == chi_math::UnknownStorageType::BLOCK)
        {
          for (int locJ=0; locJ<chi::mpi.process_count; ++locJ)
          {
            const int64_t local_id = global_id - sc_int64(locJ_block_address[locJ]);

            if (local_id < 0 or local_id >= locJ_block_size[locJ]) continue;

            address = sc_int64(locJ_block_address[locJ]*num_unknowns) +
                      sc_int64(locJ_block_size[locJ]*block_id) +
                      local_id;
            break;
          }
        }
        else if (storage == chi_math::UnknownStorageType::NODAL)
          address = global_id * sc_int64(num_unknowns) + sc_int64(block_id);

        dof_ids.push_back(address);
      }//for c
    }//for u
  }

  return dof_ids;
}

size_t chi_math::SpatialDiscretization_PWLC::
  GetCellNumNodes(const chi_mesh::Cell& cell) const
{
  return cell.vertex_ids.size();
}

std::vector<chi_mesh::Vector3> chi_math::SpatialDiscretization_PWLC::
  GetCellNodeLocations(const chi_mesh::Cell& cell) const
{
  std::vector<chi_mesh::Vector3> node_locations;
  node_locations.reserve(cell.vertex_ids.size());

  for (auto& vid : cell.vertex_ids)
    node_locations.emplace_back(ref_grid->vertices[vid]);

  return node_locations;
}

//###################################################################
/**Develops a localized view of a petsc vector.*/
void chi_math::SpatialDiscretization_PWLC::
LocalizePETScVector(Vec petsc_vector,
                    std::vector<double>& local_vector,
                    const chi_math::UnknownManager& unknown_manager) const
{
  size_t num_local_dofs = GetNumLocalDOFs(unknown_manager);

  chi_math::PETScUtils::CopyVecToSTLvector(petsc_vector,
                                           local_vector,
                                           num_local_dofs);
}

const chi_math::finite_element::UnitIntegralData&
chi_math::SpatialDiscretization_PWLC::
  GetUnitIntegrals(const chi_mesh::Cell& cell)
{
  if (ref_grid->IsCellLocal(cell.global_id))
  {
    if (integral_data_initialized)
      return fe_unit_integrals.at(cell.local_id);
    else
    {
      const auto& cell_mapping = GetCellMapping(cell);
      scratch_intgl_data.Reset();
      cell_mapping.ComputeUnitIntegrals(scratch_intgl_data);
      return scratch_intgl_data;
    }
  }
  else
  {
    if (nb_integral_data_initialized)
      return nb_fe_unit_integrals.at(cell.global_id);
    else
    {
      const auto& cell_mapping = GetCellMapping(cell);
      cell_mapping.ComputeUnitIntegrals(scratch_intgl_data);
      return scratch_intgl_data;
    }
  }
}

const chi_math::finite_element::InternalQuadraturePointData&
chi_math::SpatialDiscretization_PWLC::
  GetQPData_Volumetric(const chi_mesh::Cell& cell)
{
  if (ref_grid->IsCellLocal(cell.global_id))
  {
    if (qp_data_initialized)
      return fe_vol_qp_data.at(cell.local_id);
    else
    {
      const auto& cell_mapping = GetCellMapping(cell);
      cell_mapping.InitializeVolumeQuadraturePointData(scratch_vol_qp_data);
      return scratch_vol_qp_data;
    }
  }
  else
  {
    if (nb_qp_data_initialized)
      return nb_fe_vol_qp_data.at(cell.global_id);
    else
    {
      const auto& cell_mapping = GetCellMapping(cell);
      cell_mapping.InitializeVolumeQuadraturePointData(scratch_vol_qp_data);
      return scratch_vol_qp_data;
    }
  }
}

const chi_math::finite_element::FaceQuadraturePointData&
chi_math::SpatialDiscretization_PWLC::
  GetQPData_Surface(const chi_mesh::Cell& cell,
                  const unsigned int face)
{
  if (ref_grid->IsCellLocal(cell.global_id))
  {
    if (qp_data_initialized)
    {
      const auto& face_data = fe_srf_qp_data.at(cell.local_id);

      return face_data.at(face);
    }
    else
    {
      const auto& cell_mapping = GetCellMapping(cell);
      cell_mapping.InitializeFaceQuadraturePointData(face, scratch_face_qp_data);
      return scratch_face_qp_data;
    }
  }
  else
  {
    if (nb_qp_data_initialized)
    {
      const auto& face_data = nb_fe_srf_qp_data.at(cell.global_id);

      return face_data.at(face);
    }
    else
    {
      const auto& cell_mapping = GetCellMapping(cell);
      cell_mapping.InitializeFaceQuadraturePointData(face, scratch_face_qp_data);
      return scratch_face_qp_data;
    }
  }
}