#include "pwl.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiMath/PETScUtils/petsc_utils.h"

//###################################################################
/**Get the number of local degrees-of-freedom.*/
size_t chi_math::SpatialDiscretization_PWLD::
  GetNumLocalDOFs(const chi_math::UnknownManager& unknown_manager) const
{
  unsigned int N = unknown_manager.GetTotalUnknownStructureSize();

  return local_base_block_size*N;
}

//###################################################################
/**Get the number of global degrees-of-freedom.*/
size_t chi_math::SpatialDiscretization_PWLD::
  GetNumGlobalDOFs(const chi_math::UnknownManager& unknown_manager) const
{
  unsigned int N = unknown_manager.GetTotalUnknownStructureSize();

  return globl_base_block_size*N;
}

size_t chi_math::SpatialDiscretization_PWLD::
GetNumGhostDOFs(const UnknownManager& unknown_manager) const
{
  return 0;
}
std::vector<int64_t> chi_math::SpatialDiscretization_PWLD::
GetGhostDOFIndices(const UnknownManager& unknown_manager) const
{
  return {};
}

size_t chi_math::SpatialDiscretization_PWLD::
GetCellNumNodes(const chi_mesh::Cell& cell) const
{
  return cell.vertex_ids.size();
}

std::vector<chi_mesh::Vector3> chi_math::SpatialDiscretization_PWLD::
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
void chi_math::SpatialDiscretization_PWLD::
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
  chi_math::SpatialDiscretization_PWLD::
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
  chi_math::SpatialDiscretization_PWLD::
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
  chi_math::SpatialDiscretization_PWLD::
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