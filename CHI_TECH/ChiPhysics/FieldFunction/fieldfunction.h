#ifndef _chi_field_function_h
#define _chi_field_function_h

#include "../chi_physics_namespace.h"
#include "../../ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "ChiMath/UnknownManager/unknown_manager.h"

#include <petscksp.h>

namespace chi_physics
{
  enum class FieldFunctionType
  {
    CFEM_PWL = 1,
    DFEM_PWL = 2,
    FV = 3
  };
}

//###################################################################
/** Implementation of an abstracted field function object holding
 * a reference to a contiguous vector of data associated with a
 * mesh and a spatial discretization method. This
 * object's design has undergone numerous iterations and is still
 * evolving. The first iteration involved the structuring of groups
 * and moments associated with transport simulations which exhibited
 * significant performance gains when the groups were contiguous. The
 * mapping for each dof was then dof*num_grps*num_moms + num_groups*moment +
 * group. After this it became clear that the design can be generalized
 * even further. The deal.ii library has the concepts of unknown ``components"
 * on each degree of freedom. It also extends the concept to vector
 * valued components.
 *
 * For ChiTech we will generalize a field function to store sets of components.
 * Each set (akin to a moment) will have all its components contiguous and
 * will allow algorithms to to build accordingly.*/
class chi_physics::FieldFunction
{
public:
  std::string               text_name;
  int                       id;
  FieldFunctionType         type;
  chi_mesh::MeshContinuum*  grid;
  SpatialDiscretization*    spatial_discretization;
  chi_math::UnknownManager* unknown_manager;
  int                       num_components, num_sets;
  int                       ref_component, ref_set;

  std::vector<int>*         local_cell_dof_array_address;

  Vec*                      field_vector;
  std::vector<double>*      field_vector_local;
  bool                      using_petsc_field_vector;

private:
  std::vector<double>       temp_cell_dof_values;

public:

  /**Non-PETSc based field function constructor.*/
  FieldFunction(std::string ff_text_name,
                size_t ff_id,
                FieldFunctionType ff_type,
                chi_mesh::MeshContinuum* ff_grid,
                SpatialDiscretization* ff_sdm,
                int ff_num_components,
                int ff_num_sets,
                int ff_ref_component,
                int ff_ref_set,
                std::vector<int>* ff_dof_block_addresses,
                std::vector<double>* ff_field_vector) :
    text_name(std::move(ff_text_name)),
    id(ff_id),
    type(ff_type),
    grid(ff_grid),
    spatial_discretization(ff_sdm),
    num_components(ff_num_components),
    num_sets(ff_num_sets),
    ref_component(ff_ref_component),
    ref_set(ff_ref_set),
    local_cell_dof_array_address(ff_dof_block_addresses),
    field_vector(NULL),
    field_vector_local(ff_field_vector),
    using_petsc_field_vector(false)
  {
    auto uk_man = new chi_math::UnknownManager;

    for (int set=0; set<ff_num_sets; ++set)
      uk_man->AddUnknown(chi_math::UnknownType::VECTOR_N,ff_num_components);
  }

  /**PETSc-vector based field function constructor.*/
  FieldFunction(std::string ff_text_name,
                size_t ff_id,
                FieldFunctionType ff_type,
                chi_mesh::MeshContinuum* ff_grid,
                SpatialDiscretization* ff_sdm,
                int ff_num_components,
                int ff_num_sets,
                int ff_ref_component,
                int ff_ref_set,
                std::vector<int>* ff_dof_block_addresses,
                Vec* ff_field_vector) :
    text_name(std::move(ff_text_name)),
    id(ff_id),
    type(ff_type),
    grid(ff_grid),
    spatial_discretization(ff_sdm),
    num_components(ff_num_components),
    num_sets(ff_num_sets),
    ref_component(ff_ref_component),
    ref_set(ff_ref_set),
    local_cell_dof_array_address(ff_dof_block_addresses),
    field_vector(ff_field_vector),
    field_vector_local(nullptr),
    using_petsc_field_vector(true)
  {
  }

  /***/
  FieldFunction(std::string ff_tex_name,
                FieldFunctionType ff_type,
                SpatialDiscretization* ff_sdm,
                chi_math::UnknownManager* ff_unknown_manager,
                int ff_unknown_id,
                std::vector<double>* ff_field_vector) :
    text_name(ff_tex_name),
    id(0),
    type(ff_type),
    grid(ff_sdm->ref_grid),
    spatial_discretization(ff_sdm),
    unknown_manager(ff_unknown_manager),
    ref_set(ff_unknown_id),
    field_vector_local(ff_field_vector)
  {}

  std::vector<double>& GetCellDOFValues(size_t cell_local_id,
                                        size_t component,
                                        size_t set);

  //01
  void ExportToVTK(const std::string& base_name,
                   const std::string& field_name);
  void ExportToVTKG(const std::string& base_name,
                    const std::string& field_name);
  //01a
  void ExportToVTKFV(const std::string& base_name,
                     const std::string& field_name);
  void ExportToVTKFVG(const std::string& base_name,
                      const std::string& field_name);
  //01b
  void ExportToVTKPWLC(const std::string& base_name,
                       const std::string& field_name);
  void ExportToVTKPWLCG(const std::string& base_name,
                        const std::string& field_name);
  //01c
  void ExportToVTKPWLD(const std::string& base_name,
                       const std::string& field_name);
  void ExportToVTKPWLDG(const std::string& base_name,
                        const std::string& field_name);

  //fieldfunction_exportmultiple_fv.cc
  static void ExportMultipleFFToVTK(const std::string& file_base_name,
                                    std::vector<FieldFunction*> ff_list);

  void WritePVTU(std::string base_filename, std::string field_name, int num_grps=0);
};


#endif