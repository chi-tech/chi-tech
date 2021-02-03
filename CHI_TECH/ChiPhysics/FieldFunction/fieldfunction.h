#ifndef CHI_FIELD_FUNCTION_H
#define CHI_FIELD_FUNCTION_H

#include "../chi_physics_namespace.h"
#include "../../ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "ChiMath/UnknownManager/unknown_manager.h"

#include <petscksp.h>

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
//  int                       id;
  chi_mesh::MeshContinuum*  grid;
  SpatialDiscretization*    spatial_discretization;
  chi_math::UnknownManager  unknown_manager;
  const unsigned int        ref_component;
  const unsigned int        ref_unknown;

  Vec*                      field_vector;
  std::vector<double>*      field_vector_local;
  bool                      using_petsc_field_vector;

public:
  /**Non-PETSc type field vector*/
  FieldFunction(std::string ff_tex_name,
                SpatialDiscretization* ff_sdm,
                std::vector<double>* ff_field_vector,
                chi_math::UnknownManager& ff_unknown_manager,
                int ff_unknown_id=0,
                int ff_unknown_component_number=0) :
    text_name(std::move(ff_tex_name)),
//    id(0),
    grid(ff_sdm->ref_grid),
    spatial_discretization(ff_sdm),
    unknown_manager(ff_unknown_manager),
    ref_component(ff_unknown_component_number),
    ref_unknown(ff_unknown_id),
    field_vector(nullptr),
    field_vector_local(ff_field_vector),
    using_petsc_field_vector(false)
  {}

  /**PETSc type field vector*/
  FieldFunction(std::string ff_tex_name,
                SpatialDiscretization* ff_sdm,
                Vec* ff_field_vector,
                chi_math::UnknownManager& ff_unknown_manager,
                int ff_unknown_id=0,
                int ff_unknown_component_number=0) :
    text_name(std::move(ff_tex_name)),
//    id(0),
    grid(ff_sdm->ref_grid),
    spatial_discretization(ff_sdm),
    unknown_manager(ff_unknown_manager),
    ref_component(ff_unknown_component_number),
    ref_unknown(ff_unknown_id),
    field_vector(ff_field_vector),
    field_vector_local(nullptr),
    using_petsc_field_vector(true)
  {}

//  std::vector<double>& GetCellDOFValues(size_t cell_local_id,
//                                        size_t component,
//                                        size_t set);

  //mapping
  void
  CreateFVMappingLocal(std::vector<std::pair<uint64_t,uint>>& cell_component_pairs,
                       std::vector<uint64_t>& mapping);

  void
  CreateCFEMMappingLocal(Vec& x_mapped,
                         std::vector<std::pair<uint64_t,uint>>& node_component_pairs,
                         std::vector<uint64_t>& mapping);

  void
  CreatePWLDMappingLocal(
    std::vector<std::tuple<uint64_t,uint,uint>>& cell_node_component_tuples,
    std::vector<uint64_t>& mapping);

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