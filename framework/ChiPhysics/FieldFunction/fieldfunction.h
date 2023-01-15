#ifndef CHI_FIELD_FUNCTION_H
#define CHI_FIELD_FUNCTION_H

#include "../chi_physics_namespace.h"
#include "../../ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "ChiMath/UnknownManager/unknown_manager.h"

#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>

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
protected:
  std::string               text_name;
public:
  std::shared_ptr<chi_math::SpatialDiscretization> spatial_discretization;
  chi_math::UnknownManager  unknown_manager;
  const unsigned int        ref_component;
  const unsigned int        ref_variable;

  Vec*                      field_vector;
  std::vector<double>*      field_vector_local;

public:
  /**Non-PETSc type field vector*/
  FieldFunction(std::string ff_tex_name,
                std::shared_ptr<chi_math::SpatialDiscretization>& ff_sdm,
                std::vector<double>* ff_field_vector,
                chi_math::UnknownManager& ff_unknown_manager,
                int ff_unknown_id=0,
                int ff_unknown_component_number=0) :
    text_name(std::move(ff_tex_name)),
    spatial_discretization(ff_sdm),
    unknown_manager(ff_unknown_manager),
    ref_component(ff_unknown_component_number),
    ref_variable(ff_unknown_id),
    field_vector(nullptr),
    field_vector_local(ff_field_vector)
  {}

  /**PETSc type field vector*/
  FieldFunction(std::string ff_tex_name,
                std::shared_ptr<chi_math::SpatialDiscretization>& ff_sdm,
                Vec* ff_field_vector,
                chi_math::UnknownManager& ff_unknown_manager,
                int ff_unknown_id=0,
                int ff_unknown_component_number=0) :
    text_name(std::move(ff_tex_name)),
    spatial_discretization(ff_sdm),
    unknown_manager(ff_unknown_manager),
    ref_component(ff_unknown_component_number),
    ref_variable(ff_unknown_id),
    field_vector(ff_field_vector),
    field_vector_local(nullptr)
  {}

  //Getters
  const std::string& TextName() const {return text_name;}

  //mapping
  void
  CreateFVMappingLocal(std::vector<std::pair<uint64_t,uint>>& cell_component_pairs,
                       std::vector<uint64_t>& mapping);

  void
  CreateCFEMMappingLocal(Vec& x_mapped,
                         std::vector<std::tuple<uint64_t,uint,uint>>& cell_node_component_tuples,
                         std::vector<uint64_t>& mapping);

  void
  CreatePWLDMappingLocal(
    std::vector<std::tuple<uint64_t,uint,uint>>& cell_node_component_tuples,
    std::vector<uint64_t>& mapping);

  //01
  void ExportToVTKComponentOnly(const std::string& base_name,
                                const std::string& field_name);
  void ExportToVTK(const std::string& base_name,
                   const std::string& field_name);
  //01a
  void ExportToVTKFV(const std::string& base_name,
                     const std::string& field_name,
                     bool all_components=false);
  //01b
  void ExportToVTKPWLC(const std::string& base_name,
                       const std::string& field_name,
                       bool all_components=false);
  //01c
  void ExportToVTKPWLD(const std::string& num_nodes,
                       const std::string& field_name,
                       bool all_components=false);

  //fieldfunction_exportmultiple_fv.cc
  static void ExportMultipleFFToVTK(const std::string& file_base_name,
                                    const std::vector<std::shared_ptr<chi_physics::FieldFunction>>& ff_list);


  static void WritePVTU(const std::string& base_filename,
                        const std::string& field_name,
                        const std::vector<std::string>& component_names);

  static
  void UploadCellGeometry(const chi_mesh::MeshContinuum& grid,
                          const chi_mesh::Cell& cell,
                          int64_t& node_counter,
                          vtkNew<vtkPoints>& points,
                          vtkNew<vtkUnstructuredGrid>& ugrid);

public:
  std::vector<double> GetGhostedFieldVector() const;
};


#endif