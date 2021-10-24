#ifndef VOLUME_MESHER_EXTRUDER_H
#define VOLUME_MESHER_EXTRUDER_H

#include "../chi_volumemesher.h"
#include "ChiMesh/Cell/cell.h"



//###################################################################
/**An extruder mesher taking a flat surface and extruding it.*/
class chi_mesh::VolumeMesherExtruder : public chi_mesh::VolumeMesher
{
public:
  enum class TemplateType : int
  {
    SURFACE_MESH       = 1,
    UNPARTITIONED_MESH = 2
  };
  struct MeshLayer
  {
    std::string name;
    double height;
    int    sub_divisions;
  };
private:
  const TemplateType template_type;
  SurfaceMesh*       template_surface_mesh = nullptr;
  UnpartitionedMesh* template_unpartitioned_mesh = nullptr;
public:
  std::vector<MeshLayer> input_layers;
  std::vector<double> vertex_layers;
  size_t node_z_index_incr=0;

public:
  explicit
  VolumeMesherExtruder(chi_mesh::SurfaceMesh* in_surface_mesh) :
    template_type(TemplateType::SURFACE_MESH),
    template_surface_mesh(in_surface_mesh)
  {}
  explicit
  VolumeMesherExtruder(chi_mesh::UnpartitionedMesh* in_unpartitioned_mesh) :
    template_type(TemplateType::UNPARTITIONED_MESH),
    template_unpartitioned_mesh(in_unpartitioned_mesh)
  {}
  //02
  void Execute() override;
  //03
  //ReorderDOFs
  //04
  void SetupLayers(int default_layer_count=1);

  void CreateLocalNodes(chi_mesh::MeshContinuum& template_grid,
                        chi_mesh::MeshContinuum& grid);

  void ExtrudeCells(chi_mesh::MeshContinuum& template_grid,
                    chi_mesh::MeshContinuum& grid);

  //utils
  chi_mesh::Vector3 ComputeTemplateCell3DCentroid(
                      const chi_mesh::Cell& n_template_cell,
                      const chi_mesh::MeshContinuum& template_continuum,
                      int z_level_begin,int z_level_end);

  int GetCellPartitionIDFromCentroid(chi_mesh::Vector3& centroid);

  bool IsTemplateCellNeighborToThisPartition(
    const chi_mesh::Cell& template_cell,
    const chi_mesh::MeshContinuum& template_continuum,
    int z_level, int tc_index);




};



#endif //VOLUME_MESHER_EXTRUDER_H