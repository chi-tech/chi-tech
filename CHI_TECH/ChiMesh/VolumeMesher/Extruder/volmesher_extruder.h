#ifndef _volmesher_extruder_h
#define _volmesher_extruder_h

#include "../chi_volumemesher.h"

struct MeshLayer
{
  std::string name;
  double height;
  int    sub_divisions;
  int    major_id;
  int    minor_id;
};

//###################################################################
/**An extruder mesher taking a flat surface and extruding it.*/
class chi_mesh::VolumeMesherExtruder : public chi_mesh::VolumeMesher
{
public:
  std::vector<MeshLayer*> input_layers;
  std::vector<double> vertex_layers;
  int node_z_index_incr;

  int bot_boundary_index;
  int top_boundary_index;

public:
  //02
  void Execute();
  //03
  //ReorderDOFs
  //04
  void SetupLayers(int default_layer_count=1);
  //05
  void ExtrudeCells(chi_mesh::MeshContinuum* template_continuum,
                    chi_mesh::MeshContinuum* vol_continuum);

};



#endif