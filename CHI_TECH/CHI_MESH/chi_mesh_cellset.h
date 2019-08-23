#ifndef _chi_mesh_cellset_h
#define _chi_mesh_cellset_h



struct chi_mesh::CELL_SET
{
  int i;
  int j;
  int k;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double zmin;
  double zmax;
  std::vector<int> cells_allocated;
  chi_mesh::MeshContinuum* mesh_continuum;

  CELL_SET()
  {
    i=-1;
    j=-1;
    k=-1;
    xmin = -1.0e6;
    xmax =  1.0e6;
    ymin = -1.0e6;
    ymax =  1.0e6;
    zmin = -1.0e6;
    zmax =  1.0e6;

    mesh_continuum = nullptr;
  }
};

#endif