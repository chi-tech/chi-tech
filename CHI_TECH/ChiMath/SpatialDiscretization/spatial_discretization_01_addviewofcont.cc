#include"spatial_discretization.h"

void SpatialDiscretization::AddViewOfLocalContinuum(
  chi_mesh::MeshContinuum* vol_continuum,
  int num_cells,
  int* cell_indices)
{
  //This function is meant to be overwritten
}

void SpatialDiscretization::AddViewOfLocalContinuum(
  chi_mesh::MeshContinuum* vol_continuum)
{
  //This function is meant to be overwritten
}