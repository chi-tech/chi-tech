#include "sldfe_sq.h"

//###################################################################
/**Prints the quadrature to file.*/
void chi_math::SimplifiedLDFESQ::Quadrature::PrintQuadratureToFile()
{
  std::ofstream vert_file,cell_file,points_file;
  vert_file.open(output_filename_prefix + "verts.txt");
  {
    for (const auto& sq : deployed_SQs)
      for (const auto& vert : sq.vertices_xyz)
        vert_file << vert.x << " " << vert.y << " " << vert.z << "\n";
  }
  vert_file.close();

  cell_file.open(output_filename_prefix + "cells.txt");
  {
    int vi=0;
    for (const auto& sq : deployed_SQs)
    {
      for (const auto& vert : sq.vertices_xyz)
        cell_file << vi++ << " ";
      cell_file << "\n";
    }
  }
  cell_file.close();

  double total_weight=0.0;
  points_file.open(output_filename_prefix + "points.txt");
  {
    for (auto& sq : deployed_SQs)
    {
      int ss=-1;
      for (const auto& point : sq.sub_sqr_points)
      {
        ++ss;
        for (int i=0; i<3; ++i)
          points_file << point[i] << " ";
        points_file << sq.sub_sqr_weights[ss];
        total_weight += sq.sub_sqr_weights[ss];
        points_file << "\n";
      }
    }
  }
  points_file.close();
}