#ifndef _cell_triangle_h
#define _cell_triangle_h

#include"cell.h"
#include"../CHI_REGION/chi_region.h"

//######################################################### Class def
/***/
class chi_mesh::CellTriangle : public chi_mesh::Cell
{
public:
  int v_index[3];
  int n_index[3];
  int e_index[3][4];

public:
  CellTriangle():Cell()
  {
    for (int k=0;k<3;k++)
    {
      v_index[k]=-1;
      n_index[k]=-1;
      e_index[k][0]=-1;
      e_index[k][1]=-1;
      e_index[k][2]=-1;
      e_index[k][3]=-1;
    }
  }

  //01
  void FindBoundary2D(chi_mesh::Region* region);
  //02
  bool CheckBoundary2D();
};

#endif