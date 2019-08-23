#ifndef _pwl_cellbase_h
#define _pwl_cellbase_h

//###################################################################
/** Base class for all cell FE views.*/
class CellFEView
{
public:
  int dofs;

  CellFEView(int num_dofs)
  {
    dofs=num_dofs;
  }

};


#endif