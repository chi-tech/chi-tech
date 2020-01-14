#ifndef _fv_cellbase_h
#define _fv_cellbase_h

//######################################################### Class def
/**Base cell class for Finite Volume Method.*/
class CellFVView
{
public:
  int dofs;
  double                volume;
  std::vector<double>   face_area; ///< Actually length unity h

  CellFVView(int num_dofs)
  {
    dofs=num_dofs;
  }

};

#endif