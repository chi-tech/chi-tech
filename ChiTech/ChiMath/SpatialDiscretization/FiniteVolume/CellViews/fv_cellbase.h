#ifndef CELL_FVDATA_BASE_H
#define CELL_FVDATA_BASE_H

//######################################################### Class def
/**Base cell class for Finite Volume Method.*/
class CellFVValues
{
public:
  int dofs;
  double                volume=0.0;
  std::vector<double>   face_area = {}; ///< Actually areas

  explicit CellFVValues(int num_dofs) :
    dofs(num_dofs),
    face_area({})
    {}

};

#endif