#ifndef _pwl_cellbase_h
#define _pwl_cellbase_h

//###################################################################
/** Base class for all cell FE views.*/
class CellFEView
{
public:
  int dofs;

  std::vector<std::vector<double>>              IntV_gradShapeI_gradShapeJ;
  std::vector<std::vector<chi_mesh::Vector>>    IntV_shapeI_gradshapeJ;
  std::vector<std::vector<double>>              IntV_shapeI_shapeJ;
  std::vector<double>                           IntV_shapeI;

  std::vector<std::vector<std::vector<double>>> IntS_shapeI_shapeJ;
  std::vector<std::vector<double>>              IntS_shapeI;
  std::vector<std::vector<std::vector<chi_mesh::Vector>>> IntS_shapeI_gradshapeJ;

  std::vector<std::vector<int>> face_dof_mappings;

  CellFEView(int num_dofs)
  {
    dofs=num_dofs;
  }

};


#endif