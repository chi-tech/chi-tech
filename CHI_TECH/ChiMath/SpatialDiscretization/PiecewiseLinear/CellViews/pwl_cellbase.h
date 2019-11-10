#ifndef _pwl_cellbase_h
#define _pwl_cellbase_h

#include <ChiMesh/chi_mesh.h>

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

  /** Virtual function evaluation of the shape function. */
  virtual double ShapeValue(int i, const chi_mesh::Vector& xyz)
  {
    return 0.0;
  }

  /** Virtual function returning the all the shape function evaluations
   * at the point.*/
   virtual std::vector<double> ShapeValues(const chi_mesh::Vector& xyz)
  {
     return std::vector<double>(dofs,0.0);
  }

};


#endif