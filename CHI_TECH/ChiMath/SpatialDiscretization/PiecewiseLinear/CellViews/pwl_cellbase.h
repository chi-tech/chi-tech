#ifndef PWL_CELL_VALUES_BASE_H
#define PWL_CELL_VALUES_BASE_H

#include <ChiMesh/chi_mesh.h>

//###################################################################
/** Base class for all cell FE views.*/
class CellPWLFEValues
{
public:
  const int dofs;

  std::vector<std::vector<double>>              IntV_gradShapeI_gradShapeJ;
  std::vector<std::vector<chi_mesh::Vector3>>   IntV_shapeI_gradshapeJ;
  std::vector<std::vector<double>>              IntV_shapeI_shapeJ;
  std::vector<double>                           IntV_shapeI;

  std::vector<std::vector<std::vector<double>>> IntS_shapeI_shapeJ;
  std::vector<std::vector<double>>              IntS_shapeI;
  std::vector<std::vector<std::vector<chi_mesh::Vector3>>> IntS_shapeI_gradshapeJ;

  std::vector<std::vector<int>> face_dof_mappings;

  explicit CellPWLFEValues(int num_dofs) :
    dofs(num_dofs)
  {}

  virtual ~CellPWLFEValues() = default;

  /** Virtual function evaluation of the shape function. */
  virtual double ShapeValue(const int i, const chi_mesh::Vector3& xyz)
  {
    return 0.0;
  }

  /** Virtual function returning the all the shape function evaluations
   * at the point.*/
  virtual void ShapeValues(const chi_mesh::Vector3& xyz,
                           std::vector<double>& shape_values)
  {
    shape_values.resize(dofs,0.0);
  }

  /** Virtual function evaluation of the grad-shape function. */
  virtual chi_mesh::Vector3 GradShapeValue(const int i,
                                           const chi_mesh::Vector3& xyz)
  {
    return chi_mesh::Vector3(0.0, 0.0, 0.0);
  }

  /** Virtual function evaluation of the grad-shape function. */
  virtual void GradShapeValues(const chi_mesh::Vector3& xyz,
                               std::vector<chi_mesh::Vector3>& gradshape_values)
  {
    gradshape_values.resize(dofs,chi_mesh::Vector3());
  }

};


#endif