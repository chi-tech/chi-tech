#include "pwl_slab.h"

//###################################################################
/**Computes cell volume and surface integrals.*/
void SlabPWLFEView::
  ComputeUnitIntegrals(chi_math::finite_element::UnitIntegralData& ui_data)
{
  //gradShapeI_gradShapeJ
  ui_data.IntV_gradShapeI_gradShapeJ.emplace_back(2, 0.0);
  ui_data.IntV_gradShapeI_gradShapeJ.emplace_back(2, 0.0);

  ui_data.IntV_gradShapeI_gradShapeJ[0][0] = 1 / h;
  ui_data.IntV_gradShapeI_gradShapeJ[0][1] = -1 / h;
  ui_data.IntV_gradShapeI_gradShapeJ[1][0] = -1 / h;
  ui_data.IntV_gradShapeI_gradShapeJ[1][1] = 1 / h;

  //shapeI_gradShapeJ
  ui_data.IntV_shapeI_gradshapeJ.resize(2);
  ui_data.IntV_shapeI_gradshapeJ[0].resize(2);
  ui_data.IntV_shapeI_gradshapeJ[1].resize(2);

  ui_data.IntV_shapeI_gradshapeJ[0][0] = chi_mesh::Vector3(0.0, 0.0, -1 / 2.0);
  ui_data.IntV_shapeI_gradshapeJ[0][1] = chi_mesh::Vector3(0.0, 0.0, 1 / 2.0);
  ui_data.IntV_shapeI_gradshapeJ[1][0] = chi_mesh::Vector3(0.0, 0.0, -1 / 2.0);
  ui_data.IntV_shapeI_gradshapeJ[1][1] = chi_mesh::Vector3(0.0, 0.0, 1 / 2.0);

  //shapeI_shapeJ
  ui_data.IntV_shapeI_shapeJ.emplace_back(2, 0.0);
  ui_data.IntV_shapeI_shapeJ.emplace_back(2, 0.0);

  ui_data.IntV_shapeI_shapeJ[0][0] = h / 3;
  ui_data.IntV_shapeI_shapeJ[0][1] = h / 6;
  ui_data.IntV_shapeI_shapeJ[1][0] = h / 6;
  ui_data.IntV_shapeI_shapeJ[1][1] = h / 3;

  //shapeI
  ui_data.IntV_shapeI.push_back(h / 2);
  ui_data.IntV_shapeI.push_back(h / 2);

  //IntS_shapeI_shapeJ
  typedef std::vector<VecDbl> VecVecDbl;
  ui_data.IntS_shapeI_shapeJ.resize(2, VecVecDbl(2, VecDbl(2, 0.0)));
  ui_data.IntS_shapeI_shapeJ.resize(2, VecVecDbl(2, VecDbl(2, 0.0)));

  //Left face
  ui_data.IntS_shapeI_shapeJ[0][0][0] =  1.0;
  ui_data.IntS_shapeI_shapeJ[0][0][1] =  0.0;
  ui_data.IntS_shapeI_shapeJ[0][1][0] =  0.0;
  ui_data.IntS_shapeI_shapeJ[0][1][1] =  1.0;

  //Right face
  ui_data.IntS_shapeI_shapeJ[1][0][0] =  1.0;
  ui_data.IntS_shapeI_shapeJ[1][0][1] =  0.0;
  ui_data.IntS_shapeI_shapeJ[1][1][0] =  0.0;
  ui_data.IntS_shapeI_shapeJ[1][1][1] =  1.0;

  //intS_shapeI
  ui_data.IntS_shapeI.emplace_back(2, 0.0);
  ui_data.IntS_shapeI.emplace_back(2, 0.0);

  //f then i
  ui_data.IntS_shapeI[0][0] = 1.0;
  ui_data.IntS_shapeI[0][1] = 0.0;
  ui_data.IntS_shapeI[1][0] = 0.0;
  ui_data.IntS_shapeI[1][1] = 1.0;

  //IntS_shapeI_gradshapeJ
  ui_data.IntS_shapeI_gradshapeJ.resize(2);
  ui_data.IntS_shapeI_gradshapeJ[0].resize(2);
  ui_data.IntS_shapeI_gradshapeJ[1].resize(2);

  //Left face
  ui_data.IntS_shapeI_gradshapeJ[0][0].resize(2);
  ui_data.IntS_shapeI_gradshapeJ[0][1].resize(2);

  //Right face
  ui_data.IntS_shapeI_gradshapeJ[1][0].resize(2);
  ui_data.IntS_shapeI_gradshapeJ[1][1].resize(2);

  //Left face
  ui_data.IntS_shapeI_gradshapeJ[0][0][0] = chi_mesh::Vector3(0.0, 0.0, -1.0 / h);
  ui_data.IntS_shapeI_gradshapeJ[0][0][1] = chi_mesh::Vector3(0.0, 0.0, 1.0 / h);
  ui_data.IntS_shapeI_gradshapeJ[0][1][0] = chi_mesh::Vector3(0.0, 0.0, 0.0  );
  ui_data.IntS_shapeI_gradshapeJ[0][1][1] = chi_mesh::Vector3(0.0, 0.0, 0.0  );

  //Right face
  ui_data.IntS_shapeI_gradshapeJ[1][0][0] = chi_mesh::Vector3(0.0, 0.0, 0.0  );
  ui_data.IntS_shapeI_gradshapeJ[1][0][1] = chi_mesh::Vector3(0.0, 0.0, 0.0  );
  ui_data.IntS_shapeI_gradshapeJ[1][1][0] = chi_mesh::Vector3(0.0, 0.0, -1.0 / h);
  ui_data.IntS_shapeI_gradshapeJ[1][1][1] = chi_mesh::Vector3(0.0, 0.0, 1.0 / h);

  ui_data.face_dof_mappings = face_dof_mappings;
  ui_data.num_nodes = num_nodes;

}
