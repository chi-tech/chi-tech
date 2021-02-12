#include "pwl_slab.h"

//###################################################################
/**Computes cell volume and surface integrals.*/
void SlabPWLFEView::ComputeUnitIntegrals()
{
  //gradShapeI_gradShapeJ
  IntV_gradShapeI_gradShapeJ.emplace_back(2, 0.0);
  IntV_gradShapeI_gradShapeJ.emplace_back(2, 0.0);

  IntV_gradShapeI_gradShapeJ[0][0] = 1 / h;
  IntV_gradShapeI_gradShapeJ[0][1] = -1 / h;
  IntV_gradShapeI_gradShapeJ[1][0] = -1 / h;
  IntV_gradShapeI_gradShapeJ[1][1] = 1 / h;

  //shapeI_gradShapeJ
  IntV_shapeI_gradshapeJ.resize(2);
  IntV_shapeI_gradshapeJ[0].resize(2);
  IntV_shapeI_gradshapeJ[1].resize(2);

  IntV_shapeI_gradshapeJ[0][0] = chi_mesh::Vector3(0.0, 0.0, -1 / 2.0);
  IntV_shapeI_gradshapeJ[0][1] = chi_mesh::Vector3(0.0, 0.0, 1 / 2.0);
  IntV_shapeI_gradshapeJ[1][0] = chi_mesh::Vector3(0.0, 0.0, -1 / 2.0);
  IntV_shapeI_gradshapeJ[1][1] = chi_mesh::Vector3(0.0, 0.0, 1 / 2.0);

  //shapeI_shapeJ
  IntV_shapeI_shapeJ.emplace_back(2, 0.0);
  IntV_shapeI_shapeJ.emplace_back(2, 0.0);

  IntV_shapeI_shapeJ[0][0] = h / 3;
  IntV_shapeI_shapeJ[0][1] = h / 6;
  IntV_shapeI_shapeJ[1][0] = h / 6;
  IntV_shapeI_shapeJ[1][1] = h / 3;

  //shapeI
  IntV_shapeI.push_back(h / 2);
  IntV_shapeI.push_back(h / 2);

  //IntS_shapeI_shapeJ
  typedef std::vector<VecDbl> VecVecDbl;
  IntS_shapeI_shapeJ.resize(2, VecVecDbl(2, VecDbl(2, 0.0)));
  IntS_shapeI_shapeJ.resize(2, VecVecDbl(2, VecDbl(2, 0.0)));

  //Left face
  IntS_shapeI_shapeJ[0][0][0] =  1.0;
  IntS_shapeI_shapeJ[0][0][1] =  0.0;
  IntS_shapeI_shapeJ[0][1][0] =  0.0;
  IntS_shapeI_shapeJ[0][1][1] =  1.0;

  //Right face
  IntS_shapeI_shapeJ[1][0][0] =  1.0;
  IntS_shapeI_shapeJ[1][0][1] =  0.0;
  IntS_shapeI_shapeJ[1][1][0] =  0.0;
  IntS_shapeI_shapeJ[1][1][1] =  1.0;

  //intS_shapeI
  IntS_shapeI.emplace_back(2, 0.0);
  IntS_shapeI.emplace_back(2, 0.0);

  IntS_shapeI[0][0] = 1.0;
  IntS_shapeI[0][1] = 0.0;
  IntS_shapeI[1][0] = 0.0;
  IntS_shapeI[1][1] = 1.0;

  //IntS_shapeI_gradshapeJ
  IntS_shapeI_gradshapeJ.resize(2);
  IntS_shapeI_gradshapeJ[0].resize(2);
  IntS_shapeI_gradshapeJ[1].resize(2);

  //Left face
  IntS_shapeI_gradshapeJ[0][0].resize(2);
  IntS_shapeI_gradshapeJ[0][1].resize(2);

  //Right face
  IntS_shapeI_gradshapeJ[1][0].resize(2);
  IntS_shapeI_gradshapeJ[1][1].resize(2);

  //Left face
  IntS_shapeI_gradshapeJ[0][0][0] = chi_mesh::Vector3(0.0, 0.0, -1.0 / h);
  IntS_shapeI_gradshapeJ[0][0][1] = chi_mesh::Vector3(0.0, 0.0, 1.0 / h);
  IntS_shapeI_gradshapeJ[0][1][0] = chi_mesh::Vector3(0.0, 0.0, 0.0  );
  IntS_shapeI_gradshapeJ[0][1][1] = chi_mesh::Vector3(0.0, 0.0, 0.0  );

  //Right face
  IntS_shapeI_gradshapeJ[1][0][0] = chi_mesh::Vector3(0.0, 0.0, 0.0  );
  IntS_shapeI_gradshapeJ[1][0][1] = chi_mesh::Vector3(0.0, 0.0, 0.0  );
  IntS_shapeI_gradshapeJ[1][1][0] = chi_mesh::Vector3(0.0, 0.0, -1.0 / h);
  IntS_shapeI_gradshapeJ[1][1][1] = chi_mesh::Vector3(0.0, 0.0, 1.0 / h);

  precomputed = true;
}

void SlabPWLFEView::PreComputeValues()
{
  ComputeUnitIntegrals();
}