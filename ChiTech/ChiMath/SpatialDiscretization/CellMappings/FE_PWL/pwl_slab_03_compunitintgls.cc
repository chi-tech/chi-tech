#include "pwl_slab.h"

//###################################################################
/**Computes cell volume and surface integrals.*/
void SlabMappingFE_PWL::
  ComputeUnitIntegrals(chi_math::finite_element::UnitIntegralData& ui_data) const
{
  typedef std::vector<chi_mesh::Vector3> VecVec3;
  typedef std::vector<VecVec3> MatVec3;

  MatDbl   IntV_gradShapeI_gradShapeJ;
  MatVec3  IntV_shapeI_gradshapeJ    ;
  MatDbl   IntV_shapeI_shapeJ        ;
  VecDbl   IntV_shapeI               ;
  VecVec3  IntV_gradshapeI           ;

  std::vector<MatDbl>  IntS_shapeI_shapeJ    ;
  std::vector<VecDbl>  IntS_shapeI           ;
  std::vector<MatVec3> IntS_shapeI_gradshapeJ;

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
  IntS_shapeI_shapeJ[0][1][1] =  0.0;

  //Right face
  IntS_shapeI_shapeJ[1][0][0] =  0.0;
  IntS_shapeI_shapeJ[1][0][1] =  0.0;
  IntS_shapeI_shapeJ[1][1][0] =  0.0;
  IntS_shapeI_shapeJ[1][1][1] =  1.0;

  //intS_shapeI
  IntS_shapeI.emplace_back(2, 0.0);
  IntS_shapeI.emplace_back(2, 0.0);

  //f then i
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

  ui_data.Initialize(IntV_gradShapeI_gradShapeJ,
                     IntV_shapeI_gradshapeJ,
                     IntV_shapeI_shapeJ    ,
                     IntV_shapeI           ,
                     IntV_gradshapeI       ,
                     IntS_shapeI_shapeJ    ,
                     IntS_shapeI           ,
                     IntS_shapeI_gradshapeJ,
                     face_dof_mappings,
                     num_nodes);

}
