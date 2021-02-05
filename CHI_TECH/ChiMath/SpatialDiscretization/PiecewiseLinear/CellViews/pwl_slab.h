#ifndef _pwl_slab_h
#define _pwl_slab_h

#include "../pwl.h"
#include <vector>
#include <ChiMesh/Cell/cell_slab.h>

//###################################################################
/**Object for handling slab shaped piecewise linear shape functions.*/
class SlabFEView : public CellFEView
{
private:
  chi_mesh::MeshContinuum* grid;
  int v0i;
  int v1i;
public:
  double h;
public:

  /**Constructor for a slab view.*/
  SlabFEView(chi_mesh::CellSlab *slab_cell,
             chi_mesh::MeshContinuum* vol_continuum) :
    CellFEView(2)
  {
    grid = vol_continuum;
    v0i = slab_cell->vertex_ids[0];
    v1i = slab_cell->vertex_ids[1];
    chi_mesh::Vertex v0 = *grid->vertices[v0i];
    chi_mesh::Vertex v1 = *grid->vertices[v1i];

    chi_mesh::Vector3 v01 = v1 - v0;
    h = v01.Norm();

    IntV_shapeI.push_back(h/2);
    IntV_shapeI.push_back(h/2);

    IntV_shapeI_shapeJ.emplace_back(2, 0.0);
    IntV_shapeI_shapeJ.emplace_back(2, 0.0);

    IntV_shapeI_shapeJ[0][0] = h/3;
    IntV_shapeI_shapeJ[0][1] = h/6;
    IntV_shapeI_shapeJ[1][0] = h/6;
    IntV_shapeI_shapeJ[1][1] = h/3;

    IntV_gradShapeI_gradShapeJ.emplace_back(2, 0.0);
    IntV_gradShapeI_gradShapeJ.emplace_back(2, 0.0);

    IntV_gradShapeI_gradShapeJ[0][0] = 1/h;
    IntV_gradShapeI_gradShapeJ[0][1] = -1/h;
    IntV_gradShapeI_gradShapeJ[1][0] = -1/h;
    IntV_gradShapeI_gradShapeJ[1][1] = 1/h;

//    IntV_shapeI_gradshapeJ.push_back(new double[2]);
//    IntV_shapeI_gradshapeJ.push_back(new double[2]);

    IntV_shapeI_gradshapeJ.resize(2);
    IntV_shapeI_gradshapeJ[0].resize(2);
    IntV_shapeI_gradshapeJ[1].resize(2);

    IntV_shapeI_gradshapeJ[0][0] = chi_mesh::Vector3(0.0, 0.0, -1 / 2.0);
    IntV_shapeI_gradshapeJ[0][1] = chi_mesh::Vector3(0.0, 0.0, 1 / 2.0);
    IntV_shapeI_gradshapeJ[1][0] = chi_mesh::Vector3(0.0, 0.0, -1 / 2.0);
    IntV_shapeI_gradshapeJ[1][1] = chi_mesh::Vector3(0.0, 0.0, 1 / 2.0);

    IntS_shapeI.emplace_back(2, 0.0);
    IntS_shapeI.emplace_back(2, 0.0);

    IntS_shapeI[0][0] = 1.0;
    IntS_shapeI[0][1] = 0.0;
    IntS_shapeI[1][0] = 0.0;
    IntS_shapeI[1][1] = 1.0;

    typedef std::vector<double> VecDbl;
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

    face_dof_mappings.emplace_back(1,0);
    face_dof_mappings.emplace_back(1,1);

  }

  /**Shape function i evaluated at given point for the slab.*/
  double ShapeValue(const int i, const chi_mesh::Vector3& xyz) override
  {
    chi_mesh::Vector3& p0 = *grid->vertices[v0i];
    chi_mesh::Vector3& p1 = *grid->vertices[v1i];
    chi_mesh::Vector3 xyz_ref = xyz - p0;

    chi_mesh::Vector3 v01 = p1 - p0;

    double xi   = v01.Dot(xyz_ref)/v01.Norm()/h;

    if ((xi>=-1.0e-6) and (xi<=1.0+1.0e-6))
    {
      if (i==0)
      {
        return 1.0 - xi;
      }
      else
      {
        return xi;
      }
    }

    return 0.0;
  }

  //#################################################################
  /**Populates shape_values with the value of each shape function's
   * value evaluate at the supplied point.*/
  void ShapeValues(const chi_mesh::Vector3& xyz,
                   std::vector<double>& shape_values) override
  {
    shape_values.resize(dofs,0.0);
    chi_mesh::Vector3& p0 = *grid->vertices[v0i];
    chi_mesh::Vector3& p1 = *grid->vertices[v1i];
    chi_mesh::Vector3 xyz_ref = xyz - p0;

    chi_mesh::Vector3 v01 = p1 - p0;

    double xi   = v01.Dot(xyz_ref)/v01.Norm()/h;

    if ((xi>=-1.0e-6) and (xi<=1.0+1.0e-6))
    {
      for (int i=0; i<dofs; i++)
      {
        if (i==0)
        {
          shape_values[i] = 1.0 - xi;
        }
        else
        {
          shape_values[i] = xi;
        }
      }//for dof

      return;
    }//if in cell

  }

  //###################################################################
  /**Returns the evaluation of grad-shape function i at the supplied point.*/
  chi_mesh::Vector3 GradShapeValue(const int i, const chi_mesh::Vector3& xyz) override
  {
    if (i==0)
      return chi_mesh::Vector3(0.0, 0.0, -1.0 / h);
    else
      return chi_mesh::Vector3(0.0, 0.0, 1.0 / h);
  }

  //###################################################################
  /**Populates shape_values with the value of each shape function's
   * value evaluate at the supplied point.*/
  void GradShapeValues(const chi_mesh::Vector3& xyz,
                       std::vector<chi_mesh::Vector3>& gradshape_values) override
  {
    gradshape_values.clear();
    gradshape_values.emplace_back(GradShapeValue(0,xyz));
    gradshape_values.emplace_back(GradShapeValue(1,xyz));
  }
};
#endif
