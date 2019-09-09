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
  std::vector<double*>                          IntV_gradShapeI_gradShapeJ;
  std::vector<std::vector<chi_mesh::Vector>>    IntV_shapeI_gradshapeJ;
  std::vector<std::vector<double>>              IntV_shapeI_shapeJ;
  std::vector<double>                           IntV_shapeI;
  std::vector<double*>                          IntS_shapeI;
  std::vector<std::vector<double*>>             IntS_shapeI_shapeJ;
  std::vector<std::vector<std::vector<chi_mesh::Vector>>> IntS_shapeI_gradshapeJ;
  double h;
public:

  /**Constructor for a slab view.*/
  SlabFEView(chi_mesh::CellSlab *slab_cell,
                         chi_mesh::MeshContinuum *vol_continuum) :
                         CellFEView(2)
  {
    grid = vol_continuum;
    v0i = slab_cell->v_indices[0];
    v1i = slab_cell->v_indices[1];
    chi_mesh::Vertex v0 = *grid->nodes[v0i];
    chi_mesh::Vertex v1 = *grid->nodes[v1i];

    chi_mesh::Vector v01 = v1-v0;
    h = v01.Norm();

    IntV_shapeI.push_back(h/2);
    IntV_shapeI.push_back(h/2);

    IntV_shapeI_shapeJ.push_back(std::vector<double>(2,0.0));
    IntV_shapeI_shapeJ.push_back(std::vector<double>(2,0.0));

    IntV_shapeI_shapeJ[0][0] = h/3;
    IntV_shapeI_shapeJ[0][1] = h/6;
    IntV_shapeI_shapeJ[1][0] = h/6;
    IntV_shapeI_shapeJ[1][1] = h/3;

    IntV_gradShapeI_gradShapeJ.push_back(new double[2]);
    IntV_gradShapeI_gradShapeJ.push_back(new double[2]);

    IntV_gradShapeI_gradShapeJ[0][0] = 1/h;
    IntV_gradShapeI_gradShapeJ[0][1] = -1/h;
    IntV_gradShapeI_gradShapeJ[1][0] = -1/h;
    IntV_gradShapeI_gradShapeJ[1][1] = 1/h;

//    IntV_shapeI_gradshapeJ.push_back(new double[2]);
//    IntV_shapeI_gradshapeJ.push_back(new double[2]);

    IntV_shapeI_gradshapeJ.resize(2);
    IntV_shapeI_gradshapeJ[0].resize(2);
    IntV_shapeI_gradshapeJ[1].resize(2);

    IntV_shapeI_gradshapeJ[0][0] = chi_mesh::Vector(0.0,0.0,-1/2.0);
    IntV_shapeI_gradshapeJ[0][1] = chi_mesh::Vector(0.0,0.0, 1/2.0);
    IntV_shapeI_gradshapeJ[1][0] = chi_mesh::Vector(0.0,0.0,-1/2.0);
    IntV_shapeI_gradshapeJ[1][1] = chi_mesh::Vector(0.0,0.0, 1/2.0);

    IntS_shapeI.push_back(new double[2]);
    IntS_shapeI.push_back(new double[2]);

    IntS_shapeI[0][0] = 1.0;
    IntS_shapeI[0][1] = 0.0;
    IntS_shapeI[1][0] = 0.0;
    IntS_shapeI[1][1] = 1.0;




    IntS_shapeI_shapeJ.resize(2,std::vector<double*>(2));
    IntS_shapeI_shapeJ.resize(2,std::vector<double*>(2));

    //Left face
    IntS_shapeI_shapeJ[0][0]= new double[2];
    IntS_shapeI_shapeJ[0][1]= new double[2];

    //Right face
    IntS_shapeI_shapeJ[1][0]= new double[2];
    IntS_shapeI_shapeJ[1][1]= new double[2];

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
    IntS_shapeI_gradshapeJ[0][0][0] = chi_mesh::Vector(0.0,0.0,-1.0/h);
    IntS_shapeI_gradshapeJ[0][0][1] = chi_mesh::Vector(0.0,0.0, 1.0/h);
    IntS_shapeI_gradshapeJ[0][1][0] = chi_mesh::Vector(0.0,0.0, 0.0  );
    IntS_shapeI_gradshapeJ[0][1][1] = chi_mesh::Vector(0.0,0.0, 0.0  );

    //Right face
    IntS_shapeI_gradshapeJ[1][0][0] = chi_mesh::Vector(0.0,0.0, 0.0  );
    IntS_shapeI_gradshapeJ[1][0][1] = chi_mesh::Vector(0.0,0.0, 0.0  );
    IntS_shapeI_gradshapeJ[1][1][0] = chi_mesh::Vector(0.0,0.0,-1.0/h);
    IntS_shapeI_gradshapeJ[1][1][1] = chi_mesh::Vector(0.0,0.0, 1.0/h);


  }

  /**Shape function for the slab function.*/
  double Shape_x(int i, chi_mesh::Vector xyz)
  {
    chi_mesh::Vector p0 = *grid->nodes[v0i];
    chi_mesh::Vector p1 = *grid->nodes[v1i];
    chi_mesh::Vector xyz_ref = xyz - p0;

    chi_mesh::Vector v01 = p1 - p0;

    double xi   = v01.Dot(xyz_ref)/v01.Norm();

    if ((xi>=-1.0e-12) and (xi<=1.0+1.0e-12))
    {
      if (i==0)
      {
        return 1.0 - xi/h;
      }
      else
      {
        return xi/h;
      }
    }

    return 0.0;
  }
};
#endif
