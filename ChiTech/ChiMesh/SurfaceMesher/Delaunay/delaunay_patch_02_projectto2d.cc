#include "delaunay_mesher.h"

/**Develops Pstar, a projection of the patch's vertices to 2D.*/
void chi_mesh::SurfaceMesherDelaunay::DelaunayPatch::ProjectTo2D()
{
  //================================================== Copy all vertices to Pstar
  //                                                   and comp centroid
  std::vector<chi_mesh::Vertex>::iterator curvert;
  for (curvert = vertices.begin();
       curvert != vertices.end();
       curvert++)
  {
    Pstar.push_back(*curvert);
    centroid = centroid + (*curvert);
  }
  centroid = centroid/vertices.size();
  //printf("Centroid %+.4f %+.4f %+.4f\n",centroid.x,centroid.y,centroid.z);

  //================================================== Set hat_k
  hat_k = this->patch_normal;
  hat_k = hat_k/hat_k.Norm();
  //printf("Hat_k %+.4f %+.4f %+.4f\n",hat_k.x,hat_k.y,hat_k.z);

  //================================================== Compute cardinal direction
  EdgeLoopCollection::iterator curloop;
  double maxD = 0.0;
  for (curloop = this->simplices.begin();
       curloop != this->simplices.end();
       curloop++)
  {
    EdgeList::iterator curedge;
    for (curedge = (*curloop)->edges.begin();
         curedge != (*curloop)->edges.end();
         curedge++)
    {
      Vector3 v01 = curedge->vertices[1] - curedge->vertices[0];
      double cheatx = fabs((v01/v01.Norm()).Dot(Vector3(1, 0, 0)) * 1.0e-3);
      //printf("%f   %f (%f)",v01.Norm(),maxD,cheatx);
      //printf("%+.4f %+.4f\n", v01.x,v01.y);
      if ((v01.Norm()+cheatx)>=maxD)
      {
        hat_i = v01;
        maxD = v01.Norm()+cheatx;
      }
    }
  }
  hat_i = hat_i/hat_i.Norm();

  //================================================== Compute hat j
  hat_j = hat_k.Cross(hat_i);
  hat_j = hat_j/hat_j.Norm();

  //================================================== Calculate projected
  //                                                   vertices
  double dmax = 0.0;
  for (curvert = Pstar.begin();
       curvert != Pstar.end();
       curvert++)
  {
    Vertex v = *curvert;
    Vertex p;

    p.x = v.x - centroid.x;
    p.y = v.y - centroid.y;
    p.z = v.z - centroid.z;

    double dist = p.Norm();
    if (dist>dmax)
    {
      dmax = dist;
    }



    v.x = p.Dot(hat_i);
    v.y = p.Dot(hat_j);
    v.z = 0.0;

    *curvert = v;

//    printf("Vertex mapped from %+.4f %+.4f %+.4f to ",(*curvert).x,
//                                                      (*curvert).y,
//                                                      (*curvert).z);
//    printf("%+.4f %+.4f %+.4f to\n",v.x,
//                                    v.y,
//                                    v.z);

  }

  //================================================== Set auto minimum
  if (get_auto_min)
  {
    absoluteMinumumSize = dmax/15.0;
    printf("Absolute Minimum size set to %f\n", absoluteMinumumSize);
  }

}