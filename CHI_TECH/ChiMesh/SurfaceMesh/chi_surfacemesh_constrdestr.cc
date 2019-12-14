#include"chi_surfacemesh.h"

//#########################################################
/**Default constructor.*/
chi_mesh::SurfaceMesh::SurfaceMesh()
{

}

//#########################################################
/**Default destructor.*/
chi_mesh::SurfaceMesh::~SurfaceMesh()
{
  for (auto poly_face : poly_faces)
  {
    delete poly_face;
  }

  poly_faces.clear();
  std::cout << "Surface mesh destructor\n";
}

//#########################################################
std::ostream& operator<<(std::ostream& os,
   chi_mesh::SurfaceMesh& that)
{
  std::vector<chi_mesh::Face>::iterator curface;
  for (curface = that.faces.begin();
       curface != that.faces.end();
       curface++)
  {
    int index = std::distance(that.faces.begin(),curface);
    os << "Face " << index << " v:";
    os << curface->v_index[0] << "->";
    os << curface->v_index[1] << "->";
    os << curface->v_index[2] << "  ";

    os << "e:";
    for (int e=0;e<3;e++)
    {
      os << "[" << curface->e_index[e][0] << ",";
      os        << curface->e_index[e][1] << ",";
      os        << curface->e_index[e][2] << ",";
      os        << curface->e_index[e][3] << "]";
    }
    os << std::endl;

  }

  return os;
}