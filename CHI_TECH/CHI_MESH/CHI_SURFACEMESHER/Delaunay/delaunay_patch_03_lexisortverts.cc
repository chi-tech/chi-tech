#include "delaunay_mesher.h"

//#########################################################
/**Creates a lexicographical mapping of vertices in Pstar*/
void chi_mesh::SurfaceMesherDelaunay::DelaunayPatch::
          LexicographicallySortVertices()
{
  std::vector<Vertex>::iterator i;
  std::vector<Vertex>::iterator j;

  //================================================== Bubble sorting x
  for (i = Pstar.begin(); i!=Pstar.end(); i++)
  {
    for (j = Pstar.begin(); j!=(Pstar.end()-1);j++)
    {
      if (j->x > (j+1)->x)
      {
        std::iter_swap(j,j+1);
      }
    }
  }

  //================================================== Bubble sorting y
  for (i = Pstar.begin(); i!=Pstar.end(); i++)
  {
    for (j = Pstar.begin(); j!=(Pstar.end()-1);j++)
    {
      if ( (fabs(j->x - (j+1)->x)<tolerance ) && (j->y > (j+1)->y) )
      {
        std::iter_swap(j,j+1);
      }
    }
  }

  //================================================== InitializeAlphaElements unused vertices
  for (i = Pstar.begin(); i!=Pstar.end(); i++)
  {
    int iindex = std::distance(Pstar.begin(),i);
    unused_vs.insert(unused_vs.begin(),iindex);

  }

  std::vector<int>::iterator unused;
  for (unused = unused_vs.begin(); unused != unused_vs.end(); unused++)
  {
    int index = *unused;
    Vertex v = Pstar.at(index);
    printf("xyv(%d,1)=%+.4f;",index+1,v.x);
    printf("xyv(%d,2)=%+.4f;\n",index+1,v.y);


  }


}