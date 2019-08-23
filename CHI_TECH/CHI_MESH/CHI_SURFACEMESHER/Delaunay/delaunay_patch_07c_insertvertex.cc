#include "delaunay_mesher.h"

//###################################################################
/**Inserts a vertex into the triangulation.*/
bool chi_mesh::SurfaceMesherDelaunay::DelaunayPatch::
InsertVertex(chi_mesh::Vertex vc_pstar, unsigned seed_search_triangle_index)
{
  //================================================== Struct for connection graph
  struct graph_entity
  {
    unsigned t_index;
    int      level;

    graph_entity(unsigned t, int l)
    {
      t_index = t;
      level   = l;
    }
  };
  std::vector<graph_entity> connection_graph;


  //================================================== Add initial triangle
  //int level = 0;
  //connection_graph.push_back(graph_entity(seed_search_triangle_index,level));





  //================================================== Build a connection graph
  // by looping until no triangles found
  // This prevents us from searching over all triangles
//  bool tri_found=false;
//  bool first_cycle=true;
//  while (tri_found || first_cycle)
//  {
//    first_cycle = false;
//    tri_found = false;
//    //printf("Iteration\n");
//    //=========================================== Loop over existing graph points
//    std::vector<graph_entity>::iterator g;
//    for (unsigned gi=0; gi<connection_graph.size();gi++)
//    {
//      graph_entity g = connection_graph[gi];
//      if (g.level == level)
//      {
//        Tri tau_1 = triangles[g.t_index];
//        //printf("Level %d, triangle %d\n",level, g.t_index);
//        for (int e=0;e<3;e++)
//        {
//          if (tau_1.e_index[e][2]>=0)
//          {
//            unsigned t2 = tau_1.e_index[e][2];
//            Tri tau_2 = triangles[t2];
//            Vertex a = Pstar[tau_2.v_index[0]];
//            Vertex b = Pstar[tau_2.v_index[1]];
//            Vertex c = Pstar[tau_2.v_index[2]];
//            Vertex d = vc_pstar;
//
//            //============================== Check if inserted vertex is within
//            //                               Circumcenter
//            if (InCircle(a,b,c,d)>0.0001)
//            {
//              //============================ Check for duplication
//              bool already_there = false;
//              for (unsigned gj=0; gj<connection_graph.size();gj++)
//              {
//                if (connection_graph[gj].t_index == t2)
//                {
//                  already_there = true;
//                  break;
//                }
//              }
//              //============================ If not already there then set it
//              if (!already_there)
//              {
//                //printf("Adding Triangle %d\n",t2);
//                connection_graph.push_back(graph_entity(t2,level+1));
//                //printf("Added Triangle %d\n",connection_graph.back().t_index);
//
//                tri_found = true;
//              }
//
//            }//if in circle
//
//          }//if edge>0
//        }//for edge
//      }//if glevel == level
//
//    }//for g
//    level++;
//  }
  //================================================== Build a connection graph
  for (unsigned t2=0; t2<triangles.size(); t2++)
  {
    Tri tau_2 = triangles[t2];
    if (!tau_2.invalidated)
    {
      Vertex a = Pstar[tau_2.v_index[0]];
      Vertex b = Pstar[tau_2.v_index[1]];
      Vertex c = Pstar[tau_2.v_index[2]];
      Vertex d = vc_pstar;

      //============================== Check if inserted vertex is within
      //                               Circumcenter
      //printf("Incircle t=%d   %f\n",t2,InCircle(a,b,c,d));
      if (InCircle(a,b,c,d)>0.0)
      {
        //printf("Removing Triangle %d\n",t2);
        connection_graph.push_back(graph_entity(t2,0));
        //printf("Added Triangle %d\n",connection_graph.back().t_index);

      }
    }

  }

  //================================================== Check if vertex already added
  for (unsigned k=0; k<Pstar.size(); k++)
  {
    Vector v01 = vc_pstar - Pstar[k];
    if (v01.Norm()<0.000001)
    {
      return false;
    }
  }

  //================================================== Write verbose output
  printf("Inserting projected vertex %d at [%.2f,%.2f]\n", Pstar.size(),vc_pstar.x,vc_pstar.y);

  //================================================== Inserted projected vertex
  Pstar.push_back(vc_pstar);
  int inserted_v_index = Pstar.size()-1;

  //================================================== Insert true vertex
  Vertex true_v;
  true_v = true_v + hat_i*vc_pstar.x;
  true_v = true_v + hat_j*vc_pstar.y;
  true_v = true_v + centroid;
  vertices.push_back(true_v);





  //================================================== Create list of outer edges
  std::vector<Edge> unsorted_outer_edge_list;
  //============================================= Loop over graph entities
  for (unsigned gi=0; gi<connection_graph.size();gi++)
  {
    graph_entity g = connection_graph[gi];
    Tri tau_1 = triangles[g.t_index]; //Get associated triangle

    //====================================== Loop over edges
    for (int e=0;e<3; e++)
    {
      Edge outer_edge;
      outer_edge.vertices[0] = Pstar[  tau_1.e_index[e][0]  ];
      outer_edge.vertices[1] = Pstar[  tau_1.e_index[e][1]  ];

      outer_edge.v_index[0] = tau_1.e_index[e][0];  //Vertex index i
      outer_edge.v_index[1] = tau_1.e_index[e][1];  //Vertex index f

      outer_edge.f_index[0] = g.t_index;            //Triangle to left
      outer_edge.f_index[1] = e;                    //Left tri edge index
      outer_edge.f_index[2] = tau_1.e_index[e][2];  //Triangle to right
      outer_edge.f_index[3] = tau_1.e_index[e][3];  //Right tri edge index

      //========================= find edge association
      int tl = g.t_index;
      int tr = tau_1.e_index[e][2];
      for (int e2=0; e2<3; e2++)
      {
        if (  (triangles[tl].e_index[e2][0] == outer_edge.v_index[0]) &&
              (triangles[tl].e_index[e2][1] == outer_edge.v_index[1])  )
        {
          outer_edge.f_index[1] = e2;
        }

        if (  (triangles[tr].e_index[e2][1] == outer_edge.v_index[0]) &&
              (triangles[tr].e_index[e2][0] == outer_edge.v_index[1])  )
        {
          outer_edge.f_index[3] = e2;
        }
      }

      unsorted_outer_edge_list.push_back(outer_edge);
    }

    //====================================== Invalidate the triangle
    triangles[g.t_index].invalidated = true;
  }

  //============================================= Remove internal edges
  bool internal_edge_found = true;
  while (internal_edge_found)
  {
    internal_edge_found = false;
    for (unsigned i=0; i<unsorted_outer_edge_list.size(); i++)
    {
      int tau1_index = unsorted_outer_edge_list[i].f_index[0];
      int tau2_index = unsorted_outer_edge_list[i].f_index[2];

      bool left_triangle_in_set=false;
      bool rite_triangle_in_set=false;
      //printf("Checking edge %d->%d \n",unsorted_outer_edge_list[i].v_index[0],unsorted_outer_edge_list[i].v_index[1]);
      for (unsigned gi=0; gi<connection_graph.size();gi++)
      {
        if (connection_graph[gi].t_index == tau1_index)
        {
          left_triangle_in_set = true;
        }
        if (connection_graph[gi].t_index == tau2_index)
        {
          rite_triangle_in_set = true;
        }
        if (left_triangle_in_set && rite_triangle_in_set) break;
      }

      if (left_triangle_in_set && rite_triangle_in_set)
      {
        //printf("Edge is internal\n");
        unsorted_outer_edge_list.erase(unsorted_outer_edge_list.begin()+i);
        internal_edge_found = true;
        break;
      }
    }

  }




//================================================== Verbose output edge list
//  for (unsigned e1=0; e1<unsorted_outer_edge_list.size(); e1++)
//  {
//    printf("Unsorted Edge %d, %d->%d [%d,%d] [%d,%d]\n",e1, unsorted_outer_edge_list[e1].v_index[0],
//           unsorted_outer_edge_list[e1].v_index[1],
//           unsorted_outer_edge_list[e1].f_index[0],
//           unsorted_outer_edge_list[e1].f_index[1],
//           unsorted_outer_edge_list[e1].f_index[2],
//           unsorted_outer_edge_list[e1].f_index[3]);
//  }



  //================================================== Sort edges
  std::vector<Edge> sorted_edge_list;

  //============================================= Add first non-colinear edge
  for (unsigned e1=0; e1<unsorted_outer_edge_list.size(); e1++)
  {
    Edge unsorted_edge = unsorted_outer_edge_list[e1];

    //======================================= Check if edge is colinear with
    //                                        inserted vertex
    bool is_colinear = false;
    double orient = Orient2D(unsorted_edge.vertices[0],
                             unsorted_edge.vertices[1],
                             vc_pstar);
    //printf("Orient1 edge %d->%d =%f\n",unsorted_edge.v_index[0],unsorted_edge.v_index[1],orient);
    if (fabs(orient)<tolerance) is_colinear = true;

    //======================================= Only insert if not colinear
    if (!is_colinear) sorted_edge_list.push_back(unsorted_edge);

    unsorted_outer_edge_list.erase(unsorted_outer_edge_list.begin()+e1);
    if (!is_colinear) break;
  }


  //sorted_edge_list.push_back(unsorted_outer_edge_list.back());
  //unsorted_outer_edge_list.pop_back();
  while (unsorted_outer_edge_list.size()>0)
  {
    for (unsigned e1=0; e1<unsorted_outer_edge_list.size(); e1++)
    {
      Edge unsorted_edge = unsorted_outer_edge_list[e1];

      if (unsorted_edge.v_index[1] == sorted_edge_list.front().v_index[0])
      {
        //======================================= Check if edge is colinear with
        //                                        inserted vertex
        bool is_colinear = false;
        double orient = Orient2D(unsorted_edge.vertices[0],
                                 unsorted_edge.vertices[1],
                                 vc_pstar);
        //printf("Orient2 edge %d->%d =%f\n",unsorted_edge.v_index[0],unsorted_edge.v_index[1],orient);
        if (fabs(orient)<tolerance) is_colinear = true;

        //======================================= Only insert if not colinear
        if (!is_colinear) sorted_edge_list.insert(sorted_edge_list.begin(),unsorted_edge);

        unsorted_outer_edge_list.erase(unsorted_outer_edge_list.begin()+e1);
        break;
      }

      if (unsorted_edge.v_index[0] == sorted_edge_list.back().v_index[1])
      {
        //======================================= Check if edge is colinear with
        //                                        inserted vertex
        bool is_colinear = false;
        double orient = Orient2D(unsorted_edge.vertices[0],
                                 unsorted_edge.vertices[1],
                                 vc_pstar);
        //printf("Orient3 edge %d->%d =%f\n",unsorted_edge.v_index[0],unsorted_edge.v_index[1],orient);
        if (fabs(orient)<tolerance) is_colinear = true;

        //======================================= Only insert if not colinear
        if (!is_colinear) sorted_edge_list.push_back(unsorted_edge);

        unsorted_outer_edge_list.erase(unsorted_outer_edge_list.begin()+e1);
        break;
      }
    }
  }

  //================================================== Verbose output edge list
//  for (unsigned e1=0; e1<sorted_edge_list.size(); e1++)
//  {
//    printf("Edge %d, %d->%d [%d,%d] [%d,%d]\n",e1, sorted_edge_list[e1].v_index[0],
//                                                   sorted_edge_list[e1].v_index[1],
//                                                   sorted_edge_list[e1].f_index[0],
//                                                   sorted_edge_list[e1].f_index[1],
//                                                   sorted_edge_list[e1].f_index[2],
//                                                   sorted_edge_list[e1].f_index[3]);
//  }





  //================================================== Insert new triangles
  //                                                   edge-wise, starting with
  //                                                   first
  unsigned last_tri=triangles.size();
  unsigned first_tri_index = last_tri;
  Tri first_triangle;
  first_triangle.v_index[0] = inserted_v_index;
  first_triangle.v_index[1] = sorted_edge_list[0].v_index[0];
  first_triangle.v_index[2] = sorted_edge_list[0].v_index[1];

  //Edge0 v-indices
  first_triangle.e_index[0][0] = first_triangle.v_index[0];
  first_triangle.e_index[0][1] = first_triangle.v_index[1];
  //Edge1 v-indices
  first_triangle.e_index[1][0] = first_triangle.v_index[1];
  first_triangle.e_index[1][1] = first_triangle.v_index[2];
  //Edge2 v-indices
  first_triangle.e_index[2][0] = first_triangle.v_index[2];
  first_triangle.e_index[2][1] = first_triangle.v_index[0];

  //Edge0 unknown connection

  //Edge1 edge connection from edge
  first_triangle.e_index[1][2] = sorted_edge_list[0].f_index[2];
  first_triangle.e_index[1][3] = sorted_edge_list[0].f_index[3];

  //Edge2 unknown connection

//  printf("Triangle %d added %d->%d->%d\n",triangles.size(),
//                                       first_triangle.v_index[0],
//                                       first_triangle.v_index[1],
//                                       first_triangle.v_index[2]);
  triangles.push_back(first_triangle);

  //========= Update neighbors
  if (sorted_edge_list[0].f_index[2]>=0)
  {
    int neighbor_index      = sorted_edge_list[0].f_index[2];
    int neighbor_edge_index = sorted_edge_list[0].f_index[3];

    triangles[ neighbor_index ].e_index[neighbor_edge_index][2] = last_tri;
    triangles[ neighbor_index ].e_index[neighbor_edge_index][3] = 1;

  }

  //================================================== Insert other triangles in
  //                                                   a loop
  for (unsigned e=1; e<sorted_edge_list.size(); e++)
  {
    Tri new_triangle;
    new_triangle.v_index[0] = inserted_v_index;
    new_triangle.v_index[1] = sorted_edge_list[e].v_index[0];
    new_triangle.v_index[2] = sorted_edge_list[e].v_index[1];

    //Edge0 v-indices
    new_triangle.e_index[0][0] = new_triangle.v_index[0];
    new_triangle.e_index[0][1] = new_triangle.v_index[1];
    //Edge1 v-indices
    new_triangle.e_index[1][0] = new_triangle.v_index[1];
    new_triangle.e_index[1][1] = new_triangle.v_index[2];
    //Edge2 v-indices
    new_triangle.e_index[2][0] = new_triangle.v_index[2];
    new_triangle.e_index[2][1] = new_triangle.v_index[0];

    //Edge0 from previous triangle
    new_triangle.e_index[0][2] = last_tri;
    new_triangle.e_index[0][3] = 2;

    //Edge1 edge connection from edge
    new_triangle.e_index[1][2] = sorted_edge_list[e].f_index[2];
    new_triangle.e_index[1][3] = sorted_edge_list[e].f_index[3];

    //Edge2 unknown connection

//    printf("Triangle added %d->%d->%d\n",new_triangle.v_index[0],
//           new_triangle.v_index[1],
//           new_triangle.v_index[2]);
    triangles.push_back(new_triangle);

    //================================= Update last triangle
    triangles[last_tri].e_index[2][2] = triangles.size()-1;
    triangles[last_tri].e_index[2][3] = 0;

    //========= Update neighbors
    if (sorted_edge_list[e].f_index[2]>=0)
    {
      int neighbor_index      = sorted_edge_list[e].f_index[2];
      int neighbor_edge_index = sorted_edge_list[e].f_index[3];

      triangles[ neighbor_index ].e_index[neighbor_edge_index][2] = triangles.size()-1;
      triangles[ neighbor_index ].e_index[neighbor_edge_index][3] = 1;

    }

    last_tri=triangles.size()-1;
  }


  //================================================== Making first and last
  //                                                   triangle talk when not
  //                                                   edge splitting
  if (sorted_edge_list.front().v_index[0] == sorted_edge_list.back().v_index[1])
  {
    triangles[last_tri].e_index[2][2] = first_tri_index;
    triangles[last_tri].e_index[2][3] = 0;

    triangles[first_tri_index].e_index[0][2] = last_tri;
    triangles[first_tri_index].e_index[0][3] = 2;
  }






  //================================================== Clean up
  connection_graph.clear();

  return true;
}