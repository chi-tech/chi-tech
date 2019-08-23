#include "delaunay_mesher.h"

//#########################################################
/**Generates the unique simplices for a patch from a parent*/
void chi_mesh::SurfaceMesherDelaunay::DelaunayPatch::GenerateSimplices(bool add_vertices)
{
  //================================================== Get edge loops
  EdgeLoopCollection* edge_loops = raw_patch->GetEdgeLoops();

  //================================================== Create essential vertex
  //                                                   list and edges
  //Theoretically it would be here where a conformal interface
  //would be defined but this will have to be done later
  std::vector<int> vert_addition_history;
  EdgeLoopCollection::iterator cur_edge_loop;
  for (cur_edge_loop = edge_loops->begin();
       cur_edge_loop != edge_loops->end();
       cur_edge_loop++)
  {
    //=========================================== Split loop by angle
    EdgeLoopCollection* split_loop = SplitEdgeLoopByAngle(*cur_edge_loop);

    //=========================================== Loop over splitloops
    EdgeLoopCollection::iterator split_edge_loop;
    for (split_edge_loop = split_loop->begin();
         split_edge_loop != split_loop->end();
         split_edge_loop++)
    {
      Edge first_edge = (*split_edge_loop)->edges.front();
      Edge last_edge  = (*split_edge_loop)->edges.back();

      int first_v_index = first_edge.v_index[0];
      int last_v_index = last_edge.v_index[1];

      int num1_index = -1;
      int num2_index = -1;

      if (add_vertices)
      {
        //==================================== Check if vertices are already added
        bool num1_already_there = false;
        bool num2_already_there = false;

        std::vector<int>::iterator curvert_index;
        for (curvert_index = vert_addition_history.begin();
             curvert_index != vert_addition_history.end();
             curvert_index++)
        {
          if (*curvert_index == first_v_index)
          {
            num1_already_there = true;
            num1_index = std::distance(vert_addition_history.begin(),curvert_index);
            break;
          }
        }
        for (curvert_index = vert_addition_history.begin();
             curvert_index != vert_addition_history.end();
             curvert_index++)
        {
          if (*curvert_index == last_v_index)
          {
            num2_already_there = true;
            num2_index = std::distance(vert_addition_history.begin(),curvert_index);
            break;
          }
        }

        //==================================== Add the vertices
        if (!num1_already_there)
        {
          this->vertices.push_back(first_edge.vertices[0]);
          vert_addition_history.push_back(first_v_index);
          num1_index = this->vertices.size()-1;
        }
        if (!num2_already_there)
        {
          this->vertices.push_back(last_edge.vertices[1]);
          vert_addition_history.push_back(last_v_index);
          num2_index = this->vertices.size()-1;
        }
      }
      else
      {
        num1_index = first_v_index;
        num2_index = last_v_index;
      }



      EdgeLoop* new_loop = new EdgeLoop;
      Edge* new_edge = new Edge;

      new_edge->vertices[0] = first_edge.vertices[0];
      new_edge->vertices[1] = last_edge.vertices[1];

      new_edge->v_index[0] = num1_index;
      new_edge->v_index[1] = num2_index;

      //printf("edge %d->%d\n",num1_index,num2_index);

      new_loop->edges.push_back(*new_edge);

      this->simplices.push_back(new_loop);
    }

  }
}