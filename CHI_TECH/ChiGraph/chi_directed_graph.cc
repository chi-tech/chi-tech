#include "chi_directed_graph.h"

#include <chi_log.h>

//###################################################################
/** Adds a vertex to the graph.*/
void chi_graph::DirectedGraph::AddVertex(void* context)
{
  vertices.push_back(new GraphVertex(context));
  vertices.back()->id = vertices.size()-1;
}

//###################################################################
/** Adds an edge to the graph.*/
void chi_graph::DirectedGraph::AddEdge(int from, int to)
{
  if ((from<0) || (from>=vertices.size()))
  {
    chi_log.Log(LOG_ALL)
      << "chi_graph::DirectedGraph: Error added edge (from).";
    exit(EXIT_FAILURE);
  }

  if ((to<0) || (to>=vertices.size()))
  {
    chi_log.Log(LOG_ALLERROR)
      << "chi_graph::DirectedGraph: Error added edge (to).";
    exit(EXIT_FAILURE);
  }

  vertices[from]->ds_edge.insert(to);
  vertices[to]->us_edge.insert(from);


}

//###################################################################
/** Generates a topological sort.*/
std::vector<int> chi_graph::DirectedGraph::GenerateTopologicalSort()
{
  std::set<int> unsorted_verts;
  //============================================= Find vertices with no
  //                                              upstream connections
  verts_rank_r.push_back(std::vector<GraphVertex*>());

  int num_verts = vertices.size();
  for (int v=0; v<num_verts; v++)
  {
    if (vertices[v]->us_edge.size() == 0)
    {
      verts_rank_r[0].push_back(vertices[v]);
      vertices[v]->rank = 0;
    }
    else
      unsorted_verts.insert(v);
  }//for v

  //============================================= Building ranked levels
  int cur_rank = 0;
  while (unsorted_verts.size()>0)
  {
    verts_rank_r.push_back(std::vector<GraphVertex*>());
    int new_rank = cur_rank+1;

    bool contrib_made = false;
    int num_r_verts = verts_rank_r[cur_rank].size();
    for (int vr=0; vr<num_r_verts; vr++)
    {
      GraphVertex* cur_vert = verts_rank_r[cur_rank][vr];
      std::set<int>::iterator ds_id;
      for (ds_id  = cur_vert->ds_edge.begin();
           ds_id != cur_vert->ds_edge.end();
           ds_id++)
      {
        int v_ds = *ds_id;
        if (vertices[v_ds]->rank<0)
        {
          verts_rank_r[new_rank].push_back(vertices[v_ds]);
          vertices[v_ds]->rank = new_rank;
          unsorted_verts.erase(v_ds);
          contrib_made = true;
        }
        else
        {
          verts_rank_r[new_rank].push_back(vertices[v_ds]);
          vertices[v_ds]->rank = new_rank;
          unsorted_verts.erase(v_ds);
          contrib_made = true;
        }

      }//for ds_id
    }//for vr

    if ((contrib_made == false) && (unsorted_verts.size()>0))
    {
      chi_log.Log(LOG_ALLERROR)
        << "DirectedGraph::GenerateTopologicalSort: Sorting stalled.";
      exit(EXIT_FAILURE);
    }
    cur_rank++;
  }//while

  //============================================= De-levelize into sort
  for (int r=0; r<verts_rank_r.size(); r++)
  {
    int num_r_verts = verts_rank_r[r].size();
    for (int vr=0; vr<num_r_verts; vr++)
      topological_sort.push_back(verts_rank_r[r][vr]->id);
  }


  return topological_sort;
}

//###################################################################
/**Gets the underlying ranked sorting used for generating
 * the topological sort.*/
std::vector<std::vector<chi_graph::GraphVertex*>>
  chi_graph::DirectedGraph::GetGraphRanks()
{
  return verts_rank_r;
}

//###################################################################
/**Clears all the data structures associated with the graph.*/
void chi_graph::DirectedGraph::Clear()
{

  int num_verts = vertices.size();
  for (int v=0; v<num_verts; v++)
  {
    delete vertices[v];
  }
  verts_rank_r.clear();
  vertices.clear();
}