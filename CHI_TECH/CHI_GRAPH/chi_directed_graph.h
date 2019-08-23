#ifndef _chi_directed_graph_h
#define _chi_directed_graph_h

#include "chi_graph.h"

struct chi_graph::GraphVertex
{
  int id;
  void* context;
  int rank;

  std::set<int> us_edge;
  std::set<int> ds_edge;

  GraphVertex()
  {
    id = -1;
    rank = -1;
    context = nullptr;
  }

  GraphVertex(void* in_context)
  {
    id = -1;
    rank = -1;
    context = in_context;
  }

  GraphVertex(const GraphVertex& in_v)
  {
    this->id = in_v.id;
    this->context = in_v.context;
    this->rank = in_v.rank;
  }
};

//###################################################################
/**Simple implementation of a directed graph. This implementation was
 * considered to serve more versatile strategies with regards to grid
 * parallel partitioning.*/
class chi_graph::DirectedGraph
{
private:
  std::vector<std::vector<GraphVertex*>> verts_rank_r;
  std::vector<int> topological_sort;

public:
  std::vector<GraphVertex*> vertices;

  void AddVertex(void* context = nullptr);
  void AddEdge(int from, int to);

  std::vector<int> GenerateTopologicalSort();
  std::vector<std::vector<GraphVertex*>> GetGraphRanks();

  void Clear();
};

#endif