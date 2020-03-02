#ifndef _chi_directed_graph_h
#define _chi_directed_graph_h

#include "chi_graph.h"

struct chi_graph::GraphVertex
{
  int id;
  void* context;

  std::set<int> us_edge;
  std::set<int> ds_edge;

  GraphVertex(int in_id, void* in_context) :
    id(in_id),
    context(in_context)
  {}

  GraphVertex(int in_id) :
    id(in_id),
    context(nullptr)
  {}

  GraphVertex(const GraphVertex& in_v)
  {
    this->id = in_v.id;
    this->context = in_v.context;

    us_edge = in_v.us_edge;
    ds_edge = in_v.ds_edge;
  }

  GraphVertex(GraphVertex&& in_v)
  {
    this->id = in_v.id;
    this->context = in_v.context;

    us_edge = in_v.us_edge;
    ds_edge = in_v.ds_edge;

    in_v.context = nullptr;
  }

  GraphVertex& operator=(const GraphVertex& in_v)
  {
    this->id = in_v.id;
    this->context = in_v.context;

    us_edge = in_v.us_edge;
    ds_edge = in_v.ds_edge;

    return *this;
  }

  GraphVertex& operator=(GraphVertex&& in_v)
  {
    this->id = in_v.id;
    this->context = in_v.context;

    us_edge = in_v.us_edge;
    ds_edge = in_v.ds_edge;

    in_v.context = nullptr;

    return *this;
  }

  bool operator==(const GraphVertex& other)
  {
    return other.id == this->id;
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
  std::vector<GraphVertex> vertices;

  void AddVertex(void* context = nullptr);
  bool AddEdge(int from, int to, bool allow_cycle=false);
  void RemoveEdge(int from, int to);

private:
  void DFSAlgorithm(std::vector<int>& traversal,
                    std::vector<bool>& visited,
                    int cur_vid);

  void SCCAlgorithm(int u, int& time,
                    std::vector<int>& disc,
                    std::vector<int>& low,
                    std::vector<bool>& on_stack,
                    std::stack<int>& stack,
                    std::vector<std::vector<int>>& SCCs);

public:
  std::vector<int> DepthFirstSearch(int vertex_id);
  std::vector<std::vector<int>>
    FindStronglyConnectedConnectionns();

  std::vector<int> GenerateTopologicalSort();
  std::vector<std::vector<GraphVertex*>> GetGraphRanks();

  void Clear();

  ~DirectedGraph();
};

#endif