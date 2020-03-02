#include "chi_directed_graph.h"

#include <chi_log.h>

extern ChiLog chi_log;

#include <sstream>

//###################################################################
/** Adds a vertex to the graph.*/
void chi_graph::DirectedGraph::AddVertex(void* context)
{
  vertices.emplace_back(vertices.size(),context);
}

//###################################################################
/** Adds an edge to the graph.*/
bool chi_graph::DirectedGraph::AddEdge(int from, int to, bool allow_cycle)
{
  //=================================== Check "from" is in range
  if ((from<0) || (from>=vertices.size()))
  {
    chi_log.Log(LOG_ALL)
      << "chi_graph::DirectedGraph: Error added edge (from).";
    exit(EXIT_FAILURE);
  }

  //=================================== Check "to" is in range
  if ((to<0) || (to>=vertices.size()))
  {
    chi_log.Log(LOG_ALLERROR)
      << "chi_graph::DirectedGraph: Error added edge (to).";
    exit(EXIT_FAILURE);
  }

  vertices[from].ds_edge.insert(to);
  vertices[to].us_edge.insert(from);

  return true;
}

//###################################################################
/**Remove an edge from the graph.*/
void chi_graph::DirectedGraph::RemoveEdge(int from, int to)
{
  vertices[from].ds_edge.erase(to);
  vertices[to].us_edge.erase(from);
}

//###################################################################
/** Depth-First-Search main recursive algorithm.*/
void chi_graph::DirectedGraph::
  DFSAlgorithm(std::vector<int> &traversal,
               std::vector<bool> &visited,
               int cur_vid)
{
  traversal.push_back(cur_vid);
  visited[cur_vid] = true;

  for (auto v : vertices[cur_vid].ds_edge)
    if (not visited[v])
      DFSAlgorithm(traversal,visited,v);

}

//###################################################################
/** Depth-First-Search or traversal from specified vertex.*/
std::vector<int> chi_graph::DirectedGraph::DepthFirstSearch(int vertex_id)
{
  std::vector<int>  traversal;
  std::vector<bool> visited(vertices.size(),false);

  DFSAlgorithm(traversal,visited,vertex_id);

  return traversal;
}

//###################################################################
/**SCC main recursive algorithm.*/
void chi_graph::DirectedGraph::SCCAlgorithm(
  int u, int& time,
  std::vector<int>& disc,
  std::vector<int>& low,
  std::vector<bool>& on_stack,
  std::stack<int>& stack,
  std::vector<std::vector<int>>& SCCs)
{
  // Init discovery time and low value
  disc[u] = low[u] = ++time;
  stack.push(u);
  on_stack[u] = true;

  for (auto v : vertices[u].ds_edge)
  {
    if (disc[v] == -1)
    {
      SCCAlgorithm(v,time,disc,low,on_stack,stack,SCCs);
      low[u] = std::min(low[u],low[v]);
    }
    else if (on_stack[v] == true)
      low[u] = std::min(low[u],disc[v]);
  }

  int w=0;
  if (low[u] == disc[u])
  {
    std::vector<int> sub_SCC;
    while (stack.top() != u)
    {
      w = stack.top();
      sub_SCC.push_back(w);
      on_stack[w] = false;
      stack.pop();
    }
    w = stack.top();
    sub_SCC.push_back(w);
    if (sub_SCC.size() > 1) SCCs.push_back(sub_SCC);
    on_stack[w] = false;
    stack.pop();
  }

}

//###################################################################
/**Find strongly connected components.*/
std::vector<std::vector<int>> chi_graph::DirectedGraph::
  FindStronglyConnectedConnectionns()
{
  size_t V = vertices.size();

  std::vector<int>  disc(V,-1);        // Discovery times
  std::vector<int>  low(V,-1);         // Earliest visited vertex
  std::vector<bool> on_stack(V,false); // On stack flags
  std::stack<int>   stack;             // Stack

  std::vector<std::vector<int>> SCCs;  // Collection of SCCs

  int time = 0;

  for (int v=0; v<V; ++v)
    if (disc[v]==-1)
      SCCAlgorithm(v,time,disc,low,on_stack,stack,SCCs);


  return SCCs;
}


//###################################################################
/** Generates a topological sort. This method is the implementation
 * of Kahn's algorithm [1].
 *
 * [1] Kahn, Arthur B. (1962), "Topological sorting of large networks",
 *     Communications of the ACM, 5 (11): 558–562
 *
 * \return Returns the vertex ids sorted topologically. If this
 *         vector is empty the algorithm failed because it detected
 *         cyclic dependencies.*/
std::vector<int> chi_graph::DirectedGraph::GenerateTopologicalSort()
{
  bool has_cycles=false;
  std::vector<int> L;
  std::vector<GraphVertex*> S;

  L.reserve(vertices.size());
  S.reserve(vertices.size());

  //======================================== Make a copy of the graph
  auto cur_vertices = vertices;

  //======================================== Identify vertices that
  //                                         have no incoming edge
  for (auto& vertex : cur_vertices)
    if (vertex.us_edge.empty())
      S.push_back(&vertex);

  if (S.empty())
  {
    has_cycles = true;
    goto endofalgo;
  }

  //======================================== Repeatedly remove
  //                                         vertices
  while (not S.empty())
  {
    GraphVertex* node_n = S.back(); int n=node_n->id;
    S.erase(S.end()-1);

    L.push_back(n);
    auto nodes_m = node_n->ds_edge;
    for (int m : nodes_m)
    {
      GraphVertex* node_m = &cur_vertices[m];

      //remove edge from graph
      node_n->ds_edge.erase(m);
      node_m->us_edge.erase(n);

      if (node_m->us_edge.empty())
        S.push_back(node_m);
    }
  }

  endofalgo:
  {
    if (has_cycles)
      return std::vector<int>();
    if (L.size() != vertices.size())
      return std::vector<int>();
  }

  return L;
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
  verts_rank_r.clear();
  vertices.clear();
}

//###################################################################
/**Destructor.*/
chi_graph::DirectedGraph::~DirectedGraph()
{
  Clear();
}
