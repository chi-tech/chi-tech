#include "chi_directed_graph.h"

#include "chi_log.h"
#include "chi_mpi.h"

extern ChiLog&  chi_log;
extern ChiMPI& chi_mpi;

#include <sstream>
#include <algorithm>

//###################################################################
/** Adds a vertex to the graph.*/
void chi_graph::DirectedGraph::
  VertexAccessor::AddVertex(void* context, int id)
{
  if (id<0)
    vertices.emplace_back(vertices.size(), context);
  else
    vertices.emplace_back(id, context);

  vertex_valid_flags.push_back(true);
}

//###################################################################
/** Removes a vertex from the graph.*/
void chi_graph::DirectedGraph::
  VertexAccessor::RemoveVertex(int v)
{
  //=================================== Check "from" is in range
  if ((v<0) || (v>=vertices.size()))
  {
    chi_log.Log(LOG_ALL)
      << "chi_graph::DirectedGraph::VertexAccessor: "
      << "Error removing vertex " << v;
    exit(EXIT_FAILURE);
  }

  auto& vertex = vertices[v];

  //=================================== Get adjacent vertices
  auto num_us = vertex.us_edge.size();
  auto num_ds = vertex.ds_edge.size();
  std::vector<int> adj_verts;
  adj_verts.reserve(num_us+num_ds);

  for (int u : vertex.us_edge)
    adj_verts.push_back(u);

  for (int u : vertex.ds_edge)
    adj_verts.push_back(u);

  //=================================== Remove v from all u
  for (int u : adj_verts)
  {
    vertices[u].us_edge.erase(v);
    vertices[u].ds_edge.erase(v);
  }

  vertex_valid_flags[v] = false;
}

//###################################################################
/** Accesses a vertex from the graph.*/
chi_graph::GraphVertex& chi_graph::DirectedGraph::
  VertexAccessor::operator[](int v)
{
  if (not vertex_valid_flags[v])
    chi_log.Log(LOG_ALLERROR)
      << "chi_graph::DirectedGraph::VertexAccessor: "
         "Invalid vertex accessed. Vertex may have been removed.";
  return vertices[v];
}

//###################################################################
/** Adds a vertex to the graph. By default <I>context</I> is
 * assumed to be nullptr and <I>id</I> is assumed to be -1. In
 * the latter case the vertex id will be the same as the order
 * in which it was added (0,1,2,3,etc ... will have id's 0,1,2,3,etc)*/
void chi_graph::DirectedGraph::AddVertex(void* context, int id)
{
  vertices.AddVertex(context,id);
}

//###################################################################
/** Removes a vertex from the graph. This method does not
 * free any context related data.*/
void chi_graph::DirectedGraph::
  RemoveVertex(int v)
{
  vertices.RemoveVertex(v);
}

//###################################################################
/** Adds an edge to the graph. Range checks are supplied by the
 * vertex accessor.*/
bool chi_graph::DirectedGraph::AddEdge(int from, int to, double weight)
{
  vertices[from].ds_edge.insert(to);
  vertices[to].us_edge.insert(from);

  vertices[from].ds_weights[to] = weight;
  vertices[to].us_weights[from] = weight;

  return true;
}

//###################################################################
/**Remove an edge from the graph. Range checks are supplied by the
 * vertex accessor.*/
void chi_graph::DirectedGraph::RemoveEdge(int from, int to)
{
  vertices[from].ds_edge.erase(to);
  vertices[to].us_edge.erase(from);
}

//###################################################################
/** Depth-First-Search main recursive algorithm. This is the recursive
 * portion of the method below this one
 * (chi_graph::DirectedGraph::DepthFirstSearch).*/
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
/** Depth-First-Search or traversal from specified vertex. Returns
 * the order in which vertices will be traversed in a depth first sense.
 * This algorithm will return the sequence of DFS traversal.*/
std::vector<int> chi_graph::DirectedGraph::DepthFirstSearch(int vertex_id)
{
  std::vector<int>  traversal;
  std::vector<bool> visited(vertices.size(),false);

  DFSAlgorithm(traversal,visited,vertex_id);

  return traversal;
}

//###################################################################
/**SCC main recursive algorithm. This is the recursive call for the
 * method defined below this one
 * (chi_graph::DirectedGraph::FindStronglyConnectedConnections).*/
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
/**Find strongly connected components. This method is the implementation
 * of Tarjan's algorithm [1].
 *
 * [1] Tarjan R.E. "Depth-first search and linear graph algorithms",
 *     SIAM Journal on Computing, 1972.
 *
 * It returns collections of vertices that form strongly connected
 * components excluding singletons.*/
std::vector<std::vector<int>> chi_graph::DirectedGraph::
  FindStronglyConnectedComponents()
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
/**Finds a sequence that minimizes the Feedback Arc Set (FAS). This
 * algorithm implements the algorithm depicted in [1].
 *
 * [1] Eades P., Lin X., Smyth W.F., "Fast & Effective heuristic for
 *     the feedback arc set problem", Information Processing Letters,
 *     Volume 47. 1993.*/
std::vector<int> chi_graph::DirectedGraph::
  FindApproxMinimumFAS()
{
  auto GetVertexDelta = [](chi_graph::GraphVertex& vertex)
  {
//    double delta = vertex.ds_edge.size() - vertex.us_edge.size();
//    return delta;

    double delta = 0.0;
    for (int ds : vertex.ds_edge)
      delta += 1.0*vertex.ds_weights[ds];

    for (int us : vertex.us_edge)
      delta -= 1.0*vertex.us_weights[us];

    return delta;
  };

  auto& TG = *this;

  //==================================== Execute GR-algorithm
  std::vector<int> s1,s2,s;
  while (TG.vertices.GetNumValid()>0)
  {
    //======================== Remove sinks
    while (TG.GetNumSinks()>0)
    {
      for (auto& u : TG.vertices)
        if (u.ds_edge.empty())
        {
          TG.RemoveVertex(u.id);
          s2.push_back(u.id);
          break;
        }
    }//G contains sinks

    //======================== Remove sources
    while (TG.GetNumSources()>0)
    {
      for (auto& u : TG.vertices)
        if (u.us_edge.empty())
        {
          TG.RemoveVertex(u.id);
          s1.push_back(u.id);
          break;
        }
    }//G contains sinks

    //======================== Get max delta
    std::pair<int,double> max_delta(-1,-100.0);
    for (auto& u : TG.vertices)
    {
      double delta = GetVertexDelta(u);
      if (delta > max_delta.second)
        max_delta = std::make_pair(u.id,delta);
    }

    //======================== Remove max delta
    TG.RemoveVertex(max_delta.first);
    s1.push_back(max_delta.first);
  }


  //========================== Make appr. minimum FAS sequence
  s.reserve(s1.size() + s2.size());
  for (int u : s1) s.push_back(u);
  for (int u : s2) s.push_back(u);

  return s;
}


//###################################################################
/**Prints the graph in Graphviz format.*/
void chi_graph::DirectedGraph::PrintGraphviz(int location_mask)
{
  if (chi_mpi.location_id != location_mask) return;

  std::stringstream o;
  std::string offset("    ");
  o << "Printing directed graph:\n";
  o << "digraph DG {\n";

  o << offset << "splines=\"FALSE\";\n";
  o << offset << "rankdir=\"LR\";\n\n";


  o << offset << "/* Vertices */\n";
  for (auto& v : vertices)
    o << offset << v.id << " [shape=\"circle\"]\n";

  o << "\n" << offset << "/* Edges */\n";
  for (auto& v : vertices)
    for (int w : v.ds_edge)
      o << offset << v.id << " -> " << w << "\n";

  o << "}\n";

  std::cout << o.str();
}

//###################################################################
/**Prints a sub-graph in Graphviz format.*/
void chi_graph::DirectedGraph::
  PrintSubGraphviz(const std::vector<int>& verts_to_print,
                   int location_mask)
{
  if (chi_mpi.location_id != location_mask) return;

  std::stringstream o;
  std::string offset("    ");
  o << "Printing directed graph:\n";
  o << "digraph DG {\n";

  o << offset << "splines=\"FALSE\";\n";
  o << offset << "rankdir=\"LR\";\n\n";


  o << offset << "/* Vertices */\n";
  for (int v : verts_to_print)
    o << offset << v << " [shape=\"circle\"]\n";

  o << "\n" << offset << "/* Edges */\n";
  for (int v : verts_to_print)
    for (int w : vertices[v].ds_edge)
    {
      if (std::find(verts_to_print.begin(),
                    verts_to_print.end(),
                    w) != verts_to_print.end())
        o << offset << v << " -> " << w << "\n";
    }

  o << "}\n";

  std::cout << o.str();
}

//###################################################################
/**Clears all the data structures associated with the graph.*/
void chi_graph::DirectedGraph::Clear()
{
  vertices.clear();
}


//###################################################################
/**Destructor.*/
chi_graph::DirectedGraph::~DirectedGraph()
{
  Clear();
}
