#ifndef _chi_directed_graph_h
#define _chi_directed_graph_h

#include "chi_directed_graph_vertex.h"
#include <stack>

//###################################################################
/**Simple implementation of a directed graph. This implementation was
 * considered to serve more versatile strategies with regards to grid
 * parallel partitioning.*/
class chi_graph::DirectedGraph
{
public:

  //============================================= Vertex accessor definition
  /**Allows semi-sane access to vertices even if
   * they are removed from the graph.*/
  class VertexAccessor
  {
  private:
    std::vector<GraphVertex> vertices;
    std::vector<bool>        vertex_valid_flags;
  public:
    void AddVertex(void* context, int id=-1);
    void RemoveVertex(int v);

    GraphVertex& operator[](int v);

    //############################ iterator Class Definition
    /**Internal iterator class for vertex accessor.*/
    class iterator
    {
    public:
      VertexAccessor& ref_block;
      size_t      ref_element;

      iterator(VertexAccessor& in_block, size_t i) :
        ref_block(in_block),
        ref_element(i) {}

      iterator operator++()
      {
        iterator i = *this;
        ++ref_element;
        while (not ref_block.vertex_valid_flags[ref_element] and
               ref_element<ref_block.vertices.size())
          ++ref_element;
        return i;
      }
      iterator operator++(int junk)
      {
        ++ref_element;
        while (not ref_block.vertex_valid_flags[ref_element] and
               ref_element<ref_block.vertices.size())
          ++ref_element;
        return *this;
      }
      GraphVertex& operator*()
      { return ref_block.vertices[ref_element]; }
      GraphVertex* operator->()
      { return &(ref_block.vertices[ref_element]); }
      bool operator==(const iterator& rhs)
      { return ref_element == rhs.ref_element; }
      bool operator!=(const iterator& rhs)
      { return ref_element != rhs.ref_element; }
    };
    //############################ End of iterator Class Definition

    iterator begin()
    {
      size_t count=0;
      while (not vertex_valid_flags[count] and
             count<vertices.size())
        ++count;
      return {*this,count};
    }

    iterator end(){return {*this,vertices.size()};}

    size_t size() {return vertices.size();}

    int GetNumValid()
    {
      int count=0;
      for (bool val : vertex_valid_flags)
        if (val) ++count;

      return count;
    }

    void clear() {vertices.clear(); vertex_valid_flags.clear();}
  };
  //============================================= End of Vertex accessor def

  VertexAccessor vertices;

  void AddVertex(void* context = nullptr, int id=-1);
  void RemoveVertex(int v);
  bool AddEdge(int from, int to, double weight=1.0);
  void RemoveEdge(int from, int to);

  int GetNumSinks()
  {
    int count=0;
    for (auto& v : vertices)
      if (v.ds_edge.empty() and not v.us_edge.empty())
        ++count;
    return count;
  }

  int GetNumSources()
  {
    int count=0;
    for (auto& v : vertices)
      if (v.us_edge.empty() and not v.ds_edge.empty())
        ++count;
    return count;
  }

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
    FindStronglyConnectedComponents();

  std::vector<int> GenerateTopologicalSort();

  std::vector<int> FindApproxMinimumFAS();

  void PrintGraphviz(int location_mask=0);

  void PrintSubGraphviz(const std::vector<int>& verts_to_print,
                        int location_mask=0);

  void Clear();

  ~DirectedGraph();
};

#endif