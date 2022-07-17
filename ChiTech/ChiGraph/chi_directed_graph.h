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
    void AddVertex(size_t id, void* context);
    void AddVertex(void* context);
    void RemoveVertex(size_t v);

    GraphVertex& operator[](size_t v);

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
      bool operator==(const iterator& rhs) const
      { return ref_element == rhs.ref_element; }
      bool operator!=(const iterator& rhs) const
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

    size_t GetNumValid()
    {
      size_t count=0;
      for (bool val : vertex_valid_flags)
        if (val) ++count;

      return count;
    }

    void clear() {vertices.clear(); vertex_valid_flags.clear();}
  };
  //============================================= End of Vertex accessor def

  VertexAccessor vertices;

  void AddVertex(size_t id, void* context = nullptr);
  void AddVertex(void* context = nullptr);
  void RemoveVertex(size_t v);
  bool AddEdge(size_t from, size_t to, double weight=1.0);
  void RemoveEdge(size_t from, size_t to);

  size_t GetNumSinks()
  {
    size_t count=0;
    for (auto& v : vertices)
      if (v.ds_edge.empty() and not v.us_edge.empty())
        ++count;
    return count;
  }

  size_t GetNumSources()
  {
    size_t count=0;
    for (auto& v : vertices)
      if (v.us_edge.empty() and not v.ds_edge.empty())
        ++count;
    return count;
  }

private:
  void DFSAlgorithm(std::vector<size_t>& traversal,
                    std::vector<bool>& visited,
                    size_t cur_vid);

  void SCCAlgorithm(size_t u, int& time,
                    std::vector<int>& disc,
                    std::vector<int>& low,
                    std::vector<bool>& on_stack,
                    std::stack<size_t>& stack,
                    std::vector<std::vector<size_t>>& SCCs);

public:
  std::vector<std::vector<size_t>>
    FindStronglyConnectedComponents();

  std::vector<size_t> GenerateTopologicalSort();

  std::vector<size_t> FindApproxMinimumFAS();

  void PrintGraphviz(int location_mask=0);

  void PrintSubGraphviz(const std::vector<int>& verts_to_print,
                        int location_mask=0);

  std::vector<std::pair<size_t,size_t>> RemoveCyclicDependencies();

  void Clear();

  ~DirectedGraph();
};

#endif