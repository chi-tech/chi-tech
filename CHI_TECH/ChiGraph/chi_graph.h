#ifndef _chi_graph_h
#define _chi_graph_h

#include "../CHI_MESH/chi_mesh.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/index/indexable.hpp>


/**Cell pointer struct.*/
struct GRAPH_CELL_INFO
{
  chi_mesh::Cell* cell_ptr;

  GRAPH_CELL_INFO()
  {
    cell_ptr = nullptr;
  }

  GRAPH_CELL_INFO& operator=(const GRAPH_CELL_INFO& other)
  {
    cell_ptr = other.cell_ptr;
    return *this;
  }
};

typedef boost::property<boost::vertex_degree_t,
  int,
  GRAPH_CELL_INFO> vertex_deg_prop;
typedef boost::property<boost::vertex_color_t,
                        boost::default_color_type,
                        vertex_deg_prop> vertexProperties;
typedef boost::adjacency_list<
        boost::vecS, //Use std::vector for EdgeList
        boost::vecS, //Use std::vector for VertexList
        boost::undirectedS, //Graph type selector
        vertexProperties> CHI_UD_GRAPH;

typedef boost::adjacency_list<
        boost::vecS, //Use std::vector for EdgeList
        boost::vecS, //Use std::vector for VertexList
        boost::directedS, //Graph type selector
        vertexProperties> CHI_D_GRAPH;

namespace chi_graph
{
  void CuthillMckee(CHI_UD_GRAPH& in_graph,
                    std::vector<int>* mapping);
  struct GraphVertex;
  class DirectedGraph;
}


#endif