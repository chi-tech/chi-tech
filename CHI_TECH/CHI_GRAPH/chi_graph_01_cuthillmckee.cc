#include "chi_graph.h"

#include <boost/graph/bandwidth.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/graphviz.hpp>
#include <chi_log.h>

extern ChiLog chi_log;

//###################################################################
/**This method takes an undirected graph and performs
 * a Cuthill-Mckee reordering to reduce the bandwidth.*/
void chi_graph::
CuthillMckee(CHI_UD_GRAPH &in_graph,
             std::vector<int>* mapping)
{
  typedef boost::graph_traits<CHI_UD_GRAPH>::vertex_descriptor Vertex;
  typedef boost::graph_traits<CHI_UD_GRAPH>::vertices_size_type size_type;


  boost::graph_traits<CHI_UD_GRAPH>::vertex_iterator  ui, ui_end;

  boost::property_map<CHI_UD_GRAPH,boost::vertex_degree_t>::type deg =
    get(boost::vertex_degree,in_graph);

  //============================================= Compute vertex degrees
  for (boost::tie(ui, ui_end) = vertices(in_graph); ui != ui_end; ++ui)
    deg[*ui] = degree(*ui, in_graph);

  chi_log.Log(LOG_ALLVERBOSE_2) << "original bandwidth: " <<
            boost::bandwidth(in_graph) << std::endl;

  //=============================================
  boost::property_map<CHI_UD_GRAPH, boost::vertex_index_t>::type
    index_map = get(boost::vertex_index, in_graph);
  std::vector<Vertex>    inv_perm(num_vertices(in_graph));
  std::vector<size_type> perm(num_vertices(in_graph));

  //============================================= Compute Reverse
  //                                              Cuthill-Mckee_ordering
  Vertex s = vertex(0, in_graph);
  boost::cuthill_mckee_ordering(
    in_graph, s, inv_perm.rbegin(),
    get(boost::vertex_color, in_graph),
    get(boost::vertex_degree, in_graph));

  //============================================= Assign ordering
  //                                              to mapping
  for (std::vector<Vertex>::const_iterator i = inv_perm.begin();
       i != inv_perm.end(); ++i)
  {
    (*mapping)[*i] = index_map[*i];
  }


  for (size_type c = 0; c != inv_perm.size(); ++c)
  {
    perm[index_map[inv_perm[c]]] = c;
  }

  chi_log.Log(LOG_ALLVERBOSE_2)  << "  bandwidth: "
            << bandwidth(in_graph, make_iterator_property_map(&perm[0], index_map, perm[0]))
            << std::endl;

  //write_graphviz(std::cout, in_graph);
}