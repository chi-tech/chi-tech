#include "chi_meshcontinuum.h"
#include "../Cell/cell_polygon.h"
#include "../Cell/cell_polyhedron.h"

#include <boost/graph/bandwidth.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/graphviz.hpp>

#include <chi_mpi.h>
#include <chi_log.h>

extern ChiMPI chi_mpi;
extern ChiLog chi_log;

//###################################################################
/**Connects the nodes of the mesh to each other.*/
void chi_mesh::MeshContinuum::ConnectGrid()
{
  //================================================== Loop over all item_id
//  for (int c=0; c<this->cells.size(); c++)
//  {
//    auto cell = (chi_mesh::CellPolyhedron*)this->grid_graph[c].cell_ptr;
//
//    for (int f=0; f<cell->faces.size(); f++)
//    {
//      if (cell->faces[f]->face_indices[0]>=0)
//      {
//        int cd = cell->faces[f]->face_indices[0];
//        boost::add_edge(c,cd,this->grid_graph);
//      }
//    }
//  }
/*
  //================================================== Iterate
  typedef boost::graph_traits<UndirectedGraph>::vertex_descriptor Vertex;
  typedef boost::graph_traits<UndirectedGraph>::vertices_size_type size_type;


  boost::graph_traits<UndirectedGraph>::vertex_iterator  ui, ui_end;

  boost::property_map<UndirectedGraph,boost::vertex_degree_t>::type deg =
    get(boost::vertex_degree,grid_graph);

  for (boost::tie(ui, ui_end) = vertices(grid_graph); ui != ui_end; ++ui)
    deg[*ui] = degree(*ui, grid_graph);


  boost::property_map<UndirectedGraph, boost::vertex_index_t>::type
    index_map = get(boost::vertex_index, grid_graph);

  std::cout << "original bandwidth: " <<
  boost::bandwidth(grid_graph) << std::endl;



  std::vector<Vertex> inv_perm(num_vertices(grid_graph));
  std::vector<size_type> perm(num_vertices(grid_graph));

  //reverse cuthill_mckee_ordering
  Vertex s = vertex(6, grid_graph);
  boost::cuthill_mckee_ordering(grid_graph, s, inv_perm.rbegin(),
    get(boost::vertex_color, grid_graph),
                         get(boost::vertex_degree, grid_graph));
  std::cout << "Reverse Cuthill-McKee ordering starting at: " << s <<
  std::endl;
  std::cout << "  ";

  for (std::vector<Vertex>::const_iterator i = inv_perm.begin();
       i != inv_perm.end(); ++i)
    std::cout << index_map[*i] << " ";
  std::cout << std::endl;

  for (size_type c = 0; c != inv_perm.size(); ++c)
    perm[index_map[inv_perm[c]]] = c;
  std::cout << "  bandwidth: "
            << bandwidth(grid_graph, make_iterator_property_map(&perm[0], index_map, perm[0]))
            << std::endl;

  write_graphviz(std::cout, grid_graph);
  */
}


/**Check whether a cell is local*/
bool chi_mesh::MeshContinuum::IsCellLocal(int cell_global_index)
{
  if (cell_global_index<0)
  {
    return false;
  }
  else
  {
    auto cell = cells[cell_global_index];
    if (cell->partition_id == chi_mpi.location_id)
    {
      return true;
    } else
    {
      return false;
    }

  }
  return false;
}

/**Check whether a cell is a boundary*/
bool chi_mesh::MeshContinuum::IsCellBndry(int cell_global_index)
{
  if (cell_global_index<0)
    return true;

  return false;
}

//###################################################################
/**Polyhedron find associated cell face*/
int chi_mesh::MeshContinuum::FindAssociatedFace(chi_mesh::PolyFace *cur_face,
                                                int adj_cell_g_index,bool verbose)
{
  //======================================== Check index validity
  if (IsCellBndry(adj_cell_g_index) || (!IsCellLocal(adj_cell_g_index)))
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid cell index encountered in call to "
      << "MeshContinuum::FindAssociatedFace. Index points to either a boundary"
      << "or a non-local cell.";
    exit(EXIT_FAILURE);
  }

  //======================================== Check cell validity by index
  if (adj_cell_g_index >= cells.size())
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid cell index encountered in call to "
      << "MeshContinuum::FindAssociatedFace. Index is out of cell index bounds.";
    exit(EXIT_FAILURE);
  }

  chi_mesh::Cell* cell = cells[adj_cell_g_index];

  if (cell->Type() != chi_mesh::CellType::POLYHEDRON)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid cell type encountered in "
      << "MeshContinuum::FindAssociatedFace. Adjacent cell is expected to"
         "be polyhedron but is found not to be. Index given "
      << adj_cell_g_index << " " << cur_face->face_centroid.PrintS();
    exit(EXIT_FAILURE);
  }

  chi_mesh::CellPolyhedron* adj_polyhcell = (chi_mesh::CellPolyhedron*)cell;

  int associated_face = -1;

  //======================================== Loop over adj cell faces
  for (int af=0; af<adj_polyhcell->faces.size(); af++)
  {
    //Assume face matches
    bool face_matches = true; //Now disprove it
    //================================= Loop over adj cell face verts
    for (int afv=0; afv<adj_polyhcell->faces[af]->v_indices.size(); afv++)
    {
      //========================== Try and find them in the reference face
      bool found = false;
      for (int cfv=0; cfv<cur_face->v_indices.size(); cfv++)
      {
        if (cur_face->v_indices[cfv] == adj_polyhcell->faces[af]->v_indices[afv])
        {
          found = true;
          break;
        }
      }//for cfv

      if (!found) {face_matches = false; break;}
    }//for afv

    if (face_matches) {associated_face = af; break;}
  }


  //======================================== Check associated face validity
  if (associated_face<0)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Could not find associated face in call to "
      << "MeshContinuum::FindAssociatedFace. Reference face with centroid at \n"
      << cur_face->face_centroid.PrintS();
    for (int af=0; af<adj_polyhcell->faces.size(); af++)
    {
      chi_log.Log(LOG_ALLERROR)
        << "Adjacent cell face " << af << " centroid "
        << adj_polyhcell->faces[af]->face_centroid.PrintS();
    }
    exit(EXIT_FAILURE);
  }

  //======================================== Verbose output
  if (verbose)
  {
    std::stringstream out_string;

    out_string
    << "Adj cell " << adj_cell_g_index << ":\n"
    << "face " << associated_face << "\n";
    for (int v=0; v<adj_polyhcell->faces[associated_face]->v_indices.size(); v++)
    {
      out_string
      << "vertex " << v << " "
      << adj_polyhcell->faces[associated_face]->v_indices[v]
      << "\n";
    }
    out_string
      << "Cur cell "
      << adj_polyhcell->faces[associated_face]->face_indices[NEIGHBOR] << ":\n";
    for (int v=0; v<cur_face->v_indices.size(); v++)
    {
      out_string
      << "vertex " << v << " "
      << cur_face->v_indices[v]
      << "\n";
    }
    chi_log.Log(LOG_ALL) << out_string.str();
  }

  return associated_face;
}

//###################################################################
/**Polyhedron find associated cell face*/
int chi_mesh::MeshContinuum::FindAssociatedEdge(int* edgeinfo,
                                                int adj_cell_g_index,bool verbose)
{
  //======================================== Check index validity
  if (IsCellBndry(adj_cell_g_index) || (!IsCellLocal(adj_cell_g_index)))
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid cell index encountered in call to "
      << "MeshContinuum::FindAssociatedEdge. Index points to either a boundary"
      << "or a non-local cell.";
    exit(EXIT_FAILURE);
  }

  //======================================== Check cell validity by index
  if (adj_cell_g_index >= cells.size())
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid cell index encountered in call to "
      << "MeshContinuum::FindAssociatedEdge. Index is out of cell index bounds.";
    exit(EXIT_FAILURE);
  }

  chi_mesh::Cell* cell = cells[adj_cell_g_index];

  if (cell->Type() != chi_mesh::CellType::POLYGON)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid cell type encountered in "
      << "MeshContinuum::FindAssociatedEdge. Adjacent cell is expected to"
         "be polygon but is found not to be. Index given "
      << adj_cell_g_index << " " << edgeinfo[EDGE_NEIGHBOR];
    exit(EXIT_FAILURE);
  }

  chi_mesh::CellPolygon* adj_polycell = (chi_mesh::CellPolygon*)cell;

  int associated_face = -1;

  //======================================== Loop over adj cell faces
  for (int af=0; af<adj_polycell->edges.size(); af++)
  {
    //Assume face matches
    bool face_matches = false; //Now disprove it
    if ((adj_polycell->edges[af][0] == edgeinfo[1]) and
        (adj_polycell->edges[af][1] == edgeinfo[0]))
    {
      face_matches = true;
    }

    if (face_matches) {associated_face = af; break;}
  }


  //======================================== Check associated face validity
  if (associated_face<0)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Could not find associated face in call to "
      << "MeshContinuum::FindAssociatedEdge.";
    exit(EXIT_FAILURE);
  }

  //======================================== Verbose output
//  if (verbose)
//  {
//    std::stringstream out_string;
//
//    out_string
//      << "Adj cell " << adj_cell_g_index << ":\n"
//      << "face " << associated_face << "\n";
//    for (int v=0; v<adj_polyhcell->faces[associated_face]->v_indices.size(); v++)
//    {
//      out_string
//        << "vertex " << v << " "
//        << adj_polyhcell->faces[associated_face]->v_indices[v]
//        << "\n";
//    }
//    out_string
//      << "Cur cell "
//      << adj_polyhcell->faces[associated_face]->face_indices[NEIGHBOR] << ":\n";
//    for (int v=0; v<cur_face->v_indices.size(); v++)
//    {
//      out_string
//        << "vertex " << v << " "
//        << cur_face->v_indices[v]
//        << "\n";
//    }
//    chi_log.Log(LOG_ALL) << out_string.str();
//  }

  return associated_face;
}

//###################################################################
/**Polyhedron map vertices*/
void chi_mesh::MeshContinuum::
       FindAssociatedVertices(chi_mesh::PolyFace* cur_face,
                              int adj_cell_g_index, int associated_face,
                              std::vector<int>& dof_mapping)
{
  //======================================== Check index validity
  if (IsCellBndry(adj_cell_g_index) || (!IsCellLocal(adj_cell_g_index)))
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid cell index encountered in call to "
      << "MeshContinuum::FindAssociatedVertices. Index points to either a boundary"
      << "or a non-local cell.";
    exit(EXIT_FAILURE);
  }

  //======================================== Check cell validity by index
  if (adj_cell_g_index >= cells.size())
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid cell index encountered in call to "
      << "MeshContinuum::FindAssociatedVertices. Index is out of cell index bounds.";
    exit(EXIT_FAILURE);
  }

  chi_mesh::Cell* cell = cells[adj_cell_g_index];

  if (cell->Type() != chi_mesh::CellType::POLYHEDRON)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid cell type encountered in "
      << "MeshContinuum::FindAssociatedVertices. Adjacent cell is expected to"
         "be polyhedron but is found not to be. Index given "
      << adj_cell_g_index << " " << cur_face->face_centroid.PrintS();
    exit(EXIT_FAILURE);
  }

  chi_mesh::CellPolyhedron* adj_polyhcell = (chi_mesh::CellPolyhedron*)cell;


  for (int cfv=0; cfv<cur_face->v_indices.size(); cfv++)
  {
    bool found = false;
    for (int afv=0;
         afv<adj_polyhcell->faces[associated_face]->v_indices.size(); afv++)
    {
      if (cur_face->v_indices[cfv] ==
          adj_polyhcell->faces[associated_face]->v_indices[afv])
      {
        dof_mapping.push_back(afv);
        found = true;
        break;
      }
    }

    if (!found)
    {
      chi_log.Log(LOG_ALLERROR)
        << "Face DOF mapping failed in call to "
        << "MeshContinuum::FindAssociatedVertices. Could not find a matching"
           "node."
        << adj_cell_g_index << " " << cur_face->face_centroid.PrintS();
      exit(EXIT_FAILURE);
    }

  }//for cfv

}