#include "chi_meshcontinuum.h"
#include "../Cell/cell_slab.h"
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
/**Populates a face histogram.
 *
 * \param master_tolerance Multiple histograms will only be attempted
 * if the ratio of the maximum dofs-per-face to the average dofs-per-face
 * is greater than this value. Default 1.2.
 *
 * \param slave_tolerance While traversing a sorted list of dofs-per-face,
 * a new bin will only be generated when the ratio of the listed dofs-per-face
 * to a running bin average exceeds this value. Defualt 1.1.
 *
 * The function populates face_categories which is a structure containing
 * pairs. Pair.first is the max dofs-per-face for the category and Pair.second
 * is the number of faces in this category.
 *
 * */
void chi_mesh::MeshContinuum::
  BuildFaceHistogramInfo(double master_tolerance, double slave_tolerance)
{
  if (face_histogram_available) return;

  //================================================== Fill histogram
  std::vector<size_t> face_size_histogram;
  for (auto c : local_cell_glob_indices)
  {
    auto cell = cells[c];

    if (cell->Type() == chi_mesh::CellType::SLAB)
    {
      CellSlab* slab_cell = static_cast<CellSlab*>(cell);

      size_t num_faces = 2;
      for (int f=0; f<num_faces; f++)
        face_size_histogram.push_back(1);

    }//if slab
    else if (cell->Type() == chi_mesh::CellType::POLYGON)
    {
      CellPolygon* poly_cell = static_cast<CellPolygon*>(cell);

      for (auto face : poly_cell->edges)
        face_size_histogram.push_back(2);
    }//if polygon
    else if (cell->Type() == chi_mesh::CellType::POLYHEDRON)
    {
      CellPolyhedron* polyh_cell = static_cast<CellPolyhedron*>(cell);

      for (auto face : polyh_cell->faces)
        face_size_histogram.push_back(face->v_indices.size());
    }//if polyhedron
  }
  std::stable_sort(face_size_histogram.begin(), face_size_histogram.end());

  //================================================== Determine total face dofs
  size_t total_face_dofs_count = 0;
  for (auto face_size : face_size_histogram)
    total_face_dofs_count += face_size;

  //================================================== Compute average and ratio
  size_t smallest_face = face_size_histogram.front();
  size_t largest_face = face_size_histogram.back();
  size_t total_num_faces = face_size_histogram.size();
  double average_dofs_per_face = (double)total_face_dofs_count/total_num_faces;

  std::stringstream outstr;
  outstr << "\nSmallest face = " << smallest_face;
  outstr << "\nLargest face = " << largest_face;
  outstr << "\nTotal face dofs = " << total_face_dofs_count;
  outstr << "\nTotal faces = " << face_size_histogram.size();
  outstr << "\nAverage dofs/face = " << average_dofs_per_face;
  outstr << "\nMax to avg ratio = " << largest_face/average_dofs_per_face;
  chi_log.Log(LOG_ALLVERBOSE_1) << outstr.str();

  //================================================== Determine number of bins
  size_t last_bin_num_faces = total_num_faces;
  if ((largest_face/average_dofs_per_face) > master_tolerance)
  {
    chi_log.Log(LOG_ALLVERBOSE_1)
    << "The ratio of max face dofs to average face dofs "
    << "is larger than " << master_tolerance
    << ", therefore a binned histogram "
    << "will be constructed.";

    //====================================== Build categories
    size_t running_total_face_dofs = 0;
    size_t running_face_count = 0;
    size_t running_face_size = face_size_histogram[0];

    double running_average = face_size_histogram[0];

    for (size_t f=0; f<total_num_faces; ++f)
    {
      if ((face_size_histogram[f]/running_average) > slave_tolerance)
      {
        face_categories.emplace_back(running_face_size,running_face_count);
        running_total_face_dofs = 0;
        running_face_count = 0;
      }

      running_face_size = face_size_histogram[f];
      running_total_face_dofs += face_size_histogram[f];
      running_face_count++;
      running_average = (double)running_total_face_dofs/running_face_count;
      last_bin_num_faces = running_face_count;
    }
  }
  face_categories.emplace_back(largest_face,last_bin_num_faces);

  //================================================== Verbose print bins
  outstr.str(std::string());
  outstr
  << "A total of " << face_categories.size()
  << " bins were created:\n";

  size_t bin_counter = -1;
  for (auto bins : face_categories)
  {
    outstr
    << "Bin " << ++bin_counter << ": "
    << bins.second << " faces with max face dofs " << bins.first << "\n";
  }

  chi_log.Log(LOG_ALLVERBOSE_1) << outstr.str();

  face_histogram_available = true;
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

//###################################################################
/**Gets the number of face-histogram categories.*/
size_t chi_mesh::MeshContinuum::NumberOfFaceHistogramBins()
{
  if (!face_histogram_available) BuildFaceHistogramInfo();

  return face_categories.size();
}

//###################################################################
/**Maps the face-histogram category number for a given face size.*/
size_t chi_mesh::MeshContinuum::MapFaceHistogramBins(size_t num_face_dofs)
{
  if (!face_histogram_available) BuildFaceHistogramInfo();

  size_t category_counter = -1;
  for (auto category : face_categories)
  {
    category_counter++;
    if (num_face_dofs <= category.first)
      return category_counter;
  }

  return 0;
}

//###################################################################
/**Maps the face-histogram category number for a given face size.*/
size_t chi_mesh::MeshContinuum::GetFaceHistogramBinDOFSize(size_t category)
{
  if (!face_histogram_available) BuildFaceHistogramInfo();

  size_t face_dof_size = 0;

  try {
    face_dof_size = face_categories.at(category).first;
  }
  catch (std::out_of_range o){
    chi_log.Log(LOG_ALLWARNING)
    << "Fault detected in chi_mesh::MeshContinuum::"
    << "GetFaceHistogramBinDOFSize.";
    return 0;
  }

  return face_dof_size;
}

//###################################################################
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