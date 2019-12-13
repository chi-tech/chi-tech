#include "chi_surfacemesh.h"

#include <chi_log.h>
extern ChiLog chi_log;

#include <ChiTimer/chi_timer.h>
extern ChiTimer    chi_program_timer;

#include <algorithm>

//#########################################################
/** Runs over the faces of the surfacemesh and determines
 * neighbors. */
void chi_mesh::SurfaceMesh::UpdateInternalConnectivityOld()
{
  chi_log.Log(LOG_0)
    << chi_program_timer.GetTimeString()
    << " Updating surface mesh internal connectivity.";
  std::vector<chi_mesh::Face>::iterator curFace;
  for (curFace = this->faces.begin(); curFace!=this->faces.end(); curFace++)
  {
    int outer_index = std::distance(this->faces.begin(),curFace);

    //Reset face connections
    for (int e=0;e<3;e++) //Loop over current face's edges
    {
      curFace->e_index[e][2] = -1; //tri index
      curFace->e_index[e][3] = -1;
    }
    //printf("Outer index %d\n",outer_index);
    //=========================================== Determine face connectivity
    //If the head v_index of the current triangle is equal to the tail
    //of the other triangle and vice versa on the tail of the current
    //triangle, then they are connected.
    std::vector<chi_mesh::Face>::iterator other_face;
    for (other_face = faces.begin(); other_face!=faces.end(); other_face++)
    {
      int inner_index = std::distance(this->faces.begin(),other_face);
      //if (outer_index!=inner_index)
      {
        for (int e=0;e<3;e++) //Loop over current face's edges
        {

          for (int e2=0;e2<3;e2++)  //Loop over other face's edges
          {
            if ( (curFace->e_index[e][0]==other_face->e_index[e2][1]) &&
                 (curFace->e_index[e][1]==other_face->e_index[e2][0]) )
            {
              curFace->e_index[e][2] = inner_index; //tri index
              curFace->e_index[e][3] = e2;                       //edge index
            }
          }
        }
      }

    }
  }

  for (int f=0; f<this->poly_faces.size(); f++)
  {
    chi_mesh::PolyFace* curFace = this->poly_faces[f];

    for (int f2=0; f2<this->poly_faces.size(); f2++)
    {
      if (f!=f2)
      {
        chi_mesh::PolyFace* other_face = this->poly_faces[f2];
        for (int e=0;e<curFace->edges.size();e++) //Loop over current face's edges
        {

          for (int e2=0;e2<other_face->edges.size();e2++)  //Loop over other face's edges
          {
            if ( (curFace->edges[e][0]==other_face->edges[e2][1]) &&
                 (curFace->edges[e][1]==other_face->edges[e2][0]) )
            {
              curFace->edges[e][2] = f2; //tri index
              curFace->edges[e][3] = e2;                       //edge index
            }
          }
        }
      }
    }
  }

  chi_log.Log(LOG_0)
    << chi_program_timer.GetTimeString()
    << " Completed surface mesh internal connectivity.";
}

//#########################################################
/** Runs over the faces of the surfacemesh and determines
 * neighbors. */
void chi_mesh::SurfaceMesh::UpdateInternalConnectivity()
{
  //======================================== Get number of cells
  size_t number_of_cells = poly_faces.size() + faces.size();

  if ((number_of_cells <= 1000) or (faces.size()>0))
  { UpdateInternalConnectivityOld(); return;}

  chi_log.Log(LOG_0)
    << chi_program_timer.GetTimeString()
    << " Updating surface mesh internal connectivity.";

  //======================================== Determine bound box
  double x_min = 0.0;
  double x_max = 0.0;
  double y_min = 0.0;
  double y_max = 0.0;
  bool first_vertex = true;
  for (auto& vertex : vertices)
  {
    if (first_vertex)
    {
      first_vertex=false;
      x_min = vertex.x; x_max = vertex.x;
      y_min = vertex.y; y_max = vertex.y;
    }

    if (vertex.x > x_max) x_max = vertex.x;
    if (vertex.y > y_max) y_max = vertex.y;
    if (vertex.x < x_min) x_min = vertex.x;
    if (vertex.y < y_min) y_min = vertex.y;
  }//for vertex

  //======================================== Make abbreviated cell info list
  //Utility typedef
  typedef std::tuple<int,double,double> CellInfo; //cell_index,centroid x, y
  typedef std::vector<CellInfo> CellList;         //vector of CellInfo
  typedef std::vector<std::pair<double,double>> DimLimits; //min,max

  CellList cells;
  cells.reserve(number_of_cells);

  int cell_id = 0;
  for (auto& poly_face : poly_faces)
    cells.emplace_back(cell_id++,
                       poly_face->face_centroid.x,
                       poly_face->face_centroid.y);
  for (auto& face : faces)
    cells.emplace_back(cell_id++,
                       face.face_centroid.x,
                       face.face_centroid.y);

  //======================================== Build Dimension limits
  int N = std::ceil(sqrt(number_of_cells));
  double dx = (x_max-x_min)/N;
  double dy = (y_max-y_min)/N;

  DimLimits x_dims(N);
  DimLimits y_dims(N);

  for (int i=0; i<N; i++)
  {
    x_dims[i].first = x_min + i*dx;
    x_dims[i].second= x_min + i*dx + dx;

    y_dims[i].first = y_min + i*dy;
    y_dims[i].second= y_min + i*dy + dy;
  }

  //======================================== Lambda for finding cell
  auto FindCellIJ = [](CellInfo& cell_info,
                       DimLimits& x_dims,
                       DimLimits& y_dims)
  {
    double cell_x = std::get<1>(cell_info);
    double cell_y = std::get<2>(cell_info);

    //Find I
    int I=-1,i=-1;
    for (auto& dim : x_dims){ ++i;
      if ( (cell_x >= dim.first) and (cell_x < dim.second) )
      { I = i; break;}
    }

    //Find J
    int J=-1,j=-1;
    for (auto& dim : y_dims){ ++j;
      if ( (cell_y >= dim.first) and (cell_y < dim.second) )
      { J = j; break;}
    }

    return std::pair<int,int>(I,J);
  };

  //======================================== Test to find all cells
  std::vector<std::vector<CellList>> cell_xy_list(N);
  for (auto& cell_y_list : cell_xy_list) cell_y_list.resize(N);

  bool all_cells_found = true;
  int num_unfound = 0;
  for (auto& cell : cells)
  {
    auto ij = FindCellIJ(cell,x_dims,y_dims);

    if ( (ij.first < 0) or (ij.second < 0) )
    {
      all_cells_found = false;
      ++num_unfound;
    }
    else
    {
      cell_xy_list[ij.first][ij.second].push_back(cell);
    }
  }

  if (!all_cells_found)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Not all cells found in call to "
      << "chi_mesh::SurfaceMesh::UpdateInternalConnectivity. "
      << "Number of cells not found: " << num_unfound;
    exit(EXIT_FAILURE);
  }

  std::vector<std::pair<int,int>> search_stencil =
    {{0,0},{1,0},{-1,0},{0,1},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}};

  int ci = 0;
  for (auto curFace : poly_faces)
  {
    auto ij = FindCellIJ(cells[ci],x_dims,y_dims);

    int ir = ij.first;
    int jr = ij.second;
    for (auto& stencil : search_stencil)
    {
      int i = ir + stencil.first;
      int j = jr + stencil.second;

      if ((i<0 or i>=N) or
          (j<0 or j>=N))
      { continue; }

      for (auto& cell : cell_xy_list[i][j])
      {
        auto other_face = poly_faces[std::get<0>(cell)];

        //Loop over current face's edges
        for (size_t e=0;e<curFace->edges.size();e++)
        {
          //Loop over other face's edges
          for (size_t e2=0;e2<other_face->edges.size();e2++)
          {
            if ( (curFace->edges[e][0]==other_face->edges[e2][1]) &&
                 (curFace->edges[e][1]==other_face->edges[e2][0]) )
            {
              curFace->edges[e][2] = std::get<0>(cell); //tri index
              curFace->edges[e][3] = e2;                //edge index
            }
          }//for e2
        }//for e
      }//for cell in stencil
    }

    ++ci;
  }

  chi_log.Log(LOG_0)
    << chi_program_timer.GetTimeString()
    << " Completed surface mesh internal connectivity.";
}