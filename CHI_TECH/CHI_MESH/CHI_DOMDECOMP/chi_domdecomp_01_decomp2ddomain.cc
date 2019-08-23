#include "chi_domdecomp.h"
#include "../CHI_MESHHANDLER/chi_meshhandler.h"
#include "../CHI_MESHCONTINUUM/chi_meshcontinuum.h"
#include "../CHI_REGION/chi_region.h"
#include "../CHI_CELL/cell.h"
#include "../CHI_CELL/cell_triangle.h"

#include <typeinfo>

//###################################################################
/**Decomposes a 2D domain*/
void CHI_DOMDECOMP::Decompose2DDomain(int Px, int Py,
  chi_mesh::SurfaceMesherTriangle* mesher, chi_mesh::Region* region)
{
  //================================================== Get the current handler
  chi_mesh::MeshHandler* mesh_handler = chi_mesh::GetCurrentHandler();

  //================================================== Get continuum
  chi_mesh::MeshContinuum* mesh_continuum = region->volume_mesh_continua.back();


  //================================================== Find minimums & maximums
  double max_x = -1.0e10;
  double max_y = -1.0e10;

  double min_x = 1.0e10;
  double min_y = 1.0e10;
  for (unsigned v=0; v<mesh_continuum->nodes.size(); v++)
  {
    if (mesh_continuum->nodes[v]->x > max_x)
    {max_x = mesh_continuum->nodes[v]->x;}

    if (mesh_continuum->nodes[v]->y > max_y)
    {max_y = mesh_continuum->nodes[v]->y;}

    if (mesh_continuum->nodes[v]->x < min_x)
    {min_x = mesh_continuum->nodes[v]->x;}

    if (mesh_continuum->nodes[v]->y < min_y)
    {min_y = mesh_continuum->nodes[v]->y;}
  }

  //================================================== Build auxiliary cutlist
  std::vector<double> aux_xcuts;
  std::vector<double> aux_ycuts;

  for (unsigned i=0; i<mesher->xcuts.size(); i++)
  {
    double cut = mesher->xcuts[i];
    aux_xcuts.push_back(cut);
  }
  for (unsigned i=0; i<mesher->ycuts.size(); i++)
  {
    double cut = mesher->ycuts[i];
    aux_ycuts.push_back(cut);
  }

  aux_xcuts.push_back(max_x);
  aux_ycuts.push_back(max_y);



  //================================================== Build domain bounds
  //std::vector<chi_mesh::CELL_SET*> cell_sets;
  int dpx = aux_xcuts.size()/Px;
  int dpy = aux_ycuts.size()/Py;

  for (int j=0;j<Py;j++)
  {
    for (int i=0; i<Px;i++)
    {
      chi_mesh::CELL_SET* new_cell_set = new chi_mesh::CELL_SET;
      if (i==0) {new_cell_set->xmin = min_x;}
      if (j==0) {new_cell_set->ymin = min_y;}

      if (i!=0) {new_cell_set->xmin = aux_xcuts[dpx + (i-1)*dpx-1];}
      if (j!=0) {new_cell_set->ymin = aux_ycuts[dpy + (j-1)*dpy-1];}

      new_cell_set->xmax = aux_xcuts[dpx + i*dpx-1];
      new_cell_set->ymax = aux_ycuts[dpy + j*dpy-1];

      new_cell_set->i = i;
      new_cell_set->j = j;

      new_cell_set->mesh_continuum = mesh_continuum;

      mesh_handler->cell_sets.push_back(new_cell_set);
    }
  }

  //================================================== Allocate item_id
  for (unsigned c=0; c<mesh_continuum->cells.size(); c++)
  {
    chi_mesh::Cell* cell = mesh_continuum->cells[c];
    for (unsigned d=0; d<mesh_handler->cell_sets.size(); d++)
    {
      if ( (cell->centroid.x >= mesh_handler->cell_sets[d]->xmin) &&
           (cell->centroid.x <= mesh_handler->cell_sets[d]->xmax) &&
           (cell->centroid.y >= mesh_handler->cell_sets[d]->ymin) &&
           (cell->centroid.y <= mesh_handler->cell_sets[d]->ymax) )
      {
        mesh_handler->cell_sets[d]->cells_allocated.push_back(c);
        break;
      }
    }//for d
  }//for c

  //================================================== Verbose output
  for (unsigned d=0; d<mesh_handler->cell_sets.size(); d++)
  {
    printf("Domain %d, [ %d,%d ], %f<=x<=%f, %f<=y<=%f, cellcount=%d\n", d,
           mesh_handler->cell_sets[d]->i,
           mesh_handler->cell_sets[d]->j,
           mesh_handler->cell_sets[d]->xmin,
           mesh_handler->cell_sets[d]->xmax,
           mesh_handler->cell_sets[d]->ymin,
           mesh_handler->cell_sets[d]->ymax,
           mesh_handler->cell_sets[d]->cells_allocated.size());
  }


}