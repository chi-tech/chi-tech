#include "chi_volumemesher.h"
#include <iostream>
#include <vector>
#include "../MeshHandler/chi_meshhandler.h"
#include "../Region/chi_region.h"
#include "../Cell/cell_triangle.h"
#include "../Cell/cell_polygon.h"
#include "../Cell/cell_polyhedron.h"


#define CELL_TYPE_TRIANGLE   1
#define CELL_TYPE_POLYGON    2
#define CELL_TYPE_POLYHEDRON 3

//#########################################################
/***/
void chi_mesh::VolumeMesher::
  ReOrderDOF(chi_mesh::MeshContinuum* vol_continuum)
{
  this->ReOrderDOF_CuthillMckee(vol_continuum);
}

//#########################################################
/***/
void chi_mesh::VolumeMesher::
  ReOrderDOF_CuthillMckee(chi_mesh::MeshContinuum* vol_continuum)
{
  std::vector<chi_mesh::CellIndexMap*> unused_cells;
  std::vector<chi_mesh::NodeIndexMap*> unused_nodes;
  std::vector<chi_mesh::NodeIndexMap*> used_nodes;

  //============================================= Create a mapping for each cell
  std::vector<chi_mesh::Cell*>::iterator cell;
  for (cell = vol_continuum->cells.begin();
       cell != vol_continuum->cells.end();
       cell++)
  {
    chi_mesh::CellIndexMap* new_map = new chi_mesh::CellIndexMap;
    int from_index = std::distance(vol_continuum->cells.begin(),cell);
    new_map->mapped_from = from_index;

    this->cell_ordering.push_back(new_map);
    unused_cells.push_back(new_map);
  }

  //============================================= Create a mapping for each node
  std::vector<chi_mesh::Node*>::iterator node;
  for (node = vol_continuum->nodes.begin();
       node != vol_continuum->nodes.end();
       node++)
  {
    chi_mesh::NodeIndexMap* new_map = new chi_mesh::NodeIndexMap;
    int from_index = std::distance(vol_continuum->nodes.begin(),node);
    new_map->mapped_from = from_index;

    this->node_ordering.push_back(new_map);
    unused_nodes.push_back(new_map);
  }

  //============================================= Start with seed node
  int seed_index = 0;
  int map_level = 1;
  //======================================== Remove it from the unused_stack
  chi_mesh::NodeIndexMap* seed = unused_nodes.back();
  unused_nodes.pop_back();

  //======================================== Add it to the used stack
  seed->mapped_to    = seed_index;
  seed->mapped_level = map_level;
  used_nodes.push_back(seed);

  int num_unused_prev = unused_cells.size();
  while (unused_cells.size()>0)
  {
    //printf("Mapping Level %d Unused Cells =%d\n",map_level,unused_cells.size());
    seed_index = this->MapLevel(&used_nodes,
                                &unused_cells,
                                vol_continuum,
                                seed_index,map_level);
    int num_unused = unused_cells.size();

    if (num_unused==num_unused_prev)
    {
      break;
    } else
    {
      num_unused_prev = num_unused;
    }
    map_level++;
  }

//  std::vector<chi_mesh::NodeIndexMap *>::iterator used_node;
//  for (used_node = used_nodes.begin();
//       used_node != used_nodes.end();
//       used_node++)
//  {
//    chi_mesh::NodeIndexMap *actual_node = (*used_node);
//    int index = std::distance(used_nodes.begin(),used_node);
//    printf("Used node %3d, f=%3d, t=%3d\n",index,
//      actual_node->mapped_from,
//      actual_node->mapped_to);
//  }
//  printf("Remapped nodes= %d\n",used_nodes.size());

}

//#########################################################
/***/
bool chi_mesh::VolumeMesher::
  IsInList(std::vector<chi_mesh::NodeIndexMap*>* list, int index)
{
  std::vector<chi_mesh::NodeIndexMap*>::iterator cur_node;
  for (cur_node = list->begin();
       cur_node != list->end();
       cur_node++)
  {
    if ((*cur_node)->mapped_from == index)
    {
      return true;
    }
  }

  return false;
}


//#########################################################
/***/
int chi_mesh::VolumeMesher::
  MapLevel(std::vector<chi_mesh::NodeIndexMap *> *used_list,
                     std::vector<chi_mesh::NodeIndexMap*>* unused_cells,
                     chi_mesh::MeshContinuum* vol_continuum,
                     int seeded_index,
                     int level)
{
  int seed_index = seeded_index;
  std::vector<chi_mesh::NodeIndexMap*> nodes_to_be_added;

  //============================================= Run through all the used nodes
  std::vector<chi_mesh::NodeIndexMap *>::iterator used_node;
  for (used_node = used_list->begin();
       used_node != used_list->end();
       used_node++)
  {
    chi_mesh::NodeIndexMap* actual_node = (*used_node);
    //====================================== If map level equals required
    if (actual_node->mapped_level == level)
    {

      //============================================= Find all the connected item_id
      std::vector<chi_mesh::CellIndexMap*>::iterator cell_notyetused;
      for (cell_notyetused = unused_cells->begin();
           cell_notyetused != unused_cells->end();
           cell_notyetused++)
      {
        //====================================== Get the actual cell
        int old_index = (*cell_notyetused)->mapped_from;
        chi_mesh::CellTriangle* tri_cell;
        chi_mesh::CellPolygon* poly_cell;
        chi_mesh::CellPolyhedron* polyh_cell;
        int cell_type;

        if (dynamic_cast<chi_mesh::CellTriangle*>
            (vol_continuum->cells.at(old_index)))
        //if (typeid(*vol_continuum->cells.at(old_index)) ==
        //    typeid(chi_mesh::CellTriangle))
        {
          tri_cell = (chi_mesh::CellTriangle*)vol_continuum->cells.at(
                               old_index);
          cell_type = CELL_TYPE_TRIANGLE;
        }
        else if (dynamic_cast<chi_mesh::CellPolygon*>
                 (vol_continuum->cells.at(old_index)))
        {
          poly_cell = (chi_mesh::CellPolygon*)vol_continuum->cells.at(
                               old_index);
          cell_type = CELL_TYPE_POLYGON;
        }
        else if (dynamic_cast<chi_mesh::CellPolyhedron*>
                 (vol_continuum->cells.at(old_index)))
        {
          polyh_cell = (chi_mesh::CellPolyhedron*)vol_continuum->cells.at(
                               old_index);
          cell_type = CELL_TYPE_POLYHEDRON;
        }
        else
        {
          fprintf(stderr, "ERROR: Unrecognized cell type encountered "
                          "in DOF renumbering algorithm.\n");
          exit(EXIT_FAILURE);
        }




        //====================================== Check if it is connected to the seed
        bool is_connected = false;
        if (cell_type == CELL_TYPE_TRIANGLE)
        {
          for (int e=0;e<3;e++)
          {
            if (tri_cell->e_index[e][0]==actual_node->mapped_from)
            {
              is_connected = true;
              break;
            }
          }
        }
        else if (cell_type == CELL_TYPE_POLYGON)
        {
          for (int e=0;e<poly_cell->edges.size();e++)
          {
            if (poly_cell->edges[e][0]==actual_node->mapped_from)
            {
              is_connected = true;
              break;
            }
          }
        }
        else if (cell_type == CELL_TYPE_POLYHEDRON)
        {
          for (int f=0;f<polyh_cell->faces.size();f++)
          {
            for (int e=0; e<polyh_cell->faces[f]->edges.size(); e++)
            {
              if (polyh_cell->faces[f]->edges[e][0]==actual_node->mapped_from)
              {
                is_connected = true;
                break;
              }
            }
          }
        }



        //====================================== If it is connected
        //                                       add its nodes to the list and
        //                                       pull the cell from the unused list
        if (is_connected)
        {
          //======================== Push up all its nodes
          if (cell_type == CELL_TYPE_TRIANGLE)
          {
            for (int n=0;n<3;n++)
            {
              if (tri_cell->e_index[n][0] != actual_node->mapped_from)
              {
                //================= Check if its not already pushed
                if ((!this->IsInList(used_list,
                                     tri_cell->e_index[n][0])) &&
                    (!this->IsInList(&nodes_to_be_added,
                                     tri_cell->e_index[n][0])))
                {
                  chi_mesh::NodeIndexMap* map
                    = this->node_ordering.at(tri_cell->e_index[n][0]);

                  //=============== Push it to the used stack
                  seed_index++;
                  map->mapped_to = seed_index;
                  map->mapped_level = level+1;
                  nodes_to_be_added.push_back(map);
                }
              }
            }
          }
          else if (cell_type == CELL_TYPE_POLYGON)
          {
            for (int n=0;n<poly_cell->edges.size();n++)
            {
              if (poly_cell->edges[n][0] != actual_node->mapped_from)
              {
                //================= Check if its not already pushed
                if ((!this->IsInList(used_list,
                                     poly_cell->edges[n][0])) &&
                    (!this->IsInList(&nodes_to_be_added,
                                     poly_cell->edges[n][0])))
                {
                  chi_mesh::NodeIndexMap* map
                    = this->node_ordering.at(poly_cell->edges[n][0]);

                  //=============== Push it to the used stack
                  seed_index++;
                  map->mapped_to = seed_index;
                  map->mapped_level = level+1;
                  nodes_to_be_added.push_back(map);
                }
              }
            }
          }
          else if (cell_type == CELL_TYPE_POLYHEDRON)
          {
            for (int f=0; f<polyh_cell->faces.size(); f++)
            {
              for (int n=0;n<polyh_cell->faces[f]->edges.size();n++)
              {
                if (polyh_cell->faces[f]->edges[n][0] !=
                    actual_node->mapped_from)
                {
                  //================= Check if its not already pushed
                  if ((!this->IsInList(used_list,
                               polyh_cell->faces[f]->edges[n][0])) &&
                      (!this->IsInList(&nodes_to_be_added,
                               polyh_cell->faces[f]->edges[n][0])))
                  {
                    chi_mesh::NodeIndexMap* map
                      = this->node_ordering.at(
                                 polyh_cell->faces[f]->edges[n][0]);

                    //=============== Push it to the used stack
                    seed_index++;
                    map->mapped_to = seed_index;
                    map->mapped_level = level+1;
                    nodes_to_be_added.push_back(map);
                  }
                }
              }
            }

          }

          //======================== Remove the cell from unused list
          unused_cells->erase(cell_notyetused);
        }

        //====================================== Double check loop end
        if (cell_notyetused == unused_cells->end())
        {
          break;
        }
      }
    }
  }

  //============================================= Add all the new nodes
  for (used_node = nodes_to_be_added.begin();
       used_node != nodes_to_be_added.end();
       used_node++)
  {
    chi_mesh::NodeIndexMap *actual_node = (*used_node);
    used_list->push_back(actual_node);
  }


  return seed_index;
}