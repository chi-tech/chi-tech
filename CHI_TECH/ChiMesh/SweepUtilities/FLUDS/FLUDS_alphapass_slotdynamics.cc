#include "FLUDS.h"
#include "ChiMesh/SweepUtilities/SPDS/SPDS.h"

#include <ChiMesh/Cell/cell.h>

typedef std::vector<std::pair<int,short>> LockBox;

#include <chi_log.h>

extern ChiLog chi_log;

//###################################################################
/**Performs slot dynamics for Polyhedron cell.*/
void chi_mesh::sweep_management::PRIMARY_FLUDS::
  SlotDynamics(chi_mesh::Cell *cell,
               chi_mesh::sweep_management::SPDS* spds,
               std::vector<std::vector<std::pair<int,short>>>& lock_boxes,
               std::vector<std::pair<int,short>>& delayed_lock_box,
               std::set<int>& location_boundary_dependency_set)
{
  chi_mesh::MeshContinuum* grid = spds->grid;

  chi_mesh::Vector ihat(1.0,0.0,0.0);
  chi_mesh::Vector jhat(0.0,1.0,0.0);
  chi_mesh::Vector khat(0.0,0.0,1.0);

  //=================================================== Loop over faces
  //           INCIDENT                                 but process
  //                                                    only incident faces
  std::vector<int> inco_face_face_category;
  int bndry_face_counter = 0;
  for (short f=0; f < cell->faces.size(); f++)
  {
    CellFace&  face = cell->faces[f];
    double     mu   = spds->omega.Dot(face.normal);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Incident face
    if (mu<(0.0-1.0e-16))
    {
      int neighbor = face.neighbor;

      //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ LOCAL CELL DEPENDENCE
      if (grid->IsCellLocal(neighbor))
      {
        size_t num_face_dofs = face.vertex_ids.size();
        size_t face_categ = grid->MapFaceHistogramBins(num_face_dofs);

        inco_face_face_category.push_back(face_categ);

        LockBox& lock_box = lock_boxes[face_categ];

        //========================================== Check if part of cyclic
        //                                           dependency
        auto adj_cell = grid->cells[neighbor];
        bool is_cyclic = false;
        for (auto cyclic_dependency : spds->local_cyclic_dependencies)
        {
          if ((cyclic_dependency.first == cell->cell_local_id) &&
              (cyclic_dependency.second == adj_cell->cell_local_id) )
          {
            is_cyclic = true;
            inco_face_face_category.back() *= -1;
            inco_face_face_category.back() -= 1;
          }
        }
        if (is_cyclic) continue;

        //======================================== Find associated face for
        //                                         dof mapping and lock box
        short ass_face = grid->FindAssociatedFace(face, neighbor);

        //Now find the cell (index,face) pair in the lock box and empty slot
        bool found = false;
        for (int k=0; k<lock_box.size(); k++)
        {
          if ((lock_box[k].first == neighbor) &&
              (lock_box[k].second== ass_face))
          {
            lock_box[k].first = -1;
            lock_box[k].second= -1;
            found = true;
            break;
          }
        }
        if (!found)
        {
          chi_log.Log(LOG_ALLERROR)
            << "Lock-box location not found in call to "
            << "InitializeAlphaElements. Local Cell "
            << cell->cell_local_id
            << " face " << f
            << " looking for cell "
            << adj_cell->cell_local_id
            << " face " << ass_face
            << " cat: " << face_categ
            << " omg=" << spds->omega.PrintS()
            << " lbsize=" << lock_box.size();
          exit(EXIT_FAILURE);
        }

      }//if local
      //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ BOUNDARY DEPENDENCE
      else if (grid->IsCellBndry(neighbor))
      {
        chi_mesh::Vector& face_norm = cell->faces[f].normal;

        if (face_norm.Dot(ihat)>0.999)
          location_boundary_dependency_set.insert(0);
        else if (face_norm.Dot(ihat)<-0.999)
          location_boundary_dependency_set.insert(1);
        else if (face_norm.Dot(jhat)>0.999)
          location_boundary_dependency_set.insert(2);
        else if (face_norm.Dot(jhat)<-0.999)
          location_boundary_dependency_set.insert(3);
        else if (face_norm.Dot(khat)>0.999)
          location_boundary_dependency_set.insert(4);
        else if (face_norm.Dot(khat)<-0.999)
          location_boundary_dependency_set.insert(5);
      }
    }//if incident

  }//for f

  so_cell_inco_face_face_category.push_back(inco_face_face_category);

  //=================================================== Loop over faces
  //                OUTGOING                            but process
  //                                                    only outgoing faces
  std::vector<int>                outb_face_slot_indices;
  std::vector<int>                outb_face_face_category;
  for (short f=0; f < cell->faces.size(); f++)
  {
    CellFace&  face         = cell->faces[f];
    double     mu           = spds->omega.Dot(face.normal);
    int        neighbor     = face.neighbor;
    int        cell_g_index = cell->cell_global_id;

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Outgoing face
    if (mu>=(0.0+1.0e-16))
    {
      size_t num_face_dofs = face.vertex_ids.size();
      size_t face_categ = grid->MapFaceHistogramBins(num_face_dofs);

      outb_face_face_category.push_back(face_categ);

      LockBox* temp_lock_box = &lock_boxes[face_categ];

      //========================================== Check if part of cyclic
      //                                           dependency
      if (grid->IsCellLocal(neighbor))
      {
        auto adj_cell = grid->cells[neighbor];
        for (auto cyclic_dependency : spds->local_cyclic_dependencies)
        {
          if ((cyclic_dependency.first == cell->cell_local_id) &&
              (cyclic_dependency.second == adj_cell->cell_local_id) )
          {
            temp_lock_box = &delayed_lock_box;
            outb_face_face_category.back() *= -1;
            outb_face_face_category.back() -= 1;
          }
        }
      }

      LockBox& lock_box = *temp_lock_box;

      //========================================== Check if this face is
      //                                           the max size
      if (num_face_dofs>largest_face)
        largest_face = num_face_dofs;

      //========================================== Find a open slot
      bool slot_found = false;
      for (int k=0; k<lock_box.size(); k++)
      {
        if (lock_box[k].first < 0)
        {
          outb_face_slot_indices.push_back(k);
          lock_box[k].first = cell_g_index;
          lock_box[k].second= f;
          slot_found = true;
          break;
        }
      }
      //========================================= If an open slot was not found
      //                                          push a new one
      if (!slot_found)
      {
        outb_face_slot_indices.push_back(lock_box.size());
        lock_box.push_back(std::pair<int,short>(cell_g_index,f));
      }

      //========================================== Non-local outgoing
      if ((!grid->IsCellLocal(neighbor)) && (!grid->IsCellBndry(neighbor)))
      {
        auto adj_cell     = grid->cells[neighbor];
        int  locJ         = adj_cell->partition_id;
        int  deplocI      = spds->MapLocJToDeplocI(locJ);
        int  face_slot    = deplocI_face_dof_count[deplocI];

        deplocI_face_dof_count[deplocI]+= face.vertex_ids.size();

        nonlocal_outb_face_deplocI_slot.emplace_back(deplocI,face_slot);

        AddFaceViewToDepLocI(deplocI, cell_g_index,
                             face_slot, face);

      }//non-local neighbor
    }//if outgoing

  }//for f

  so_cell_outb_face_slot_indices.push_back(outb_face_slot_indices);
  so_cell_outb_face_face_category.push_back(outb_face_face_category);
}


