#include "chi_FLUDS.h"
#include "chi_SPDS.h"

#include "../../ChiMesh/Cell/cell_slab.h"

typedef std::vector<std::pair<int,short>> LockBox;

#include <chi_log.h>

extern ChiLog chi_log;

//###################################################################
/**Performs slot dynamics for Polyhedron cell.*/
void chi_mesh::SweepManagement::FLUDS::
SlotDynamics(TSlab *slab_cell,
             chi_mesh::SweepManagement::SPDS* spds,
             std::vector<std::vector<std::pair<int,short>>>& lock_boxes,
             std::vector<std::pair<int,short>>& delayed_lock_box,
             std::set<int>& location_boundary_dependency_set)
{
  chi_mesh::MeshContinuum* grid = spds->grid;

  chi_mesh::Vector ihat(1.0,0.0,0.0);
  chi_mesh::Vector jhat(0.0,1.0,0.0);
  chi_mesh::Vector khat(0.0,0.0,1.0);

  short        outgoing_face_count=0;

  LockBox& lock_box = lock_boxes.front();

  //=================================================== Loop over faces
  //           INCIDENT                                 but process
  //                                                    only incident faces
  std::vector<int> inco_face_face_category;
  int bndry_face_counter = 0;
  int num_faces = 2;
  for (short f=0; f<num_faces; f++)
  {
    double mu = spds->omega.Dot(slab_cell->face_normals[f]);

    int lock_box_bound = lock_box.size();

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Incident face
    if (mu<0.0)
    {
      int neighbor = slab_cell->edges[f];

      //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ LOCAL CELL DEPENDENCE
      if (grid->IsCellLocal(neighbor))
      {
        inco_face_face_category.push_back(0);

        //======================================== Find associated face for
        //                                         dof mapping and lock box
        short ass_face = 0;
        if (f==0)
          ass_face = 1;
        else
          ass_face = 0;

        //Now find the cell (index,face) pair in the lock box and empty slot
        //if possible
        bool found = false;
        for (int k=0; k<lock_box_bound; k++)
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
        //If it is was not found
        if (!found)
        {
          chi_log.Log(LOG_ALLERROR)
            << "Lock-box location not found in call to "
            << "InitializeAlphaElements:Slab";
          exit(EXIT_FAILURE);
        }

      }//if local
      else if (grid->IsCellBndry(neighbor))
      {
        chi_mesh::Vector& face_norm = slab_cell->face_normals[f];

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
  for (short f=0; f<num_faces; f++)
  {
    double mu        = spds->omega.Dot(slab_cell->face_normals[f]);
    int    neighbor  = slab_cell->edges[f];
    int    cell_g_index = slab_cell->cell_global_id;

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Outgoing face
    if (mu>=0.0)
    {
      outb_face_face_category.push_back(0);

      //========================================== Check if this face is
      //                                           the max size
      largest_face = 1;

      outgoing_face_count++;

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

        deplocI_face_dof_count[deplocI]+= 1;

        nonlocal_outb_face_deplocI_slot.
          push_back(std::pair<int,int>(deplocI,face_slot));

        AddFaceViewToDepLocI(deplocI,cell_g_index,
                             face_slot,slab_cell->v_indices[f]);

      }//if non-local neighbor
    }//if outgoing

  }//for f

  so_cell_outb_face_slot_indices.push_back(outb_face_slot_indices);
  so_cell_outb_face_face_category.push_back(outb_face_face_category);
}


//###################################################################
/**Performs Incident mapping for slab cell.*/
void chi_mesh::SweepManagement::FLUDS::
IncidentMapping(TSlab *slab_cell,
                chi_mesh::SweepManagement::SPDS* spds,
                std::vector<int>&  local_so_cell_mapping)
{
  chi_mesh::MeshContinuum* grid = spds->grid;
  std::vector<std::pair<int,std::vector<int>>> inco_face_dof_mapping;

  short        incoming_face_count=-1;

  //=================================================== Loop over faces
  //           INCIDENT                                 but process
  //                                                    only incident faces
  int num_faces = 2;
  for (short f=0; f<num_faces; f++)
  {
    double mu = slab_cell->face_normals[f].Dot(spds->omega);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Incident face
    if (mu<0.0)
    {
      int neighbor = slab_cell->edges[f];

      if (grid->IsCellLocal(neighbor))
      {
        incoming_face_count++;
        //======================================== Find associated face for
        //                                         dof mapping
        int ass_face = 0;
        if (f == 0)
          ass_face = 1;
        else
          ass_face = 0;

//        grid->FindAssociatedFace(poly_face,neighbor);

        std::pair<int,std::vector<int>> dof_mapping;
        dof_mapping.second.push_back(0);

//        grid->FindAssociatedVertices(poly_face,
//                                     neighbor,
//                                     ass_face,
//                                     dof_mapping.second);

        //======================================== Find associated face
        //                                         counter for slot lookup
        auto adj_cell     = grid->cells[neighbor];
        int  adj_so_index = local_so_cell_mapping[adj_cell->cell_local_id];
        int  ass_f_counter=-1;

        auto adj_slab_cell = (chi_mesh::CellSlab*)adj_cell;
        int out_f = -1;
        for (short af=0; af<2; af++)
        {
          double mur = adj_slab_cell->face_normals[af].Dot(spds->omega);

          if (mur>=0.0) {out_f++;}
          if (af == ass_face)
          {
            ass_f_counter = out_f;
            break;
          }
        }
        if (ass_f_counter<0)
        {
          chi_log.Log(LOG_ALLERROR)
            << "Slab. Associated face counter not found "
            << ass_face << " " << neighbor;
          exit(EXIT_FAILURE);
        }

        dof_mapping.first = /*local_psi_stride*G**/
          so_cell_outb_face_slot_indices[adj_so_index][ass_f_counter];

        dof_mapping.second.shrink_to_fit();
        inco_face_dof_mapping.push_back(dof_mapping);
      }//if local
    }//if incident
  }//for incindent f

  so_cell_inco_face_dof_indices.push_back(inco_face_dof_mapping);
}