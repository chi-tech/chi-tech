#include "chi_FLUDS.h"
#include "chi_SPDS.h"

#include "../../ChiMesh/Cell/cell_polyhedron.h"

//###################################################################
/**Performs slot dynamics for Polyhedron cell.*/
void chi_mesh::SweepManagement::FLUDS::
  SlotDynamics(TPolyhedron *polyh_cell,
               chi_mesh::SweepManagement::SPDS* spds,
               std::vector<std::pair<int,short>>& lock_box,
               std::set<int>& location_boundary_dependency_set)
{
  chi_mesh::MeshContinuum* grid = spds->grid;

  chi_mesh::Vector ihat(1.0,0.0,0.0);
  chi_mesh::Vector jhat(0.0,1.0,0.0);
  chi_mesh::Vector khat(0.0,0.0,1.0);

  short        outgoing_face_count=0;

  //=================================================== Loop over faces
  //           INCIDENT                                 but process
  //                                                    only incident faces
  int bndry_face_counter = 0;
  for (short f=0; f<polyh_cell->faces.size(); f++)
  {
    TPolyFace* poly_face = polyh_cell->faces[f];
    double     mu        = spds->omega.Dot(poly_face->geometric_normal);

    int lock_box_bound = lock_box.size();

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Incident face
    if (mu<0.0)
    {
      int neighbor = poly_face->face_indices[NEIGHBOR];

      //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ LOCAL CELL DEPENDENCE
      if (grid->IsCellLocal(neighbor))
      {
        //======================================== Find associated face for
        //                                         dof mapping and lock box
        short ass_face = grid->FindAssociatedFace(poly_face,neighbor);

        //Now find the cell (index,face) pair in the lock box and empty slot
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
        if (!found)
        {
          chi_log.Log(LOG_ALLERROR)
            << "Lock-box location not found in call to "
            << "InitializeAlphaElements";
          exit(EXIT_FAILURE);
        }

      }//if local
      else if (grid->IsCellBndry(neighbor))
      {
        chi_mesh::Vector& face_norm = polyh_cell->faces[f]->geometric_normal;

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

  //=================================================== Loop over faces
  //                OUTGOING                            but process
  //                                                    only outgoing faces
  std::vector<int>                outb_face_slot_indices;
  for (short f=0; f<polyh_cell->faces.size(); f++)
  {
    TPolyFace* poly_face = polyh_cell->faces[f];
    double     mu        = spds->omega.Dot(poly_face->geometric_normal);
    int        neighbor  = poly_face->face_indices[NEIGHBOR];
    int        cell_g_index = polyh_cell->cell_global_id;

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Outgoing face
    if (mu>=0.0)
    {
      //========================================== Check if this face is
      //                                           the max size
      if (poly_face->edges.size()>largest_face)
        largest_face = poly_face->edges.size();

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

        deplocI_face_dof_count[deplocI]+= poly_face->v_indices.size();

        nonlocal_outb_face_deplocI_slot.
          push_back(std::pair<int,int>(deplocI,face_slot));

        AddFaceViewToDepLocI(deplocI,cell_g_index,
                             face_slot,poly_face);

      }//non-local neighbor
    }//if outgoing

  }//for f

  so_cell_outb_face_slot_indices.push_back(outb_face_slot_indices);
}


//###################################################################
/**Performs Incident mapping for Polyhedron cell.*/
void chi_mesh::SweepManagement::FLUDS::
  IncidentMapping(TPolyhedron *polyh_cell,
                  chi_mesh::SweepManagement::SPDS* spds,
                  std::vector<int>&  local_so_cell_mapping)
{
  chi_mesh::MeshContinuum* grid = spds->grid;
  std::vector<std::pair<int,std::vector<int>>> inco_face_dof_mapping;

  short        incoming_face_count=-1;

  //=================================================== Loop over faces
  //           INCIDENT                                 but process
  //                                                    only incident faces
  for (short f=0; f<polyh_cell->faces.size(); f++)
  {
    TPolyFace* poly_face = polyh_cell->faces[f];
    double     mu        = poly_face->geometric_normal.Dot(spds->omega);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Incident face
    if (mu<0.0)
    {
      int neighbor = poly_face->face_indices[NEIGHBOR];

      if (grid->IsCellLocal(neighbor))
      {
        incoming_face_count++;
        //======================================== Find associated face for
        //                                         dof mapping
        int ass_face = grid->FindAssociatedFace(poly_face,neighbor);

        std::pair<int,std::vector<int>> dof_mapping;
        grid->FindAssociatedVertices(poly_face,
                                     neighbor,
                                     ass_face,
                                     dof_mapping.second);

        //======================================== Find associated face
        //                                         counter for slot lookup
        auto adj_cell     = grid->cells[neighbor];
        int  adj_so_index = local_so_cell_mapping[adj_cell->cell_local_id];
        int  ass_f_counter=-1;

        auto adj_polyh_cell = (chi_mesh::CellPolyhedron*)adj_cell;
        int out_f = -1;
        for (short af=0; af<adj_polyh_cell->faces.size(); af++)
        {
          double mur = adj_polyh_cell->faces[af]->
            geometric_normal.Dot(spds->omega);

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
            << "Associated face counter not found"
            << ass_face << " " << neighbor;
          grid->FindAssociatedFace(poly_face,neighbor,true);
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