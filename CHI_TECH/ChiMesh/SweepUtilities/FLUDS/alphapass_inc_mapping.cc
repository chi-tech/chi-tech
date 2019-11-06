#include "FLUDS.h"
#include "ChiMesh/SweepUtilities/SPDS/SPDS.h"

#include <ChiMesh/Cell/cell_newbase.h>

typedef std::vector<std::pair<int,short>> LockBox;

#include <chi_log.h>

extern ChiLog chi_log;

//###################################################################
/**Performs Incident mapping for Polyhedron cell.*/
void chi_mesh::sweep_management::FLUDS::
LocalIncidentMapping(chi_mesh::CellBase *cell_base,
                     chi_mesh::sweep_management::SPDS* spds,
                     std::vector<int>&  local_so_cell_mapping)
{
  chi_mesh::MeshContinuum* grid = spds->grid;
  std::vector<std::pair<int,std::vector<int>>> inco_face_dof_mapping;

  short        incoming_face_count=-1;

  //=================================================== Loop over faces
  //           INCIDENT                                 but process
  //                                                    only incident faces
  for (short f=0; f < cell_base->faces.size(); f++)
  {
    CellFace& face = cell_base->faces[f];
    double     mu  = face.normal.Dot(spds->omega);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Incident face
    if (mu<0.0)
    {
      int neighbor = face.neighbor;

      if (grid->IsCellLocal(neighbor))
      {
        incoming_face_count++;
        //======================================== Find associated face for
        //                                         dof mapping
        int ass_face = grid->FindAssociatedFace(face, neighbor);

        std::pair<int,std::vector<int>> dof_mapping;
        grid->FindAssociatedVertices(face,
                                     neighbor,
                                     ass_face,
                                     dof_mapping.second);

        //======================================== Find associated face
        //                                         counter for slot lookup
        auto adj_cell     = grid->cells[neighbor];
        int  adj_so_index = local_so_cell_mapping[adj_cell->cell_local_id];
        int  ass_f_counter=-1;

        auto adj_cell_base = (chi_mesh::CellBase*)adj_cell;
        int out_f = -1;
        for (short af=0; af < adj_cell_base->faces.size(); af++)
        {
          double mur = adj_cell_base->faces[af].normal.Dot(spds->omega);

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
          grid->FindAssociatedFace(face, neighbor, true);
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