#ifndef _npt_sweepchunk_pwl_h
#define _npt_sweepchunk_pwl_h


#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiMesh/SweepUtilities/sweep_namespace.h"
#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "ChiMath/SpatialDiscretization/PiecewiseLinear/pwl.h"
#include "ChiMath/SpatialDiscretization/PiecewiseLinear/CellViews/pwl_polyhedron.h"

#include "ChiMesh/Cell/cell.h"
#include <ChiPhysics/chi_physics.h>

#include "ChiMath/chi_math.h"
#include "../GroupSet/lbs_groupset.h"
#include "../lbs_linear_boltzman_solver.h"
#include "ChiMath/Quadratures/product_quadrature.h"

#include "ChiMesh/SweepUtilities/SPDS/SPDS.h"
#include "ChiMesh/SweepUtilities/AngleAggregation/angleaggregation.h"

#include "ChiTimer/chi_timer.h"

#include <chi_mpi.h>
#include <chi_log.h>

extern ChiMath&     chi_math_handler;
extern ChiMPI&      chi_mpi;
extern ChiLog&     chi_log;

typedef std::vector<chi_physics::TransportCrossSections*> TCrossSections;

//###################################################################
/**Sweep chunk to compute the fixed source.*/
class LBSSweepChunkPWL : public chi_mesh::sweep_management::SweepChunk
{
protected:
  chi_mesh::MeshContinuum*    grid_view;
  SpatialDiscretization_PWL*     grid_fe_view;
  std::vector<LinearBoltzman::CellViewBase*>* grid_transport_view;
  LinearBoltzman::Solver&     ref_solver;
//std::vector<double>*        x;                   BASE CLASS
  std::vector<double>*        q_moments;
  LBSGroupset*               groupset;
  TCrossSections*             xsections;
  int                         num_moms;

  int                         G;
  int                         g;
  chi_mesh::Vector3            omega;
  double                      wn;

  int                         max_cell_dofs;

  bool                        a_and_b_initialized;

//bool                        suppress_surface_src; BASE CLASS

  std::vector<std::vector<double>> Amat;
  std::vector<std::vector<double>> Atemp;
  std::vector<std::vector<double>> b;
  std::vector<double>              source;

  int LOCAL;
  double test_source;
  std::vector<double> test_mg_src;
  std::vector<double> zero_mg_src;

public:
  //################################################## Constructor
  LBSSweepChunkPWL(chi_mesh::MeshContinuum* vol_continuum,
                   SpatialDiscretization_PWL* discretization,
                   std::vector<LinearBoltzman::CellViewBase*>* cell_transport_views,
                   LinearBoltzman::Solver& in_ref_solver,
                   std::vector<double>* destination_phi,
                   std::vector<double>* source_moments,
                   LBSGroupset* in_groupset,
                   TCrossSections* in_xsections,
                   int in_num_moms,
                   int in_max_cell_dofs) :
                   ref_solver(in_ref_solver)
  {
    grid_view           = vol_continuum;
    grid_fe_view        = discretization;
    grid_transport_view = cell_transport_views;
    x                   = destination_phi;
    q_moments           = source_moments;
    groupset            = in_groupset;
    xsections           = in_xsections;
    num_moms            = in_num_moms;
    max_cell_dofs       = in_max_cell_dofs;

    G                   = in_groupset->groups.size();

    a_and_b_initialized = false;
    suppress_surface_src= false;

    LOCAL = chi_mpi.location_id;

    test_source = 100.0/4.0/M_PI;
    test_mg_src.resize(G,test_source);
    test_mg_src[0] = test_source;
    zero_mg_src.resize(G,0.0);
  }


  //############################################################ Actual chunk
  virtual void Sweep(chi_mesh::sweep_management::AngleSet* angle_set)
  {
    int outface_master_counter=0;

    if (!a_and_b_initialized)
    {
      Amat.resize(max_cell_dofs,std::vector<double>(max_cell_dofs));
      Atemp.resize(max_cell_dofs,std::vector<double>(max_cell_dofs));
      b.resize(G,std::vector<double>(max_cell_dofs,0.0));
      source.resize(max_cell_dofs,0.0);

      a_and_b_initialized = true;
    }

    chi_mesh::sweep_management::SPDS* spds = angle_set->GetSPDS();
    chi_mesh::sweep_management::FLUDS* fluds = angle_set->fluds;

    GsSubSet& subset = groupset->grp_subsets[angle_set->ref_subset];
    int gs_ss_size  = groupset->grp_subset_sizes[angle_set->ref_subset];

    int gs_ss_begin = subset.first;
    int gs_ss_end   = subset.second;

    //Groupset subset first group number
    int gs_gi = groupset->groups[gs_ss_begin]->id;


    int deploc_face_counter = -1;
    int preloc_face_counter = -1;
    int bndry_face_counter  = -1;

    double* phi        = x->data();
    double* psi        = zero_mg_src.data();
    double* q_mom      = q_moments->data();


    //========================================================== Loop over each cell
    size_t num_loc_cells = spds->spls->item_id.size();
    for (int cr_i=0; cr_i<num_loc_cells; cr_i++)
    {
      int  cell_local_id = spds->spls->item_id[cr_i];
      auto cell          = &grid_view->local_cells[cell_local_id];
      int  cell_g_index  = cell->global_id;

      auto cell_fe_view   = (CellFEView*)grid_fe_view->MapFeViewL(cell->local_id);
      auto transport_view =
        (LinearBoltzman::CellViewFull*)(*grid_transport_view)[cell->local_id];

      int     cell_dofs    = cell_fe_view->dofs;
      int     xs_id        = transport_view->xs_id;
      double* sigma_tg = (*xsections)[xs_id]->sigma_tg.data();

      std::vector<bool> face_incident_flags(cell->faces.size(),false);


      //=================================================== Get Cell matrices
      const std::vector<std::vector<chi_mesh::Vector3>>& L =
        cell_fe_view->IntV_shapeI_gradshapeJ;

      const std::vector<std::vector<double>>& M =
        cell_fe_view->IntV_shapeI_shapeJ;

      const std::vector<std::vector<std::vector<double>>>& N =
        cell_fe_view->IntS_shapeI_shapeJ;

      //=================================================== Loop over angles in set
      int ni_deploc_face_counter = deploc_face_counter;
      int ni_preloc_face_counter = preloc_face_counter;
      int ni_bndry_face_counter  = bndry_face_counter;
      for (int n=0; n<angle_set->angles.size(); n++)
      {
        deploc_face_counter = ni_deploc_face_counter;
        preloc_face_counter = ni_preloc_face_counter;
        bndry_face_counter  = ni_bndry_face_counter;

        angle_num = angle_set->angles[n];
        omega = groupset->quadrature->omegas[angle_num];
        wn    = groupset->quadrature->weights[angle_num];

        //============================================ Gradient matrix
        for (int i=0; i<cell_dofs; i++)
        {
          for (int j=0; j<cell_dofs; j++)
          {
            Amat[i][j] = omega.Dot(L[i][j]);
          }//for j
        }//for i

        for (int gsg=0; gsg<gs_ss_size; gsg++)
          b[gsg].assign(cell_dofs,0.0);


        //============================================ Surface integrals
        int num_faces = cell->faces.size();
        int in_face_counter=-1;
        int internal_face_bndry_counter = -1;
        for (int f=0; f<num_faces; f++)
        {

          double mu              = omega.Dot(cell->faces[f].normal);
          auto& face             = cell->faces[f];
          bool  face_on_boundary = grid_view->IsCellBndry(face.neighbor);
          bool neighbor_is_local = (transport_view->face_local[f]);

          //============================= Set flags
          if (mu>=0.0) face_incident_flags[f] = false;
          else         face_incident_flags[f] = true;

          //This counter update-logic is for mapping an incident boundary
          //condition. Because it is cheap, the cell faces was mapped to a
          //corresponding boundary during initialization and is
          //independent of angle. Accessing things like reflective boundary
          //angular fluxes (and complex boundary conditions), requires the
          //more general bndry_face_counter.
          int bndry_map = -1;
          if (face.neighbor<0)
          {
            internal_face_bndry_counter++;
            bndry_map = -(face.neighbor+1);
          }

          if (mu < 0.0) //UPWIND
          {
            //============================== Increment face counters
            if (neighbor_is_local)
              in_face_counter++;

            if ((not neighbor_is_local) && (not face_on_boundary))
              preloc_face_counter++;

            if (face_on_boundary)
              bndry_face_counter++;


            //============================== Loop over face vertices
            int num_face_indices = cell->faces[f].vertex_ids.size();
            for (int fi=0; fi<num_face_indices; fi++)
            {
              int i = cell_fe_view->face_dof_mappings[f][fi];

              //=========== Loop over face unknowns
              for (int fj=0; fj<num_face_indices; fj++)
              {
                int j = cell_fe_view->face_dof_mappings[f][fj];

                // %%%%% LOCAL CELL DEPENDENCY %%%%%
                if (neighbor_is_local)
                {psi = fluds->UpwindPsi(cr_i,in_face_counter,fj,0,n);}
                  // %%%%% NON-LOCAL CELL DEPENDENCY %%%%%
                else if (not face_on_boundary)
                {psi = fluds->NLUpwindPsi(preloc_face_counter,fj,0,n);}
                  // %%%%% BOUNDARY CELL DEPENDENCY %%%%%
                else
                {psi = angle_set->PsiBndry(bndry_map,
                                           angle_num,
                                           cell->local_id,
                                           f,fj,gs_gi,gs_ss_begin,
                                           suppress_surface_src);
                }


                double mu_Nij = -mu*N[f][i][j];

                Amat[i][j] += mu_Nij;

                for (int gsg=0; gsg<gs_ss_size; gsg++)
                  b[gsg][i] += psi[gsg]*mu_Nij;
              }
            };

          }//if mu<0.0

        }//for f

        //========================================== Looping over groups
        double sigma_tgr = 0.0;
        double temp_src = 0.0;
        int gi_deploc_face_counter = deploc_face_counter;
        int gi_preloc_face_counter = preloc_face_counter;
        for (int gsg=0; gsg<gs_ss_size; gsg++)
        {
          deploc_face_counter = gi_deploc_face_counter;
          preloc_face_counter = gi_preloc_face_counter;

          g = gs_gi+gsg;

          //============================= Contribute source moments
          double m2d = 0.0;
          for (int i=0; i<cell_fe_view->dofs; i++)
          {
            temp_src = 0.0;
            for (int m=0; m<num_moms; m++)
            {
              m2d = groupset->m2d_op[m][angle_num];

              int ir = transport_view->MapDOF(i,m,g);
              temp_src += m2d*q_mom[ir];
            }
            source[i] = temp_src;
          }

          //============================= Mass Matrix and Source
          sigma_tgr = sigma_tg[g];
          for (int i=0; i<cell_fe_view->dofs; i++)
          {
            double temp = 0.0;
            for (int j=0; j<cell_fe_view->dofs; j++)
            {
              double Mij = M[i][j];
              Atemp[i][j] = Amat[i][j] + Mij*sigma_tgr;
              temp += Mij*source[j];
            }//for j
            b[gsg][i] += temp;
          }//for i


          //============================= Solve system
          chi_math::GaussElimination(Atemp,b[gsg],cell_fe_view->dofs);


        }//for g




        //============================= Accumulate flux
        double wn_d2m = 0.0;
        for (int m=0; m<num_moms; m++)
        {
          wn_d2m = groupset->d2m_op[m][angle_num];
          for (int i=0; i<cell_fe_view->dofs; i++)
          {
            int ir = transport_view->MapDOF(i,m,gs_gi);

            for (int gsg=0; gsg<gs_ss_size; gsg++)
              phi[ir+gsg] += wn_d2m*b[gsg][i];

            for (auto callback : groupset->moment_callbacks)
              if(callback)
                for (int gsg=0; gsg<gs_ss_size; gsg++)
                  callback(ir+gsg, m, angle_num, b[gsg][i]);
          }
        }





        //============================================= Outgoing fluxes
        int out_face_counter=-1;
        internal_face_bndry_counter = -1;
        for (int f=0; f<cell->faces.size(); f++)
        {
          if (face_incident_flags[f]) continue;

          //============================= Set flags and counters
          out_face_counter++;

          auto& face = cell->faces[f];
          bool  face_on_boundary = grid_view->IsCellBndry(face.neighbor);

          int bndry_index = -1;
          if (grid_view->IsCellBndry(face.neighbor))
          {
            internal_face_bndry_counter++;
            bndry_index = -(face.neighbor + 1);
          }



          //============================= Store outgoing Psi Locally
          if (transport_view->face_local[f])
          {
            for (int fi=0; fi<cell->faces[f].vertex_ids.size(); fi++)
            {
              int i = cell_fe_view->face_dof_mappings[f][fi];
              psi = fluds->OutgoingPsi(cr_i,out_face_counter,fi,n);

              for (int gsg=0; gsg<gs_ss_size; gsg++)
                psi[gsg] = b[gsg][i];
            }
          }//
          //============================= Store outgoing Psi Non-Locally
          else if (not face_on_boundary)
          {
            deploc_face_counter++;
            for (int fi=0; fi<cell->faces[f].vertex_ids.size(); fi++)
            {
              int i = cell_fe_view->face_dof_mappings[f][fi];
              psi = fluds->NLOutgoingPsi(deploc_face_counter,fi,n);

              for (int gsg=0; gsg<gs_ss_size; gsg++)
                psi[gsg] = b[gsg][i];
            }//for fdof
          }//if non-local
            //============================= Store outgoing reflecting Psi
          else if (angle_set->ref_boundaries[bndry_index]->IsReflecting())
          {
            for (int fi=0; fi<cell->faces[f].vertex_ids.size(); fi++)
            {
              int i = cell_fe_view->face_dof_mappings[f][fi];
              psi = angle_set->ReflectingPsiOutBoundBndry(bndry_index, angle_num,
                                                          cell->local_id, f,
                                                          fi, gs_ss_begin);

              for (int gsg=0; gsg<gs_ss_size; gsg++)
                psi[gsg] = b[gsg][i];
            }//for fdof
          }//reflecting
        }//for f


      }//for n

    }// for cell

  }//Sweep function
};//class def

#endif
