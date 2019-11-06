#include "mc_bndry_source.h"

#include "../../RandomNumberGenerator/montecarlon_rng.h"

#include <ChiMesh/MeshContinuum/chi_meshcontinuum.h>
#include <ChiMesh/MeshHandler/chi_meshhandler.h>
#include <ChiMesh/VolumeMesher/chi_volumemesher.h>
#include <ChiMesh/VolumeMesher/Linemesh1D/volmesher_linemesh1d.h>
#include <ChiMesh/Cell/cell_slabv2.h>
#include <ChiMesh/Cell/cell_polygonv2.h>
#include <ChiMesh/Cell/cell_polyhedronv2.h>

#include <FiniteVolume/fv.h>
#include <FiniteVolume/CellViews/fv_slab.h>
#include <FiniteVolume/CellViews/fv_polygon.h>

#include <ChiMath/Statistics/cdfsampler.h>

#include <chi_log.h>

extern ChiLog chi_log;


//###################################################################
/**Default constructor.*/
chi_montecarlon::BoundarySource::BoundarySource()
{
  type_index = MC_BNDRY_SRC;
  ref_bndry = MC_ALL_BOUNDARIES;
}

//###################################################################
/**Initializes references to necessary mesh elements.
 *
 * Since this is a boundary source the major difficulty here is
 * to get the sampling direction right. Sampled directions need
 * to be rotation to the opposite of the outward pointing boundary
 * normals. For this we assemble a rotation matrix for each face
 * that needs to be sampled.*/
void chi_montecarlon::BoundarySource::
  Initialize(chi_mesh::MeshContinuum* ref_grid,
             SpatialDiscretization_FV*   ref_fv_sdm)
{
  grid = ref_grid;
  fv_sdm = ref_fv_sdm;

  int num_local_cells = grid->local_cell_glob_indices.size();
  for (int lc=0; lc<num_local_cells; lc++)
  {
    int cell_glob_index = grid->local_cell_glob_indices[lc];
    auto cell = grid->cells[cell_glob_index];

    if (cell->Type() == chi_mesh::CellType::CELL_NEWBASE)
    {
      auto cell_base = (chi_mesh::CellBase*)cell;

      //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ SLAB
      if (cell_base->Type2() == chi_mesh::CellType::SLABV2)
      {
        auto slab_cell = (chi_mesh::CellSlabV2*)cell_base;

        int num_faces = 2;
        for (int f=0; f<num_faces; f++)
        {
          int neighbor = slab_cell->faces[f].neighbor;

          //================================== If all boundaries
          if ((neighbor < 0) and (ref_bndry == MC_ALL_BOUNDARIES))
          {
            FACE_REF* new_face_ref = new FACE_REF(cell_glob_index,f);
            chi_mesh::Matrix3x3& R = new_face_ref->RotMatrix;

            R.SetRowIVec(0,chi_mesh::Vector(1.0,0.0,0.0));
            R.SetRowIVec(1,chi_mesh::Vector(0.0,1.0,0.0));
            R.SetRowIVec(2,slab_cell->faces[f].normal*-1.0);

            ref_cell_faces.push_back(new_face_ref);
          }
            //================================== If specific bndries
          else if ((neighbor < 0) and (ref_bndry == abs(neighbor)))
          {
            FACE_REF* new_face_ref = new FACE_REF(cell_glob_index,f);
            chi_mesh::Matrix3x3& R = new_face_ref->RotMatrix;

            R.SetRowIVec(0,chi_mesh::Vector(1.0,0.0,0.0));
            R.SetRowIVec(1,chi_mesh::Vector(0.0,1.0,0.0));
            R.SetRowIVec(2,slab_cell->faces[f].normal*-1.0);

            ref_cell_faces.push_back(new_face_ref);
          }

        }//for faces
      }//if slab
        // $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ POLYGON
      else if (cell_base->Type2() == chi_mesh::CellType::POLYGONV2)
      {
        auto poly_cell = (chi_mesh::CellPolygonV2*)cell_base;
        auto cell_fv_view = (PolygonFVView*)fv_sdm->MapFeView(cell_glob_index);

        int num_faces = poly_cell->faces.size();
        for (int f=0; f<num_faces; f++)
        {
          int neighbor = poly_cell->faces[f].neighbor;

          //================================== If all boundaries
          if ((neighbor < 0) and (ref_bndry == MC_ALL_BOUNDARIES))
          {
            FACE_REF* new_face_ref = new FACE_REF(cell_glob_index,f);
            chi_mesh::Matrix3x3& R = new_face_ref->RotMatrix;

            chi_mesh::Vector n = poly_cell->faces[f].normal*-1.0;
            chi_mesh::Vector khat(0.0,0.0,1.0);

            if (n.Dot(khat) > 0.9999)
              R.SetDiagonalVec(1.0,1.0,1.0);
            else
            {
              chi_mesh::Vector binorm = khat.Cross(n);
              binorm = binorm/binorm.Norm();

              chi_mesh::Vector tangent = binorm.Cross(n);
              tangent = tangent/tangent.Norm();

              R.SetColJVec(0,tangent);
              R.SetColJVec(1,binorm);
              R.SetColJVec(2,n);
            }

            new_face_ref->area = cell_fv_view->face_area[f];

            ref_cell_faces.push_back(new_face_ref);
          }
            //================================== If specific bndries
          else if ((neighbor < 0) and (ref_bndry == abs(neighbor)))
          {
            FACE_REF* new_face_ref = new FACE_REF(cell_glob_index,f);
            chi_mesh::Matrix3x3& R = new_face_ref->RotMatrix;

            chi_mesh::Vector n = poly_cell->faces[f].normal*-1.0;
            chi_mesh::Vector khat(0.0,0.0,1.0);

            int v0i = poly_cell->faces[f].vertex_ids[0];
            int v1i = poly_cell->faces[f].vertex_ids[1];

            chi_mesh::Node v0 = *grid->nodes[v0i];
            chi_mesh::Node v1 = *grid->nodes[v1i];

//          chi_log.Log(LOG_0) << v0.PrintS();
//          chi_log.Log(LOG_0) << v1.PrintS();
//
//          chi_log.Log(LOG_0) << n.PrintS();

            if (n.Dot(khat) > 0.9999)
              R.SetDiagonalVec(1.0,1.0,1.0);
            else
            {
              chi_mesh::Vector binorm = khat.Cross(n);
              binorm = binorm/binorm.Norm();

              chi_mesh::Vector tangent = binorm.Cross(n);
              tangent = tangent/tangent.Norm();

              R.SetColJVec(0,tangent);
              R.SetColJVec(1,binorm);
              R.SetColJVec(2,n);
            }
//
//          chi_log.Log(LOG_0) << R.PrintS();

            new_face_ref->area = cell_fv_view->face_area[f];

            ref_cell_faces.push_back(new_face_ref);
          }
        }
      }//for polygon
      else
      {
        chi_log.Log(LOG_ALLERROR)
          << "chi_montecarlon: Unsupported cell type encountered in "
          << "call to BoundarySource::Initialize.";
        exit(EXIT_FAILURE);
      }
    }//new cell base


  }//for local cell

  //============================================= Build cdf for surface
  //                                              sampling
  double total_area = 0.0;
  int num_ref_faces = ref_cell_faces.size();
  for (int rf=0; rf<num_ref_faces; rf++)
    total_area += ref_cell_faces[rf]->area;

  double intgl = 0.0;
  face_cdf.resize(num_ref_faces);
  for (int rf=0; rf<num_ref_faces; rf++)
  {
    intgl += ref_cell_faces[rf]->area;
    face_cdf[rf] = intgl/total_area;
  }

  surface_sampler = new chi_math::CDFSampler(face_cdf);
}

//###################################################################
/**Source routine for boundary source.*/
chi_montecarlon::Particle chi_montecarlon::BoundarySource::
  CreateParticle(chi_montecarlon::RandomNumberGenerator *rng)
{
  chi_mesh::MeshHandler*    mesh_handler = chi_mesh::GetCurrentHandler();
  chi_mesh::VolumeMesher*         mesher = mesh_handler->volume_mesher;

  chi_montecarlon::Particle new_particle;

  chi_mesh::Cell* first_cell = grid->cells[0];

  if (first_cell->Type() == chi_mesh::CellType::CELL_NEWBASE)
  {
    auto cell_base = (chi_mesh::CellBase*)first_cell;

    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ SLAB
    // Categorally slab surface sampling does not require
    // a cdf and therefore we have a specific routine for it.
    if (cell_base->Type2() == chi_mesh::CellType::SLABV2)
    {
      //====================================== Sample ref face
      int num_ref_faces = ref_cell_faces.size();
      int f = std::round(rng->Rand()*(num_ref_faces-1));
      FACE_REF* face_ref = ref_cell_faces[f];

      auto slab_cell =
        (chi_mesh::CellSlabV2*)grid->cells[face_ref->cell_glob_index];
      int ref_vert_index = slab_cell->vertex_ids[face_ref->face_num];
      chi_mesh::Vertex* ref_vert = grid->nodes[ref_vert_index];

      //====================================== Sample direction
      double costheta = rng->Rand();     //Sample half-range only
      double theta    = acos(costheta);
      double varphi   = rng->Rand()*2.0*M_PI;

      chi_mesh::Vector ref_dir;
      ref_dir.x = sin(theta)*cos(varphi);
      ref_dir.y = sin(theta)*sin(varphi);
      ref_dir.z = cos(theta);

      //====================================== Set quantities
      new_particle.pos = *ref_vert;
      new_particle.dir = face_ref->RotMatrix*ref_dir;

      new_particle.egrp = 0;
      new_particle.w = costheta;
      new_particle.cur_cell_ind = face_ref->cell_glob_index;

      return new_particle;
    }//Slab cells
    else if (cell_base->Type2() == chi_mesh::CellType::POLYGONV2)
    {
      int f = surface_sampler->Sample(rng->Rand());

      FACE_REF* face_ref = ref_cell_faces[f];
      auto poly_cell =
        (chi_mesh::CellPolygonV2*)grid->cells[face_ref->cell_glob_index];

      chi_mesh::CellFace& ref_face = poly_cell->faces[face_ref->face_num];

      chi_mesh::Vertex v0 = *grid->nodes[ref_face.vertex_ids[0]];
      chi_mesh::Vertex v1 = *grid->nodes[ref_face.vertex_ids[1]];

      //====================================== Sample direction
      double costheta = rng->Rand();     //Sample half-range only
      double theta    = acos(costheta);
      double varphi   = rng->Rand()*2.0*M_PI;

      chi_mesh::Vector ref_dir;
      ref_dir.x = sin(theta)*cos(varphi);
      ref_dir.y = sin(theta)*sin(varphi);
      ref_dir.z = cos(theta);

//    chi_log.Log(LOG_0) << v0.PrintS();
//    chi_log.Log(LOG_0) << v1.PrintS();

      //====================================== Set quantities
      double w = rng->Rand();
      new_particle.pos = v0*w + v1*(1.0-w);
      new_particle.dir = face_ref->RotMatrix*ref_dir;

      new_particle.egrp = 0;
      new_particle.w = costheta;
      new_particle.cur_cell_ind = face_ref->cell_glob_index;

      return new_particle;
    }
    else
    {
      chi_log.Log(LOG_ALLERROR)
        << "chi_montecarlon: Unsupported cell type encountered in "
        << "call to BoundarySource::CreateParticle.";
      exit(EXIT_FAILURE);
    }
  }//new cell base





  return new_particle;
}