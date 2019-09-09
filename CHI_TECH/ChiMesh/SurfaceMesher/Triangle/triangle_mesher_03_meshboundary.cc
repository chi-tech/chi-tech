#include "triangle_mesher.h"
#include "../../../ChiMPI/chi_mpi.h"

extern ChiMPI chi_mpi;

#include <chi_log.h>

extern ChiLog chi_log;

extern "C"
{
#define REAL double
#define ANSI_DECLARATORS
#define VOID int
#include "triangle.h"
#undef REAL
#undef ANSI_DECLARATORS
#undef VOID
}
#define REAL double
#include "../../Boundary/chi_boundary.h"

/*****************************************************************************/
/*                                                                           */
/*  report()   Print the input or output.                                    */
/*                                                                           */
/*****************************************************************************/

void report(struct triangulateio *io,
            int  markers,int  reporttriangles,
            int  reportneighbors,int  reportsegments,
            int reportedges, int reportnorms)
{
  int i, j;

  for (i = 0; i < io->numberofpoints; i++) {
    printf("Point %4d:", i);
    for (j = 0; j < 2; j++) {
      printf("  %.6g", io->pointlist[i * 2 + j]);
    }
    if (io->numberofpointattributes > 0) {
      printf("   attributes");
    }
    for (j = 0; j < io->numberofpointattributes; j++) {
      printf("  %.6g",
             io->pointattributelist[i * io->numberofpointattributes + j]);
    }
    if (markers) {
      printf("   marker %d\n", io->pointmarkerlist[i]);
    } else {
      printf("\n");
    }
  }
  printf("\n");

  if (reporttriangles || reportneighbors) {
    for (i = 0; i < io->numberoftriangles; i++) {
      if (reporttriangles) {
        printf("Triangle %4d points:", i);
        for (j = 0; j < io->numberofcorners; j++) {
          printf("  %4d", io->trianglelist[i * io->numberofcorners + j]);
        }
        if (io->numberoftriangleattributes > 0) {
          printf("   attributes");
        }
        for (j = 0; j < io->numberoftriangleattributes; j++) {
          printf("  %.6g", io->triangleattributelist[i *
                                io->numberoftriangleattributes + j]);
        }
        printf("\n");
      }
      if (reportneighbors) {
        printf("Triangle %4d neighbors:", i);
        for (j = 0; j < 3; j++) {
          printf("  %4d", io->neighborlist[i * 3 + j]);
        }
        printf("\n");
      }
    }
    printf("\n");
  }

  if (reportsegments) {
    for (i = 0; i < io->numberofsegments; i++) {
      printf("Segment %4d points:", i);
      for (j = 0; j < 2; j++) {
        printf("  %4d", io->segmentlist[i * 2 + j]);
      }
      if (markers) {
        printf("   marker %d\n", io->segmentmarkerlist[i]);
      } else {
        printf("\n");
      }
    }
    printf("\n");
  }

  if (reportedges) {
    for (i = 0; i < io->numberofedges; i++) {
      printf("Edge %4d points:", i);
      for (j = 0; j < 2; j++) {
        printf("  %4d", io->edgelist[i * 2 + j]);
      }
      if (reportnorms && (io->edgelist[i * 2 + 1] == -1)) {
        for (j = 0; j < 2; j++) {
          printf("  %.6g", io->normlist[i * 2 + j]);
        }
      }
      if (markers) {
        printf("   marker %d\n", io->edgemarkerlist[i]);
      } else {
        printf("\n");
      }
    }
    printf("\n");
  }
}





void chi_mesh::SurfaceMesherTriangle::MeshBoundary(chi_mesh::Boundary *boundary)
{
  chi_mesh::SurfaceMesh* new_surface_mesh;
  chi_mesh::SurfaceMesh* temp_surface_mesh;

  //####################################################### Master process
  if (chi_mpi.location_id == 0)
  {
    struct triangulateio in, mid, out, vorout;
    struct triangulateio* output;

    chi_mesh::SurfaceMesh* surf_mesh =
      boundary->initial_mesh_continuum.surface_mesh;

    //================================================== Parse vertices
    in.numberofpoints = surf_mesh->vertices.size();
    in.numberofpointattributes = 0;
    in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));
    in.pointmarkerlist = (int *) malloc(in.numberofpoints * sizeof(int));


    int k=-1;
    for (unsigned v=0; v< surf_mesh->vertices.size();v++)
    {
      k++; in.pointlist[k] = surf_mesh->vertices[v].x;
      k++; in.pointlist[k] = surf_mesh->vertices[v].y;

      in.pointmarkerlist[v]=0;
    }

    //================================================== Parse segments
    in.numberofsegments = surf_mesh->lines.size();
    //in.numberofsegments = 0;
    in.numberofholes = 0;
    in.numberofregions = 0;
    in.segmentlist = (int *) malloc(in.numberofsegments * 2 * sizeof(int));
    k=-1;
    for (unsigned l=0; l< surf_mesh->lines.size();l++)
    {
      k++; in.segmentlist[k] = surf_mesh->lines[l].v_index[0];
      k++; in.segmentlist[k] = surf_mesh->lines[l].v_index[1];
    }

    /* Triangulate the points.  Switches are chosen to read and write a  */
    /*   PSLG (p), preserve the convex hull (c), number everything from  */
    /*   zero (z), assign a regional attribute to each element (A), and  */
    /*   produce an edge list (e), a Voronoi diagram (v), and a triangle */
    /*   neighbor list (n).                                              */

    mid.pointlist = (REAL *) NULL;
    mid.pointattributelist = (REAL *) NULL;
    mid.pointmarkerlist = (int *) NULL;

    mid.trianglelist = (int *) NULL;
    mid.triangleattributelist = (REAL *) NULL;

    mid.segmentlist = (int *) NULL;
    mid.segmentmarkerlist = (int *) NULL;

    mid.edgelist = (int *) NULL;



    //report(&in, 0, 0, 0, 0, 0, 0);
    chi_log.Log(LOG_0VERBOSE_1) << "Generating delaunay triangulation";
    triangulate("zcDpQ", &in, &mid, NULL);
    chi_log.Log(LOG_0VERBOSE_1) << "Completed delaunay triangulation";
    //report(&mid, 0, 1, 0, 1, 0, 0);


    mid.trianglearealist = (REAL *) malloc(mid.numberoftriangles * sizeof(REAL));
    for (int t=0;t<mid.numberoftriangles;t++)
    {
      mid.trianglearealist[t] = area_constraint;
    }




    out.pointlist = (REAL *) NULL;
    out.pointmarkerlist = (int *) NULL;

    out.trianglelist = (int *) NULL;
    out.triangleattributelist = (REAL *) NULL;

    out.segmentlist = (int *) NULL;
    out.segmentmarkerlist = (int *) NULL;
    chi_log.Log(LOG_0VERBOSE_1) << "Generating triangulation refinement";
    triangulate("zcarsDpQq", &mid, &out, NULL);
    chi_log.Log(LOG_0VERBOSE_1) << "Completed triangulation refinement";
    //report(&out, 0, 1, 0, 1, 0, 0);

    output = &out;

    //================================================== Create remeshed surface
    new_surface_mesh = new chi_mesh::SurfaceMesh;
    temp_surface_mesh = new chi_mesh::SurfaceMesh;

    //============================================= Create vertices
    for (int v=0; v<output->numberofpoints; v++)
    {
//    printf("Vertex %d   ", v);
//    printf("%f  ",output->pointlist[v * 2 + 0]);
//    printf("%f  ",output->pointlist[v * 2 + 1]);
//    printf("\n");

      Vertex new_vert;
      new_vert.x = output->pointlist[v * 2 + 0];
      new_vert.y = output->pointlist[v * 2 + 1];
      new_vert.z = 0.0;

      temp_surface_mesh->vertices.push_back(new_vert);
    }

    //============================================= Create triangles
    for (int t=0; t<output->numberoftriangles; t++)
    {
//    printf("Triangle %d   ", t);
//    printf("%d  ",output->trianglelist[t * 3 + 0]);
//    printf("%d  ",output->trianglelist[t * 3 + 1]);
//    printf("%d  ",output->trianglelist[t * 3 + 2]);
//    printf("\n");

      chi_mesh::Face* new_face = new chi_mesh::Face;
      new_face->SetIndices(output->trianglelist[t * 3 + 0],
                           output->trianglelist[t * 3 + 1],
                           output->trianglelist[t * 3 + 2]);
      new_face->geometric_normal = chi_mesh::Normal(0.0,0.0,1.0);

      temp_surface_mesh->faces.push_back(*new_face);
    }
    chi_log.Log(LOG_0) << "SurfaceMesherTriangle: Number of triangles="
                                        << temp_surface_mesh->faces.size();
    chi_log.Log(LOG_0) << "SurfaceMesherTriangle: Number of vertices="
                       << temp_surface_mesh->vertices.size();

    //============================================= Update connectivity
    temp_surface_mesh->UpdateInternalConnectivity();

    free(in.pointlist);
//  free(in.pointattributelist);
    free(in.pointmarkerlist);
//  free(in.regionlist);
//  free(mid.pointlist);
//  free(mid.pointattributelist);
//  free(mid.pointmarkerlist);
    free(mid.trianglelist);
    free(mid.triangleattributelist);
    free(mid.trianglearealist);
//  free(mid.neighborlist);
    free(mid.segmentlist);
    free(mid.segmentmarkerlist);
    free(mid.edgelist);
//  free(mid.edgemarkerlist);
//  free(vorout.pointlist);
//  free(vorout.pointattributelist);
//  free(vorout.edgelist);
//  free(vorout.normlist);
    free(out.pointlist);
//  free(out.pointattributelist);
    free(out.trianglelist);
    free(out.triangleattributelist);

    ReorderCells(temp_surface_mesh,
                 new_surface_mesh);
//    new_surface_mesh = temp_surface_mesh;
  }
    //##################################################### CHILD PROCESSES
  else
  {
    //================================================== Create remeshed surface
    new_surface_mesh = new chi_mesh::SurfaceMesh;
  }

  if (chi_mpi.location_id == 0)
  {
    chi_log.Log(LOG_ALLVERBOSE_2) << "Sending nodes";
    chi_mpi.BroadcastNodes(&(new_surface_mesh->vertices));
    chi_log.Log(LOG_ALLVERBOSE_2) << "Sending faces";
    chi_mpi.BroadcastTriFaces(&(new_surface_mesh->faces));
  }
  else
  {
    chi_log.Log(LOG_ALLVERBOSE_2)
      << "Process "
      << chi_mpi.location_id
      << "receiving nodes";
    chi_mpi.ReceiveNodes(&(new_surface_mesh->vertices));
    chi_log.Log(LOG_ALLVERBOSE_2)
      << "Process "
      << chi_mpi.location_id
      << "receiving faces";
    chi_mpi.ReceiveTriFaces(&(new_surface_mesh->faces));
    chi_log.Log(LOG_ALLVERBOSE_2)
      << "Process "
      << chi_mpi.location_id
      << "done receiving items";
  }



  chi_mesh::MeshContinuum* new_cont = new chi_mesh::MeshContinuum;
  new_cont->surface_mesh = new_surface_mesh;

  boundary->mesh_continua.clear();
  boundary->mesh_continua.push_back(new_cont);

  if (chi_log.GetVerbosity() > LOG_0VERBOSE_1)
  {
    std::string filename = std::string("TestP") +
                           std::to_string(chi_mpi.location_id) +
                           std::string(".obj");
    new_surface_mesh->ExportToOBJFile(filename.c_str());
  }

}