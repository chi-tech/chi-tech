#include "delaunay_mesher.h"
#include "../../SurfaceMesh/chi_surfacemesh.h"

/**Meshes a region-context
 *
*      - For each surface mesh
 *          - Split into into patches (Co-planar and contigous)
 *          - For each patch
 *              - Perform Lexicographical meshing
 *                  - Get edge loops and essential vertices
 *                  - If not first patch
 *                      - Add vertices from in-region interfaces
 *                  - Project to 2D
 *
 *              */
void chi_mesh::SurfaceMesherDelaunay::MeshRegion(
  RegionDelaunyContexts& region_context)
{
  //================================================== Loop over all sub contexts
  std::vector<DelaunayMeshContext*>::iterator delaunay_context;
  for (delaunay_context = region_context.begin();
       delaunay_context != region_context.end();
       delaunay_context++)
  {
    if ((*delaunay_context)->context_surface_mesh!= nullptr)
    {
      //========================================= Split by patch
      (*delaunay_context)->context_surface_mesh->SplitByPatch((*delaunay_context)->patches);
      printf("Number of surfaces =%d\n",(*delaunay_context)->patches.size());

      //========================================= Mesh each patch
      std::vector<chi_mesh::SurfaceMesh*>::iterator patch;
      for (patch = (*delaunay_context)->patches.begin();
           patch != (*delaunay_context)->patches.end();
           patch++)
      {
        DelaunayPatch* new_del_patch = new DelaunayPatch(tolerance);
        //MeshLexically((*patch),*new_del_patch);
        MeshDelaunay((*patch),*new_del_patch);
        (*delaunay_context)->delaunay_patches.push_back(new_del_patch);

        //================================== Apply interfaces interior to region

      }
    }
  }//for delaunay contexts

  //================================================== Create surface for each context
  for (unsigned c=0; c<region_context.size(); c++)
  {
    if (region_context[c]->context_surface_mesh!= nullptr)
    {
      chi_mesh::SurfaceMesh* new_surf = new chi_mesh::SurfaceMesh;
      region_context[c]->remeshed_surface_mesh = new_surf;

      for (unsigned p=0; p<region_context[c]->delaunay_patches.size();p++)
      {
        DelaunayPatch* curPatch = (DelaunayPatch*)
          region_context[c]->delaunay_patches[p];

        //======================================= Create vertex add history
        struct vhist
        {
          int original_index;
          int new_index;
        };
        std::vector<vhist*> vertex_add_hist;
        printf("Number of triangles in patch %d\n",curPatch->triangles.size());

        //======================================= Loop over all triangles
        for (unsigned t=0; t<curPatch->triangles.size(); t++)
        {
          if (!(curPatch->triangles[t].invalidated))
          {
            //===================================== Copy vertices
            for (int e=0; e<3; e++)
            {
              //========================= Check if vertex is already added
              bool vertex_added_already=false;
              for (unsigned vh=0; vh<vertex_add_hist.size(); vh++)
              {
                if (curPatch->triangles[t].v_index[e] ==
                    vertex_add_hist[vh]->original_index)
                {
                  vertex_added_already = true; break;
                }
              }

              if (!vertex_added_already)
              {
                int orig_index = curPatch->triangles[t].v_index[e];
                new_surf->vertices.push_back(curPatch->vertices[orig_index]);

                vhist* new_hist = new vhist;
                new_hist->original_index = orig_index;
                new_hist->new_index = new_surf->vertices.size()-1;
                vertex_add_hist.push_back(new_hist);

//                printf("Vertex[%3d] copied:",orig_index);
//                curPatch->vertices[orig_index].Print();
//                printf(", new index=%d\n",new_hist->new_index);
              }

            }

            //===================================== Copy Triangle to face
            Tri* cur_tri = &curPatch->triangles[t];
            Face* new_face = new Face;
            int v_orig[3];
            v_orig[0] = cur_tri->v_index[0];
            v_orig[1] = cur_tri->v_index[1];
            v_orig[2] = cur_tri->v_index[2];

            int v_new[3]={-1,-1,-1};
            bool done[]={false,false,false};
            for (unsigned vh=0; vh<vertex_add_hist.size(); vh++)
            {
              if (vertex_add_hist[vh]->original_index == v_orig[0])
              {v_new[0] = vertex_add_hist[vh]->new_index; done[0]=true;}

              if (vertex_add_hist[vh]->original_index == v_orig[1])
              {v_new[1] = vertex_add_hist[vh]->new_index; done[1]=true;}

              if (vertex_add_hist[vh]->original_index == v_orig[2])
              {v_new[2] = vertex_add_hist[vh]->new_index; done[2]=true;}


              if (done[0] && done[1] && done[2]) break;
            }

            new_face->SetIndices(v_new[0],v_new[1],v_new[2]);
            new_surf->faces.push_back(*new_face);

          }//if valid triangle

        }//for triangle
        vertex_add_hist.clear();
      }//for patch
      //===================================== Re-establish connections
      //new_surf->UpdateInternalConnectivity();
    }//if surfacemesh

  }//for boundary context
}