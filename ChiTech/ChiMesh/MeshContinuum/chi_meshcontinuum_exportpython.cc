#include "chi_meshcontinuum.h"
#include <cinttypes>
#include <fstream>
#include "ChiMesh/Cell/cell.h"

#include "chi_mpi.h"
#include "chi_log.h"

extern ChiMPI& chi_mpi;
extern ChiLog& chi_log;

//###################################################################
/**Export cells to python.*/
void chi_mesh::MeshContinuum::
ExportCellsToPython(const char* fileName, bool surface_only,
                    std::vector<int>* cell_flags,
                    int options) const
{
  FILE* of = fopen(fileName,"w");

  if (of==NULL)
  {
    fprintf(stdout,"Could not open file %s\n",fileName);
    chi_log.Log(LOG_ALLWARNING) << "Could not open file: "
                                << std::string(fileName);
    return;
  }

  //fprintf(of,"Template\n");

  //============================================= Import plot functions
  fprintf(of,"import matplotlib.pyplot as plt\n"
             "from mpl_toolkits.mplot3d import Axes3D\n"
             "import numpy as np\n"
             "import mpl_toolkits.mplot3d as a3\n"
             "from matplotlib import colors\n"
             "import scipy as sp\n");

  //============================================= Parse nodes
  fprintf(of, "xyz=np.zeros((%llu,%d))\n", GetGlobalVertexCount(), 3);
  for (int n=0; n < GetGlobalVertexCount(); n++)
  {
    fprintf(of, "xyz[%d][%d]=%f;  ", n, 0, vertices[n].x);
    fprintf(of, "xyz[%d][%d]=%f;  ", n, 1, vertices[n].y);
    fprintf(of, "xyz[%d][%d]=%f;\n", n, 2, vertices[n].z);
  }

  //============================================= Find amount of faces
  int num_faces=0;
  for (int lc=0; lc<local_cell_glob_indices.size(); lc++)
  {
    uint64_t cell_g_index = local_cell_glob_indices[lc];
    const auto& cell = cells[cell_g_index];

    if (cell.Type() == chi_mesh::CellType::POLYGON)
    {
      num_faces++;
    }

    if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
    {
      if (surface_only)
      {
        for (int f=0; f<cell.faces.size(); f++)
        {
          if (not cell.faces[f].has_neighbor)
          {
            num_faces++;
          }
          else
          {
            auto& adj_cell = cells[cell.faces[f].neighbor_id];
            if (adj_cell.partition_id != chi_mpi.location_id)
            {
              num_faces++;
            }
          }
        }
      }
      else
        num_faces+=cell.faces.size();

    }
  }

  chi_log.Log(LOG_ALL) << "Number of faces to be exported = "
                            << num_faces;

  //============================================= Parse faces
  fprintf(of,"face_numverts=np.zeros((%d,2))   \n",num_faces);
  fprintf(of,"face_vertindi=np.zeros((%d,10))\n",num_faces);

  int f=-1;
  for (int lc=0; lc<local_cell_glob_indices.size(); lc++)
  {
    int cell_g_index = local_cell_glob_indices[lc];
    auto& cell = cells[cell_g_index];

    //################################################### POLYGON
    if (cell.Type() == chi_mesh::CellType::POLYGON)
    {
      fprintf(of,"face_numverts[%d,0]=%lu   \n",cell_g_index,
              cell.vertex_ids.size());

      bool flagged = false;
      if (cell_flags != nullptr)
      {
        for (int i=0; i<cell_flags->size(); i++)
        {
          if (cell_g_index==(*cell_flags)[i])
          {flagged = true;}

//          if (!poly_cell->CheckBoundary2D())
//          {flagged = true;}
        }
        if (flagged)
        {
          fprintf(of,"face_numverts[%d,1]=%d   \n",cell_g_index,1);
        }
      }

      for (size_t v=0; v<cell.vertex_ids.size(); v++)
      {
        fprintf(of,"face_vertindi[%d][%lu]=%" PRIu64 "\n",cell_g_index,v,cell.vertex_ids[v]);
      }
    }//polygon

    if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
    {
      if (surface_only)
      {
        for (int s=0; s< cell.faces.size(); s++)
        {
          bool export_face = false;
          if (not cell.faces[s].has_neighbor)
          {
            export_face = true;
          }
          else
          {
            auto& adj_cell = cells[cell.faces[s].neighbor_id];
            if (adj_cell.partition_id != chi_mpi.location_id)
            {
              export_face = true;
            }
          }

          if (export_face)
          {
            f++;
            fprintf(of,"face_numverts[%d,0]=%lu   \n",f,
                    cell.faces[s].vertex_ids.size());

            bool flagged = false;
            if (cell_flags != nullptr)
            {
              for (int i=0; i<cell_flags->size(); i++)
              {
                if (cell_g_index==(*cell_flags)[i])
                {flagged = true;}
              }
              if (flagged)
              {
                fprintf(of,"face_numverts[%d,1]=%d   \n",f,1);
              }
            }


            const chi_mesh::CellFace& face = cell.faces[s];

            for (size_t v=0; v<face.vertex_ids.size(); v++)
            {
              fprintf(of,"face_vertindi[%d][%lu]=%" PRIu64 "\n",f,v,face.vertex_ids[v]);
            }
          }


        }
      }
      else
      {
        for (int s=0; s< cell.faces.size(); s++)
        {
          f++;
          fprintf(of,"face_numverts[%d,0]=%lu   \n",f,
                  cell.faces[s].vertex_ids.size());

          bool flagged = false;
          if (cell_flags != nullptr)
          {
            for (int i=0; i<cell_flags->size(); i++)
            {
              if (cell_g_index==(*cell_flags)[i])
              {flagged = true;}
            }
            if (flagged)
            {
              fprintf(of,"face_numverts[%d,1]=%d   \n",f,1);
            }
          }

          const chi_mesh::CellFace& face = cell.faces[s];

          for (size_t v=0; v<face.vertex_ids.size(); v++)
          {
            fprintf(of,"face_vertindi[%d][%lu]=%" PRIu64 "\n",f,v,face.vertex_ids[v]);
          }

        }
      }
    }//polyhedron
  }

  //============================================= Initialize plot
  fprintf(of,"fig = plt.figure(1)\n"
             "ax = a3.Axes3D(plt.figure(1))\n");

  //============================================= Plot each face
  fprintf(of,"num_faces=%d\n",num_faces);
  fprintf(of,"\n"
             "for f in range(0,num_faces):\n"
             "    vcount = int(face_numverts[f,0])\n"
             "    verts=np.zeros((vcount+1,3))\n"
             "    \n"
             "    for v in range(0,vcount):\n"
             "        index = int(face_vertindi[f][v])\n"
             "        verts[v][0] = xyz[index][0]\n"
             "        verts[v][1] = xyz[index][1]\n"
             "        verts[v][2] = xyz[index][2]\n"
             "    \n"
             "    index = int(face_vertindi[f][0])\n"
             "    verts[vcount][0] = xyz[index][0]\n"
             "    verts[vcount][1] = xyz[index][1]\n"
             "    verts[vcount][2] = xyz[index][2]\n"
             "    \n"
             "    #ax.plot(verts[:,0],verts[:,1],verts[:,2],'k-',linewidth=0.5)\n"
             "    pol = a3.art3d.Poly3DCollection([verts],linewidths=0.5)\n"
             "    pol.set_edgecolor('k')\n"
             "    pol.set_facecolor(colors.rgb2hex([0.8,0.8,0.8]))\n"
             "    if (int(face_numverts[f,1])==1):\n"
             "        pol.set_facecolor(colors.rgb2hex([1.0,0.0,0.0]))\n"
             "    ax.add_collection3d(pol)\n"
             "\n"
             "xmax = np.max(xyz[:,0])\n"
             "xmin = np.min(xyz[:,0])\n"
             "ymax = np.max(xyz[:,1])\n"
             "ymin = np.min(xyz[:,1])\n"
             "zmax = np.max(xyz[:,2])\n"
             "zmin = np.min(xyz[:,2])\n"
             "\n"
             "ax.set_xlim([xmin, xmax])\n"
             "ax.set_ylim([ymin, ymax])\n"
             "ax.set_zlim([zmin, zmax])\n"
             "\n"
             "ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))\n"
             "ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))\n"
             "ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))\n"
             "ax.view_init(elev=90.0,azim=0.0)\n");

  if (options == 1)
  {
    fprintf(of,"plt.savefig('");
    fprintf(of,"%s",fileName);
    fprintf(of,".png')\n");
  }
  else
  {
    fprintf(of,"plt.show()\n");
  }


  fclose(of);

}