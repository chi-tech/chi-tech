#include "delaunay_mesher.h"
#include "ChiMath/chi_math.h"

extern ChiMath chi_math_handler;

//###################################################################
/**Determine the orientation of point c relative to point a and b.
 *
 * Synonomous with
\f[
	(ab {\times} bc) {\bullet} [0 0 1]
\f]

but gives an indication of where $c$ lays relative to the line $ab$.

\return >0 if $c$ is left of $ab$, =0 if $abc$ are co-linear and <0 if $c$ is
 right of $ab$.

\image html Orient2D.png

\author Jan*/
double chi_mesh::SurfaceMesherDelaunay::DelaunayPatch::Orient2D(
  chi_mesh::Vertex a, chi_mesh::Vertex b, chi_mesh::Vertex c)
{
  Vector ab = b-a;
  Vector bc = c-b;
  ab = ab/ab.Norm();
  bc = bc/bc.Norm();
  Vector khat(0.0,0.0,1.0);

  Vector ab_x_bc = ab.Cross(bc);

  return ab_x_bc.Dot(khat);
}


//###################################################################
/** Checks whether the line, given by the two vertices v0 and v1,
 * any of the simplices defined.
 *
 * This method calculates the intersection of two lines

\f[
 v = v_{0} + t_v \hat{n}_v
\f]
and
\f[
 w = w_0 + t_w \hat{n}_w
\f]

 * */
bool chi_mesh::SurfaceMesherDelaunay::DelaunayPatch::CheckCrossSimplices(
  chi_mesh::Vertex v0, chi_mesh::Vertex v1)
{
  ChiMath& math = chi_math_handler;
  Vector v01 = v1 - v0;
  Vector n_v = v01;

//  printf("Line from (%+.2f,%+.2f)  ",v0.x,v0.y);
//  printf("to (%+.2f,%+.2f)\n",v1.x,v1.y);


  EdgeLoopCollection::iterator cursimplexcoll;
  for (cursimplexcoll = simplices.begin();
       cursimplexcoll != simplices.end();
       cursimplexcoll++)
  {
    EdgeList::iterator cursimplex;
    for (cursimplex = (*cursimplexcoll)->edges.begin();
         cursimplex != (*cursimplexcoll)->edges.end();
         cursimplex++)
    {
      //========================================= Get unprojected vertices
      Vertex w0_up = (*cursimplex).vertices[0] - centroid;
      Vertex w1_up = (*cursimplex).vertices[1] - centroid;

      //========================================= Project the vertices to 2D
      Vertex w0(w0_up.Dot(hat_i), w0_up.Dot(hat_j), 0.0);
      Vertex w1(w1_up.Dot(hat_i), w1_up.Dot(hat_j), 0.0);

      Vector w01 = w1 - w0;
      Vector n_w = w01;

      if (fabs((w01.Cross(v01)).z)>tolerance)
      {
        VecDbl rhs(2);
        rhs[0] = w0.x-v0.x;
        rhs[1] = w0.y-v0.y;

        MatDbl A(2,VecDbl(2,0.0));
        A[0][0] = n_v.x; A[0][1] = -n_w.x;
        A[1][0] = n_v.y; A[1][1] = -n_w.y;

        math.GaussElimination(A,rhs,2);
        VecDbl t = rhs;

//        printf("Simplex from (%+.2f,%+.2f)  ",w0.x,w0.y);
//        printf("to (%+.2f,%+.2f)  t(%+.4f,%+.4f\n",w1.x,w1.y,t(0),t(1));


        if (   (t[0]>(0.0+1.0e-4)) && (t[0]<(1.0-1.0e-4)) &&
               (t[1]>(0.0+1.0e-4)) && (t[1]<(1.0-1.0e-4)) )
        {
          return true;
        }
      }
    }
  }

  return false;
}

//###################################################################
/**Determines if point d is within the circumcircle of triangle abc.*/
double chi_mesh::SurfaceMesherDelaunay::
  DelaunayPatch::InCircle(chi_mesh::Vertex a,
                          chi_mesh::Vertex b,
                          chi_mesh::Vertex c,
                          chi_mesh::Vertex d)
{
  MatDbl m(4,VecDbl(4,0.0));
  m[0][0] = a.x; m[0][1] = a.y; m[0][2] = a.x*a.x + a.y*a.y; m[0][3] = 1.0;
  m[1][0] = b.x; m[1][1] = b.y; m[1][2] = b.x*b.x + b.y*b.y; m[1][3] = 1.0;
  m[2][0] = c.x; m[2][1] = c.y; m[2][2] = c.x*c.x + c.y*c.y; m[2][3] = 1.0;
  m[3][0] = d.x; m[3][1] = d.y; m[3][2] = d.x*d.x + d.y*d.y; m[3][3] = 1.0;

  return chi_math_handler.Determinant(m);
}


//###################################################################
/**Accumulates non-locally Delaunay edges.*/
void chi_mesh::SurfaceMesherDelaunay::
  DelaunayPatch::FindNonLocDelEdge(std::vector<Edge>& non_loc_del_edges)
{
  //======================================================= Run through each triangle
  for (unsigned t1=0; t1 < triangles.size(); t1++)
  //for (int i=0;i<in_triangles.size();i++)
  {
    //unsigned t1 = in_triangles.at(i);
    //printf("Checking tri %d\n",t1);
    Tri* tau_1 = &triangles[t1];
    //================================================ For each edge of the tri
    for (int e=0; e<3; e++)
    {
      if (tau_1->e_index[e][2]>-1) //Check if its an interior edge
      {
        //======================================= Assign neighboring triangle
        unsigned t2 = tau_1->e_index[e][2];
        Tri* tau_2 = &triangles[t2]; //This has to be on the master stack

//        printf("Tau1[%d] %d->%d->%d\n",t1,tau_1->v_index[0],
//                                          tau_1->v_index[1],
//                                          tau_1->v_index[2]);
//        printf("Tau2[%d] %d->%d->%d\n",t2,tau_2->v_index[0],
//                                          tau_2->v_index[1],
//                                          tau_2->v_index[2]);

        //=============================== Find indexes for abc
        int aindex = tau_1->v_index[0];
        int bindex = tau_1->v_index[1];
        int cindex = tau_1->v_index[2];
        int dindex = 0;

        //=============================== Find the index of d (the one that is not a b or c)
        for (int v=0;v<3;v++)
        {
          if ((tau_2->v_index[v]!=tau_1->e_index[e][0]) &&
              (tau_2->v_index[v]!=tau_1->e_index[e][1]))
          {
            dindex = tau_2->v_index[v]; break;
          }
        }

        //=============================== Get the vertex values of abcd
        Vertex a = Pstar[aindex]; //a.Print(); std::cout<<std::endl;
        Vertex b = Pstar[bindex]; //b.Print(); std::cout<<std::endl;
        Vertex c = Pstar[cindex]; //c.Print(); std::cout<<std::endl;
        Vertex d = Pstar[dindex]; //d.Print(); std::cout<<std::endl;


        //=============================== Check if locally delaunay
        double localDelaunay = InCircle(a,b,c,d);
        //printf("%e\n",localDelaunay);std::cout<<std::endl;
        if (localDelaunay>0.0001)
        {
          //====================== Check for duplicate
          bool duplicateFound = false;
          for (unsigned ie=0; ie<non_loc_del_edges.size(); ie++)
          {
            Edge* curEdge = &non_loc_del_edges[ie];
            /*printf("Dup check %d->%d vs %d->%d\n",
                                  curEdge->vi,curEdge->vf,
                                  tau_1->edges[e][1],tau_1->edges[e][0]);*/
            if ( ((curEdge->v_index[0]==tau_1->e_index[e][1]) &&
                  (curEdge->v_index[1]==tau_1->e_index[e][0])) ||
                 ((curEdge->v_index[1]==tau_1->e_index[e][1]) &&
                  (curEdge->v_index[0]==tau_1->e_index[e][0])) )
            {
              duplicateFound = true; break;
            }
          }

          //====================== If no duplicate found then add edge to list
          if (!duplicateFound)
          {
            Edge* newEdge = new Edge;
            newEdge->v_index[0] = tau_1->e_index[e][0];
            newEdge->v_index[1] = tau_1->e_index[e][1];

            newEdge->f_index[0]=t1;
            newEdge->f_index[1]=e;

            newEdge->f_index[2]=tau_1->e_index[e][2];

            //Find which edge of Tau2
            for (int e2=0;e2<3; e2++)
            {
              if ((tau_2->e_index[e2][0]==newEdge->v_index[1]) &&
                  (tau_2->e_index[e2][1]==newEdge->v_index[0]))
              {
                newEdge->f_index[3] = e2; break;
              }
            }
            if (newEdge->f_index[3]<0)
            {
              printf("Error finding tau_2_edgenumber %d\n", t2);
            }
            //printf("Non-delaunay interior edge from vertex %3d to %3d, tau1=%3d(%3d), tau2=%3d(%3d) %.9f\n",
            //        newEdge->vi, newEdge->vf,
            //        newEdge->tau_1, newEdge->tau_1_edgeNumber,
            //        newEdge->tau_2, newEdge->tau_2_edgeNumber,localDelaunay);
            non_loc_del_edges.push_back(*newEdge);
            return;
          }

        }

      }//if interior edge
    }//for e
  }//for t1
}



//###################################################################
/** Dumps the current triangle collection to scilab format.
 */
void chi_mesh::SurfaceMesherDelaunay::
  DelaunayPatch::DumpToScilab(const char *file_name, bool verbose)
{
  if (this->triangles.size()==0)
  {
    std::cout << "Cannot export empty SurfaceMesh to scilab format\n";
    return;
  }
  FILE* of = fopen(file_name,"w");
  if (of==NULL)
  {
    printf("Error creating file %s!\n",file_name);
    return;
  }

  fprintf(of,"clear\n");
  fprintf(of,"clc\n");

  fprintf(of,"verts=zeros(%lu,2);\n",Pstar.size());

  for (int v=0; v<Pstar.size(); v++)
  {
    fprintf(of,"verts(%3d,1)=%+.5f; verts(%3d,2)=%+.5f;\n",
            v+1,Pstar[v].x,
            v+1,Pstar[v].y);
  }

  fprintf(of,"\n");
  fprintf(of,"tris=zeros(%lu,6);\n",triangles.size());

  for (int f=0; f<triangles.size(); f++)
  {
    if (!triangles[f].invalidated)
    {
      fprintf(of,"tris(%3d,1)=%3d;  ",f+1,triangles[f].v_index[0]);
      fprintf(of,"tris(%3d,2)=%3d;  ",f+1,triangles[f].v_index[1]);
      fprintf(of,"tris(%3d,3)=%3d;\n",f+1,triangles[f].v_index[2]);
      fprintf(of,"tris(%3d,4)=%3d;  ",f+1,triangles[f].e_index[0][2]);
      fprintf(of,"tris(%3d,5)=%3d;  ",f+1,triangles[f].e_index[1][2]);
      fprintf(of,"tris(%3d,6)=%3d;\n",f+1,triangles[f].e_index[2][2]);
    }
  }

  if (verbose)
  {
    fprintf(of,"\nscf(0)\n"
               "clf(0)\n"
               "\n"
               "scatter(verts(:,1),verts(:,2));\n\n");

    fprintf(of,"[v_count,cols]=size(verts)\n"
               "for v=1:v_count do\n"
               "    dx=-0.025\n"
               "    dy=-0.005\n"
               "    xstring(verts(v,1)+dx,verts(v,2)+dy,string(v-1))\n"
               "    id=color(\"red\")\n"
               "    ta=get(\"hdl\")\n"
               "    ta.font_foreground=id\n"
               "end\n\n");

    fprintf(of,"[tri_count,cols]=size(tris)\n"
               "for t=1:tri_count do\n"
               "    v0 = verts(tris(t,1)+1,:);\n"
               "    v1 = verts(tris(t,2)+1,:);\n"
               "    v2 = verts(tris(t,3)+1,:);\n"
               "    points=[v0; v1; v2; v0]\n"
               "    plot2d(points(:,1),points(:,2))\n"
               "    \n"
               "    xmean = mean(points(1:3,1));\n"
               "    ymean = mean(points(1:3,2));\n"
               "    e0mean=[0.5*(v0(1,1)+v1(1,1)) 0.5*(v0(1,2)+v1(1,2)) ]\n"
               "    e1mean=[0.5*(v1(1,1)+v2(1,1)) 0.5*(v1(1,2)+v2(1,2)) ]\n"
               "    e2mean=[0.5*(v2(1,1)+v0(1,1)) 0.5*(v2(1,2)+v0(1,2)) ]\n"
               "    \n"
               "    dx=-0.005\n"
               "    dy=-0.015\n"
               "    xstring(xmean+dx,ymean+dy,string(t-1))\n"
               "    \n"
               "    \n"
               "    d=0.4\n"
               "    \n"
               "    dx=d*(e0mean(1)-xmean)-0.005\n"
               "    dy=d*(e0mean(2)-ymean)-0.015\n"
               "    xstring(xmean+dx,ymean+dy,string(tris(t,4)))\n"
               "    id=color(\"gray\")\n"
               "    ta=get(\"hdl\")\n"
               "    ta.font_foreground=id\n"
               "    \n"
               "    dx=d*(e1mean(1)-xmean)-0.005\n"
               "    dy=d*(e1mean(2)-ymean)-0.015\n"
               "    xstring(xmean+dx,ymean+dy,string(tris(t,5)))\n"
               "    id=color(\"gray\")\n"
               "    ta=get(\"hdl\")\n"
               "    ta.font_foreground=id\n"
               "\n"
               "    dx=d*(e2mean(1)-xmean)-0.005\n"
               "    dy=d*(e2mean(2)-ymean)-0.015\n"
               "    xstring(xmean+dx,ymean+dy,string(tris(t,6)))\n"
               "    id=color(\"gray\")\n"
               "    ta=get(\"hdl\")\n"
               "    ta.font_foreground=id\n"
               "end\n"
               "    \n"
               "a=gca();\n");
  } else
  {
    fprintf(of,"\nscf(0)\n"
               "clf(0)\n\n");

    fprintf(of,"[tri_count,cols]=size(tris)\n"
               "for t=1:tri_count do\n"
               "    v0 = verts(tris(t,1)+1,:);\n"
               "    v1 = verts(tris(t,2)+1,:);\n"
               "    v2 = verts(tris(t,3)+1,:);\n"
               "    points=[v0; v1; v2; v0]\n"
               "    plot2d(points(:,1),points(:,2))\n"
               "    \n"
               "end\n"
               "    \n"
               "a=gca();\n");
  }


  fclose(of);
  printf("Exported mesh to %s\n",file_name);
}


//###################################################################
/** Checks if an edge can be split.*/
bool chi_mesh::SurfaceMesherDelaunay::DelaunayPatch::
  CheckEdgeCanBeSplit(unsigned i_vindex,unsigned f_vindex)
{
  Vertex a = Pstar[i_vindex];
  Vertex b = Pstar[f_vindex];

  Vector ab = b - a;

  if ((ab.Norm()/2.0)>absoluteMinumumSize)
  {
    return true;
  } else
  {
    return false;
  }
}

//###################################################################
/** Checks if an edge is encroached.
\param i_vindex int Index of the vertex where the edge start.
\param f_vindex int Index of the vertex where the edge ends.
\param other_vindex int Index of the other vertex of the triangle.

\return true if the edge is encroached.*/
bool chi_mesh::SurfaceMesherDelaunay::DelaunayPatch::
  CheckEdgeEncroached(unsigned i_vindex,
                      unsigned f_vindex,
                      unsigned other_vindex)
{
  Vertex v0 = Pstar[i_vindex];
  Vertex v1 = Pstar[f_vindex];
  Vertex v2 = Pstar[other_vindex];

  Vertex vc = v0*0.5 + v1*0.5;

  Vector vc2i = vc - v0;
  Vector vc2o = vc - v2;

  double r = vc2i.Norm();
  double d = vc2o.Norm();

  //printf("r,d=%f,%f for edge %d->%d#%d\n",r,d,i_vindex,f_vindex,other_vindex);

  if ((d<(r+0.0000001)) )
  {

    return true;
  }

  return false;
}

bool chi_mesh::SurfaceMesherDelaunay::DelaunayPatch::
CheckEdgeEncroached(unsigned i_vindex,
                    unsigned f_vindex,
                    Vertex v2)
{
  Vertex v0 = Pstar[i_vindex];
  Vertex v1 = Pstar[f_vindex];

  Vertex vc = v0*0.5 + v1*0.5;

  Vector vc2i = vc - v0;
  Vector vc2o = vc - v2;

  double r = vc2i.Norm();
  double d = vc2o.Norm();

  //printf("r,d=%f,%f for edge %d->%d#%d\n",r,d,i_vindex,f_vindex,-1);

  if ((d<(r+0.0000001)) )
  {

    return true;
  }

  return false;
}

//###################################################################
/**Finds the circum circle (center and radius) of a triangle.*/
void chi_mesh::SurfaceMesherDelaunay::DelaunayPatch::
FindCircumCircle(chi_mesh::Vertex a, chi_mesh::Vertex b, chi_mesh::Vertex c,
  chi_mesh::Vertex &center,
  double &radius)
{
  bool line1Vertical=false;
  bool line2Vertical=false;

  //===================================================== Slope 1
  float m1 = 0.0;
  if (fabs(b.x-a.x)>tolerance)
  {
    float mBA = (b.y-a.y)/(b.x-a.x);
    if (fabs(mBA)>tolerance)
    {
      m1  = -1.0/mBA;
    } else
    {
      line1Vertical = true;
    }

  }

  float x1 = (b.x+a.x)/2.0;
  float y1 = (b.y+a.y)/2.0;

  //y = m*x +c
  //c1 = y1-m1*x1
  //yc = m1*xc +y1-m1*x1

  //===================================================== Slope 2
  float m2 = 0.0;
  if (fabs(c.x-b.x)>tolerance)
  {
    float mCB = (c.y-b.y)/(c.x-b.x);
    if (fabs(mCB)>tolerance)
    {
      m2  = -1.0/mCB;
    }
    else
    {
      line2Vertical = true;
    }
  }

  float x2 = (c.x+b.x)/2.0;
  float y2 = (c.y+b.y)/2.0;

  //y = m*x +c
  //c2 = y2-m2*x2
  //yc = m2*xc +y2-m2*x2

  //m1*xc +y1-m1*x1 = m2*xc +y2-m2*x2
  //m1*xc -m2*xc  =  m1*x1-y1 - (m2*x2 - y2)

  if (line1Vertical)
  {
    center.x = x1;
    center.y = m2*(center.x-x2)+y2;
    center.z = 0.0;
  }
  else if (line2Vertical)
  {
    center.x = x2;
    center.y = m1*(center.x-x1)+y1;
    center.z = 0.0;
  }
  else
  {
    center.x = (  (m1*x1 - y1) - (m2*x2-y2)  )/(m1-m2);
    center.y = m2*(center.x-x2)+y2;
    center.z = 0.0;
  }


  radius = sqrt(    (center.x-b.x)*(center.x-b.x)  + (center.y-b.y)*(center.y-b.y)  );
}
