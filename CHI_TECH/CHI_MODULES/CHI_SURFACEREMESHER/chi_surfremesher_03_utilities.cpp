#include"chi_surfaceremesher.h"


//################################################################### Get vertex
/**Gets the actual xyz value of a vertex from an index from the given mesh.

\author Jan*/
bool CHI_SURFACEREMESHER::GetVertexFromIndex(CHI_SURFACE* mesh,int index, float* vertex)
{
	float* vertexOnStack = mesh->vertexStack.GetItem(index);

	if (vertexOnStack==NULL)
	{
		return false;
	}
	else
	{
		vertex[0] = vertexOnStack[0];
		vertex[1] = vertexOnStack[1];
		vertex[2] = vertexOnStack[2];
	}


	return true;
}


//################################################################### Get Face Vertices
/**Populates a structure with non-pointer values of the vertices for a face.
 *
 \author Jan*/
void CHI_SURFACEREMESHER::GetFaceVertices(CHI_SURFACE* mesh,CST_FACE* face,CHI_TRIPLET_VERTEX* vertices)
{
    GLfloat* rawVertex;
	
    rawVertex = mesh->vertexStack.GetItem(face->vertex[0]);
    vertices->v0(0) = rawVertex[0];
    vertices->v0(1) = rawVertex[1];
    vertices->v0(2) = rawVertex[2];

    rawVertex = mesh->vertexStack.GetItem(face->vertex[1]);
    vertices->v1(0) = rawVertex[0];
    vertices->v1(1) = rawVertex[1];
    vertices->v1(2) = rawVertex[2];

    rawVertex = mesh->vertexStack.GetItem(face->vertex[2]);
    vertices->v2(0) = rawVertex[0];
    vertices->v2(1) = rawVertex[1];
    vertices->v2(2) = rawVertex[2];
	
    vertices->i0 = face->vertex[0];
    vertices->i1 = face->vertex[1];
    vertices->i2 = face->vertex[2];

    return;
}

//################################################################### Get Edge Connection
/**Determines if the edge, given by vertex indices, is connected to
 * a face (which also only hold vertex index information).
 *
 \author Jan*/
bool CHI_SURFACEREMESHER::CheckEdgeConnection(int* edge,CST_FACE* face)
{

    bool found1 = false;
    bool found2 = false;

    for (int k=0;k<3;k++)
    {
        if (edge[0]==face->vertex[k]) {found1=true;}
        if (edge[1]==face->vertex[k]) {found2=true;}
    }

    if (found1 && found2) {return true;}

    return false;
}



//################################################################### Check edges coliniar
/**Determines if two edges are colinear by comparing their slopes.

\param edgeA The first edge
\param edgeB The second edge
\param mesh This function actually requires information about the vertex x,y,z values
 * and therefore requires a reference to which mesh the vertex indices refer.
 *
 \author Jan*/
bool CHI_SURFACEREMESHER::CheckEdgesColinear(int* edgeA, int* edgeB, CHI_SURFACE* mesh)
{
    int lineA_v0_index = edgeA[0];
    int lineA_v1_index = edgeA[1];
    
    int lineB_v0_index = edgeB[0];
    int lineB_v1_index = edgeB[1];
    
    float* vertexA0 = mesh->vertexStack.GetItem(lineA_v0_index);
	float* vertexA1 = mesh->vertexStack.GetItem(lineA_v1_index);
	float* vertexB0 = mesh->vertexStack.GetItem(lineB_v0_index);
	float* vertexB1 = mesh->vertexStack.GetItem(lineB_v1_index);
 
	if ((vertexA0==NULL) || (vertexA1==NULL) ||
	    (vertexB0==NULL) || (vertexB1==NULL))
	{
		printf("Error during CheckEdgesColinear\n");
		return false;
	}
	
	Eigen::Vector3f vertA0(vertexA0[0],vertexA0[1],vertexA0[2]);
	Eigen::Vector3f vertA1(vertexA1[0],vertexA1[1],vertexA1[2]);
	Eigen::Vector3f vertB0(vertexB0[0],vertexB0[1],vertexB0[2]);
	Eigen::Vector3f vertB1(vertexB1[0],vertexB1[1],vertexB1[2]);
	
	Eigen::Vector3f vecA01 = vertA1 - vertA0;
	Eigen::Vector3f vecB01 = vertB1 - vertB0;
	Eigen::Vector3f vecAB0 = vertB0 - vertA0;
	
	vecA01.normalize();
	vecB01.normalize();
	vecAB0.normalize();
	
	float dotProductSlopes = fabs(vecA01.dot(vecB01));
	float dotProductInters = fabs(vecA01.dot(vecAB0));
	
	if ((dotProductSlopes>(1.0-this->precision)) && (dotProductInters>(1.0-this->precision)))
	{
		return true;
	}
    return false;
}



//################################################################### Resolve edge extent
/**Given a collection of colinear edges, determines the outer two vertices of
 * the collection.
 *
\param colinearEdges The collection for which the extent needs to be determined.
\param mesh          The mesh to which the indices point.
 *
 \author Jan*/
void CHI_SURFACEREMESHER::ResolveEdgeExtent(CHI_VECTOR<CHI_COLINEDGELIST>* colinearEdges,CHI_SURFACE* mesh)
{
	for (int colin=0;colin<colinearEdges->itemCount;colin++)
	{
		CHI_COLINEDGELIST* curList = colinearEdges->GetItem(colin);
		
		curList->extentLength = 0.0;
		
		for (int e=0;e<curList->edges.itemCount;e++)
		{
			int* mastrEdge = curList->edges.GetItem(e);
			for (int e2=0;e2<curList->edges.itemCount;e2++)
			{
				int* slaveEdge = curList->edges.GetItem(e2);
				float* vertexA0 = mesh->vertexStack.GetItem(mastrEdge[0]);
				float* vertexA1 = mesh->vertexStack.GetItem(mastrEdge[1]);
				float* vertexB0 = mesh->vertexStack.GetItem(slaveEdge[0]);
				float* vertexB1 = mesh->vertexStack.GetItem(slaveEdge[1]);
				
				Eigen::Vector3f vertA0(vertexA0[0],vertexA0[1],vertexA0[2]);
				Eigen::Vector3f vertA1(vertexA1[0],vertexA1[1],vertexA1[2]);
				Eigen::Vector3f vertB0(vertexB0[0],vertexB0[1],vertexB0[2]);
				Eigen::Vector3f vertB1(vertexB1[0],vertexB1[1],vertexB1[2]);
				
				Eigen::Vector3f combination1;
				Eigen::Vector3f combination2;
				Eigen::Vector3f combination3;
				float combination1Size;
				float combination2Size;
				float combination3Size;
				
				combination1 = vertA0-vertA1; combination1Size=combination1.norm();
				combination2 = vertA0-vertB0; combination2Size=combination2.norm();
				combination3 = vertA0-vertB1; combination3Size=combination3.norm();
				
				if (combination1Size>curList->extentLength)
				{
					curList->extentLength=combination1Size;
					curList->extentVertices[0]=mastrEdge[0];
					curList->extentVertices[1]=mastrEdge[1];
				}
				if (combination2Size>curList->extentLength)
				{
					curList->extentLength=combination2Size;
					curList->extentVertices[0]=mastrEdge[0];
					curList->extentVertices[1]=slaveEdge[0];
				}
				if (combination3Size>curList->extentLength)
				{
					curList->extentLength=combination3Size;
					curList->extentVertices[0]=mastrEdge[0];
					curList->extentVertices[1]=slaveEdge[1];
				}
				
				combination1 = vertA1-vertB0; combination1Size=combination1.norm();
				combination2 = vertA1-vertB1; combination2Size=combination2.norm();

				
				if (combination1Size>curList->extentLength)
				{
					curList->extentLength=combination1Size;
					curList->extentVertices[0]=mastrEdge[1];
					curList->extentVertices[1]=slaveEdge[0];
				}
				if (combination2Size>curList->extentLength)
				{
					curList->extentLength=combination2Size;
					curList->extentVertices[0]=mastrEdge[1];
					curList->extentVertices[1]=slaveEdge[1];
				}
				
				
			}
			
		}
	}
	return;
}




//############################################################################# Find open edges
/** Find open edges of a patch.
 *
\author Jan*/
void CHI_SURFACEREMESHER::FindPatchOpenEdges(CHI_PATCH *patch, CHI_VECTOR<CHI_OPENEDGE> *openEdges)
{
	openEdges->ClearVector();
	
	//===================================================== Create temp list of all open edges
	CHI_VECTOR<CHI_OPENEDGE> temp_openEdges;
	for (int t=0; t<patch->faceList.itemCount; t++)
	{

		CST_FACE* curFace = patch->faceList.GetItem(t);
		if (curFace!=NULL)
		{

			for (int e=0; e<3;e++)
			{
				if (curFace->edges[e][2]<0)
				{
					CHI_OPENEDGE* newOpenEdge = new CHI_OPENEDGE;


					newOpenEdge->vi = curFace->edges[e][0];
					newOpenEdge->vf = curFace->edges[e][1];
					newOpenEdge->ownerTriangle = t;
					newOpenEdge->ownerEdgeNumber = e;

					int newInd = temp_openEdges.PushItem(newOpenEdge);
					//printf("Temp Open-edge %3d, from %3d to %3d at t=%3d(%3d)\n",newInd,newOpenEdge->vi,newOpenEdge->vf, newOpenEdge->ownerTriangle, newOpenEdge->ownerEdgeNumber);
				}
			}
		}

	}
	
	//printf("New open edge loop:\n");
	//===================================================== Seed the open edge list
	CHI_OPENEDGE* seedEdge = temp_openEdges.GetItem(0);
	openEdges->PushItem(seedEdge);
	int lastVertex = seedEdge->vf;
	//printf("   From %3d to %3d\n", seedEdge->vi,seedEdge->vf);
	
	bool forceStop = false;
	while ((openEdges->itemCount < temp_openEdges.itemCount) && (!forceStop))
	{
		bool connectionFound = false;
		for (int e=0; e<temp_openEdges.itemCount; e++)
		{
			CHI_OPENEDGE* curEdge = temp_openEdges.GetItem(e);
			if (curEdge->vi==lastVertex)
			{

				openEdges->PushItem(curEdge);
				lastVertex=curEdge->vf;
				connectionFound=true;
				//printf("   From %3d to %3d\n", curEdge->vi,curEdge->vf);
			}
		}
		
		if (!connectionFound)
		{
			if (openEdges->itemCount < temp_openEdges.itemCount)
			{
				printf("Error a broken edge has been found! Count=%3d Temp=%3d\n",openEdges->itemCount,temp_openEdges.itemCount);
				forceStop = true;
			}
			
		}
	}
	
	
	
	//===================================================== Clean up
	temp_openEdges.EmptyVector();
}


//############################################################################# Find circumcircle
/** Finds the circumcenter and radius of a triangle.
 *
\author Jan*/
void CHI_SURFACEREMESHER::FindCircumCircle(CHI_PATCH *patch, CST_FACE *tri, float *center, float *radius)
{
	//===================================================== Get vertices
	float* a = patch->Pstar.GetItem(tri->vertex[0]);
	float* b = patch->Pstar.GetItem(tri->vertex[1]);
	float* c = patch->Pstar.GetItem(tri->vertex[2]);
	bool line1Vertical=false;
	bool line2Vertical=false;
	
	//===================================================== Slope 1
	float m1 = 0.0;
	if (fabs(b[0]-a[0])>this->precision)
	{
		float mBA = (b[1]-a[1])/(b[0]-a[0]);
		if (fabs(mBA)>this->precision)
		{
			m1  = -1.0/mBA;
		} else
		{
			line1Vertical = true;
		}
		
	}
	
	float x1 = (b[0]+a[0])/2.0;
	float y1 = (b[1]+a[1])/2.0;
	
	//y = m*x +c
	//c1 = y1-m1*x1
	//yc = m1*xc +y1-m1*x1
	
	//===================================================== Slope 2
	float m2 = 0.0;
	if (fabs(c[0]-b[0])>this->precision)
	{
		float mCB = (c[1]-b[1])/(c[0]-b[0]);
		if (fabs(mCB)>this->precision)
		{
			m2  = -1.0/mCB;
		}
		else
		{
			line2Vertical = true;
		}
	}
	
	float x2 = (c[0]+b[0])/2.0;
	float y2 = (c[1]+b[1])/2.0;
	
	//y = m*x +c
	//c2 = y2-m2*x2
	//yc = m2*xc +y2-m2*x2
	
	//m1*xc +y1-m1*x1 = m2*xc +y2-m2*x2
	//m1*xc -m2*xc  =  m1*x1-y1 - (m2*x2 - y2)
	
	if (line1Vertical)
	{
		center[0] = x1;
		center[1] = m2*(center[0]-x2)+y2;
		center[2] = 0.0;
	}
	else if (line2Vertical)
	{
		center[0] = x2;
		center[1] = m1*(center[0]-x1)+y1;
		center[2] = 0.0;
	}
	else
	{
		center[0] = (  (m1*x1 - y1) - (m2*x2-y2)  )/(m1-m2);
		center[1] = m2*(center[0]-x2)+y2;
		center[2] = 0.0;
	}
	
	
	*radius = sqrt(    (center[0]-b[0])*(center[0]-b[0])  + (center[1]-b[1])*(center[1]-b[1])  );
}

void CHI_SURFACEREMESHER::DumpPatchToScilab(CHI_PATCH *patch)
{/*
	printf("==================================================== SCILAB SCRIPT\n");
	printf("scf(2)\n"
			       "clf(2)\n");
	printf("ponts = zeros(%d,3);\n",patch->Pstar.itemCount);
	for (int v=0;v<patch->Pstar.itemCount;v++)
	{
		//float* curV = patch->Pstar.GetItem(v);
		float* curV = remeshedSurface->vertexStack.GetItem(v);
		printf("ponts(%d,1)=%6.3f; ponts(%d,2)=%6.3f; ponts(%d,3)=%d; ponts(%d,4)=%6.3f;\n",v+1,curV[0],v+1,curV[1],v+1,v+1,v+1,curV[2]);
	}
	printf("cHull = [ ]\n");
	for (int v=0;v<patch->faceList.itemCount;v++)
	{
		CST_FACE* curV = patch->faceList.GetItem(v);
		if (curV!=NULL)
		{
			printf("cHull=[cHull; %d %d %d %d %d %d ]\n",curV->vertex[0],
			       curV->vertex[1],
			       curV->vertex[2],
			       curV->edges[0][2],
			       curV->edges[1][2],
			       curV->edges[2][2]);
		}
		else
		{
			printf("cHull=[cHull; 1 1 1 1 1 1 ]\n");
		}
		
	}

	for (int e=0;e<patch->edges.itemCount;e++)
	{
		auto curEdge = patch->edges.GetItem(e);
		printf("edge(%d,1)=%d; edge(%d,2)=%d\n",e,curEdge[0],e,curEdge[1]);
	}

	for (int el=0;el<patch->edgeLoops.itemCount;el++)
	{
		auto curEdgeLoop = patch->edgeLoops.GetItem(el);
		if (curEdgeLoop!=NULL)
		{
			printf("//Edge loop %d\n",el);
			for (int e=0;e<curEdgeLoop->edges.itemCount;e++)
			{
				auto curEdge = curEdgeLoop->edges.GetItem(e);
				if (curEdge!=NULL)
				{
					printf("edge(%d,1)=%d; edge(%d,2)=%d\n",e,curEdge[0],e,curEdge[1]);
				}

			}
		}

	}
*/
	
	/*
	printf("loud = \%%T\n"
			       "for t=1:(size(cHull)(1))\n"
			       "    firstTri=[\n"
			       "    ponts(cHull(t,1)+1,:)\n"
			       "    ponts(cHull(t,2)+1,:)\n"
			       "    ponts(cHull(t,3)+1,:)\n"
			       "    ponts(cHull(t,1)+1,:)\n"
			       "    ]\n"
			       "    //id=color(\"green\")\n"
			       "    plot2d(firstTri(:,1),firstTri(:,2))\n"
			       "\n"
			       "    if (loud) then\n"
			       "        point1x = mean(firstTri(1:3,1)) +0.2*(mean(firstTri(1:2,1)) - mean(firstTri(1:3,1)));\n"
			       "        point1y = mean(firstTri(1:3,2)) +0.2*(mean(firstTri(1:2,2)) - mean(firstTri(1:3,2)));\n"
			       "        \n"
			       "        point2x = mean(firstTri(1:3,1)) +0.2*(mean(firstTri(2:3,1)) - mean(firstTri(1:3,1)));\n"
			       "        point2y = mean(firstTri(1:3,2)) +0.2*(mean(firstTri(2:3,2)) - mean(firstTri(1:3,2)));\n"
			       "        \n"
			       "        point3x = mean(firstTri(1:3,1)) +0.2*(mean(firstTri(3:4,1)) - mean(firstTri(1:3,1)));\n"
			       "        point3y = mean(firstTri(1:3,2)) +0.2*(mean(firstTri(3:4,2)) - mean(firstTri(1:3,2)));\n"
			       "        xstring(point1x-0.02,point1y-0.02,string(cHull(t,3+1)));\n"
			       "        xstring(point2x-0.02,point2y-0.02,string(cHull(t,3+2)));\n"
			       "        xstring(point3x-0.02,point3y-0.02,string(cHull(t,3+3)));\n"
			       "    \n"
			       "        dx=-0.02\n"
			       "        dy=-0.02\n"
			       "        xstring(mean(firstTri(1:3,1))+dx,mean(firstTri(1:3,2))+dy,string(t-1))\n"
			       "    \n"
			       "        id=color(\"blue\")\n"
			       "        t=get(\"current_entity\")\n"
			       "        t.font_foreground=id\n"
			       "    end\n"
			       "end\n"
			       "scatter(ponts(:,1),ponts(:,2),,\"black\",\".\")\n"
			       "dx=0.005\n"
			       "dy=-0.01\n"
			       "if (loud) then\n"
			       "    //xstring(ponts(:,1)+dx,ponts(:,2)+dy,string(ponts(:,3)-1))\n"
			       "end\n"
			       "\n"
			       "a=gca();\n"
			       "//a.axes_visible = [\"off\" \"off\" \"off\"];\n"
			       "a.box = \"on\"\n"
			       "a.data_bounds = [-0.1,-0.1;1.1,1.1]\n");*/
	//printf("==================================================== SCILAB SCRIPT END\n");
}









/**Finds all triangles in the given patch that have vc in their circumdisk.
\author Jan*/
void CHI_SURFACEREMESHER::CollectTrianglesWithVertexInCircumdisc(CHI_PATCH *patch, CST_FACE* master, CHI_PATCH *bucket, Eigen::Vector3f vc)
{
	for (int e=0; e<3; e++)
	{
		if (master->edges[e][2]>=0)
		{
			CST_FACE* tau_1 = patch->faceList.GetItem(master->edges[e][2]);
			if (tau_1!=NULL)
			{
				float v_0[3];
				float r;
				
				FindCircumCircle(patch, tau_1,&v_0[0],&r);
				
				//printf("Tri %d, center=[%6.3f,%6.3f] r=%6.3f..",master->edges[e][2],v_0[0],v_0[1],r);
				
				Eigen::Vector3f v0(v_0[0],v_0[1],v_0[2]);
				
				if ((vc-v0).norm() < r)
				{
					//printf("collected\n");
					patch->faceList.ClearItem(master->edges[e][2]);                    //Remove from master list
					bucket->faceList.PushItem(tau_1);             //Add to temp list
					//printf("tau %d, %d->%d->%d\n",bucket->faceList.itemCount,tau_1[0],tau_1[1],tau_1[2]);
					CollectTrianglesWithVertexInCircumdisc(patch,tau_1,bucket,vc);
				}
				else
				{
					//printf("skipped\n");
				}
			}
			
		}
	}
}


/** Checks whether the simplex defined by a,b,c intersects any of the prescribed simplices of the patch.
\author Jan*/
bool CHI_SURFACEREMESHER::FindIntersectingSimplices(CHI_PATCH *patch, int va, int vb, int vc)
{

    float* v1 = patch->Pstar.GetItem(va); if (v1==NULL){return false;}
    float* v2 = patch->Pstar.GetItem(vb); if (v2==NULL){return false;}
    float* v3 = patch->Pstar.GetItem(vc); if (v3==NULL){return false;}

    for (int e=0;e<patch->edgeLoops.itemCount;e++)
    {
        auto curLoop = patch->edgeLoops.GetItem(e);
        for (int k=0;k<curLoop->edges.itemCount;k++)
        {
            auto curEdge = curLoop->edges.GetItem(k);
            if (curEdge!=NULL)
            {
                //printf("Checking edge %d %d against triangle %d %d %d ...",curEdge[0],curEdge[1],va,vb,vc);
                float* b = patch->Pstar.GetItem(curEdge[0]);
                float* c = patch->Pstar.GetItem(curEdge[1]);

                float ori1 = Orient2D(v1,b,c)*Orient2D(v2,b,c);
                float ori2 = Orient2D(v2,b,c)*Orient2D(v3,b,c);
                float ori3 = Orient2D(v3,b,c)*Orient2D(v1,b,c);

                Eigen::Vector3f V1(v1[0],v1[1],0.0);
                Eigen::Vector3f V2(v2[0],v2[1],0.0);
                Eigen::Vector3f V3(v3[0],v3[1],0.0);

                Eigen::Vector3f B(b[0],b[1],0.0);
                Eigen::Vector3f C(c[0],c[1],0.0);

                Eigen::Vector3f BC = C-B;
                Eigen::Vector3f V12 = V2-V1;
                Eigen::Vector3f V23 = V3-V2;
                Eigen::Vector3f V31 = V1-V3;

                Eigen::Matrix<float,2,2> A;
                Eigen::Vector2f rhs;
                Eigen::Vector2f x;

                Eigen::Vector3f A12 = V12;

                //========================================= Triangle edge 1
                A12 = V12;
                if (      fabs(BC.cross(A12)(2)) > this->precision)
                {
                    A(0,0) = BC(0); A(0,1) = -1*A12(0);
                    A(1,0) = BC(1); A(1,1) = -1*A12(1);
                    rhs(0) = V1(0)-B(0);
                    rhs(1) = V1(1)-B(1);

                    x = A.colPivHouseholderQr().solve(rhs);
                    if (  (  (x(0)>(0.0+this->precision)) && (x(0)<(1.0-this->precision)) ) &&
                          (  (x(1)>(0.0+this->precision)) && (x(1)<(1.0-this->precision)) )  )
                    {

                        //printf("intersecting V12, %6.3f %6.3f %.8f\n",x(0),x(1),fabs(BC.cross(A12)(2)));
                        return true;
                    }
                }

                //========================================= Triangle edge 2
                A12 = V23;
                if (      fabs(BC.cross(A12)(2)) > this->precision)
                {
                    A(0,0) = BC(0); A(0,1) = -1*A12(0);
                    A(1,0) = BC(1); A(1,1) = -1*A12(1);
                    rhs(0) = V2(0)-B(0);
                    rhs(1) = V2(1)-B(1);

                    x = A.colPivHouseholderQr().solve(rhs);
                    if (  (  (x(0)>(0.0+this->precision)) && (x(0)<(1.0-this->precision)) ) &&
                          (  (x(1)>(0.0+this->precision)) && (x(1)<(1.0-this->precision)) )  )
                    {
                        //printf("intersecting V23, %6.3f %6.3f %.8f\n",x(0),x(1),fabs(BC.cross(A12)(2)));
                        return true;
                    }
                }

                //========================================= Triangle edge 3
                A12 = V31;
                if (      fabs(BC.cross(A12)(2)) > this->precision)
                {
                    A(0,0) = BC(0); A(0,1) = -1*A12(0);
                    A(1,0) = BC(1); A(1,1) = -1*A12(1);
                    rhs(0) = V3(0)-B(0);
                    rhs(1) = V3(1)-B(1);

                    x = A.colPivHouseholderQr().solve(rhs);
                    if (  (  (x(0)>(0.0+this->precision)) && (x(0)<(1.0-this->precision)) ) &&
                          (  (x(1)>(0.0+this->precision)) && (x(1)<(1.0-this->precision)) )  )
                    {
                        //printf("intersecting V31, %6.3f %6.3f %.8f\n",x(0),x(1),fabs(BC.cross(A12)(2)));
                        return true;
                    }
                }
                //printf("clear\n");





            }
        }

    }

    return false;
}



/**Removes filler triangles from a mesh which used ConvexifyHull.
During lexicographical meshing the algorithm convexifies the hull after each vertex attachment. This
function removes triangles that might have been inserted for this purpose but no lay outside of the
original path domain.

The algorithm runs through each triangle of the remeshed surface. It calculates its centroid and then runs
through all the triangles in the initial mesh. For each triangle in the original mesh it first project the
 three vertices to the same 2D plane as the projected vertices used for triangulation of the new mesh.
 It then uses Orient2D on each of the edges of the triangle to see if the calculated centroid is within it.
\author Jan*/
void CHI_SURFACEREMESHER::RemoveFillerTriangles(CHI_PATCH *patch)
{
	CHI_PATCH* curPatch = patch;
	CHI_PATCH* oldPatch = patch->parent;

	for (int t1=0;t1<curPatch->faceList.itemCount;t1++)
	{

		CST_FACE* tau_1 = curPatch->faceList.GetItem(t1);
		if (tau_1 != NULL)
		{
			//======================================================= Calculating centroid
			int vi[] = {tau_1->vertex[0], tau_1->vertex[1],tau_1->vertex[2]};

			float* vf0 = curPatch->Pstar.GetItem(vi[0]);
			float* vf1 = curPatch->Pstar.GetItem(vi[1]);
			float* vf2 = curPatch->Pstar.GetItem(vi[2]);

			float centroid[3];
			int i;
			i=0;centroid[i] = (vf0[i]+vf1[i]+vf2[i])/3.0;
			i=1;centroid[i] = (vf0[i]+vf1[i]+vf2[i])/3.0;
			i=2;centroid[i] = (vf0[i]+vf1[i]+vf2[i])/3.0;

			printf("Face %d, centroid %6.3f %6.3f %6.3f\n",t1,centroid[0],centroid[1],centroid[2]);

			//======================================================= Find triangles with this in there
			bool homeFound = false;
			for (int t2=0;t2<oldPatch->faceList.itemCount;t2++)
			{
				CST_FACE* tau_2 = oldPatch->faceList.GetItem(t2);
				if (tau_2 != NULL)
				{
					int vbi[] = {tau_2->vertex[0], tau_2->vertex[1],tau_2->vertex[2]};
					float* vbf0 = initialMesh->vertexStack.GetItem(vbi[0]);
					float* vbf1 = initialMesh->vertexStack.GetItem(vbi[1]);
					float* vbf2 = initialMesh->vertexStack.GetItem(vbi[2]);

					Eigen::Vector3f V0(vbf0[0],vbf0[1],vbf0[2]);
					Eigen::Vector3f V1(vbf1[0],vbf1[1],vbf1[2]);
					Eigen::Vector3f V2(vbf2[0],vbf2[1],vbf2[2]);

					Eigen::Vector3f P;

					P(0) = V0[0]-curPatch->C_0(0);
					P(1) = V0[1]-curPatch->C_0(1);
					P(2) = V0[2]-curPatch->C_0(2);

					float vpf0[] = {P.dot(curPatch->hat_i),P.dot(curPatch->hat_j),0.0};

					P(0) = V1[0]-curPatch->C_0(0);
					P(1) = V1[1]-curPatch->C_0(1);
					P(2) = V1[2]-curPatch->C_0(2);

					float vpf1[] = {P.dot(curPatch->hat_i),P.dot(curPatch->hat_j),0.0};

					P(0) = V2[0]-curPatch->C_0(0);
					P(1) = V2[1]-curPatch->C_0(1);
					P(2) = V2[2]-curPatch->C_0(2);

					float vpf2[] = {P.dot(curPatch->hat_i),P.dot(curPatch->hat_j),0.0};

					//printf("Vertex 0 %6.3f %6.3f %6.3f, ori=%.10f\n",t1,vpf0[0],vpf0[1],vpf0[2],Orient2D(centroid,vpf0,vpf1));
					//printf("Vertex 1 %6.3f %6.3f %6.3f, ori=%.10f\n",t1,vpf1[0],vpf1[1],vpf1[2],Orient2D(centroid,vpf1,vpf2));
					//printf("Vertex 2 %6.3f %6.3f %6.3f, ori=%.10f\n",t1,vpf2[0],vpf2[1],vpf2[2],Orient2D(centroid,vpf2,vpf0));

					//===================================== Checking orientations
					bool to_the_left0=false; if (Orient2D(centroid,vpf0,vpf1)>(0.0-this->precision)){to_the_left0=true;}
					bool to_the_left1=false; if (Orient2D(centroid,vpf1,vpf2)>(0.0-this->precision)){to_the_left1=true;}
					bool to_the_left2=false; if (Orient2D(centroid,vpf2,vpf0)>(0.0-this->precision)){to_the_left2=true;}

					if (to_the_left0 && to_the_left1 && to_the_left2)
					{
						//printf("Home found for triangle %d %d %d\n",vi[0],vi[1],vi[2]);
						homeFound = true;
						break;
					}
				}
			}

			//======================================================= If no matching triangle, remove triangle
			if (!homeFound)
			{
				printf("Home not found for triangle %d %d %d\n",vi[0],vi[1],vi[2]);
				curPatch->faceList.SetItem(t1,NULL);

				for (int e=0;e<3;e++)
				{
					if (tau_1->edges[e][2]>=0)
					{
						CST_FACE* tau_2 = curPatch->faceList.GetItem(tau_1->edges[e][2]);
						if (tau_2!=NULL)
						{
							for (int e2=0;e2<3;e2++)
							{
								if ( (tau_2->edges[e2][0] == tau_1->edges[e][1]) &&
									 (tau_2->edges[e2][1] == tau_1->edges[e][0]) )
								{
									tau_2->edges[e2][2] = -1;
								}
							}
						}
					}

				}
			}



		}
	}
}