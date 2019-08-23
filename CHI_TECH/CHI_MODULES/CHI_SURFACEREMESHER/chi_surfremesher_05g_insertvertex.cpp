#include "chi_surfaceremesher.h"

//############################################################################# Insert vertex
/**Inserts a vertex into the triangulation.
\author Jan*/
void CHI_SURFACEREMESHER::InsertVertex(CHI_PATCH *patch, float *vertex, CHI_INTERIOREDGE *encrSeg)
{
	//printf("================= Inserting vertex algorithm: %d\n");
	//=============================================================== Insert the vertex into the triangulation
	float* newVert = new float[3];
	newVert[0] = vertex[0];
	newVert[1] = vertex[1];
	newVert[2] = vertex[2];
	this->remeshedSurface->vertexStack.PushItem(newVert);
	//printf("Inserting real Vertex %8.5f %8.5f %8.5f\n",vertex[0],vertex[1],vertex[2]);
	
	//=============================================================== Calculate projected vertex
	Eigen::Vector3f P;
	Eigen::Vector3f vc;
	
	P(0) = vertex[0]-patch->C_0(0);
	P(1) = vertex[1]-patch->C_0(1);
	P(2) = vertex[2]-patch->C_0(2);
	vc(0) =P.dot(patch->hat_i);
	vc(1) =P.dot(patch->hat_j);
	vc(2) =0.0;
	
	float* projectedVertex = new float[3];
	projectedVertex[0] = vc(0);
	projectedVertex[1] = vc(1);
	projectedVertex[2] = vc(2);
	
	//printf("Inserting Projected Vertex %8.5f %8.5f %8.5f\n",projectedVertex[0],projectedVertex[1],projectedVertex[2]);


	//=============================================================== Insert the projected vertex globally
	int insertedVertexIndex = -1;
	for (int p=0;p<remeshedSurface->patches.itemCount;p++)
	{
		CHI_PATCH* curPatch = remeshedSurface->patches.GetItem(p);
		insertedVertexIndex = curPatch->Pstar.PushItem(projectedVertex);
	}

	
	
	//=============================================================== Find triangles with the vertex in their circumcircles
	//                                                                Add them to a temporary patch
	CHI_PATCH tempPatch;
	for (int t=patch->faceList.itemCount-1; t>=0; t--)
	{
		CST_FACE* tau_1 = patch->faceList.GetItem(t);
		if (tau_1!=NULL)
		{
			float v_0[3];
			float r;
			
			FindCircumCircle(patch, tau_1,&v_0[0],&r);

			Eigen::Vector3f v0(v_0[0],v_0[1],v_0[2]);
			
			if ((vc-v0).norm() < r)
			{
				patch->faceList.ClearItem(t);                    //Remove from master list
				tempPatch.faceList.PushItem(tau_1);             //Add to temp list
				//printf("tau %d, %d->%d->%d\n",tempPatch.faceList.itemCount,tau_1[0],tau_1[1],tau_1[2]);
				//printf("Recursing down\n");
				CollectTrianglesWithVertexInCircumdisc(patch,tau_1,&tempPatch,vc);
				//printf("Done recursing\n");
				break;
			}
		}
	}
	
	//=============================================================== Find the open edges and edge loop of the temp patch
	FindOpenEdges(&tempPatch);
	FindEdgeLoops(&tempPatch);
	CHI_EDGELOOP* curLoop = tempPatch.edgeLoops.GetItem(0);
	
	
	
	//=============================================================== Define new triangles per edge
	int countBeforeAddition=patch->faceList.itemCount-1;
	int lastTri=-1;
	int firstTri=-1;
	int prevTri=-1;
	//printf("Edgeloop:\n");
	for (int e=0;e<curLoop->edges.itemCount;e++)
	{
		int* curEdge = curLoop->edges.GetItem(e);
		float* a = projectedVertex;
		float* b = patch->Pstar.GetItem(curEdge[0]);
		float* c = patch->Pstar.GetItem(curEdge[1]);
		
		
		
		double orientation = Orient2D(a,b,c);
		//printf("  -edge %d->%d ori %6.3f\n",curEdge[0],curEdge[1], orientation);
		//if (!((curEdge[0]==encrSeg->vi) && (curEdge[1]==encrSeg->vf)))
		if (fabs(orientation)>this->precision)
		{
			//printf("Edge %3d, from %3d->%3d\n",e,curEdge[0],curEdge[1]);
			//============================================= Instantiate a new triangle
			CST_FACE* newTri = new CST_FACE;
			newTri->vertex[0] = curEdge[0];
			newTri->vertex[1] = curEdge[1];
			newTri->vertex[2] = insertedVertexIndex;
            newTri->faceNormal[0] = patch->normal(0);
            newTri->faceNormal[1] = patch->normal(1);
            newTri->faceNormal[2] = patch->normal(2);
			
			//printf("Tri %d->%d->%d\n",newTri->vertex[0],newTri->vertex[1],newTri->vertex[2]);
			//============================================= Initialize its edges (no neighbor)
			for (int s=0; s<3; s++)
			{
				if (s<2)
				{
					newTri->edges[s][0] = newTri->vertex[s];
					newTri->edges[s][1] = newTri->vertex[s+1];
				} else
				{
					newTri->edges[s][0] = newTri->vertex[s];
					newTri->edges[s][1] = newTri->vertex[0];
				}
			}

			//============================================= Set the edge neighbor for the adjacent triangle
			//												in the original patch
			newTri->edges[0][2] = curEdge[2];  //Only the neighbor in the original patch is known

			//============================================= If you are at the end of the edge loop
			//                                              the outgoing edge is connected to the first triangle
			if (e==(curLoop->edges.itemCount-1))
			{
				newTri->edges[1][2] = firstTri;
			}

			//============================================= If this vertex was inserted upon an encroached edge
			//												and this triangle's base is on that edge
			//												Set its neighbor to empty
			if (encrSeg!=NULL)
			{
				if (  (newTri->edges[1][0] == encrSeg->vi)         && (newTri->edges[1][1] == insertedVertexIndex)  )
				{
					//printf("encroached edge tau_2=%d\n",encrSeg->tau_2);
					if (encrSeg->tau_2<0) {newTri->edges[1][2]=-1;}
				}
			}

			//============================================= The last triangle's can now be connected to the
			//												current triangle
			if (lastTri>=0)
			{
				newTri->edges[2][2] = lastTri;
			}
			//============================================= If this vertex was inserted upon an encroached edge
			//												and this triangle's base is on that edge
			//												Set its neighbor to empty
			if (encrSeg!=NULL)
			{
				if (  (newTri->edges[2][0] == insertedVertexIndex)         && (newTri->edges[2][1] == encrSeg->vf)  )
				{
					if (encrSeg->tau_2<0) {newTri->edges[2][2]=-1;lastTri=-1;}
				}
			}
			
			//printf("New Tri %d %d %d\n",newTri->vertex[0],
			//newTri->vertex[1],
			//newTri->vertex[2]);
			
			
			//============================================= Add the triangle
			//tempPatch2.faceList.PushItem(newTri);
			prevTri = lastTri;
			lastTri = patch->faceList.AddItem(newTri);
			if (e==0)
			{
				firstTri = lastTri;
			}
			
			//============================================= Notify neighbors
			CST_FACE* last_tau = patch->faceList.GetItem(prevTri);
			if (last_tau!=NULL)
			{
				last_tau->edges[1][2] = lastTri;
			}
			
			if (e==(curLoop->edges.itemCount-1))
			{
				CST_FACE* first_tau = patch->faceList.GetItem(firstTri);
				if (first_tau!=NULL)
				{
					first_tau->edges[2][2] = lastTri;
				}
			}
			CST_FACE* neighbour = patch->faceList.GetItem(newTri->edges[0][2]);
			if (neighbour!=NULL)
			{
				for (int s=0;s<3;s++)
				{
					if ((neighbour->edges[s][1]==newTri->edges[0][0]) && (neighbour->edges[s][0]==newTri->edges[0][1]))
					{
						neighbour->edges[s][2] = lastTri;
					}
				}
			}
			
			
		}
		
		
	}
	
}