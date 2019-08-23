#include "chi_surfaceremesher.h"

//############################################################################# Delete a vertex
/**Delete a vertex from a patch.
 *
\author Jan*/
void CHI_SURFACEREMESHER::DeleteVertex(CHI_PATCH *patch, int vertexIndex)
{
	printf("================= Deleting vertex algorithm: %d\n",vertexIndex);
	CHI_PATCH tempPatch;
	
	float* deletingVertex = patch->Pstar.GetItem(vertexIndex);
	Eigen::Vector3f vc(deletingVertex[0],deletingVertex[1],deletingVertex[2]);
	
	printf("  -Extracting connected triangles\n");
	//=============================================================== Extracting connected triangles
	for (int t=patch->faceList.itemCount-1; t>=0; t--)
	{
		CST_FACE* curFace = patch->faceList.GetItem(t);
		if (curFace!=NULL)
		{
			for (int v=0; v<3; v++)
			{
				if (curFace->vertex[v]==vertexIndex)
				{
					printf("    -Pulling face %d\n",t);
					patch->faceList.ClearItem(t);
					tempPatch.faceList.PushItem(curFace);
				}
			}
		}
	}
	
	printf("  -Extracting triangles within the circumcircle of removed vertex\n");
	//=============================================================== Extracting Triangles within the circumcircle of removed vertex
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
				//printf("Recursing down\n");
				CollectTrianglesWithVertexInCircumdisc(patch,tau_1,&tempPatch,vc);
				//printf("Done recursing\n");
				break;
			}
		}
	}
	
	printf("  -Finding edge loops\n");
	//=============================================================== Find edgeloop
	FindOpenEdges(&tempPatch);
	FindEdgeLoops(&tempPatch);
	
	
	
	//=============================================================== Form first triangle
	CHI_EDGELOOP* firstLoop = tempPatch.edgeLoops.GetItem(0);
	
	//===================================================== Check if edgeloop contains the removed vertex
	bool containsRemoved = false;
	printf("  -Edge loop:\n");
	for (int e=0;e<firstLoop->edges.itemCount;e++)
	{
		int* curEdge = firstLoop->edges.GetItem(e);
		printf("    -edge %d, %d->%d\n",e,curEdge[0],curEdge[1]);
		
		if (curEdge[0]==vertexIndex)
		{
			int* prevEdge = NULL;
			if (e>0) {prevEdge = firstLoop->edges.GetItem(e-1);}
			else     {prevEdge = firstLoop->edges.GetItem(firstLoop->edges.itemCount-1);}
			prevEdge[1] = curEdge[1];
			firstLoop->edges.PullItem(e);
		}
	}
	
	printf("  -Form first triangle\n");
	int* edgeNm1 = firstLoop->edges.GetItem(firstLoop->edges.itemCount-1);
	int* edgeNm2 = firstLoop->edges.GetItem(firstLoop->edges.itemCount-2);
	CST_FACE* firstFace = new CST_FACE;
	firstFace->vertex[0] = edgeNm1[1];
	firstFace->vertex[1] = edgeNm2[0];
	firstFace->vertex[2] = edgeNm1[0]; //or edgeNm2[1]
	
	firstFace->edges[0][0] = firstFace->vertex[0]; firstFace->edges[0][1] = firstFace->vertex[1];
	firstFace->edges[1][0] = firstFace->vertex[1]; firstFace->edges[1][1] = firstFace->vertex[2];
	firstFace->edges[2][0] = firstFace->vertex[2]; firstFace->edges[2][1] = firstFace->vertex[0];
	
	//firstFace->edges[0][2] = unknown;
	firstFace->edges[1][2] = edgeNm2[2];
	firstFace->edges[2][2] = edgeNm1[2];

	firstFace->faceNormal[0] = patch->normal(0);
	firstFace->faceNormal[1] = patch->normal(1);
	firstFace->faceNormal[2] = patch->normal(2);
	
	
	int lastTri = patch->faceList.AddItem(firstFace);
	printf("  -First Face %d,  %d->%d->%d neighbors %d,%d,%d\n",lastTri,firstFace->vertex[0],firstFace->vertex[1],firstFace->vertex[2],
	       firstFace->edges[0][2],
	       firstFace->edges[1][2],
	       firstFace->edges[2][2]);
	
	printf("  -Update first neighbors\n");
	//=============================================================== Update first neighbors
	for (int e=1;e<3;e++)
	{
		if (firstFace->edges[e][2]>=0)
		{
			CST_FACE* tau_2 = patch->faceList.GetItem(firstFace->edges[e][2]);
			
			//================================================= Check which side is adjacent
			for (int s=0; s<3; s++)
			{
				if (  (tau_2->edges[s][0]==firstFace->edges[e][1]) && (tau_2->edges[s][1]==firstFace->edges[e][0])  )
				{
					tau_2->edges[s][2] = lastTri;
					break;
				}
			}
		}
		
	}
	
	printf("  -Running over other edges\n");
	//=============================================================== Running over rest of edges
	for (int e=firstLoop->edges.itemCount-3; e>=1; e--) //must be greater than 1
	{
		int* curEdge = firstLoop->edges.GetItem(e);
		int vB_index = curEdge[0];
		int vf_index = curEdge[1];
		int vi_index = edgeNm1[1];
		int vA_index = -1;
		printf("    -Edge %d->%d neighbor %d\n",vB_index,vf_index,curEdge[2]);
		
		printf("    -Create new face\n");
		//================================================= Create new face
		CST_FACE* tau_2 = new CST_FACE;
		tau_2->vertex[0] = vB_index;
		tau_2->vertex[1] = vf_index;
		tau_2->vertex[2] = vi_index;
		
		tau_2->edges[0][0] = tau_2->vertex[0]; tau_2->edges[0][1] = tau_2->vertex[1];
		tau_2->edges[1][0] = tau_2->vertex[1]; tau_2->edges[1][1] = tau_2->vertex[2];
		tau_2->edges[2][0] = tau_2->vertex[2]; tau_2->edges[2][1] = tau_2->vertex[0];
		
		tau_2->edges[0][2] = curEdge[2];
		tau_2->edges[1][2] = lastTri;
		//newFace->edges[2][2] = unknown;

		tau_2->faceNormal[0] = patch->normal(0);
		tau_2->faceNormal[1] = patch->normal(1);
		tau_2->faceNormal[2] = patch->normal(2);

		if (e==1)
		{
			int* lastEdge = firstLoop->edges.GetItem(0);
			tau_2->edges[2][2] = lastEdge[2];
		}
		
		lastTri = patch->faceList.AddItem(tau_2);
		printf("    -New Face %d,  %d->%d->%d neighbors %d,%d,%d\n",lastTri,tau_2->vertex[0],tau_2->vertex[1],tau_2->vertex[2],
		       tau_2->edges[0][2],
		       tau_2->edges[1][2],
		       tau_2->edges[2][2]);
		
		printf("    -Update neighbors\n");
		//================================================= Update neighbors
		int tau_1_edgenum=-1;
		for (int e2=0;e2<3;e2++)
		{
			if (tau_2->edges[e2][2]>=0)
			{
				CST_FACE* tau_1 = patch->faceList.GetItem(tau_2->edges[e2][2]);
				
				//=============================== Check which side is adjacent
				for (int s=0; s<3; s++)
				{
					if (  (tau_1->edges[s][0]==tau_2->edges[e2][1]) && (tau_1->edges[s][1]==tau_2->edges[e2][0])  )
					{
						tau_1->edges[s][2] = lastTri;
						if (e2==1)
						{
							if (s==2) {vA_index = tau_1->edges[0  ][1];}
							else      {vA_index = tau_1->edges[s+1][1];}
							tau_1_edgenum = s;
						}
						break;
					}
				}
			}
		}
		printf("    -Check if edge-flipping is needed\n");
		//================================================= Check if edge-flipping is needed
		//=============================== Find indexes for abc
		int aindex = vB_index;
		int bindex = vf_index;
		int cindex = vi_index;
		int dindex = vA_index;
		printf("    -a,b,c,d %d,%d,%d,%d\n",aindex,	bindex,	cindex,	dindex);
		
		//=============================== Get the vertex values of abcd
		float* a = patch->Pstar.GetItem(aindex);
		float* b = patch->Pstar.GetItem(bindex);
		float* c = patch->Pstar.GetItem(cindex);
		float* d = patch->Pstar.GetItem(dindex);
		
		double localDelaunay = InCircle(a,b,c,d);
		if (localDelaunay>0.000001)
		{
			CHI_INTERIOREDGE newEdge;
			newEdge.vi = vi_index;
			newEdge.vf = vf_index;
			newEdge.tau_1 = tau_2->edges[1][2];
			newEdge.tau_1_edgeNumber = tau_1_edgenum;
			newEdge.tau_2 = lastTri;
			newEdge.tau_2_edgeNumber = 1;
			
			EdgeFlip(patch,&newEdge);
		}
	}
	
	
}