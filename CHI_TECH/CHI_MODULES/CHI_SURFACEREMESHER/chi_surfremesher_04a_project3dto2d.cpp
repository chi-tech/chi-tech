#include "chi_surfaceremesher.h"

//############################################################################# Project 3D points onto 2D plane
/**Projects 3D points onto a 2D plane using the patch normal.

\author Jan*/
void CHI_SURFACEREMESHER::Project3DVerticesTo2D(CHI_REMESHEDSURFACE *mesh, CHI_PATCH *patch)
{
	//Eigen::Vector3f hati;
	//Eigen::Vector3f hatj;
	//Eigen::Vector3f hatk;
	
	//===================================================================== Getting hatk from patch normal
	patch->hat_k = patch->normal;
	patch->hat_k.normalize();
	
	//===================================================================== Copy all vertices to the temporary stack
	patch->Pstar.ClearVector();
	for (int i=0;i<mesh->vertexStack.itemCount;i++)
	{
		float* curVert = mesh->vertexStack.GetItem(i);
		float* pstari = new float[3];
		
		pstari[0] = curVert[0];
		pstari[1] = curVert[1];
		pstari[2] = curVert[2];
		
		patch->Pstar.PushItem(pstari);
	}
	
	//========================================================================= Calculate centroid of essential vertices
	//Eigen::Vector3f C_0; //Centroid
	if (!keep2Dorientation)
	{
		float sumx=0.0;
		float sumy=0.0;
		float sumz=0.0;
		int counter=0;
		for (int i=0;i<patch->essentialVertices.itemCount;i++)
		{
			int* vertIndex = patch->essentialVertices.GetItem(i);
			float* curVert = mesh->vertexStack.GetItem(*vertIndex);
			
			sumx += curVert[0];
			sumy += curVert[1];
			sumz += curVert[2];
			counter++;
		}
		patch->C_0(0) = sumx/counter;
		patch->C_0(1) = sumy/counter;
		patch->C_0(2) = sumz/counter;
	}
	
	
	//========================================================================= Calculate maximum D and hati
	//This is one method to calculate \hat{i}, one other method could be to
	//look for cardinal directions and align the x-axis to that. WIP
	Eigen::Vector3f P; //Point of interest
	float maxD=0.0;

	for (int el=0;el<patch->edgeLoops.itemCount;el++)
	{
		CHI_EDGELOOP* curLoop=patch->edgeLoops.GetItem(el);
		for (int e = 0; e < curLoop->edges.itemCount; e++)
		{
			int *curEdge = curLoop->edges.GetItem(e);
			float *v0 = mesh->vertexStack.GetItem(curEdge[0]);
			float *v1 = mesh->vertexStack.GetItem(curEdge[1]);
			
			P = Eigen::Vector3f(v1[0] - v0[0],
			                    v1[1] - v0[1],
			                    v1[2] - v0[2]);
			if (P.norm() > maxD)
			{
				maxD = P.norm();
				patch->hat_i = P;
			}
		}
	}
	patch->hat_i.normalize();
	
	
	//========================================================================= Calculate hatj
	patch->hat_j = patch->hat_k.cross(patch->hat_i);
	//printf("Hat j = [%6.3f,%6.3f,%6.3f]\n",hatj(0),hatj(1),hatj(2));
	//printf("Hat k = [%6.3f,%6.3f,%6.3f]\n",hatk(0),hatk(1),hatk(2));
	
	if (keep2Dorientation)
	{
		patch->hat_i = Eigen::Vector3f(1,0,0);
		patch->hat_j = Eigen::Vector3f(0,1,0);
		patch->hat_k = Eigen::Vector3f(0,0,1);
		patch->C_0 = Eigen::Vector3f(0,0,0);
	}

	//printf("Hat i = [%6.3f,%6.3f,%6.3f]\n",patch->hat_i(0),patch->hat_i(1),patch->hat_i(2));
	//printf("Hat j = [%6.3f,%6.3f,%6.3f]\n",patch->hat_j(0),patch->hat_j(1),patch->hat_j(2));
	//printf("Hat k = [%6.3f,%6.3f,%6.3f]\n",patch->hat_k(0),patch->hat_k(1),patch->hat_k(2));
	//printf("Cntrd = [%6.3f,%6.3f,%6.3f]\n",patch->C_0(0),patch->C_0(1),patch->C_0(2));
	
	//========================================================================= Calculate transformed vertices
	for (int i=0;i<patch->Pstar.itemCount;i++)
	{
		float* curVert = patch->Pstar.GetItem(i);
		
		//printf("Transformed vert %3d, from [%6.3f,%6.3f,%6.3f]", i, curVert[0], curVert[1], curVert[2]);
		
		
		P(0) = curVert[0]-patch->C_0(0);
		P(1) = curVert[1]-patch->C_0(1);
		P(2) = curVert[2]-patch->C_0(2);
		
		curVert[0] = P.dot(patch->hat_i);
		curVert[1] = P.dot(patch->hat_j);
		curVert[2] = 0.0;
		
		//printf(" to [%6.3f,%6.3f,%6.3f]\n", curVert[0], curVert[1], curVert[2]);
	}
}