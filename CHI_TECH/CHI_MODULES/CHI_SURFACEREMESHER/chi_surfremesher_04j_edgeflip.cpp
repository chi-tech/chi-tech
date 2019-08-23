#include"chi_surfaceremesher.h"

//############################################################################# Edge flip algorithm
/**Flips the edge a triangle pair associated with an interior edge that
 * is not locally delaunay.
 *
\author Jan*/
void CHI_SURFACEREMESHER::EdgeFlip(CHI_PATCH *patch, CHI_INTERIOREDGE* non_loc_del_edge)
{
	if (non_loc_del_edge==NULL){return;}
	CHI_INTERIOREDGE* curEdge = non_loc_del_edge;

	printf("Flipping Edge %d->%d\n",curEdge->vi,curEdge->vf);
	
	//=============================================================== Get edges triangles
	CST_FACE* tau_1 = patch->faceList.GetItem(curEdge->tau_1);
	CST_FACE* tau_2 = patch->faceList.GetItem(curEdge->tau_2);
	
	//=============================================================== Find neighbour A and B
	int neighborA = -1;
	int neighborB = -1;
	     if (curEdge->tau_1_edgeNumber==0)
	{
		neighborA = tau_1->edges[1][2];
		neighborB = tau_1->edges[2][2];
	}
	else if (curEdge->tau_1_edgeNumber==1)
	{
		neighborA = tau_1->edges[2][2];
		neighborB = tau_1->edges[0][2];
	}
	else
	{
		neighborA = tau_1->edges[0][2];
		neighborB = tau_1->edges[1][2];
	}
	
	
	//=============================================================== Find neighbour C and D
	int neighborC = -1;
	int neighborD = -1;
	if (curEdge->tau_2_edgeNumber==0)
	{
		neighborC = tau_2->edges[1][2];
		neighborD = tau_2->edges[2][2];
	}
	else if (curEdge->tau_2_edgeNumber==1)
	{
		neighborC = tau_2->edges[2][2];
		neighborD = tau_2->edges[0][2];
	}
	else
	{
		neighborC = tau_2->edges[0][2];
		neighborD = tau_2->edges[1][2];
	}
	
	//=============================================================== Find vertex A and B
	int vertA = -1;
	for (int v=0; v<3; v++)
	{
		if ( (tau_1->vertex[v]!=curEdge->vi) && (tau_1->vertex[v]!=curEdge->vf) )
		{
			vertA = tau_1->vertex[v];
			break;
		}
	}
	int vertB = -1;
	for (int v=0; v<3; v++)
	{
		if ( (tau_2->vertex[v]!=curEdge->vi) && (tau_2->vertex[v]!=curEdge->vf) )
		{
			vertB = tau_2->vertex[v];
			break;
		}
	}
	/*printf("Tau1 %3d,%3d,%3d -- %3d,%3d,%3d   Tau2 %3d,%3d,%3d -- %3d,%3d,%3d\n",
	        tau_1->vertex[0],tau_1->vertex[1],tau_1->vertex[2],tau_1->edges[0][2],tau_1->edges[1][2],tau_1->edges[2][2],
	        tau_2->vertex[0],tau_2->vertex[1],tau_2->vertex[2],tau_2->edges[0][2],tau_2->edges[1][2],tau_2->edges[2][2]);*/
	//=============================================================== Modify tau_1 and tau_2
	tau_1->vertex[0] = vertA;       tau_1->edges[0][2]=curEdge->tau_2;
	tau_1->vertex[1] = vertB;       tau_1->edges[1][2]=neighborD;
	tau_1->vertex[2] = curEdge->vf; tau_1->edges[2][2]=neighborA;
	
	tau_1->edges[0][0] = tau_1->vertex[0]; tau_1->edges[0][1] = tau_1->vertex[1];
	tau_1->edges[1][0] = tau_1->vertex[1]; tau_1->edges[1][1] = tau_1->vertex[2];
	tau_1->edges[2][0] = tau_1->vertex[2]; tau_1->edges[2][1] = tau_1->vertex[0];
	
	tau_2->vertex[0] = vertB;       tau_2->edges[0][2]=curEdge->tau_1;
	tau_2->vertex[1] = vertA;       tau_2->edges[1][2]=neighborB;
	tau_2->vertex[2] = curEdge->vi; tau_2->edges[2][2]=neighborC;
	
	tau_2->edges[0][0] = tau_2->vertex[0]; tau_2->edges[0][1] = tau_2->vertex[1];
	tau_2->edges[1][0] = tau_2->vertex[1]; tau_2->edges[1][1] = tau_2->vertex[2];
	tau_2->edges[2][0] = tau_2->vertex[2]; tau_2->edges[2][1] = tau_2->vertex[0];
	
	/*printf("Tau1 %3d,%3d,%3d -- %3d,%3d,%3d   Tau2 %3d,%3d,%3d -- %3d,%3d,%3d\n",
	       tau_1->vertex[0],tau_1->vertex[1],tau_1->vertex[2],tau_1->edges[0][2],tau_1->edges[1][2],tau_1->edges[2][2],
	       tau_2->vertex[0],tau_2->vertex[1],tau_2->vertex[2],tau_2->edges[0][2],tau_2->edges[1][2],tau_2->edges[2][2]);*/
	
	//=============================================================== Update neighbors B and D

	for (int e=0; e<3; e++)
	{
		if (neighborB>0)
		{
			CST_FACE* tau_neighborB = patch->faceList.GetItem(neighborB);
			if (tau_neighborB->edges[e][2]==curEdge->tau_1)
			{
				tau_neighborB->edges[e][2]=curEdge->tau_2;
			}
		}
		
		if (neighborD>0)
		{
			CST_FACE* tau_neighborD = patch->faceList.GetItem(neighborD);
			if (tau_neighborD->edges[e][2]==curEdge->tau_2)
			{
				tau_neighborD->edges[e][2]=curEdge->tau_1;
			}
		}
	}
	
	
	
}