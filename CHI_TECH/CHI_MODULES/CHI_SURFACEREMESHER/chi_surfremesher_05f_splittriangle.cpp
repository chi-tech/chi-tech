#include "chi_surfaceremesher.h"

//############################################################################# Split triangle
/**Splits a triangle.

\author Jan*/
void CHI_SURFACEREMESHER::SplitTriangle(CHI_PATCH *patch, CST_FACE *tri)
{
	if (tri==NULL){return;}
	//=============================================================== Find the circumcenter of the triangle
	float vcenter[3];
	float r;
	FindCircumCircle(patch,tri,vcenter,&r);
	//printf("Tri being split %d->%d->%d\n",tri->vertex[0],tri->vertex[1],tri->vertex[2]);
	//printf("vc %6.3f %6.3f %6.3f\n",vcenter[0],vcenter[1],vcenter[2]);
	//=============================================================== Find the open edges and edge loops of the current patch
	FindOpenEdges(patch);
	//printf("Find edgeloops\n");
	FindEdgeLoops(patch);
	
	//=============================================================== Check if the new vertex would encroach
	//printf("Check encroachment\n");
	bool encroachedFound = false;
	CHI_INTERIOREDGE* encrEdge = NULL;
	for (int el=0; el<patch->edgeLoops.itemCount; el++)
	{
		CHI_EDGELOOP* curLoop = patch->edgeLoops.GetItem(el);
		for (int e=0; e<curLoop->edges.itemCount; e++)
		{
			int* curEdge = curLoop->edges.GetItem(e);
			if (curEdge!=NULL)
			{
				float* v[3];
				v[0] = patch->Pstar.GetItem(curEdge[0]);
				v[1] = patch->Pstar.GetItem(curEdge[1]);
				v[2] = &vcenter[0];
				
				Eigen::Vector3f v_0(v[0][0],v[0][1],v[0][2]);
				Eigen::Vector3f v_1(v[1][0],v[1][1],v[1][2]);
				Eigen::Vector3f v_2(v[2][0],v[2][1],v[2][2]);
				
				Eigen::Vector3f v_c = 0.5*(v_0+v_1);
				float r = (v_1-v_c).norm();
				float d = (v_2-v_c).norm();
				
				
				if (d<(r+0.000001))
				{
					encroachedFound = true;
					encrEdge = new CHI_INTERIOREDGE;
					encrEdge->vi    = curEdge[0];
					encrEdge->vf    = curEdge[1];
					encrEdge->tau_1 = curEdge[2]; //Assign triangle
					break;
				}
			}
		}
		if (encroachedFound) {break;}
	}
	
	//=============================================================== If an encroached edge was found
	if (encroachedFound)
	{
		SplitSubSegment(patch,encrEdge);
	}
	//=============================================================== ELSE split triangle
	else
	{
		Eigen::Vector3f vc(vcenter[0],vcenter[1],vcenter[2]);
		float newVert[3];
		Eigen::Vector3f Reproject(0,0,0);
		Reproject = Reproject + vc(0)*patch->hat_i;
		Reproject = Reproject + vc(1)*patch->hat_j;
		//Reproject += vc(2)*patch->hat_k;
		Reproject = Reproject + patch->C_0;
		newVert[0] = Reproject(0);
		newVert[1] = Reproject(1);
		newVert[2] = Reproject(2);
		InsertVertex(patch,newVert);
	}
	
}