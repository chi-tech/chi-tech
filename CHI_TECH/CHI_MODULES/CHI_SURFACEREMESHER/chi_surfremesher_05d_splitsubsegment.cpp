#include"chi_surfaceremesher.h"

//############################################################################# Find Encroached edge
/** Finds an encroached subsegment.
\author Jan*/
bool CHI_SURFACEREMESHER::SplitSubSegment(CHI_PATCH *patch, CHI_INTERIOREDGE *encrSeg)
{
	if (encrSeg==NULL){return false;}
	//=============================================================== Get edge center v_c
	float* v_i = patch->Pstar.GetItem(encrSeg->vi);
	float* v_f = patch->Pstar.GetItem(encrSeg->vf);
	
	Eigen::Vector3f vi(v_i[0],v_i[1],v_i[2]);
	Eigen::Vector3f vf(v_f[0],v_f[1],v_f[2]);
	
	//printf("v_i = %8.5f %8.5f %8.5f \n",v_i[0],v_i[1],v_i[2]);
	//printf("v_f = %8.5f %8.5f %8.5f \n",v_f[0],v_f[1],v_f[2]);

	Eigen::Vector3f vc = 0.5*(vi+vf);

	//=============================================================== Insert the vertex into the mesh
	//float* newVert = new float[3];
	float newVert[3];
	Eigen::Vector3f Reproject(0,0,0);
	Reproject = Reproject + vc(0)*patch->hat_i;
	Reproject = Reproject + vc(1)*patch->hat_j;
	//Reproject += vc(2)*patch->hat_k;
	Reproject = Reproject + patch->C_0;
	newVert[0] = Reproject(0);
	newVert[1] = Reproject(1);
	newVert[2] = Reproject(2);

	//printf("Hat i = [%6.3f,%6.3f,%6.3f]\n",patch->hat_i(0),patch->hat_i(1),patch->hat_i(2));
	//printf("Hat j = [%6.3f,%6.3f,%6.3f]\n",patch->hat_j(0),patch->hat_j(1),patch->hat_j(2));
	//printf("Hat k = [%6.3f,%6.3f,%6.3f]\n",patch->hat_k(0),patch->hat_k(1),patch->hat_k(2));

	//printf("Inserting Vertex %8.5f %8.5f %8.5f\n",newVert[0],newVert[1],newVert[2]);

	//printf("Splitting Edge %d->%d %6.3f %6.3f %6.3f\n",encrSeg->vi,encrSeg->vf,vc(0),vc(1),vc(2));
	InsertVertex(patch, &newVert[0],encrSeg);

	return false;

	
}