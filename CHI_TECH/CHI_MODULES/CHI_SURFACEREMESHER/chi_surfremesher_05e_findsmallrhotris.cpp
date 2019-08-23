#include "chi_surfaceremesher.h"

//############################################################################# Find small rho triangles
/** Finds triangles with a low edge to circumradius ratio.
 *
\author Jan*/
void CHI_SURFACEREMESHER::FindSmallRhoTriangles(CHI_PATCH *patch, CST_FACE*& smallRhoTri,float criteria)
{
	smallRhoTri=NULL;
	//========================================================================= For each triangle in the patch
	for (int t=0; t<patch->faceList.itemCount; t++)
	{
		CST_FACE* tau_1 = patch->faceList.GetItem(t);
		if (tau_1!=NULL)
		{
			//============================================================ Calculate circumcircle
			float v_0[3];
			float r;
			
			FindCircumCircle(patch, tau_1,&v_0[0],&r);
			
			//============================================================ Find smallest edge and largest edge
			float* a = patch->Pstar.GetItem(tau_1->vertex[0]);
			float* b = patch->Pstar.GetItem(tau_1->vertex[1]);
			float* c = patch->Pstar.GetItem(tau_1->vertex[2]);
			
			Eigen::Vector3f va(a[0],a[1],a[2]);
			Eigen::Vector3f vb(b[0],b[1],b[2]);
			Eigen::Vector3f vc(c[0],c[1],c[2]);
			
			float length_ab = (vb-va).norm();
			float length_bc = (vc-vb).norm();
			float length_ca = (va-vc).norm();
			
			float smallestEdge = length_ab;
			float largestEdge = length_ab;

			if (length_bc<smallestEdge) {smallestEdge = length_bc;}
			if (length_ca<smallestEdge) {smallestEdge = length_ca;}

			if (length_bc>largestEdge) {largestEdge = length_bc;}
			if (length_ca>largestEdge) {largestEdge = length_ca;}

			//============================================================ Check if it meets the criteria
			if ((r/smallestEdge)>criteria)
			{
				//=================================================== If it does check whether the edge was seditious
				bool hasSeditiousEdge = false;
				bool isTooSmallToSplit = false;
				for (int e=0;e<3;e++)
				{
					if (tau_1->edges[e][2]==-1)
					{
						//====================================== Obtain the vertices of the edge
						float* v[3];
						v[0] = patch->Pstar.GetItem(tau_1->edges[e][0]);
						v[1] = patch->Pstar.GetItem(tau_1->edges[e][1]);

						//====================================== Find the index of the vertex not belonging to the edge
						int encrInd = -1;
						if (e<2){v[2] = patch->Pstar.GetItem(tau_1->edges[e+1][1]);encrInd=tau_1->edges[e+1][1];}
						else    {v[2] = patch->Pstar.GetItem(tau_1->edges[0][1]);encrInd=tau_1->edges[0][1];}

						//====================================== Calculate distances
						Eigen::Vector3f v_0(v[0][0],v[0][1],v[0][2]);
						Eigen::Vector3f v_1(v[1][0],v[1][1],v[1][2]);
						Eigen::Vector3f v_2(v[2][0],v[2][1],v[2][2]);

						Eigen::Vector3f v_c = 0.5*(v_0+v_1);
						float r = (v_1-v_c).norm();
						float d = (v_2-v_c).norm();

						//====================================== Check if it is encroached
						if (d<(r+this->precision))
						{
							//============================= If it is encroached check if it is splittable
							float splitLength = (v_1-v_0).norm()/2.0;
							if (splitLength<this->absoluteMinumumSize)
							{
								hasSeditiousEdge = true;
							}
						}
					}
					if (hasSeditiousEdge){break;}
				}

				//============================================= Check if too small
				if ((largestEdge/2)<this->absoluteMinumumSize)
				{
					isTooSmallToSplit = true;
				}

				//=================================================== If not then assign the triangle
				if ((!hasSeditiousEdge) && (!isTooSmallToSplit))
				{
					smallRhoTri = tau_1;
					//printf("Small rho tri %d->%d->%d\n",tau_1->vertex[0],tau_1->vertex[1],tau_1->vertex[2]);
					break;
				}


			}
		}
	}
}