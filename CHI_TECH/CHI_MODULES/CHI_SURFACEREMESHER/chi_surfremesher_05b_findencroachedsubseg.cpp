#include"chi_surfaceremesher.h"

//############################################################################# Find Encroached edge
/** Finds an encroached subsegment.
\author Jan*/
void CHI_SURFACEREMESHER::FindEncroachedSubSegments(CHI_PATCH *patch, CHI_INTERIOREDGE *&encrSeg, bool excludeNonsplittable)
{
	encrSeg=NULL;
	//========================================================================= For each triangle of the patch
	for (int t=0; t<patch->faceList.itemCount; t++)
	{
		CST_FACE* tau_1 = patch->faceList.GetItem(t);
		if (tau_1!=NULL)
		{
			bool encroachedFound = false;
			//=========================================================== Loop through all the edges
			for (int e=0; e<3; e++)
			{

				//============================================= If the edge is open
				if (tau_1->edges[e][2]==-1)
				{
					//=============================== Obtain the vertices of the edge
					float* v[3];
					v[0] = patch->Pstar.GetItem(tau_1->edges[e][0]);
					v[1] = patch->Pstar.GetItem(tau_1->edges[e][1]);

					//=============================== Find the index of the vertex not belonging to the edge
					int encrInd = -1;
					if (e<2){v[2] = patch->Pstar.GetItem(tau_1->edges[e+1][1]);encrInd=tau_1->edges[e+1][1];}
					else    {v[2] = patch->Pstar.GetItem(tau_1->edges[0][1]);encrInd=tau_1->edges[0][1];}

					//=============================== Calculate distances
					Eigen::Vector3f v_0(v[0][0],v[0][1],v[0][2]);
					Eigen::Vector3f v_1(v[1][0],v[1][1],v[1][2]);
					Eigen::Vector3f v_2(v[2][0],v[2][1],v[2][2]);

					Eigen::Vector3f v_c = 0.5*(v_0+v_1);
					float r = (v_1-v_c).norm();
					float d = (v_2-v_c).norm();

					//=============================== Check if minimum size needs to be enforced
					float splitLength = (v_1-v_0).norm()/2.0;
					bool continueFlag = true;

					if ((excludeNonsplittable) && (splitLength<this->absoluteMinumumSize))
					{
						continueFlag = false;
					}

					//=============================== Determine if the edge is encroached
					if ((d<(r+this->precision)) && (continueFlag))
					{
						//CST_FACE* tau_2 = patch->faceList.GetItem(tau_1->edges[e][2]);

						//================= Fill the edge structure
						encrSeg = new CHI_INTERIOREDGE;
						encrSeg->vi = tau_1->edges[e][0];
						encrSeg->vf = tau_1->edges[e][1];
						encrSeg->tau_1 = t;
						encrSeg->tau_1_edgeNumber = e;
						encrSeg->encroachingVertex = encrInd;

						/*
                        if (tau_2!=NULL)
                        {
                            encrSeg->tau_2 = tau_1->edges[e][2];

                            for (int e2=0;e2<3; e2++)
                            {
                                if ((tau_2->edges[e2][0]==encrSeg->vf) && (tau_2->edges[e2][1]==encrSeg->vi))
                                {
                                    encrSeg->tau_2_edgeNumber = e2; break;
                                }
                            }
                        }
                        else
                        {
                            encrSeg->tau_2 = -1;
                            encrSeg->tau_2_edgeNumber = -1;
                        }*/
						//printf("Encroached Edge Found, %d %d c[%6.3f,%6.3f]\n",encrSeg->vi,encrSeg->vf,v_c(0),v_c(1));
						encroachedFound=true;
						break;
					}
				}





			}
			if (encroachedFound)
			{
				break;
			}
		}

	}
}