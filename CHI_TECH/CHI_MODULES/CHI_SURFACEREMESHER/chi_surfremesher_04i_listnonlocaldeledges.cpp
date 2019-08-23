#include"chi_surfaceremesher.h"

void CHI_SURFACEREMESHER::ListNonLocallyDelaunayEdges(CHI_PATCH *patch,
                                                      CHI_VECTOR<CHI_INTERIOREDGE> *non_loc_del_edges)
{
	CHI_VECTOR<float> *projectedVertices = &patch->Pstar;
	non_loc_del_edges->ClearVector();
	
	//===================================================== Run through each triangle
	for (int t=0; t<patch->faceList.itemCount; t++)
	{
		//============================================ For each edge of the triangle
		CST_FACE* tau_1 = patch->faceList.GetItem(t);
		if (tau_1!=NULL)
		{
			for (int e=0;e<3;e++)
			{
				if (tau_1->edges[e][2]>-1) //Check if its an interior edge
				{
					//=============================== Assign neighboring triangle
					CST_FACE* tau_2 = patch->faceList.GetItem(tau_1->edges[e][2]);
					if (tau_2==NULL)
					{
						//printf("Error: Triangle2 not found ListNonLocallyDelaunayEdges\n");
						return;
					}

					//=============================== Find indexes for abc
					int aindex = tau_1->vertex[0];
					int bindex = tau_1->vertex[1];
					int cindex = tau_1->vertex[2];
					int dindex;

					//=============================== Find the index of d (the one that is not a b or c)
					for (int v=0;v<3;v++)
					{
						if ((tau_2->vertex[v]!=tau_1->edges[e][0]) && (tau_2->vertex[v]!=tau_1->edges[e][1]))
						{
							dindex = tau_2->vertex[v]; break;
						}
					}

					//=============================== Get the vertex values of abcd
					float* a = projectedVertices->GetItem(aindex);
					float* b = projectedVertices->GetItem(bindex);
					float* c = projectedVertices->GetItem(cindex);
					float* d = projectedVertices->GetItem(dindex);

					//=============================== Check if locally delaunay
					double localDelaunay = InCircle(a,b,c,d);
					if (localDelaunay>0.0001)
					{
						//====================== Check for duplicate
						bool duplicateFound = false;
						for (int ie=0; ie<non_loc_del_edges->itemCount; ie++)
						{
							CHI_INTERIOREDGE* curEdge = non_loc_del_edges->GetItem(ie);
							/*printf("Dup check %d->%d vs %d->%d\n",
                                    curEdge->vi,curEdge->vf,
                                    tau_1->edges[e][1],tau_1->edges[e][0]);*/
							if ( ((curEdge->vi==tau_1->edges[e][1]) && (curEdge->vf==tau_1->edges[e][0])) ||
								 ((curEdge->vf==tau_1->edges[e][1]) && (curEdge->vi==tau_1->edges[e][0])) )
							{
								duplicateFound = true; break;
							}
						}

						//====================== If no duplicate found then add edge to list
						if (!duplicateFound)
						{
							CHI_INTERIOREDGE* newEdge = new CHI_INTERIOREDGE;
							newEdge->vi = tau_1->edges[e][0];
							newEdge->vf = tau_1->edges[e][1];
							newEdge->tau_1 = t;
							newEdge->tau_1_edgeNumber = e;
							newEdge->tau_2 = tau_1->edges[e][2];

							for (int e2=0;e2<3; e2++)
							{
								if ((tau_2->edges[e2][0]==newEdge->vf) && (tau_2->edges[e2][1]==newEdge->vi))
								{
									newEdge->tau_2_edgeNumber = e2; break;
								}
							}
							if (newEdge->tau_2_edgeNumber<0)
							{
								printf("Error finding tau_2_edgenumber\n");
							}
							//printf("Non-delaunay interior edge from vertex %3d to %3d, tau1=%3d(%3d), tau2=%3d(%3d) %.9f\n",
                            //        newEdge->vi, newEdge->vf,
                            //        newEdge->tau_1, newEdge->tau_1_edgeNumber,
                            //        newEdge->tau_2, newEdge->tau_2_edgeNumber,localDelaunay);
							non_loc_del_edges->PushItem(newEdge);
						}
					}
				}
			}
		}

	}
}