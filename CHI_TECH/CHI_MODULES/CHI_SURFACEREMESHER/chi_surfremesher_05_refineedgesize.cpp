#include "chi_surfaceremesher.h"

//############################################################################# Refine edge sizes
/**Refines interior and/or outer edges sizes.

\return Returns true if an edge was refined.

\author Jan*/
bool CHI_SURFACEREMESHER::RefineEdgeSize(CHI_PATCH *patch)
{
    //========================================================================= For each triangle of the patch
    for (int t=0; t<patch->faceList.itemCount; t++)
    {
        CST_FACE* tau_1 = patch->faceList.GetItem(t);

        //=========================================================== Loop through all the edges
        for (int e=0; e<3; e++)
        {
            //=============================== Obtain the vertices of the edge
            float* v[2];
            v[0] = patch->Pstar.GetItem(tau_1->edges[e][0]);
            v[1] = patch->Pstar.GetItem(tau_1->edges[e][1]);

            //=============================== Calculate distances
            Eigen::Vector3f v_0(v[0][0],v[0][1],v[0][2]);
            Eigen::Vector3f v_1(v[1][0],v[1][1],v[1][2]);

            Eigen::Vector3f v01 = v_1 - v_0;
            float edgeLength = v01.norm();

            if (edgeLength>this->baseSize)
            {
                if ((edgeLength/2.0) > this->absoluteMinumumSize)
                {
                    Eigen::Vector3f vc = 0.5*(v_0+v_1);

                    //======================= Created project vertex
                    float newVert[3];
                    Eigen::Vector3f Reproject(0,0,0);
                    Reproject += vc(0)*patch->hat_i;
                    Reproject += vc(1)*patch->hat_j;
                    //Reproject += vc(2)*patch->hat_k;
                    Reproject += patch->C_0;
                    newVert[0] = Reproject(0);
                    newVert[1] = Reproject(1);
                    newVert[2] = Reproject(2);

                    InsertVertex(patch,newVert);
                    return true;
                }
            }


        }

    }

    return false;
}
