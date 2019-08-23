#include"chi_surfaceremesher.h"

//############################################################################# Convexify hull
/** Checks  a hull edges for convexity and adds a triangle if needed.
\author Jan*/
void CHI_SURFACEREMESHER::ConvexifyHull(CHI_PATCH *patch)
{
	CHI_VECTOR<int> *lexlist = &patch->lexicalListOfVertices;
	CHI_VECTOR<float> *projectedVertices = &patch->Pstar;
	
	//===================================================== Get all open edge
	CHI_VECTOR<CHI_OPENEDGE> hullEdges;
	FindPatchOpenEdges(patch,&hullEdges);
	
	bool concaveEdgesFound = true;
	while (concaveEdgesFound)
	{
		concaveEdgesFound = false;
		for (int e = 0; e < hullEdges.itemCount; e++)
		{
			CHI_OPENEDGE *edge1 = hullEdges.GetItem(e);
			CHI_OPENEDGE *edge2;
			if (e < (hullEdges.itemCount - 1))
			{
				edge2 = hullEdges.GetItem(e + 1);
			} else
			{
				edge2 = hullEdges.GetItem(0);
			}
			
			int aind = edge2->vf;
			int bind = edge1->vi;
			int cind = edge1->vf;
			
			float *a = projectedVertices->GetItem(aind);
			float *b = projectedVertices->GetItem(bind);
			float *c = projectedVertices->GetItem(cind);
			
			double orientation = Orient2D(a,b,c);
			//printf("Convexify ori=%10.7f\n",orientation);
			if (orientation<(-1*this->precision))
			{
                if (!FindIntersectingSimplices(patch,edge1->vi,edge2->vf,edge1->vf))
                {
                    //printf("Found Concave Edges\n");
                    concaveEdgesFound=true;
                    //=================================== Tell the owners of the edges that they now have a neighbor
                    CST_FACE* ownerTriangle = patch->faceList.GetItem(edge1->ownerTriangle);
                    ownerTriangle->edges[edge1->ownerEdgeNumber][2] = patch->faceList.itemCount;

                    ownerTriangle = patch->faceList.GetItem(edge2->ownerTriangle);
                    ownerTriangle->edges[edge2->ownerEdgeNumber][2] = patch->faceList.itemCount;

                    //=================================== Tell the new triangle which triangle its attaching to
                    CST_FACE* newTriangle = new CST_FACE;
                    newTriangle->normal[0] = ownerTriangle->normal[0];
                    newTriangle->normal[1] = ownerTriangle->normal[1];
                    newTriangle->normal[2] = ownerTriangle->normal[2];
                    newTriangle->faceNormal[0] = ownerTriangle->faceNormal[0];
                    newTriangle->faceNormal[1] = ownerTriangle->faceNormal[1];
                    newTriangle->faceNormal[2] = ownerTriangle->faceNormal[2];

                    newTriangle->vertex[0] = edge1->vi;
                    newTriangle->vertex[1] = edge2->vf;
                    newTriangle->vertex[2] = edge1->vf;

                    //printf("New triangle from convex %d->%d->%d\n",newTriangle->vertex[0],newTriangle->vertex[1],newTriangle->vertex[2]);

                    newTriangle->edges[0][0] = newTriangle->vertex[0];  newTriangle->edges[0][1] = newTriangle->vertex[1];
                    newTriangle->edges[1][0] = newTriangle->vertex[1];  newTriangle->edges[1][1] = newTriangle->vertex[2];
                    newTriangle->edges[2][0] = newTriangle->vertex[2];  newTriangle->edges[2][1] = newTriangle->vertex[0];

                    newTriangle->edges[1][2] = edge2->ownerTriangle;
                    newTriangle->edges[2][2] = edge1->ownerTriangle;
                    //remeshedSurface->faceStack.PushItem(newTriangle);
                    patch->faceList.PushItem(newTriangle);
                    break;
                }

			}
			
		}
		if (concaveEdgesFound)
		{
			FindPatchOpenEdges(patch,&hullEdges);
		}
	}
}