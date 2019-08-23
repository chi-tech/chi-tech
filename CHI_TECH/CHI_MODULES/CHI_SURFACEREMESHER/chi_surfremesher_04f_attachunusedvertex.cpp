#include"chi_surfaceremesher.h"

//############################################################################# Attach unused vertex
/** Takes the first unused vertex and attempts to attach it to the existing faces.
 *
\author Jan*/
void CHI_SURFACEREMESHER::AttachUnusedVertex(CHI_PATCH *patch, CHI_VECTOR<int> *unusedVertices)
{
	CHI_VECTOR<int> *lexlist = &patch->lexicalListOfVertices;
	CHI_VECTOR<float> *projectedVertices = &patch->Pstar;
	
	//===================================================== Get lexical index of unused vertex
	int* unusedVertLexIndex = unusedVertices->GetItem(0);
	if (unusedVertLexIndex==NULL){return;}
	
	//===================================================== Get physical index of unused vertex
	int* vertexIndex = lexlist->GetItem(*unusedVertLexIndex);
	if (vertexIndex==NULL){return;}
	
	//===================================================== Get vertex xyz-value
	float* unusedVertex = projectedVertices->GetItem(*vertexIndex);
	
	//===================================================== Get all open edge
	CHI_VECTOR<CHI_OPENEDGE> hullEdges;
	FindPatchOpenEdges(patch,&hullEdges);
	
	//===================================================== Attempt to attach unused vertex to any edge
	//printf("Checkpoint 1 - find edge for lexVert %d\n", *vertexIndex);
	for (int e=0; e<hullEdges.itemCount; e++)
	{
		CHI_OPENEDGE* curEdge = hullEdges.GetItem(e);
		
		float* a = unusedVertex;
		float* b = projectedVertices->GetItem(curEdge->vi);
		float* c = projectedVertices->GetItem(curEdge->vf);
		
		double orientation = Orient2D(a,b,c);
		
		if (orientation<(-1*this->precision))
		{
			//=================================== Check if we will cross simplices
			if (!FindIntersectingSimplices(patch,curEdge->vi,*vertexIndex,curEdge->vf))
            {
                //=================================== Tell the owner of the edge that he now has a neighbor
                CST_FACE* ownerTriangle = patch->faceList.GetItem(curEdge->ownerTriangle);
                ownerTriangle->edges[curEdge->ownerEdgeNumber][2] = patch->faceList.itemCount;

                //=================================== Tell the new triangle which triangle its attaching to
                CST_FACE* newTriangle = new CST_FACE;
                newTriangle->normal[0] = ownerTriangle->normal[0];
                newTriangle->normal[1] = ownerTriangle->normal[1];
                newTriangle->normal[2] = ownerTriangle->normal[2];
                newTriangle->faceNormal[0] = ownerTriangle->faceNormal[0];
                newTriangle->faceNormal[1] = ownerTriangle->faceNormal[1];
                newTriangle->faceNormal[2] = ownerTriangle->faceNormal[2];

                newTriangle->vertex[0] = curEdge->vi;
                newTriangle->vertex[1] = *vertexIndex;
                newTriangle->vertex[2] = curEdge->vf;

                //printf("New triangle from lexica %d->%d->%d\n",newTriangle->vertex[0],newTriangle->vertex[1],newTriangle->vertex[2]);

                newTriangle->edges[0][0] = newTriangle->vertex[0];  newTriangle->edges[0][1] = newTriangle->vertex[1];
                newTriangle->edges[1][0] = newTriangle->vertex[1];  newTriangle->edges[1][1] = newTriangle->vertex[2];
                newTriangle->edges[2][0] = newTriangle->vertex[2];  newTriangle->edges[2][1] = newTriangle->vertex[0];

                newTriangle->edges[2][2] = curEdge->ownerTriangle;
                //remeshedSurface->faceStack.PushItem(newTriangle);
                patch->faceList.PushItem(newTriangle);
                break;
            }
			//printf("Checkpoint 2\n");

		}
	}
}