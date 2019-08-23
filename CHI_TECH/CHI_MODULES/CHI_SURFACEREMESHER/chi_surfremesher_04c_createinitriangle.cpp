#include"chi_surfaceremesher.h"

//############################################################################# Create initial triangle
/**Creates an initial right-hand rule triangle.

\author Jan*/
void CHI_SURFACEREMESHER::CreateInitialTriangle(CHI_PATCH *patch)
{
	CHI_VECTOR<int> *lexlist = &patch->lexicalListOfVertices;
	CHI_VECTOR<float> *projectedVertices = &patch->Pstar;
	//========================================================================= Initial hull point
	int* hullPoint1ref = lexlist->GetItem(0);
	int* hullPoint2ref = NULL;
	int* hullPoint3ref = NULL;
	if (hullPoint1ref==NULL)
	{
		printf("Error on hull point 1 in CreateInitialTriangle\n");
		return;
	}
	int    hullPoint1ind = *hullPoint1ref;
	float* hullPoint1Vrt = projectedVertices->GetItem(hullPoint1ind);
	//printf("First hull point=%d\n",hullPoint1ind);
	
	
	//========================================================================= Second hull point
	for (int i=1;i<lexlist->itemCount;i++)
	{
		int*   tempRef = lexlist->GetItem(i);
		float* tempVrt = projectedVertices->GetItem(*tempRef);
		//printf("tempRef hull point=%d %6.3f %6.3f %6.3f\n",*tempRef,tempVrt[0],tempVrt[1],tempVrt[2]);
		//printf("        hull point=%d %6.3f %6.3f %6.3f\n",*hullPoint1ref,hullPoint1Vrt[0],hullPoint1Vrt[1],hullPoint1Vrt[0] - tempVrt[0]);
		
		if (fabs(hullPoint1Vrt[0] - tempVrt[0])>this->precision)   //x-values not equal
		{
			hullPoint2ref = lexlist->GetItem(i);
			break;
		}
	}
	int    hullPoint2ind = *hullPoint2ref;
	//printf("Second hull point=%d\n",hullPoint2ind);
	float* hullPoint2Vrt = projectedVertices->GetItem(hullPoint2ind);

	
	//========================================================================= Third hull point
	for (int i=1;i<lexlist->itemCount;i++)
	{
		int*   tempRef = lexlist->GetItem(i);
		float* tempVrt = projectedVertices->GetItem(*tempRef);
		
		if (*tempRef!=hullPoint2ind)   //not hullpoint2
		{
			hullPoint3ref = lexlist->GetItem(i);
			break;
		}
	}
	int    hullPoint3ind = *hullPoint3ref;
	float* hullPoint3Vrt = projectedVertices->GetItem(hullPoint3ind);
	//printf("Third hull point=%d\n",hullPoint3ind);
	
	//========================================================================= Create new face
	CST_FACE* newFace = new CST_FACE;
	float* newNormal = new float[3];
	newNormal[0] = patch->normal(0);
	newNormal[1] = patch->normal(1);
	newNormal[2] = patch->normal(2);
	int newNormalIndex = remeshedSurface->normalStack.PushItem(newNormal);
	newFace->normal[0]=newNormalIndex;
	newFace->normal[1]=newNormalIndex;
	newFace->normal[2]=newNormalIndex;
	
	//========================================================================= Check right-hand rule
	double orientation = Orient2D(hullPoint1Vrt,hullPoint2Vrt,hullPoint3Vrt);
	if (orientation>0)
	{
		newFace->vertex[0] = hullPoint1ind;
		newFace->vertex[1] = hullPoint2ind;
		newFace->vertex[2] = hullPoint3ind;
	}
	else
	{
		newFace->vertex[0] = hullPoint1ind;
		newFace->vertex[1] = hullPoint3ind;
		newFace->vertex[2] = hullPoint2ind;
	}
	newFace->edges[0][0] = newFace->vertex[0]; newFace->edges[0][1] = newFace->vertex[1];
	newFace->edges[1][0] = newFace->vertex[1]; newFace->edges[1][1] = newFace->vertex[2];
	newFace->edges[2][0] = newFace->vertex[2]; newFace->edges[2][1] = newFace->vertex[0];
	//remeshedSurface->faceStack.PushItem(newFace);
	patch->faceList.PushItem(newFace);
	
	//printf("Initial Triangle %d->%d->%d\n",newFace->vertex[0],newFace->vertex[1],newFace->vertex[2]);
	
}