#include "chi_surfaceremesher.h"
using namespace chiSurfaceMeshing;

//###################################################################
/**
\author Jan*/
bool CHI_SURFACEREMESHER::ExecuteMeshing(CHI_SURFACE* initialMesh)
{
    printf("################################################ Executing Surface Remesher\n");
	bool success = true;

	this->initialMesh = initialMesh;
    CHI_VECTOR<CHI_FACELIST> coPlanarFaceCollections;
    CHI_VECTOR<CHI_PATCH> patchCollections;
	
	printf("Collecting co-planar faces\n");
	CollectCoPlanarFaces(&coPlanarFaceCollections);
	printf("Splitting by patch\n");
    CreatePatchList(&coPlanarFaceCollections,&patchCollections);
	printf("Finding open edges\n");
    FindOpenEdges(&patchCollections);
	printf("Find co-linear edges\n");
	FindCoLinearEdges(&patchCollections);
	printf("Find edge loops\n");
	FindEdgeLoops(&patchCollections);
	printf("Find essential vertices\n");
    ListEssentialVertices(&patchCollections);
	
	//printf("################################################ Starting Delaunay triangulation\n");
	this->remeshedSurface = new CHI_REMESHEDSURFACE;
	CopyInformationToRemeshedSurface(&patchCollections, remeshedSurface);


	Create2DDelaunayTriangulation(remeshedSurface);


	for (int p=0; p<remeshedSurface->patches.itemCount; p++)
	//for (int p=0; p<2; p++)
	{
		CHI_PATCH* curPatch = remeshedSurface->patches.GetItem(p);
		RefinePatchMesh(curPatch);
	}


    printf("Assign faces to remeshed surface\n");
	for (int p=0; p<remeshedSurface->patches.itemCount; p++)
	{
		CHI_PATCH* curPatch = remeshedSurface->patches.GetItem(p);
		for (int f=0; f<curPatch->faceList.itemCount; f++)
		{
			CST_FACE* curFace = curPatch->faceList.GetItem(f);
			if (curFace!=NULL)
			{
				int newIndex = remeshedSurface->faceStack.PushItem(curFace);
				GLfloat* newNormal = new GLfloat[3];
				newNormal[0] = curPatch->normal(0);
				newNormal[1] = curPatch->normal(1);
				newNormal[2] = curPatch->normal(2);
				int normIndex = remeshedSurface->normalStack.PushItem(newNormal);
				curFace->normal[0] = normIndex;
				curFace->normal[1] = normIndex;
				curFace->normal[2] = normIndex;
				//printf("Face %3d, %3d %3d %3d\n",newIndex,curFace->vertex[0],curFace->vertex[1],curFace->vertex[2]);
			}

		}
	}
	for (int p=0;p<remeshedSurface->vertexStack.itemCount;p++)
	{
		float* v = remeshedSurface->vertexStack.GetItem(p);

		//printf("Vertex %3d, %6.3f %6.3f %6.3f\n",p,v[0],v[1],v[2]);
	}
	

    printf("################################################ Surface Remesher Stopped\n");
	return false;
}