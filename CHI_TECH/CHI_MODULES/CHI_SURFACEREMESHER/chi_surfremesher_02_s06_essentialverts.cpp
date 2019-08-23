#include "chi_surfaceremesher.h"

using namespace chiSurfaceMeshing;

//###################################################################List Essential Vertices
/** Populates the essential vertices for each patch.
\author Jan*/
bool CHI_SURFACEREMESHER::ListEssentialVertices(CHI_VECTOR<CHI_PATCH>* patchCollections)
{

	for (int p=0;p<patchCollections->itemCount;p++)
	{
		CHI_PATCH* curPatch = patchCollections->GetItem(p);

		//================================================= Loop through edge loops
		for (int loop=0;loop<curPatch->edgeLoops.itemCount;loop++)
		{
			CHI_EDGELOOP* curLoop = curPatch->edgeLoops.GetItem(loop);
            //=================================== Compare each edge to existing essential vertex list
			for (int e=0;e<curLoop->edges.itemCount;e++)
			{
				int* curEdge = curLoop->edges.GetItem(e);

                //===================== Compare the first vertex of every edge
                bool matchFound=false;
                for (int ev=0;ev<curPatch->essentialVertices.itemCount; ev++)
                {
                    int essVertIndex = *curPatch->essentialVertices.GetItem(ev);
                    if (curEdge[0]==essVertIndex)
                    {
                        matchFound=true;
                    }
                }
                if (!matchFound)
                {
                    int* newEssVertex = new int;
                    *newEssVertex = curEdge[0];
                    curPatch->essentialVertices.PushItem(newEssVertex);
                }

                //===================== If on the last edge then compare second vertex also
                if (e==(curLoop->edges.itemCount-1))
                {
                    bool matchFound=false;
                    for (int ev=0;ev<curPatch->essentialVertices.itemCount; ev++)
                    {
                        int essVertIndex = *curPatch->essentialVertices.GetItem(ev);
                        if (curEdge[1]==essVertIndex)
                        {
                            matchFound=true;
                        }
                    }
                    if (!matchFound)
                    {
                        int* newEssVertex = new int;
                        *newEssVertex = curEdge[1];
                        curPatch->essentialVertices.PushItem(newEssVertex);
                    }
                }


			}


		}
        printf("Patch %3d, has %3d essential vertices: ",p,curPatch->essentialVertices.itemCount);
        for (int ev=0;ev<curPatch->essentialVertices.itemCount; ev++)
        {
            printf(" %3d",*curPatch->essentialVertices.GetItem(ev));
        }
        printf("\n");
	}
	return false;
}