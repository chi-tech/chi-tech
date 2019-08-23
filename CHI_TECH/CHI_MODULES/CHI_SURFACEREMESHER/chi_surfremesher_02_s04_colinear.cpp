#include "chi_surfaceremesher.h"

using namespace chiSurfaceMeshing;

//################################################################### Find colinear edges
/** Builds a list of colinear edges from a list of open edges.
\author Jan*/
bool CHI_SURFACEREMESHER::FindCoLinearEdges(CHI_VECTOR<CHI_PATCH>* patchCollections)
{
	for (int p=0; p<patchCollections->itemCount; p++)
	{
		CHI_PATCH *curPatch = patchCollections->GetItem(p);
		
		//================================================= Add all open edges to unused cache
		CHI_VECTOR<int> unusedEdges;
		for (int e=0;e<curPatch->edges.itemCount;e++)
		{
			int* curEdge = curPatch->edges.GetItem(e);
			unusedEdges.PushItem(curEdge);
		}
		
		//================================================= Seed the first edge
		if (curPatch->colinearEdges.itemCount==0)
		{
			int* firstEdge = unusedEdges.PullItem(unusedEdges.itemCount-1);
			CHI_COLINEDGELIST* newColinEdge = new CHI_COLINEDGELIST;
			newColinEdge->edges.PushItem(firstEdge);
			curPatch->colinearEdges.PushItem(newColinEdge);
		}
		
		//================================================= Recursively use the unused edges
		bool stopRecursive = false;
		int iter = 0;
		while (!stopRecursive)
		{
			//=================================== Run through each list
			bool matchFound = false;
			for (int c=0;c<curPatch->colinearEdges.itemCount;c++)
			{
				CHI_COLINEDGELIST* curList = curPatch->colinearEdges.GetItem(c);
				
				//===================== Run through each edge of the list as master
				for (int e=0;e<curList->edges.itemCount;e++)
				{
					int* masterEdge = curList->edges.GetItem(e);
					for (int e2=unusedEdges.itemCount-1;e2>=0;e2--)
					{
						int* slaveEdge = unusedEdges.GetItem(e2);
						if (CheckEdgesColinear(masterEdge,slaveEdge,initialMesh))
						{
							curList->edges.PushItem(slaveEdge);
							unusedEdges.PullItem(e2);
							matchFound = true;
						}
					}
				}
				
			}
			iter++;
			if ((unusedEdges.itemCount==0) || (iter>1000))
			{
				stopRecursive = true;
			}
			else
			{
				if (!matchFound)
				{
					int* firstEdge = unusedEdges.PullItem(unusedEdges.itemCount-1);
					CHI_COLINEDGELIST* newColinEdge = new CHI_COLINEDGELIST;
					newColinEdge->edges.PushItem(firstEdge);
					curPatch->colinearEdges.PushItem(newColinEdge);
				}
			}
		}
		
		//printf("Patch %3d has %3d co-linear edges.\n",p,curPatch->colinearEdges.itemCount);
		
		ResolveEdgeExtent(&curPatch->colinearEdges,initialMesh);
		for (int co=0;co<curPatch->colinearEdges.itemCount;co++)
		{
			CHI_COLINEDGELIST* curList = curPatch->colinearEdges.GetItem(co);
			//printf("     -> Edge %3d goes from vertex [%3d] to vertex [%3d] with length %6.3f\n", co,curList->extentVertices[0],curList->extentVertices[1],curList->extentLength);
		}
	}

	return false;
}