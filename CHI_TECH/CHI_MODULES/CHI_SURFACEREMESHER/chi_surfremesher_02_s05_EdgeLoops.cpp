#include"chi_surfaceremesher.h"

using namespace chiSurfaceMeshing;


//################################################################### Find edge loops
/** Builds a list of edge loops from the colinear edges.
\author Jan*/
bool CHI_SURFACEREMESHER::FindEdgeLoops(CHI_VECTOR<CHI_PATCH> *patchCollections)
{
	//=============================================================== Loop through all patches
	for (int p=0;p<patchCollections->itemCount;p++)
	{
		CHI_PATCH* curPatch = patchCollections->GetItem(p);
		
		//================================================= Initialize unusedEdges
		CHI_VECTOR<CHI_COLINEDGELIST> unusedEdges;
		for (int e=0;e<curPatch->colinearEdges.itemCount;e++)
		{
			CHI_COLINEDGELIST* curEdge = curPatch->colinearEdges.GetItem(e);
			unusedEdges.PushItem(curEdge);
		}
		
		//================================================= Seed first edgeloop
		CHI_EDGELOOP* newEdgeLoop = new CHI_EDGELOOP;
		CHI_COLINEDGELIST* firstEdge = unusedEdges.PullItem(unusedEdges.itemCount-1);
		int* initialEdge = new int[2];
		
		initialEdge[0] = firstEdge->extentVertices[0];
		initialEdge[1] = firstEdge->extentVertices[1];
		newEdgeLoop->edges.PushItem(initialEdge);
		
		newEdgeLoop->initialVertex = firstEdge->extentVertices[0];
		newEdgeLoop->finalVertex   = firstEdge->extentVertices[1];
		curPatch->edgeLoops.PushItem(newEdgeLoop);
		//printf("Patch %3d, iVert= %3d, fVert= %3d\n",p,newEdgeLoop->initialVertex,newEdgeLoop->finalVertex);
		//================================================== Find a place for each edge
		bool stopRecursive = false;
		int iter=0;
		while (!stopRecursive)
		{
			bool updateMade=false;
			
			for (int el=0;el<curPatch->edgeLoops.itemCount;el++)
			{
				CHI_EDGELOOP* curLoop = curPatch->edgeLoops.GetItem(el);
				//printf("Loop %3d, iVert= %3d, fVert= %3d\n",el,curLoop->initialVertex,curLoop->finalVertex);
				
				for (int e=unusedEdges.itemCount-1;e>=0;e--)
				{
					CHI_COLINEDGELIST *curEdge = unusedEdges.GetItem(e);
					//printf("Unused edge %3d,  %3d -> %3d\n",e,curEdge->extentVertices[0],curEdge->extentVertices[1]);
					if (curEdge->extentVertices[0] == curLoop->finalVertex)
					{
						int *newEdge = new int[2];
						newEdge[0] = curEdge->extentVertices[0];
						newEdge[1] = curEdge->extentVertices[1];
						curLoop->finalVertex = curEdge->extentVertices[1];
						unusedEdges.PullItem(e);
						curLoop->edges.PushItem(newEdge);
						updateMade = true;
					} else if (curEdge->extentVertices[1] == curLoop->finalVertex)
					{
						int *newEdge = new int[2];
						newEdge[0] = curEdge->extentVertices[1];
						newEdge[1] = curEdge->extentVertices[0];
						curLoop->finalVertex = curEdge->extentVertices[0];
						unusedEdges.PullItem(e);
						curLoop->edges.PushItem(newEdge);
						updateMade = true;
					}
				}

			}
			
			iter++;
			//printf("Iteration %3d, loops found %3d\n", iter, curPatch->edgeLoops.itemCount);
			if ((iter>100) || (unusedEdges.itemCount==0))
			{
				stopRecursive=true;
			}
			else
			{
				if (!updateMade)
				{
					newEdgeLoop = new CHI_EDGELOOP;
					firstEdge = unusedEdges.PullItem(unusedEdges.itemCount-1);
					initialEdge = new int[2];
					
					initialEdge[0] = firstEdge->extentVertices[0];
					initialEdge[1] = firstEdge->extentVertices[1];
					newEdgeLoop->edges.PushItem(initialEdge);
					
					newEdgeLoop->initialVertex = firstEdge->extentVertices[0];
					newEdgeLoop->finalVertex   = firstEdge->extentVertices[1];
					curPatch->edgeLoops.PushItem(newEdgeLoop);
				}
			}
		}
		//printf("Patch %3d, has %3d edge loops: \n",p,curPatch->edgeLoops.itemCount);
		for (int el=0;el<curPatch->edgeLoops.itemCount;el++)
		{
			//printf("   Loop %3d:\n", el);
			CHI_EDGELOOP* curLoop = curPatch->edgeLoops.GetItem(el);
			for (int e=0;e<curLoop->edges.itemCount;e++)
			{
				int* curEdge = curLoop->edges.GetItem(e);
				//printf("%3d:     %3d -> %3d\n", e, curEdge[0],curEdge[1]);
			}
		}
	}
	return false;
}


//################################################################### Overload FindEdgeLoop
/** Finds the first edge loop from a closed hull.
\author Jan*/
bool CHI_SURFACEREMESHER::FindEdgeLoops(CHI_PATCH* patch)
{
	CHI_PATCH* curPatch = patch;
	patch->edgeLoops.ClearVector();
	
	//===================================================== Seed the first edge
	//printf("Finding Edgeloops from %d open edges\n",curPatch->edges.itemCount);
	int* firstEdge = patch->edges.GetItem(0);
	CHI_EDGELOOP* seedLoop = new CHI_EDGELOOP;
	seedLoop->edges.PushItem(firstEdge);
	seedLoop->initialVertex = firstEdge[0];
	seedLoop->finalVertex   = firstEdge[1];
	
	patch->edgeLoops.ClearVector();
	patch->edgeLoops.PushItem(seedLoop);
	//printf("Checkpoint 1 firstEdge %d->%d\n",firstEdge[0],firstEdge[1]);
	//===================================================== Populate unused edges
	CHI_VECTOR<int> unusedEdges;
	for (int e=1;e<patch->edges.itemCount;e++)
	{
		int* curEdge = patch->edges.GetItem(e);
		unusedEdges.PushItem(curEdge);
	}
	//printf("Checkpoint 2\n");
	bool stopRecursive = false;
	int iter=0;
	while (!stopRecursive)
	{
		//printf("Unused edges: %d\n",unusedEdges.itemCount);
		bool updateMade=false;
		
		for (int e=unusedEdges.itemCount-1;e>=0;e--)
		{
			int* curEdge = unusedEdges.GetItem(e);
			//printf("    curEdge %d->%d\n",curEdge[0],curEdge[1]);
			for (int el=0;el<patch->edgeLoops.itemCount;el++)
			{
				CHI_EDGELOOP* curLoop = patch->edgeLoops.GetItem(el);
				if (curLoop!=NULL)
				{
					if (curEdge[0]==curLoop->finalVertex)
					{
						curEdge = unusedEdges.PullItem(e);
						curLoop->edges.PushItem(curEdge);
						curLoop->finalVertex = curEdge[1];
						updateMade=true;
						//printf("       pushed\n");


					}
				}
			}

		}

		if ((!updateMade) && (unusedEdges.itemCount>0))
		{
			//printf("Creating new loop\n");
			CHI_EDGELOOP* newLoop = new CHI_EDGELOOP;
			int* curEdge = unusedEdges.PullItem(unusedEdges.itemCount-1);
			newLoop->edges.PushItem(curEdge);
			newLoop->finalVertex = curEdge[1];
			patch->edgeLoops.PushItem(newLoop);
		}
		
		if (unusedEdges.itemCount==0)
		{
			stopRecursive=true;
		}
	}

	return false;
}