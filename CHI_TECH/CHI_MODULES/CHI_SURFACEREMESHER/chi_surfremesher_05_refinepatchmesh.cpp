#include"chi_surfaceremesher.h"



//############################################################################# Ruppert's Delaunay refinement
/** Refines a 2D mesh according to Ruppert's Delaunay refinement algorithm.
\author Jan*/
void CHI_SURFACEREMESHER::RefinePatchMesh(CHI_PATCH *patch)
{
	//DumpPatchToScilab(patch);
	int iter=0;
	bool largeEdgeFound=true;
	while (largeEdgeFound)
	{
		iter++;
		printf("Edge refinement iteration %3d, triangle count=%5d\n",iter,patch->faceList.itemCount);
		RestoreOptimality(patch);
		largeEdgeFound = RefineEdgeSize(patch);

		//if(iter++>=100){largeEdgeFound=false;}
	}


	DumpPatchToScilab(patch);
	
}