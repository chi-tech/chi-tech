#include"chi_surfaceremesher.h"

//############################################################################# Restore 2D triangulation
/** Restores a Delaunay triangulation of a patch.
\author Jan*/
void CHI_SURFACEREMESHER::Make2DDelaunayTriangulation(CHI_PATCH* patch)
{
	
		//================================================= Project 3D vertices to 2D
		CHI_PATCH* curPatch = patch;
		
		//================================================= List non-locally delaunay edges
		//printf("        -Creating initial list of non-local Delaunay edges");
		CHI_VECTOR<CHI_INTERIOREDGE> non_local_del_edges;
		ListNonLocallyDelaunayEdges(curPatch,&non_local_del_edges);
		//printf("..Done\n");
		
		//================================================= Iterate to remove non-locally delaunay edges
		CHI_VECTOR<CHI_OPENEDGE> hullEdges;
		//FindPatchOpenEdges(curPatch,&hullEdges);
		int initialOpenEdgeCount = hullEdges.itemCount;
		int iter=0;
		bool forcestop=false;
		while ((non_local_del_edges.itemCount>0) && (!forcestop))
		{
			iter++;
			//printf("        -Edge flip iteration %d\n",iter);
			EdgeFlip(curPatch,non_local_del_edges.GetItem(0));
			
			ListNonLocallyDelaunayEdges(curPatch,&non_local_del_edges);
			
			//FindPatchOpenEdges(curPatch,&hullEdges);
			
			if (hullEdges.itemCount!= initialOpenEdgeCount)
			{
				printf("Error in Delaunay generation\n");
				forcestop=true;
			}
			//if (iter>1000) {forcestop=true;}
		}

}