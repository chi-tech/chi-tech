#include"chi_surfaceremesher.h"

//############################################################################# Create 2D triangulation
/** Creates a Delaunay triangulation of a surface mesh.
\author jan*/
void CHI_SURFACEREMESHER::Create2DDelaunayTriangulation(CHI_REMESHEDSURFACE *mesh)
{
	for (int p=0;p<mesh->patches.itemCount;p++)
	{
		//================================================= Project 3D vertices to 2D
		CHI_PATCH* curPatch = mesh->patches.GetItem(p);
		Project3DVerticesTo2D(mesh,curPatch);
		
		//================================================= Sort the vertices lexicographically and create first triangle
		SortLexicographically2D(curPatch);
		//printf("Creating initial triangle\n");
		CreateInitialTriangle(curPatch);

		//================================================= Obtain list of unknown vertices
		CHI_VECTOR<int> unusedVertices;
		GetUnusedVertices(curPatch,&unusedVertices);
		
		//================================================= Iterate until all vertices are used
		int iter=0;
		bool forcestop=false;
		while ((unusedVertices.itemCount>0) && (!forcestop))
		{
			AttachUnusedVertex(curPatch,&unusedVertices);

			ConvexifyHull(curPatch);
			GetUnusedVertices(curPatch,&unusedVertices);
			iter++;
			//if (iter>=2) {forcestop=true;}
		}
		ConvexifyHull(curPatch);

		if (curPatch->edgeLoops.itemCount>1)
		{
			RemoveFillerTriangles(curPatch);
		}

		FindOpenEdges(curPatch);
		FindEdgeLoops(curPatch);
		DumpPatchToScilab(curPatch);



		//================================================= List non-locally delaunay edges
		//printf("Creating initial list of non-local Delaunay edges\n");
		CHI_VECTOR<CHI_INTERIOREDGE> non_local_del_edges;
		ListNonLocallyDelaunayEdges(curPatch,&non_local_del_edges);
		
		//================================================= Iterate to remove non-locally delaunay edges
		//printf("Iterating to remove non-locally delaunay edges\n");
		CHI_VECTOR<CHI_OPENEDGE> hullEdges;
		//FindPatchOpenEdges(curPatch,&hullEdges);

		int initialOpenEdgeCount = hullEdges.itemCount;
		iter=0;
		forcestop=false;
		while ((non_local_del_edges.itemCount>0) && (!forcestop))
		{
			iter++;
			//printf("Edge Flip %d\n",iter);
			EdgeFlip(curPatch,non_local_del_edges.GetItem(0));
			ListNonLocallyDelaunayEdges(curPatch,&non_local_del_edges);

			//FindPatchOpenEdges(curPatch,&hullEdges);
			if (hullEdges.itemCount!= initialOpenEdgeCount)
			{
				printf("Error in Delaunay generation\n");
				forcestop=true;
			}
			//if (iter>10) {forcestop=true;}
		}

		DumpPatchToScilab(curPatch);

	}
}