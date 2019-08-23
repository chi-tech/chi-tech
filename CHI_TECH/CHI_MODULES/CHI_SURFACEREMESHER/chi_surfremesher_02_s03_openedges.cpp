#include "chi_surfaceremesher.h"

using namespace chiSurfaceMeshing;

//###################################################################
/** Searches each patch for edges that are open, meaning that they
 * are connected to only a single face.
 *
\author Jan*/
bool CHI_SURFACEREMESHER::FindOpenEdges(CHI_VECTOR<CHI_PATCH>* patchCollections)
{


	for (int p=0; p<patchCollections->itemCount; p++)
	{
		CHI_PATCH* curPatch = patchCollections->GetItem(p);
		curPatch->edges.ClearVector();

		for (int f=0; f<curPatch->faceList.itemCount; f++)
		{
			CST_FACE* curFace = curPatch->faceList.GetItem(f);

			int edge0[2];
			int edge1[2];
			int edge2[2];

			edge0[0] = curFace->vertex[0]; edge0[1] = curFace->vertex[1];
			edge1[0] = curFace->vertex[1]; edge1[1] = curFace->vertex[2];
			edge2[0] = curFace->vertex[2]; edge2[1] = curFace->vertex[0];

			bool connectionFound0=false;
			bool connectionFound1=false;
			bool connectionFound2=false;
			for (int f2=0; f2<curPatch->faceList.itemCount; f2++)
			{
				if (f2 != f)
				{
					CST_FACE *comFace = curPatch->faceList.GetItem(f2);

					if (CheckEdgeConnection(edge0,comFace)) {connectionFound0 = true;}
					if (CheckEdgeConnection(edge1,comFace)) {connectionFound1 = true;}
					if (CheckEdgeConnection(edge2,comFace)) {connectionFound2 = true;}

				}
			}
			if (connectionFound0==false)
			{
				int* newEdge = new int[2];
				newEdge[0] = edge0[0];
				newEdge[1] = edge0[1];
				curPatch->edges.PushItem(newEdge);
			}
			if (connectionFound1==false)
			{
				int* newEdge = new int[2];
				newEdge[0] = edge1[0];
				newEdge[1] = edge1[1];
				curPatch->edges.PushItem(newEdge);
			}
			if (connectionFound2==false)
			{
				int* newEdge = new int[2];
				newEdge[0] = edge2[0];
				newEdge[1] = edge2[1];
				curPatch->edges.PushItem(newEdge);
			}

		}
		//printf("Patch %3d has %3d open edges and %3d faces.\n",p,curPatch->edges.itemCount,curPatch->faceList.itemCount);
	}
	return true;
}

bool CHI_SURFACEREMESHER::FindOpenEdges(CHI_PATCH* patch)
{

	CHI_PATCH* curPatch = patch;
	patch->edges.ClearVector();
	
	for (int f=0; f<curPatch->faceList.itemCount; f++)
	{
		CST_FACE* curFace = curPatch->faceList.GetItem(f);
		if (curFace!=NULL)
		{
			int edge0[2];
			int edge1[2];
			int edge2[2];

			edge0[0] = curFace->vertex[0]; edge0[1] = curFace->vertex[1];
			edge1[0] = curFace->vertex[1]; edge1[1] = curFace->vertex[2];
			edge2[0] = curFace->vertex[2]; edge2[1] = curFace->vertex[0];

			bool connectionFound0=false;
			bool connectionFound1=false;
			bool connectionFound2=false;
			for (int f2=0; f2<curPatch->faceList.itemCount; f2++)
			{
				if (f2 != f)
				{
					CST_FACE *comFace = curPatch->faceList.GetItem(f2);
					if (comFace!=NULL)
					{
						if (CheckEdgeConnection(edge0,comFace)) {connectionFound0 = true;}
						if (CheckEdgeConnection(edge1,comFace)) {connectionFound1 = true;}
						if (CheckEdgeConnection(edge2,comFace)) {connectionFound2 = true;}
					}


				}
			}
			if (connectionFound0==false)
			{
				int* newEdge = new int[4];
				newEdge[0] = edge0[0];
				newEdge[1] = edge0[1];
				newEdge[2] = curFace->edges[0][2];
				curPatch->edges.PushItem(newEdge);
			}
			if (connectionFound1==false)
			{
				int* newEdge = new int[4];
				newEdge[0] = edge1[0];
				newEdge[1] = edge1[1];
				newEdge[2] = curFace->edges[1][2];
				curPatch->edges.PushItem(newEdge);
			}
			if (connectionFound2==false)
			{
				int* newEdge = new int[4];
				newEdge[0] = edge2[0];
				newEdge[1] = edge2[1];
				newEdge[2] = curFace->edges[2][2];
				curPatch->edges.PushItem(newEdge);
			}
		}

		
	}

	return true;
}