#include "chi_surfaceremesher.h"

using namespace chiSurfaceMeshing;

//###################################################################
/** This function copies all essential vertices and simplices of the
 * original mesh to a new mesh.
\author Jan*/
bool CHI_SURFACEREMESHER::CopyInformationToRemeshedSurface(CHI_VECTOR<CHI_PATCH> *patchCollections,
                                                           CHI_REMESHEDSURFACE *remeshedSurface)
{
	//printf("############################################### Starting surface refinement\n");
	CHI_VECTOR<int> vertexAdditionHistory;
	
	//printf("############################################### Copying patches, essential vertices and simplexes\n");
	for (int p=0;p<patchCollections->itemCount;p++)
	{
		//=========================================================== Mirror spawn patch
		CHI_PATCH* patch = patchCollections->GetItem(p);
		CHI_PATCH* newPatch = new CHI_PATCH;
        newPatch->parent = patch;
		remeshedSurface->patches.PushItem(newPatch);
		newPatch->normal = patch->normal;
		//printf("Patch %3d:\n",p);

		//================================================= If interior faces are to be removed
		//Copy only the essential vertices
		if (removeInteriorFaces)
		{
			for (int ev = 0; ev < patch->essentialVertices.itemCount; ev++)
			{
				//============================================= Get vertex index and add to history
				int* vertexIndex = new int;
				*vertexIndex = *patch->essentialVertices.GetItem(ev);
				
				//============================================= Check for duplicate
				bool duplicateFound = false;
				for (int d=0;d<vertexAdditionHistory.itemCount;d++)
				{
					if (*vertexIndex==*vertexAdditionHistory.GetItem(d))
					{
						int* index = new int;
						*index = d;
						newPatch->essentialVertices.PushItem(index);
						duplicateFound=true;
						break;
					}
				}
				
				//============================================= No duplicate add
				if (!duplicateFound)
				{
					vertexAdditionHistory.PushItem(vertexIndex);
					
					float *vertex = new float[3];
					GetVertexFromIndex(initialMesh, *vertexIndex, vertex);
					//printf("Vertex %3d, [%6.3f,%6.3f,%6.3f]", *vertexIndex, vertex[0], vertex[1], vertex[2]);
					int* index = new int;
					*index = remeshedSurface->vertexStack.PushItem(vertex);
					newPatch->essentialVertices.PushItem(index);
					//printf(" new index->%3d\n", remeshedSurface->vertexStack.itemCount - 1);
				}
			}
		}
		else
		{
			for (int f=0; f<patch->faceList.itemCount; f++)
			{
				CST_FACE* curFace = patch->faceList.GetItem(f);
				for (int v=0; v<3; v++)
				{
					//============================================= Get vertex index and add to history
					int* vertexIndex = new int;
					*vertexIndex = curFace->vertex[v];
					
					//============================================= Check for duplicate
					bool duplicateFound = false;
					for (int d=0;d<vertexAdditionHistory.itemCount;d++)
					{
						if (*vertexIndex==*vertexAdditionHistory.GetItem(d))
						{
							duplicateFound=true;
							break;
						}
					}
					
					//============================================= No duplicate add
					if (!duplicateFound)
					{
						vertexAdditionHistory.PushItem(vertexIndex);
						
						float *vertex = new float[3];
						GetVertexFromIndex(initialMesh, *vertexIndex, vertex);
						//printf("Vertex %3d, [%6.3f,%6.3f,%6.3f]", *vertexIndex, vertex[0], vertex[1], vertex[2]);
						int* index = new int;
						*index = remeshedSurface->vertexStack.PushItem(vertex);
						newPatch->essentialVertices.PushItem(index);
						//printf(" new index->%3d\n", remeshedSurface->vertexStack.itemCount - 1);
					}
				}
			}
		}
		
		//============================================= Copy edgeloops
		for (int s=0;s<patch->edgeLoops.itemCount;s++)
		{
			CHI_EDGELOOP* old_simplex = patch->edgeLoops.GetItem(s);
			CHI_EDGELOOP* new_simplex = new CHI_EDGELOOP;
			
			//printf("Copying edge-loop %3d:\n",s);
			
			//=============================== Copy edges
			for (int e=0;e<old_simplex->edges.itemCount;e++)
			{
				int* curEdge = old_simplex->edges.GetItem(e);
				int* newEdge = new int[2];
				//printf("     Edge %3d, was %3d->%3d, ",e,curEdge[0],curEdge[1]);
				
				//================= Getting new vertex index of edge[0]'s vertex
				for (int vh=0;vh<vertexAdditionHistory.itemCount;vh++)
				{
					if (curEdge[0]==*vertexAdditionHistory.GetItem(vh))
					{
						newEdge[0]=vh;
						break;
					}
				}
				//================= Getting new vertex index of edge[1]'s vertex
				for (int vh=0;vh<vertexAdditionHistory.itemCount;vh++)
				{
					if (curEdge[1]==*vertexAdditionHistory.GetItem(vh))
					{
						newEdge[1]=vh;
						break;
					}
				}
				//printf("now %3d->%3d\n",newEdge[0],newEdge[1]);
				
				new_simplex->edges.PushItem(newEdge);
			}
			//================= Getting new vertex index of old-simplex initial vertex
			for (int vh=0;vh<vertexAdditionHistory.itemCount;vh++)
			{
				if (old_simplex->initialVertex==*vertexAdditionHistory.GetItem(vh))
				{
                    new_simplex->initialVertex=vh;
					break;
				}
			}
			//================= Getting new vertex index of old-simplex  final vertex
			for (int vh=0;vh<vertexAdditionHistory.itemCount;vh++)
			{
				if (old_simplex->finalVertex==*vertexAdditionHistory.GetItem(vh))
				{
                    new_simplex->finalVertex=vh;
					break;
				}
			}
            //printf("Old Simplex %d %d, New Simplex %d %d\n", old_simplex->initialVertex,old_simplex->finalVertex,
            //                                                 new_simplex->initialVertex,new_simplex->finalVertex);

			newPatch->edgeLoops.PushItem(new_simplex);
		}
	}
	
	return false;
}

