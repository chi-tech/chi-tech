#include"chi_surfaceremesher.h"

//############################################################################# Gets unused vertices
/**Updates a list of unused vertices.

The lexlist, v_l1 to v_lN, contains indices of vertices. The unused list will be
a list of indices in the lexlist.

\author Jan*/
void CHI_SURFACEREMESHER::GetUnusedVertices(CHI_PATCH *patch, CHI_VECTOR<int> *unusedVerts)
{
	CHI_VECTOR<int> *lexlist = &patch->lexicalListOfVertices;
	unusedVerts->ClearVector();
	for (int i=0; i<lexlist->itemCount;i++)
	{
		int vertIndex = *lexlist->GetItem(i);
		bool foundVert = false;
		for (int f=0; f<patch->faceList.itemCount;f++)
		{
			CST_FACE* curFace = patch->faceList.GetItem(f);
			for (int k=0;k<3;k++)
			{
				if (curFace->vertex[k]==vertIndex)
				{
					foundVert=true;
				}
			}
			if (foundVert) {break;}
		}
		if (!foundVert)
		{
			int* newIndex = new int;
			*newIndex = i;
			unusedVerts->PushItem(newIndex);
		}
	}
}
