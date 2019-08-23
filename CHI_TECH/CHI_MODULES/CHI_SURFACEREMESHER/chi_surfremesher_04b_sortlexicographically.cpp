#include "chi_surfaceremesher.h"

//############################################################################# Sort lexicographically
/** Sorts a 2D point set lexicographically.
\author Jan*/
void CHI_SURFACEREMESHER::SortLexicographically2D(CHI_PATCH* patch)
{
	CHI_VECTOR<int> *lexlist = &patch->lexicalListOfVertices;
	lexlist->ClearVector();
	CHI_VECTOR<float> tempVertices;
	CHI_VECTOR<float> *projectedVertices = &patch->Pstar;
	
	//========================================================================= Copying vertex list
	for (int i=0; i<patch->essentialVertices.itemCount; i++)
	{
		int vertexIndex = *patch->essentialVertices.GetItem(i);
		tempVertices.PushItem(projectedVertices->GetItem(vertexIndex));
		int* newIndex = new int; *newIndex = vertexIndex;
		lexlist->PushItem(newIndex);
	}
	
	//========================================================================= Bubble sorting by x
	for (int i=0;i<tempVertices.itemCount;i++)
	{
		
		for (int j=0;j<tempVertices.itemCount-1;j++)
		{
			float* mastrVert = tempVertices.GetItem(j);
			int*   mastrIndx = lexlist->GetItem(j)    ;
			
			float* slaveVert = tempVertices.GetItem(j+1);
			int*   slaveIndx = lexlist->GetItem(j+1);
			
			if (mastrVert[0]>slaveVert[0])
			{
				float* tempVert = tempVertices.GetItem(j+1);
				int*   tempIndx = lexlist->GetItem(j+1);
				
				tempVertices.SetItem(j+1,tempVertices.GetItem(j));
				lexlist->SetItem(j+1,    lexlist->GetItem(j)    );
				
				tempVertices.SetItem(j,tempVert);
				lexlist->SetItem(j,tempIndx);
			}
		}
	}
	
	//========================================================================= Bubble sorting by y
	 for (int i=0;i<tempVertices.itemCount;i++)
	{
		
		for (int j=0;j<tempVertices.itemCount-1;j++)
		{
			float* mastrVert = tempVertices.GetItem(j);
			int*   mastrIndx = lexlist->GetItem(j)    ;
			
			float* slaveVert = tempVertices.GetItem(j+1);
			int*   slaveIndx = lexlist->GetItem(j+1);
			
			if ( (fabs(mastrVert[0]-slaveVert[0])<this->precision) && (mastrVert[1]>slaveVert[1]) )
			{
				float* tempVert = tempVertices.GetItem(j+1);
				int*   tempIndx = lexlist->GetItem(j+1);
				
				tempVertices.SetItem(j+1,tempVertices.GetItem(j));
				lexlist->SetItem(j+1,    lexlist->GetItem(j)    );
				
				tempVertices.SetItem(j,tempVert);
				lexlist->SetItem(j,tempIndx);
			}
		}
	}

	
	//========================================================================= Print sorted vertices
	for (int i=0; i<tempVertices.itemCount; i++)
	{
		float* slaveVert = tempVertices.GetItem(i);
		int*   slaveIndx = lexlist->GetItem(i);
		
		//printf("Lexical Vert %3d, Index %3d, [%6.3f,%6.3f,%6.3f]\n",i,*slaveIndx,slaveVert[0],slaveVert[1],slaveVert[2]);
	}
}