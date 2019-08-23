#include "chi_surfaceremesher.h"

using namespace chiSurfaceMeshing;

//###################################################################
/**
\author Jan*/
bool CHI_SURFACEREMESHER::CreatePatchList(CHI_VECTOR<CHI_FACELIST>* coPlanarFaceCollections,CHI_VECTOR<CHI_PATCH>* patchCollections)
{
	//===================================================== Running through each co-planar list
	for (int coPl=0;coPl<coPlanarFaceCollections->itemCount;coPl++)
	{
		CHI_FACELIST* curCollection = coPlanarFaceCollections->GetItem(coPl);

		//======================================= Create temporary patch list
		CHI_VECTOR<CHI_PATCH>  tempPatchList;

		//======================================= Create a patch for the first face
		CHI_PATCH* newPatch = new CHI_PATCH;

		CST_FACE* firstFace = curCollection->faceList.GetItem(0);
		newPatch->normal(0)=firstFace->faceNormal[0];
		newPatch->normal(1)=firstFace->faceNormal[1];
		newPatch->normal(2)=firstFace->faceNormal[2];

		newPatch->faceList.PushItem(firstFace);
		tempPatchList.PushItem(newPatch);

		//======================================= Compile list of unused faces
		CHI_FACELIST unusedFaces;

		for (int f=1;f<curCollection->faceList.itemCount;f++)
		{
			CST_FACE* curFace = curCollection->faceList.GetItem(f);
			unusedFaces.faceList.PushItem(curFace);
		}
        //printf("Unused face count at start: %d\n", unusedFaces.faceList.itemCount);

		//======================================= Recursively search patch connections
		bool stopRecursive = false;
        int iter=0;
		while (!stopRecursive)
		{
			bool aPatchWasUpdated = false;
			for (int p=0;p<tempPatchList.itemCount;p++)
			{
				CHI_PATCH* curPatch = tempPatchList.GetItem(p);
                //printf("Patch %d of %d\n",p,tempPatchList.itemCount);

				for (int f=0;f<curPatch->faceList.itemCount;f++)
				{
                    //printf("Patch %d, face %d\n",p,f);
					CST_FACE* patchFace = curPatch->faceList.GetItem(f);

                    //============= Compare to each unused face
					for (int f2=unusedFaces.faceList.itemCount-1;f2>=0;f2--)
					{
                        //printf("Comparing..");
						CST_FACE* unusedFace = unusedFaces.faceList.GetItem(f2);

						bool vertexConnected = false;
						for (int i=0;i<3;i++)
						{
							for (int j=0;j<3;j++)
							{
								if (patchFace->vertex[i]==unusedFace->vertex[j])
								{
									vertexConnected = true;
								}
							}
						}

						//============= If connected, add it to the patch and remove from unused list
						if (vertexConnected)
						{
							curPatch->faceList.PushItem(unusedFace);
							unusedFaces.faceList.PullItem(f2);
							aPatchWasUpdated = true;
                            //printf("connection found\n");
						}
                        else
                        {
                            //printf("no connection\n");
                        }
					}
				}
			}
            iter++;
            //printf("Iteration %3d Unused faces %3d tempPatches %3d PatchUpdated %d \n",iter, unusedFaces.faceList.itemCount,tempPatchList.itemCount,aPatchWasUpdated);
			if ((unusedFaces.faceList.itemCount==0) || (iter>1000))
			{
				stopRecursive=true;
			}
            else
            {
                if (!aPatchWasUpdated)
                {
                    CST_FACE* curFace = unusedFaces.faceList.PullItem(unusedFaces.faceList.itemCount-1);
                    newPatch = new CHI_PATCH;
                    newPatch->faceList.PushItem(curFace);
	                newPatch->normal(0)=curFace->faceNormal[0];
	                newPatch->normal(1)=curFace->faceNormal[1];
	                newPatch->normal(2)=curFace->faceNormal[2];
	                
                    tempPatchList.PushItem(newPatch);
                    //printf("New patch created\n");
                }
            }

		}

        for (int k=0;k<tempPatchList.itemCount;k++)
        {
	        CHI_PATCH* curPatch = tempPatchList.GetItem(k);
            patchCollections->PushItem(curPatch);
        }
	}



	printf("Number of patches found: %d\n",patchCollections->itemCount);
	return true;
}