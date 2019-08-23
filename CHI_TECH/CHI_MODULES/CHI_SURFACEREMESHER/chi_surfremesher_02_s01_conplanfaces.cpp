#include "chi_surfaceremesher.h"

using namespace chiSurfaceMeshing;

//###################################################################
/** Builds a collection of lists. Each list contains faces that are
 * co-planar.
 *
\author Jan*/
bool CHI_SURFACEREMESHER::CollectCoPlanarFaces(CHI_VECTOR<CHI_FACELIST>* coPlanarFaceCollections)
{
	//===================================================== Update face normals
    if (this->recalcNormal)
    {
        for (int f=0;f<initialMesh->faceStack.itemCount;f++)
        {
            CST_FACE*           curFace = initialMesh->faceStack.GetItem(f);
            CHI_TRIPLET_VERTEX  triplet;

            GetFaceVertices(initialMesh,curFace,&triplet);
            Eigen::Vector3f vec1 = triplet.v1 - triplet.v0;
            Eigen::Vector3f vec2 = triplet.v2 - triplet.v1;

            Eigen::Vector3f normal = vec1.cross(vec2);
            normal.normalize();

            curFace->faceNormal[0] = normal(0);
            curFace->faceNormal[1] = normal(1);
            curFace->faceNormal[2] = normal(2);
        }
    }
    else
    {
        printf("Using vertex normals\n");
        for (int f=0;f<initialMesh->faceStack.itemCount;f++)
        {
            CST_FACE*           curFace = initialMesh->faceStack.GetItem(f);
            float*           vertNormal = initialMesh->normalStack.GetItem(curFace->normal[0]);

            curFace->faceNormal[0] = vertNormal[0];
            curFace->faceNormal[1] = vertNormal[1];
            curFace->faceNormal[2] = vertNormal[2];
        }
    }

    
    
    //===================================================== Make co-planar list;
    for (int f=0;f<initialMesh->faceStack.itemCount;f++)
    {
        CST_FACE *curFace = initialMesh->faceStack.GetItem(f);
        Eigen::Vector3f normal(curFace->faceNormal[0],
                               curFace->faceNormal[1],
                               curFace->faceNormal[2]);

        bool faceFoundACollection=false;

        //======================================= If no collections, make new one
        if (coPlanarFaceCollections->itemCount==0)
        {
            CHI_FACELIST* newCollection = new CHI_FACELIST;
            newCollection->faceList.PushItem(curFace);
            coPlanarFaceCollections->PushItem(newCollection);
            faceFoundACollection = true;
        }
        //======================================= Run through each collection
        else
        {
            for (int col=0;col<coPlanarFaceCollections->itemCount;col++)
            {
                CHI_FACELIST* curCollection = coPlanarFaceCollections->GetItem(col);
                //===================== Use only first face
                CST_FACE* compareFace = curCollection->faceList.GetItem(0);
                Eigen::Vector3f comparenormal(compareFace->faceNormal[0],
                                              compareFace->faceNormal[1],
                                              compareFace->faceNormal[2]);

                CHI_TRIPLET_VERTEX  tripletcurFace;
                CHI_TRIPLET_VERTEX  tripletcompareFace;

                GetFaceVertices(initialMesh,curFace    ,&tripletcurFace);
                GetFaceVertices(initialMesh,compareFace,&tripletcompareFace);

                //Use one vertex from the currentFace and 2 from the other
                Eigen::Vector3f vec1 = tripletcompareFace.v1 - tripletcurFace.v0;
                Eigen::Vector3f vec2 = tripletcompareFace.v2 - tripletcompareFace.v1;
                Eigen::Vector3f planarNormal = vec1.cross(vec2);
                planarNormal.normalize();


                float dotproductNormal = fabs(normal.dot(comparenormal));
                float dotproductPlanar = fabs(planarNormal.dot(comparenormal));


                //if ( (dotproductNormal>0.999999) && (dotproductPlanar>0.999999) )
                if ( (dotproductNormal>(1.0-this->precision)) )
                {
                    curCollection->faceList.PushItem(curFace);
                    faceFoundACollection = true;
                    break; //Dont look for a home anymore
                }
            }
        }

        //======================================= If it didnt find a unique collection
        if (!faceFoundACollection)
        {
            CHI_FACELIST* newCollection = new CHI_FACELIST;
            newCollection->faceList.PushItem(curFace);
            coPlanarFaceCollections->PushItem(newCollection);
        }
    }
    printf("Number of co-planar collections found: %d\n",coPlanarFaceCollections->itemCount);
	return true;
}