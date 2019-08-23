#include"chi_surfaceremesher.h"

//############################################################################# Restore optimal conditions to patch
/** Restores optimal conditions to a mesh.

\author Jan*/
void CHI_SURFACEREMESHER::RestoreOptimality(CHI_PATCH *patch)
{
    //printf("#################### Removing encroached edges\n");
    RemoveEncroachedEdges(patch);

    //DumpPatchToScilab(patch);

    //printf("#################### Splitting\n");
    CST_FACE* violatingTriangle = NULL;
    FindSmallRhoTriangles(patch,violatingTriangle,1.0);

    printf("Restoring optimality, iteration %5d",0);
    int iter=0;
    bool forceStop=false;
    while ((violatingTriangle!=NULL) && (!forceStop))
    {
        iter++;
        //printf("Refinement iteration %d\n",iter);
        //printf("Splitting triangle\n");
        SplitTriangle(patch,violatingTriangle);
        //printf("Making Delaunay\n");
        Make2DDelaunayTriangulation(patch);
        RemoveEncroachedEdges(patch);
        FindSmallRhoTriangles(patch,violatingTriangle,1.0);
        //DumpPatchToScilab(patch);

        printf("\b\b\b\b\b%5d",iter);
        //if (iter>=10)
        //{
        //    forceStop=true;
        //}
    }
    printf("\n");
}
