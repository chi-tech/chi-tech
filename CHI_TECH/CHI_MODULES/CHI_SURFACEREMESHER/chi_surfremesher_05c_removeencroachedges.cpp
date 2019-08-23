#include"chi_surfaceremesher.h"

//############################################################################# Remove encroached edges
/** Removes encroached edges by splitting them.

\author Jan*/
void  CHI_SURFACEREMESHER::RemoveEncroachedEdges(CHI_PATCH *patch)
{
	CHI_INTERIOREDGE* encrSeg = NULL;
	FindEncroachedSubSegments(patch,encrSeg,false);
	
	int iter=0;
	bool forceStop=false;
	bool rejectedSegmentSplit=false;
	while ((encrSeg!=NULL) && (!forceStop))
	{
		iter++;
		//printf("  Refinement subiteration in RemoveEncroachedEdges %d\n",iter);
		//printf("   - Splitting encroached edge %d->%d\n",encrSeg->vi,encrSeg->vf);
		rejectedSegmentSplit = SplitSubSegment(patch, encrSeg);
		
		Make2DDelaunayTriangulation(patch);
		FindEncroachedSubSegments(patch,encrSeg,false);
		//DumpPatchToScilab(patch);
		if ((iter>=400) || (rejectedSegmentSplit))
		{
			forceStop=true;
		}
	}
}