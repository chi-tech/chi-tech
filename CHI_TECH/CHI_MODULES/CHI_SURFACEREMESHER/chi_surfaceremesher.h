#ifndef CHI_SURFACEREMESHER_H
#define CHI_SURFACEREMESHER_H


#include"../../CHI_GRAPHICS/CHI_SURFACE/chi_surface.h"
#include "../../CHI_VECTOR/chi_vector.h"
#include "../../../CHI_RESOURCES/Dependencies/Eigen/Eigen"

#include"chi_surfaceremesher_structs.h"
#include"CHI_REMESHEDSURFACE/chi_remeshedsurface.h"

using namespace chiSurfaceMeshing;


//######################################################### Class Def
/** Object for performing surface remeshing.
\author Jan*/
class CHI_SURFACEREMESHER
{
public:
	double      baseSize;                   ///< Base size to strive for
	double      absoluteMinumumSize;        ///< Absolute minimum size
	bool        removeInteriorFaces;
	bool        keep2Dorientation;
	bool        recalcNormal;
	double      precision;
	
	int         structureOption;            ///< 0=Unstructured, 1=Structured_Cartesian
	
	CHI_SURFACE* initialMesh;
	CHI_REMESHEDSURFACE* remeshedSurface;
private:
	CHI_VECTOR<CHI_FACELIST> coPlanarFaces;
	CHI_VECTOR<CHI_FACELIST> patchList;
	
public:
	//00
				CHI_SURFACEREMESHER();
	//01
	bool        ExecuteMeshing(CHI_SURFACE* initialMesh);
private:
	//02
	bool        CollectCoPlanarFaces(CHI_VECTOR<CHI_FACELIST>* coPlanarFaceCollections);
	bool        CreatePatchList(CHI_VECTOR<CHI_FACELIST>* coPlanarFaceCollections,CHI_VECTOR<CHI_PATCH>* patchCollections);
	bool        FindOpenEdges(CHI_VECTOR<CHI_PATCH>* patchCollections);
	bool        FindOpenEdges(CHI_PATCH* patch);
	bool        FindCoLinearEdges(CHI_VECTOR<CHI_PATCH>* patchCollections);
	bool        FindEdgeLoops(CHI_VECTOR<CHI_PATCH>* patchCollections);
	bool        FindEdgeLoops(CHI_PATCH* patch);
	bool        ListEssentialVertices(CHI_VECTOR<CHI_PATCH>* patchCollections);
	bool        CopyInformationToRemeshedSurface(CHI_VECTOR<CHI_PATCH> *patchCollections, CHI_REMESHEDSURFACE *remeshedSurface);
	
	//03 Utilities
    bool        GetVertexFromIndex(CHI_SURFACE* mesh,int index, float* vertex);
	void        GetFaceVertices(CHI_SURFACE* mesh,CST_FACE* face,CHI_TRIPLET_VERTEX* vertices);
    bool        CheckEdgeConnection(int* edge,CST_FACE* face);
	bool        CheckEdgesColinear(int* edgeA, int* edgeB, CHI_SURFACE* mesh);
	void        ResolveEdgeExtent(CHI_VECTOR<CHI_COLINEDGELIST>* colinearEdges,CHI_SURFACE* mesh);
	void        FindPatchOpenEdges(CHI_PATCH *patch, CHI_VECTOR<CHI_OPENEDGE> *openEdges);
	void        DeleteTriangle(CHI_PATCH* patch,int tauIndex);
	void        FindCircumCircle(CHI_PATCH* patch, CST_FACE* tri, float* center, float* radius);
	void        DumpPatchToScilab(CHI_PATCH* patch);
	void        CollectTrianglesWithVertexInCircumdisc(CHI_PATCH *patch, CST_FACE* master, CHI_PATCH *bucket, Eigen::Vector3f vc);
	bool        FindIntersectingSimplices(CHI_PATCH *patch,int va, int vb, int vc);
	void        RemoveFillerTriangles(CHI_PATCH* patch);
	
	//04 2D Delaunay functions
	void        Create2DDelaunayTriangulation(CHI_REMESHEDSURFACE* mesh);
	void        Project3DVerticesTo2D(CHI_REMESHEDSURFACE* mesh,CHI_PATCH* patch);
	void        SortLexicographically2D(CHI_PATCH* patch);
	void        CreateInitialTriangle(CHI_PATCH* patch);
	double      Orient2D(float* a, float* b, float* c);
	void        GetUnusedVertices(CHI_PATCH* patch, CHI_VECTOR<int>* unusedVerts);
	void        AttachUnusedVertex(CHI_PATCH* patch,CHI_VECTOR<int>* unusedVertices);
	void        ConvexifyHull(CHI_PATCH* patch);
	double      InCircle(float* a, float* b, float* c, float* d);
	void        ListNonLocallyDelaunayEdges(CHI_PATCH* patch, CHI_VECTOR<CHI_INTERIOREDGE>* non_loc_del_edges);
	void        EdgeFlip(CHI_PATCH* patch, CHI_INTERIOREDGE* non_loc_del_edge);
	
	//05 2D Delaunay refinement
	void        RefinePatchMesh(CHI_PATCH* patch);
	void        RestoreOptimality(CHI_PATCH* patch);
	bool        RefineEdgeSize(CHI_PATCH* patch);
	void        Make2DDelaunayTriangulation(CHI_PATCH* patch);
	void        FindEncroachedSubSegments(CHI_PATCH* patch, CHI_INTERIOREDGE*& encrSeg, bool excludeNonsplittable=false);
	void        RemoveEncroachedEdges(CHI_PATCH *patch);
	bool        SplitSubSegment(CHI_PATCH* patch, CHI_INTERIOREDGE* encrSeg);
	void        FindSmallRhoTriangles(CHI_PATCH* patch, CST_FACE*& smallRhoTri,float criteria);
	void        SplitTriangle(CHI_PATCH* patch, CST_FACE* tri);
	void        InsertVertex(CHI_PATCH* patch, float* vertex,CHI_INTERIOREDGE *encrSeg=NULL);
	void        DeleteVertex(CHI_PATCH* patch, int vertexIndex);
	
	
	
	
};

#endif