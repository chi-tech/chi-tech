#ifndef CHI_SURFACEREMESHER_STRUCTS_H
#define CHI_SURFACEREMESHER_STRUCTS_H

namespace chiSurfaceMeshing
{
	struct CHI_COLINEDGELIST
	{
		CHI_VECTOR<int> edges;
		int             extentVertices[2];
		float           extentLength;
	};
	
	struct CHI_EDGELOOP
	{
		CHI_VECTOR<int> edges;
		int             initialVertex;
		int             finalVertex;
		bool            closed;
		
		CHI_EDGELOOP()
		{
			initialVertex=-1;
			finalVertex  =-1;
			closed       =false;
		}
	};
	
	
	struct CHI_FACELIST
	{
		CHI_VECTOR<CST_FACE> faceList;
	};
	
	struct CHI_PATCH
	{
		Eigen::Vector3f normal;
		CHI_VECTOR<CST_FACE> faceList;
		CHI_VECTOR<int> edges;
		CHI_VECTOR<CHI_COLINEDGELIST> colinearEdges;
		CHI_VECTOR<CHI_EDGELOOP> edgeLoops;
		CHI_VECTOR<int> essentialVertices;
		CHI_VECTOR<int> lexicalListOfVertices;
		Eigen::Vector3f hat_i;
		Eigen::Vector3f hat_j;
		Eigen::Vector3f hat_k;
		Eigen::Vector3f C_0;
		CHI_VECTOR<float> Pstar;
		CHI_PATCH* parent;
	};
	
	struct CHI_OPENEDGE
	{
		int vi;
		int vf;
		int ownerTriangle;
		int ownerEdgeNumber;
	};
	
	struct CHI_INTERIOREDGE
	{
		int vi;
		int vf;
		int tau_1;
		int tau_1_edgeNumber;
		int tau_2;
		int tau_2_edgeNumber;
		int encroachingVertex;
		
		CHI_INTERIOREDGE()
		{
			vi=-1;
			vf=-1;
			tau_1=-1;
			tau_2=-1;
			tau_1_edgeNumber=-1;
			tau_2_edgeNumber=-1;
			encroachingVertex=-1;
		}
	};
	
	
	
	struct CHI_TRIPLET_VERTEX
	{
		Eigen::Vector3f v0;
		Eigen::Vector3f v1;
		Eigen::Vector3f v2;
		int             i0;
		int             i1;
		int             i2;
	};

}

#endif