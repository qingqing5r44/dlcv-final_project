#ifndef SPLOC_H
#define SPLOC_H

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <iostream>
#include <fstream>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <OpenMesh/Core/IO/MeshIO.hh>  
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>  
#include"kdtree.hpp"
#include <omp.h>
using namespace std;
using namespace Eigen;
using namespace GeometryProcess::nanoflann;

class Node
{
public:
	typedef  OpenMesh::PolyMesh_ArrayKernelT<> PGMesh;
	double geodist;
	int nodeNo;
	PGMesh::VertexHandle vh;
	bool operator<(const Node& obj) const
	{
		return geodist < obj.geodist;
	}
	bool operator()(const Node* t1, const Node* t2)
	{
		return t1->geodist < t2->geodist;
	}
};

class DeformationGraph
{
	typedef OpenMesh::PolyMesh_ArrayKernelT<> PGMesh;
public:
	DeformationGraph(PGMesh *_mesh, std::vector<PGMesh::VertexHandle> &_vertexHandle, std::vector<int> &_vertexIndex, int seqSize_);
	PGMesh *mesh_;

	OpenMesh::VPropHandleT<PGMesh::VertexHandle> vh;
	OpenMesh::VPropHandleT<int> vnum;
	/*The rotation matrix of each node of the deformation graph*/
	std::vector<std::vector<Eigen::MatrixXd>> rotation;
	/*The translation vector of each node of the deformation graph*/
	std::vector<std::vector<Eigen::VectorXd>> translation;
	std::vector<std::vector<Eigen::VectorXd>> nodePosition;
	std::vector<PGMesh::VertexHandle> vertexHandle_;
	std::vector<int> vertexIndex_;
	int sequenceSize;
};

// This is an example of a custom data set class to construct a kd-tree
template <typename T>
struct PointCloud
{
	struct Point
	{
		T x, y, z;
		OpenMesh::VertexHandle vh;
	};
	std::vector<Point>  pts;

	// Must return the number of data points
	inline size_t kdtree_get_point_count() const { return pts.size(); }

	// Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
	inline T kdtree_distance(const T *p1, const size_t idx_p2, size_t /* size*/) const
	{
		const T d0 = p1[0] - pts[idx_p2].x;
		const T d1 = p1[1] - pts[idx_p2].y;
		const T d2 = p1[2] - pts[idx_p2].z;
		return d0*d0 + d1*d1 + d2*d2;
	}

	// Returns the dim'th component of the idx'th point in the class:
	// Since this is inlined and the "dim" argument is typically an immediate value, the
	//  "if/else's" are actually solved at compile time.
	inline T kdtree_get_pt(const size_t idx, size_t dim) const
	{
		if (dim == 0) return pts[idx].x;
		else if (dim == 1) return pts[idx].y;
		else return pts[idx].z;
	}
	// Returns the vertex handle of the idx'th point in the class:
	inline OpenMesh::VertexHandle kdtree_get_vh(const int idx) const
	{
		return pts[idx].vh;
	}
	// Optional bounding-box computation: return false to default to a standard bbox computation loop.
	//   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
	//   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
	template <class BBOX>
	bool kdtree_get_bbox(BBOX & /* bb*/) const { return false; }
};

class SparseLocalDecomposition
{
	typedef OpenMesh::PolyMesh_ArrayKernelT<> PGMesh;
	typedef GeometryProcess::nanoflann::KDTreeSingleIndexAdaptor<GeometryProcess::nanoflann::L2_Simple_Adaptor<double, PointCloud<double>>, PointCloud<double>, 3>
		my_kd_tree_simple_t;
	//typedef GeometryProcess::nanoflann::KDTreeSingleIndexAdaptor<GeometryProcess::nanoflann::L2_Adaptor<double, PointCloud<double>>, PointCloud<double>, 3>
		//my_kd_tree_t;
	typedef KDTreeEigenMatrixAdaptor< Eigen::Matrix<double, Dynamic, Dynamic> >  my_kd_tree_t;
public:
	SparseLocalDecomposition(std::vector<PGMesh *> &_meshSeq);
	std::vector<Eigen::VectorXd> deformWithFeaturePoints(std::vector<int> feature_index, std::vector<std::vector<double>> feature_location, string path, std::vector<int> v_index);
	void preprocessing1();
	void preprocessing2();
	void rigidAlignment(Eigen::MatrixXd from, Eigen::MatrixXd to, Eigen::MatrixXd * R, double * s, Eigen::VectorXd *t);
	bool sortByDist(const Node  *v1, const Node  *v2);
	void findDenseCorrespondence(std::vector<Eigen::VectorXd> from, my_kd_tree_simple_t simple_t);
	int JudgeNum(string str, double& iTmp);
public:
	DeformationGraph *defgraph_;
	std::vector<PGMesh *> meshSeq_;
	PGMesh *ref_;
	PGMesh refMesh;
	OpenMesh::VPropHandleT<std::vector<PGMesh::VertexHandle>> knNodes;
	OpenMesh::VPropHandleT<std::vector<double>> knWeights;
	OpenMesh::VPropHandleT<Eigen::VectorXd> Yaxis;
	int sequenceSize;
	PointCloud<double> Points;
	my_kd_tree_t *kd_tree_t;
	double f_error;
	double* geoDistToNodes;
	std::vector<int> feature_index;
	//std::vector<int> init_index;

};

#endif