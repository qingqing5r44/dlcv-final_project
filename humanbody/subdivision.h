#pragma once
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>  
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <iostream>
#include <fstream>
#include <OpenMesh/Core/IO/MeshIO.hh>  

using namespace std;

class AdaptiveSubdivision
{
	typedef OpenMesh::PolyMesh_ArrayKernelT<> PGMesh;

public:
	AdaptiveSubdivision(PGMesh * _mesh);

	void subdivision();

	void getHeadVertices(char* h_path);

	void setSubdividedArea();

	void setSubdivisionCoeffs();

	Eigen::MatrixXd getVertexMatrix(PGMesh * _mesh);

	Eigen::VectorXd getVertexVec(PGMesh *_mesh);

	Eigen::MatrixXd getSubdivisionPosition(Eigen::MatrixXd V0);

	std::vector<int> det(PGMesh::VertexHandle v_);

	std::vector<int> adjust_seq(std::vector<int>, int);

	void outMesh(Eigen::MatrixXd cordComponent, string path);

	int isIn(std::vector<int>, int);

public:
	OpenMesh::VPropHandleT<bool> is_head_vertice;
	OpenMesh::EPropHandleT<bool> is_head_edge;
	OpenMesh::FPropHandleT<bool> is_head_face;
	OpenMesh::EPropHandleT<int> edge_vertice_idx;
	int n_affected_vertices;
	int n_affected_faces;
	int n_affected_edges;
	std::vector<int> affected_edges;
	std::vector<int> affected_faces;
	std::vector<int> affected_vertices;
	Eigen::SparseMatrix<double> subdivisionCoeff;
	PGMesh * mesh_;
	PGMesh * subdividedMesh;
	int nV;
	int nE;
	int nF;
	std::vector<int> index;
	std::vector<double> coef;
	//std::vector<Eigen::Triplet<double> > triple;
};