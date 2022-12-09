#pragma once
#include <OpenMesh/Core/IO/MeshIO.hh>  
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>  
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <unsupported/Eigen/MatrixFunctions>
#include <iostream>
#include <fstream>
#include <map>
#include "faceBodyModel.h"

class PointRegistration
{
public:
	PointRegistration();
	void readbodyFeaturePos(string path);
	void readfaceFeaturePos(string path);
	void readHandFeaturePos(string path);
public:
	FaceBodyModel * fbModel;
	std::vector<int> faceFeatureIdx;
	std::vector<Eigen::Vector3d> faceFeaturePos;
	std::vector<int> handFeatureIdx;
	std::vector<Eigen::Vector3d> handFeaturePos;
	std::vector<int> bodyFeatureIdx;
	std::vector<Eigen::Vector3d> bodyFeaturePos;

	Eigen::VectorXd shapeParm;
	Eigen::VectorXd poseParm;
	Eigen::VectorXd expressionParm;
	Eigen::VectorXd gestureParm;
	Eigen::VectorXd motionParm;
};