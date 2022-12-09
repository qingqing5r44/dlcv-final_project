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
#include "SparseLocalizedForConnectMap.h"
#include <Python.h>
#include <windows.h>

using namespace std;
using namespace Eigen;

class FaceBodyModel
{
	typedef OpenMesh::PolyMesh_ArrayKernelT<> PGMesh;
public:
	FaceBodyModel(int kS, int kP, PGMesh *_mesh);
	void loadShapeBases(char* s_path);
	void loadPoseBases(char* p_path, int n);
	void loadMeanShape(char* m_path);
	void loadExpressionBases(char* p_path, int n);
	void loadgestureBases(char* p_path, int n);
	void outputPoseBases();
	Eigen::VectorXd generateShapeWithoutPose(Eigen::VectorXd s_coeff);
	void generateShapeWithoutPose();
	Eigen::VectorXd generateShapeWithoutPose(string shapepath);
	Eigen::VectorXd getVertexMatrix(PGMesh * _mesh);
	Eigen::MatrixXd getVertexMatrix1(PGMesh * _mesh);
	void outputBlendshapes();
	void updateRef(Eigen::VectorXd x);
	std::vector<double> getAnchors(Eigen::VectorXd s_coeff, Eigen::VectorXd motion_parm);
	void getAnchorPos(Eigen::VectorXd s_coeff, Eigen::VectorXd pose_parm);
	void loadTensorflowModel();
	void loadDenseCorrespondenceModule();
	void load_pointcloud(int seqno);
	void findCoressPairs(int seqno, int ith, int jth);
	Eigen::VectorXd getPoseDE(Eigen::VectorXd s_coeff, Eigen::VectorXd pose_parm);
	Eigen::VectorXd getPoseParm(Eigen::VectorXd s_coeff, Eigen::VectorXd motion_parm);
	Eigen::MatrixXd f_x(Eigen::VectorXd s_coeff, Eigen::VectorXd motion_parm);
	Eigen::MatrixXd f_x(Eigen::VectorXd s_coeff, Eigen::VectorXd motion_parm, Eigen::VectorXd f_parm, Eigen::VectorXd g_parm);
	Eigen::MatrixXd f_x(string shape_path, string de_path);
	Eigen::MatrixXd f_x(string shape_path, Eigen::VectorXd s_coeff, Eigen::VectorXd motion_parm);
	Eigen::VectorXd getDE(string path);
	void initFace(PGMesh *pg);
	Eigen::MatrixXd f_expression(Eigen::VectorXd f_parm);
	Eigen::MatrixXd f_gesture(Eigen::VectorXd g_parm);
	Eigen::VectorXd getBodyPosfromCurrentV(Eigen::MatrixXd V, std::vector<int> idxset);
	Eigen::VectorXd getBodyPosfromCurrentV(Eigen::MatrixXd V, std::vector<std::vector<int>> idxset, std::vector<std::vector<double>> weightset);
	Eigen::VectorXd getFacePosfromCurrentV(Eigen::MatrixXd V);
	Eigen::VectorXd getBodyPosfromCurrentV1(Eigen::MatrixXd V, std::vector<int> idxset);
	void findRigidTransform(Eigen::VectorXd Vt, Eigen::Matrix3d & R, Eigen::Vector3d &t, std::vector<Eigen::Vector3d> layerpos);
	void findRigidTransform(Eigen::VectorXd Vt, Eigen::Matrix3d & R, Eigen::Vector3d &t, Eigen::MatrixXd layerpos);
	void optimizexfromde(Eigen::VectorXd &cm, Eigen::VectorXd s_coef, Eigen::VectorXd &p_coef);
	void optimizebodyfromde(Eigen::VectorXd &cm, Eigen::VectorXd &p_coef);
	void optimizeexpressionfromde(Eigen::VectorXd &cm, Eigen::VectorXd &f_coef);
	void optimizegesturefromde(Eigen::VectorXd &cm, Eigen::VectorXd &g_coef);
	std::vector<Eigen::Vector3d> getJointsPos(Eigen::VectorXd s_coeff, Eigen::VectorXd motion_parm);
	wchar_t *GetWC(const char *c)
	{
		const size_t cSize = strlen(c) + 1;
		size_t * ret = new size_t[cSize];
		wchar_t* wc = new wchar_t[cSize];
		mbstowcs_s(ret, wc, 100, c, cSize);

		return wc;
	}
	//读取body feature的idx和weights
	void readBodyFeatureidx(string p_path);
	void readBodyFeatureidx1(string p_path);

public:
	int k_Shapes;
	int k_Poses;
	int k_bodyposes;
	int k_expressions;
	int k_gestures;
	int nE;
	int nF;
	int nV;
	Eigen::MatrixXd C_expressions;
	Eigen::LLT<MatrixXd> lltExpression;
	Eigen::LLT<MatrixXd> lltGesture;
	Eigen::MatrixXd C_gestures;
	Eigen::MatrixXd C_P;
	Eigen::MatrixXd C_S;
	Eigen::VectorXd mean_shape;
	Eigen::VectorXd shape_without_pose;
	Eigen::VectorXd shape_coeff;
	Eigen::VectorXd pose_coeff;
	Eigen::VectorXd motion_param;
	Eigen::VectorXd ref_CM;
	Eigen::VectorXd CMwithpose;
	Eigen::VectorXd CMwithexpression;
	Eigen::VectorXd CMwithgesture;
	SparseLocalizedForConnectMap *spCM;
	PyObject * pModule_coresspondence;
	PyObject * pFuncLoadPointcloud;
	PyObject * pFuncFindCoress;
	PyObject * pointcloud;
	PyObject * pmaleModelMotion;
	PyObject * pFucLoadbodymotion;
	PyObject * pFuncMotion2DE;
	PyObject * pFuncPose2DE;
	PyObject * pFuncMotion2Face;
	PyObject * pFuncMotion2Pose;
	PyObject * pFuncPose2skeletonpos;
	PyObject * pFuncOptimizex;
	PyObject * pFuncFindRigid;
	PyObject * pFuncgetJointsPos;
	PyObject * pModule2;
	PyObject * pFunc2;
	PyObject * sessPose;
	PyObject * sessFace;
	PyObject * sessHand;
	PyObject * C_body;
	PyObject * m2c_model;//the motion parameter to pose coefficient model
	PyObject * ss_x;
	PyObject * ss_y;
	PyObject * bodyDEOutput;
	PyObject * sessPoseMotion;
	PyObject * tf_body_motion_x;
	PyObject * inputmotionfaceargs;
	PyObject* inputmotionargs;
	PyObject* inputmotionargs1;
	//body feature points的idx和权重weights和pos
	std::vector<std::vector<int>> body_feature_idx;
	std::vector<std::vector<double>> body_feature_weights;
	std::vector<std::vector<int>> body_feature_idx1;
	std::vector<std::vector<double>> body_feature_weights1;
	std::vector<Eigen::Vector3d> body_feature_pos;
	std::vector<Eigen::Vector3d> body_feature_pos1;
	std::vector<int> selectbodyfeature;
	//face featur points 的idx和权重weights和pos
	std::vector<std::vector<int>> face_feature_idx;
	std::vector<std::vector<double>> face_feature_weights;
	std::vector<Eigen::Vector3d> face_feature_pos;
	//hand feature points 的idx和权重weights和pos
	std::vector<std::vector<int>> hand_feature_idx;
	std::vector<std::vector<double>> hand_feature_weights;
	std::vector<Eigen::Vector3d> hand_feature_pos;
	//anchor point pos
	std::vector<Eigen::Vector3d> anchorpointpos;
};