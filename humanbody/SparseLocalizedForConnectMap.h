#ifndef SPARSE_LOCALIZED_FOR_CONNECT_MAP_H
#define SPARSE_LOCALIZED_FOR_CONNECT_MAP_H
#include <OpenMesh/Core/IO/MeshIO.hh>  
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>  
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <iostream>
#include <fstream>
#include "GeometryMath.h"
#include <math.h>
//#include "GeometryMath.h"
//#define MS_NO_COREDLL
//#include "Python.h"
#include "time.h"

//#include <unsupported/Eigen/MatrixFunctions>
using namespace std;
using namespace Eigen;


		class SparseLocalizedForConnectMap
		{
			typedef OpenMesh::PolyMesh_ArrayKernelT<> PGMesh;
			typedef PGMesh::EdgeHandle EdgeHandle;
			typedef PGMesh::VertexHandle VertexHandle;
			typedef Eigen::VectorXd VectorXd;
		protected:
			struct EdgeGraphNode
			{
				unsigned int center_node_id;
				// * 代表每条边在不同帧对应的相对于相邻边的旋转变换矩阵* //
				std::vector<std::vector<Eigen::Matrix3d>> one_ring_local_rot;
				std::vector<std::vector<Eigen::Matrix3d>> one_ring_res_rot;
				//std::vector<>
				//std::vector<OpenMesh::Vec3d> one_ring_bending;
				std::vector<unsigned int> one_ring_idx;
				std::vector<std::vector<double>> vecscaling;
				std::vector<double> edgeLength;
				double currentLength;
				std::vector<Eigen::Matrix3d> current_one_ring_local_rot;
				std::vector<std::vector<Eigen::Vector3d>> rotationaxis;
			};
			struct regreCoef
			{
				unsigned int center_node_id;
				std::vector<unsigned int> one_ring_idx;
				Eigen::VectorXd one_ring_Coef;
			};
			struct faceanchor
			{
				int frameNo;
				int idx;
				Eigen::Matrix3d localframe;
			};
			struct vertexanchor
			{
				int frameNo;
				int idx;
				vector<int> idx_list;
				vector<double> weights;
				Eigen::VectorXd cord;
			};
		public:
			VertexHandle POINT0_HANDLE(EdgeHandle edgeHandle) {return mesh_->from_vertex_handle( mesh_->halfedge_handle(edgeHandle, 0) ); }
			VertexHandle POINT1_HANDLE(EdgeHandle edgeHandle) {return mesh_->to_vertex_handle( mesh_->halfedge_handle(edgeHandle, 0) ); }
			VectorXd POINT(VertexHandle vertexhandle, int frame ) {return SpatialTemData[vertexhandle.idx()][frame];}
			VectorXd EDGE_POINT0(EdgeHandle edgeHandle, int frame) { return POINT(POINT0_HANDLE(edgeHandle), frame);}
			VectorXd EDGE_POINT1(EdgeHandle edgeHandle, int frame) { return POINT(POINT1_HANDLE(edgeHandle), frame);}

			VectorXd POINT(VertexHandle vertexhandle, Eigen::VectorXd x) { return x.segment(3 * vertexhandle.idx(),  3); }
			VectorXd EDGE_POINT0(EdgeHandle edgeHandle, Eigen::VectorXd x) { return POINT(POINT0_HANDLE(edgeHandle), x); }
			VectorXd EDGE_POINT1(EdgeHandle edgeHandle, Eigen::VectorXd x) { return POINT(POINT1_HANDLE(edgeHandle), x); }
			VectorXd VECTOR(EdgeHandle edgeHandle, Eigen::VectorXd x) { return EDGE_POINT1(edgeHandle, x) - EDGE_POINT0(edgeHandle, x); }

			VectorXd VECTOR(EdgeHandle edgeHandle, int frame) {return EDGE_POINT1( edgeHandle,  frame) - EDGE_POINT0( edgeHandle,  frame);}
			Eigen::Vector3d computeCross (Eigen::Vector3d x1, Eigen::Vector3d x2) 
			{
				Eigen::Vector3d x;
				x(0) = x1(1) * x2(2) - x1(2) * x2(1);
				x(1) = -x1(0) * x2(2) + x2(0) * x1(2);
				x(2) = x1(0) * x2(1) - x1(1) * x2(0);
				return x;
			}
			Eigen::Vector3d computeCross(PGMesh::Normal x1, PGMesh::Normal x2)
			{
				Eigen::Vector3d x;
				x(0) = x1[1] * x2[2] - x1[2] * x2[1];
				x(1) = -x1[0] * x2[2] + x2[0] * x1[2];
				x(2) = x1[0] * x2[1] - x1[1] * x2[0];
				return x;
			}
			OpenMesh::Vec3d eigenTransToVec3d (Eigen::Vector3d x){OpenMesh::Vec3d y; y[0] = x(0); y[1] = x[1]; y[2] = x[2]; return y;}
			Eigen::Vector3d Vec3dToEigen(OpenMesh::Vec3d x){Eigen::Vector3d y; y(0) = x[0]; y(1) = x[1]; y[2]= x[2]; return y;}
		public:
			/**初始化*/
			SparseLocalizedForConnectMap(PGMesh *_mesh, std::vector<std::vector<Eigen::VectorXd>> &sd);
			SparseLocalizedForConnectMap(PGMesh *_mesh);
			SparseLocalizedForConnectMap(PGMesh *_mesh1, PGMesh *_mesh2);

			void updateDeformed(PGMesh * _mesh);
			void updateDeformed1(PGMesh * _mesh);
			void updateDeformed1(string path_);
			void updateDeformed1(std::vector<double> anchorfaces);
			void updateMeshVertices(Eigen::VectorXd x, PGMesh *_mesh);
			void setMesh(PGMesh * _mesh);
			void setMesh(Eigen::VectorXd V);
			void setMesh(Eigen::MatrixXd V);
			Eigen::VectorXd getVertexMatrix(PGMesh * _mesh);
			/**保存锚点*/
			void saveanchors(vector<int> weight);
			/** 设置锚点*/
			void setAnchor(unsigned int _idx);
			void setAnchor(vector<vector<int>> idx_set, vector<vector<double>> weights, vector<Eigen::Vector3d> pos);
			/** 添加锚点*/
			void addAnchor(unsigned int _idx);
			/** 提取第i帧的网格坐标向量*/
			Eigen::VectorXd getCords (int ith);
			/** 计算各帧每个面的面积以及每个面的法向量*/
			void calculateFaceArea();
			/** 计算各帧每条边的法向量*/
			void calculateEdgesNormals();
			/** 计算各帧每个面的法向量*/
			void calculateFacesNormals();
			Eigen::MatrixXd calculateFacesNormals(Eigen::VectorXd x);
			/** 计算单帧单个面的法向量*/
			Eigen::Vector3d calculateFaceNormal(int faceNo, int frameNo);
			/** 计算单帧单个面的局部标架*/
			faceanchor calculateFaceFrame(int faceNo, int frameNo);
			/** 计算各帧各条边的局部标架*/
			void calculateLocalFrame();
			/** 计算各帧各个面的局部标架*/
			void calculateFacesFrame();
			/** 给出旋转矩阵计算欧拉旋转角*/
			void rotMat2EAngle(Eigen::Matrix3d  rotMat,Eigen::Vector3d & eulerAg);
			/** 给出局部标架的变换矩阵计算旋转轴和旋转角*/
			void rotMat2AxisAngle(Eigen::Matrix3d rotMat, Eigen::Vector3d & axis, double & angle);
			/** 通过欧拉角给出旋转矩阵*/
			void ERAngle2rotMat(Eigen::Vector3d  eulerAg, Eigen::Matrix3d  &rotMat);
			/** 计算网格序列的连接映射*/
			void buildEdgeGraph();
			/** 可视化第i帧连接映射*/
			void visualizeConnectMap(int ith);
			/** 计算IICs数据*/
			void computeIICs();
			/** 计算二面角和边长数据*/
			void computeDiEdge();
			/** 计算某网格的二面角和边长数据*/
			Eigen::VectorXd computeDiEdge(Eigen::VectorXd x);
			/** 计算平均模型*/
			void computeMeanMesh();
			/** 细分网格*/
			void subdivision();
			/** 载入连接映射文件*/
			void loadLA();
			/** 载入锚点数据*/
			void loadAnchorData();
			/** 载入锚点数据*/
			void loadAnchorData(string str1);
			/** 初始化回归*/
			void initialCongression();
			/** 估计回归系数*/
			void estimateCoeff();
			/** 载入数据到文件中*/
			void loadData();
			/** 提取细节到文件中*/
			void loadDetail();
			/** 局部稀疏分解*/
			void sparseLocalizedDecomp(int kComp);
			/** 提取形状基*/
			void getBasis(string str1, string str2);
			/** 可视化形状基的局部范围*/
			void visualizeBasisLocal (int kth);
			/** 重构网格*/
			void reconstruction(Eigen::VectorXd edgeComponent, Eigen::VectorXd & cordComponent, int anchorIthSeq, double radius, int ithComp);
			/** 利用简化线性方程求解重构网格*/
			void reconstruction1(Eigen::VectorXd edgeComponent, Eigen::VectorXd & cordComponent, int anchorIthSeq, double radius, int ithComp);
			/** 利用牛顿法重构网格*/
			void reconstructionNewton(Eigen::VectorXd edgeComponent, Eigen::VectorXd & cordComponent, int anchorIthSeq, double radius, int ithComp);
			/** 根据IIC重构网格*/
			void reconstructionFromIICs(Eigen::VectorXd edgeComponent, Eigen::VectorXd & cordComponent, int anchorIthSeq, double radius, int ithComp);
			/** 根据二面角和边长重构网格*/
			void reconstructionFromDiEdges(Eigen::VectorXd edgeComponent, Eigen::VectorXd & cordComponent, int anchorIthSeq, double radius, int ithComp);
			/** 根据点坐标导出网格*/
			void outMesh(Eigen::VectorXd cordComponent, char* path);
			void outMesh(Eigen::MatrixXd cordComponent, char* path);
			/** 根据点坐标导出网格*/
			void outMesh1(Eigen::VectorXd cordComponent, string path);
			void outMesh1(Eigen::MatrixXd cordComponent, string path);
			/** 可视化component(给定每条边的误差)*/
			void visualizeComponent(Eigen::VectorXd cordComponent, Eigen::VectorXd edgeError, char* path);
			/** 给出顶点权重，给模型上色*/
			void colorize(std::vector<double> &errorss, Eigen::VectorXd cordComponent, char* path);
			/** 可视化网格()*/
			void visualizeMesh(Eigen::VectorXd cordComponent, int ith, char* path, double radius);
			/** 利用IIC components重构第ith帧*/
			void constructFrameFromComponents(int ith, int compNum, string str);
			void constructFrameFromComponentsError(int ith, int compNum, string str, Eigen::VectorXd & errorIncrement);
			/** 设置单独的W*/
			void setW_single();
			void setW_single(string str);
			/** 设置multiple的W*/
			void setW_multiple();
			/** 根据W_single重构形状*/
			void reconstFormWsingle();
			/** 先解presolver*/
			void presolve();
			/** 先解presolver(考虑锚点由点的线性混合)*/
			void presolve1();
			/** 输出color文件*/
			void outColor(Eigen::VectorXd comp, string pPath);
			/** 输入w，计算出能量*/
			double computeEnergy(Eigen::VectorXd Weights, Eigen::VectorXd v_ij);
			double computeDevOfWeights(Eigen::VectorXd Weights);
			void outMesh2(Eigen::VectorXd cordComponent, string path);
		public:
		/**/
			PGMesh *mesh_;//参考网格
			PGMesh *mesh1;//为了方便计算法向的向量
			PGMesh *subdivMesh;//细分后的网格
			std::vector<std::vector<Eigen::VectorXd>> SpatialTemData;//序列坐标数据
			std::vector<std::vector<Eigen::VectorXd>> ConnectMapData;//连接映射数据
			std::vector<std::vector<Eigen::VectorXd>> ConnectMapData1;//连接映射数据
			std::vector<std::vector<Eigen::Vector3d>> IICsData;//Isometry-invariant instrinsic坐标数据
			Eigen::MatrixXd ConnectMapMatrix;
			Eigen::MatrixXd ConnectMapMatrix1;
			Eigen::MatrixXd IICsMatrix;
			Eigen::MatrixXd DiEdgeDataMatrix;
			Eigen::MatrixXd shapeBasis;
			Eigen::MatrixXd LBases;
			Eigen::MatrixXd DBases;
			Eigen::VectorXd LW_single;

			Eigen::MatrixXd W;
			int nSeq;//序列的帧数
			int nV;//点的个数
			int nF;//面的个数
			int nE;//边的个数
			int nSubV;
			int nSubF;
			int kCom;
			std::vector<unsigned int> edgeanchors_;//锚点向量
			std::vector<unsigned int> Vertexanchors_;//锚点向量
			std::vector<unsigned int> Faceanchors_;//锚点向量
			std::vector<EdgeGraphNode> edge_graph_;//连接映射的向量
			int anchorPointidx;//锚点的index
			int anchorEidx;//锚点的邻边的index
			std::vector<std::vector<double>> FaceAreas;
			std::vector<std::vector<Eigen::Vector3d>> FaceNormals;//各帧各个面的面法向量
			std::vector<std::vector<Eigen::Vector3d>> EdgeNormals;//各帧各条边的边法向量
			std::vector<std::vector<Eigen::Matrix3d>> EdgeLocalFrames;//各帧各条边的局部标架
			std::vector<std::vector<Eigen::Matrix3d>> FaceLocalFrames;//各帧各个面的局部标架
			std::vector<std::vector<faceanchor>> faceAnchorFrames;//保存每帧面锚点的标架
			std::vector<std::vector<vertexanchor>> vertexAnchorCords;//保存每帧点锚点的坐标
			OpenMesh::EPropHandleT<int> edgeIndex;
			OpenMesh::VPropHandleT< PGMesh::VertexHandle > vhandle;
			OpenMesh::EPropHandleT< PGMesh::VertexHandle> ehandle;
			std::vector<std::vector<Eigen::Matrix3d>> LocalTransforms;
			std::vector<regreCoef> edge_regressions;
			Eigen::VectorXd W_single;
			SimplicialCholesky<SparseMatrix<double>> presolver;
			Eigen::SparseMatrix<double> vertexCoff;
			std::vector<Eigen::Matrix3d> currentFaceFrames;
			std::vector<Eigen::Matrix3d> currentQs;
		};

#endif
