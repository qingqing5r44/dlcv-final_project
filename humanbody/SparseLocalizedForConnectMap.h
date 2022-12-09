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
				// * ����ÿ�����ڲ�ͬ֡��Ӧ����������ڱߵ���ת�任����* //
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
			/**��ʼ��*/
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
			/**����ê��*/
			void saveanchors(vector<int> weight);
			/** ����ê��*/
			void setAnchor(unsigned int _idx);
			void setAnchor(vector<vector<int>> idx_set, vector<vector<double>> weights, vector<Eigen::Vector3d> pos);
			/** ���ê��*/
			void addAnchor(unsigned int _idx);
			/** ��ȡ��i֡��������������*/
			Eigen::VectorXd getCords (int ith);
			/** �����֡ÿ���������Լ�ÿ����ķ�����*/
			void calculateFaceArea();
			/** �����֡ÿ���ߵķ�����*/
			void calculateEdgesNormals();
			/** �����֡ÿ����ķ�����*/
			void calculateFacesNormals();
			Eigen::MatrixXd calculateFacesNormals(Eigen::VectorXd x);
			/** ���㵥֡������ķ�����*/
			Eigen::Vector3d calculateFaceNormal(int faceNo, int frameNo);
			/** ���㵥֡������ľֲ����*/
			faceanchor calculateFaceFrame(int faceNo, int frameNo);
			/** �����֡�����ߵľֲ����*/
			void calculateLocalFrame();
			/** �����֡������ľֲ����*/
			void calculateFacesFrame();
			/** ������ת�������ŷ����ת��*/
			void rotMat2EAngle(Eigen::Matrix3d  rotMat,Eigen::Vector3d & eulerAg);
			/** �����ֲ���ܵı任���������ת�����ת��*/
			void rotMat2AxisAngle(Eigen::Matrix3d rotMat, Eigen::Vector3d & axis, double & angle);
			/** ͨ��ŷ���Ǹ�����ת����*/
			void ERAngle2rotMat(Eigen::Vector3d  eulerAg, Eigen::Matrix3d  &rotMat);
			/** �����������е�����ӳ��*/
			void buildEdgeGraph();
			/** ���ӻ���i֡����ӳ��*/
			void visualizeConnectMap(int ith);
			/** ����IICs����*/
			void computeIICs();
			/** �������Ǻͱ߳�����*/
			void computeDiEdge();
			/** ����ĳ����Ķ���Ǻͱ߳�����*/
			Eigen::VectorXd computeDiEdge(Eigen::VectorXd x);
			/** ����ƽ��ģ��*/
			void computeMeanMesh();
			/** ϸ������*/
			void subdivision();
			/** ��������ӳ���ļ�*/
			void loadLA();
			/** ����ê������*/
			void loadAnchorData();
			/** ����ê������*/
			void loadAnchorData(string str1);
			/** ��ʼ���ع�*/
			void initialCongression();
			/** ���ƻع�ϵ��*/
			void estimateCoeff();
			/** �������ݵ��ļ���*/
			void loadData();
			/** ��ȡϸ�ڵ��ļ���*/
			void loadDetail();
			/** �ֲ�ϡ��ֽ�*/
			void sparseLocalizedDecomp(int kComp);
			/** ��ȡ��״��*/
			void getBasis(string str1, string str2);
			/** ���ӻ���״���ľֲ���Χ*/
			void visualizeBasisLocal (int kth);
			/** �ع�����*/
			void reconstruction(Eigen::VectorXd edgeComponent, Eigen::VectorXd & cordComponent, int anchorIthSeq, double radius, int ithComp);
			/** ���ü����Է�������ع�����*/
			void reconstruction1(Eigen::VectorXd edgeComponent, Eigen::VectorXd & cordComponent, int anchorIthSeq, double radius, int ithComp);
			/** ����ţ�ٷ��ع�����*/
			void reconstructionNewton(Eigen::VectorXd edgeComponent, Eigen::VectorXd & cordComponent, int anchorIthSeq, double radius, int ithComp);
			/** ����IIC�ع�����*/
			void reconstructionFromIICs(Eigen::VectorXd edgeComponent, Eigen::VectorXd & cordComponent, int anchorIthSeq, double radius, int ithComp);
			/** ���ݶ���Ǻͱ߳��ع�����*/
			void reconstructionFromDiEdges(Eigen::VectorXd edgeComponent, Eigen::VectorXd & cordComponent, int anchorIthSeq, double radius, int ithComp);
			/** ���ݵ����굼������*/
			void outMesh(Eigen::VectorXd cordComponent, char* path);
			void outMesh(Eigen::MatrixXd cordComponent, char* path);
			/** ���ݵ����굼������*/
			void outMesh1(Eigen::VectorXd cordComponent, string path);
			void outMesh1(Eigen::MatrixXd cordComponent, string path);
			/** ���ӻ�component(����ÿ���ߵ����)*/
			void visualizeComponent(Eigen::VectorXd cordComponent, Eigen::VectorXd edgeError, char* path);
			/** ��������Ȩ�أ���ģ����ɫ*/
			void colorize(std::vector<double> &errorss, Eigen::VectorXd cordComponent, char* path);
			/** ���ӻ�����()*/
			void visualizeMesh(Eigen::VectorXd cordComponent, int ith, char* path, double radius);
			/** ����IIC components�ع���ith֡*/
			void constructFrameFromComponents(int ith, int compNum, string str);
			void constructFrameFromComponentsError(int ith, int compNum, string str, Eigen::VectorXd & errorIncrement);
			/** ���õ�����W*/
			void setW_single();
			void setW_single(string str);
			/** ����multiple��W*/
			void setW_multiple();
			/** ����W_single�ع���״*/
			void reconstFormWsingle();
			/** �Ƚ�presolver*/
			void presolve();
			/** �Ƚ�presolver(����ê���ɵ�����Ի��)*/
			void presolve1();
			/** ���color�ļ�*/
			void outColor(Eigen::VectorXd comp, string pPath);
			/** ����w�����������*/
			double computeEnergy(Eigen::VectorXd Weights, Eigen::VectorXd v_ij);
			double computeDevOfWeights(Eigen::VectorXd Weights);
			void outMesh2(Eigen::VectorXd cordComponent, string path);
		public:
		/**/
			PGMesh *mesh_;//�ο�����
			PGMesh *mesh1;//Ϊ�˷�����㷨�������
			PGMesh *subdivMesh;//ϸ�ֺ������
			std::vector<std::vector<Eigen::VectorXd>> SpatialTemData;//������������
			std::vector<std::vector<Eigen::VectorXd>> ConnectMapData;//����ӳ������
			std::vector<std::vector<Eigen::VectorXd>> ConnectMapData1;//����ӳ������
			std::vector<std::vector<Eigen::Vector3d>> IICsData;//Isometry-invariant instrinsic��������
			Eigen::MatrixXd ConnectMapMatrix;
			Eigen::MatrixXd ConnectMapMatrix1;
			Eigen::MatrixXd IICsMatrix;
			Eigen::MatrixXd DiEdgeDataMatrix;
			Eigen::MatrixXd shapeBasis;
			Eigen::MatrixXd LBases;
			Eigen::MatrixXd DBases;
			Eigen::VectorXd LW_single;

			Eigen::MatrixXd W;
			int nSeq;//���е�֡��
			int nV;//��ĸ���
			int nF;//��ĸ���
			int nE;//�ߵĸ���
			int nSubV;
			int nSubF;
			int kCom;
			std::vector<unsigned int> edgeanchors_;//ê������
			std::vector<unsigned int> Vertexanchors_;//ê������
			std::vector<unsigned int> Faceanchors_;//ê������
			std::vector<EdgeGraphNode> edge_graph_;//����ӳ�������
			int anchorPointidx;//ê���index
			int anchorEidx;//ê����ڱߵ�index
			std::vector<std::vector<double>> FaceAreas;
			std::vector<std::vector<Eigen::Vector3d>> FaceNormals;//��֡��������淨����
			std::vector<std::vector<Eigen::Vector3d>> EdgeNormals;//��֡�����ߵı߷�����
			std::vector<std::vector<Eigen::Matrix3d>> EdgeLocalFrames;//��֡�����ߵľֲ����
			std::vector<std::vector<Eigen::Matrix3d>> FaceLocalFrames;//��֡������ľֲ����
			std::vector<std::vector<faceanchor>> faceAnchorFrames;//����ÿ֡��ê��ı��
			std::vector<std::vector<vertexanchor>> vertexAnchorCords;//����ÿ֡��ê�������
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
