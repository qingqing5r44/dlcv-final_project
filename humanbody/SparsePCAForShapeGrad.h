#include <OpenMesh/Core/IO/MeshIO.hh>  
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>  
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <unsupported/Eigen/MatrixFunctions>
#include <iostream>
#include <fstream>
#include <map>
//#include <QDebug>
#include <math.h>
#include <iomanip>


using namespace std;
using namespace Eigen;


		class SparsePCAForShapeGrad
		{
			typedef OpenMesh::PolyMesh_ArrayKernelT<> PGMesh;
		public:
			SparsePCAForShapeGrad (PGMesh * _mesh, std::vector<std::vector<Eigen::VectorXd>> & sd);
			SparsePCAForShapeGrad(PGMesh *_mesh);
			SparsePCAForShapeGrad(PGMesh *_mesh1, PGMesh *_mesh2);
			bool setDeformedMesh(std::vector<Eigen::VectorXd> & sd);
			void setDeformedMesh(PGMesh *_mesh);
			void computeDG();
			void computeCM();
			void writeCM(char* path);
			void build_matrix();
			void computeS_matrix();
			void compute_V_inv();
			void vecToMat(Eigen::VectorXd CMComponent, std::vector<Eigen::MatrixXd> &R_ij, std::vector<Eigen::MatrixXd> &S_i);
			Eigen::VectorXd reconstructFromRS(PGMesh * _mesh, std::vector<Eigen::MatrixXd> R_ij, std::vector<Eigen::MatrixXd> S_i);
			Eigen::VectorXd reconstrtuction(PGMesh * _mesh, Eigen::VectorXd gradientComponent, int frameNumber);		//后面是帧数以及对应的mesh
			bool reconstrtuction();
			Eigen::VectorXd getOneComponent( int i );
			Eigen::VectorXd getComponent();
			bool writeDeformGrad( string& filepath );
			bool writeSubdivMesh_vertex( string & filePath );			//output the vertex
			bool writeSubdivMesh_face( string & filePath );			//output the faces by listing their vertexs
			bool writeMeshProperty( string & filePath );					//output the mesh property
			bool write_nSeq_nF( string & filePath );
			bool writeVisualizeComponent_color( string & filePath );			//calculate the value of color base on components
			bool writeError( string & filePath );											// calculate the error of every face base on Xmean
			bool writeFixFacePoints();														//output the fix points
			bool readFixFacePoints();														// input the fix points;
			bool readKComponent( string & filePath );
			bool readW( string & filePath );
			bool readXmean( string & filePath );
			bool readDeformGrad();
			bool visualizeComponent( int i_th );			//visualize the i th component on the mesh
			void getColor( double color, double& r, double& g, double& b );
			bool doSPLOCS_py( string & deformGradPath, string & subdiviVertexPath, 
				string & subdiviFacePath , string & meshPropertyPath, int ndim );
			bool setK( int k );
			void setWs( std::vector< double >::iterator begin, std::vector< double >::iterator end );
			void setW_single();
			bool set_nSeq( int nSeq );
			/** 提取第i帧的网格坐标向量*/
			Eigen::VectorXd getCords (int ith);
			/** 可视化component(给定每条边的误差)*/
			void visualizeComp(Eigen::VectorXd cordComponent, Eigen::VectorXd edgeError, char* path);
			void subdivision();
			void writeFaceErrorOfComp(char * filePath, Eigen::VectorXd vertexCords, Eigen::VectorXd face_error);
			void outMesh( Eigen::VectorXd cordComponent, string path);
			void outInverFaceMesh();

			void reconstFormWsingle();
			void setAnchorPoints(std::vector<int> & _fixPointIdx, std::vector<Eigen::VectorXd> & _fixPointPos);
			void setAnchorR(std::vector<int>& _fixFaceIdx, std::vector<Eigen::MatrixXd>& _fixFaceR);
		public:
			/**/
			PGMesh *mesh_;
			PGMesh subdivMesh;
			OpenMesh::FPropHandleT< int > cogs;										//mesh property for map a face to the center of it's subdivide faces
			std::vector<std::vector<Eigen::VectorXd>> SpatialTemData;		//vertexIndex sequence vertex
			std::vector<std::vector<Eigen::VectorXd>> DeformGrad;			//face sequence vector
			std::vector<Eigen::MatrixXd> Rs;
			std::vector<Eigen::MatrixXd> R_ijs;
			std::vector<Eigen::MatrixXd> Ss;
			std::vector<Eigen::MatrixXd> DGs;
			std::vector<Eigen::VectorXd> CMs;
			Eigen::MatrixXd KComp;														//k components  k * 9nF
			std::vector< Eigen::VectorXd > W;										// nSeq * k
			Eigen::VectorXd W_single;													// single W
			Eigen::VectorXd Xmean;														//the mean of all deformation gradients
			Eigen::VectorXd settedW;														// W's value is setted by interface
			Eigen::MatrixXd DFMatrix;
			int nSeq;		//number of frames
			int nE;   // number of edges
			int nV;   //number of vertices
			int nF;			//number of faces
			int k;				// the number of component
			static int  fixFacesNum;				//the number of fix faces
			std::vector< Eigen::VectorXd > fixFacePoints;
			std::vector<int> fixPointIdx;
			std::vector<Eigen::VectorXd> fixPointPos;
			std::vector<int> fixFaceIdx;
			std::vector<Eigen::MatrixXd> fixFaceR;
			std::vector<Eigen::MatrixXd> V_inv;
			SimplicialCholesky<SparseMatrix<double>> Vsolver;
			SimplicialCholesky<SparseMatrix<double>> Ssolver;
			SimplicialCholesky<SparseMatrix<double>> Rsolver;
			Eigen::SparseMatrix<double> DGCoef;
			Eigen::SparseMatrix<double> SSCoef;
			Eigen::VectorXd VecCM;
		};
