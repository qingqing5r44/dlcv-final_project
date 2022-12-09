#pragma once
#include "global.h"
#include <iostream>
#include <fstream>
using namespace std;

namespace mns
{
	struct MyPoint
	{
		integer idx;
		decimal x;
		decimal y;
		decimal z;
		decimal nx;
		decimal ny;
		decimal nz;
		int_vector neighbor;
	};

	struct MyPointSet
	{
		std::vector<MyPoint> set;
	};

	struct MyEdge
	{
		decimal length;
		int_vector points;
	};

	struct MyEdgeSet
	{
		std::vector<MyEdge> set;
	};

	struct MyHalfEdge
	{
		decimal length;
		integer from_vertex;
		integer to_vertex;
	};

	struct MyHalfEdgeSet
	{
		std::vector<MyHalfEdge> set;
	};

	struct MyFace
	{
		decimal area;
		decimal nx;
		decimal ny;
		decimal nz;
		int_vector points;
	};

	struct MyFaceSet
	{
		std::vector<MyFace> set;
	};

	struct MyEigen
	{
		std::complex<decimal> val;
		Eigen::VectorXd vec;
	};

	class MyDeformationHandler
	{
	public:		          //(mesh,getSimMesh(),wei_c or not,sim_or_not, max_it,          stp,     bw_or_not,     eg_or_not,   ba_or_not,       e_str,       n_eg,        ep1,             ep2,              b_m_num, b_col,              contri,          a_c,          a_s,           a_b,          a_v           , a_w);
		MyDeformationHandler(std::vector<std::vector<int>>, std::vector<std::vector<double>>, MyMesh&, MyMesh&, const bool, const bool, const integer, const decimal, const bool, const bool, const bool, const std::string, const integer, const decimal, const decimal, const integer, const integer, const decimal, const decimal, const decimal, const decimal, const decimal, const decimal);

		MyDeformationHandler(MyMesh&, MyMesh&, const bool, const bool, const integer, const decimal, const bool, const bool, const bool, const std::string, const integer, const decimal, const decimal, const integer, const integer, const decimal, const decimal, const decimal, const decimal, const decimal, const decimal);
		~MyDeformationHandler();

		MyPointSet GetConSetByInput(string );
		MyPointSet GetConSetByInput();
		MyPointSet GetConSetByDif();
		MyPointSet GetConSetByV(const std::vector<integer>&, const Eigen::MatrixXd&);
		MyPointSet GetConSetByMesh();
		MyPointSet GetConSetByMesh(std::string, decimal);

		void CreateConstraintPointSet(const MyPointSet&);
		void CreateConstraintPointSet1(const MyPointSet&);
		//void AdaptConSetToSim();
		bool CheckMovedPt(MyMesh::Point, Eigen::RowVector3d);
		//void PrintPtControl();

		decimal GetConTerm(integer);
		decimal GetConTerm(integer, decimal);
		decimal GetConTermNew(integer id, decimal fac);
		decimal GetFixConTermNew(integer id, decimal fac);
		decimal Ec();
		decimal EcNew();
		void GetFeaturePos(std::vector<Eigen::Vector3d>);
		void setFixed(std::vector<std::vector<int>>, std::vector<std::vector<double>>, std::vector<Eigen::Vector3d>);
		void setconset(std::vector<std::vector<int>> idxset, std::vector<std::vector<double>> weightset, std::vector<Eigen::Vector3d> posset)
		{
			decimal wc = alpha_c / idxset.size();
			swc = std::sqrt(wc);
			fixed_feature_idx.resize(0);
			fixed_feature_pos.resize(0);
			fixed_feature_weights.resize(0);
			body_feature_idx = idxset;
			body_feature_weights = weightset;
			body_feature_pos = posset;
		}
		void addconset(std::vector<std::vector<int>> idxset, std::vector<std::vector<double>> weightset, std::vector<Eigen::Vector3d> posset)
		{
			decimal wc = alpha_c / (idxset.size() + body_feature_idx.size()+ fixed_feature_idx.size());
			swc = std::sqrt(wc);
			int beforecount = fixed_feature_idx.size();
			fixed_feature_idx.resize(body_feature_idx.size() + fixed_feature_idx.size());
			fixed_feature_weights.resize(body_feature_idx.size() + fixed_feature_idx.size());
			fixed_feature_pos.resize(body_feature_idx.size() + fixed_feature_idx.size());
			for (int i = 0; i < body_feature_idx.size(); i++)
			{
				fixed_feature_idx[beforecount +i]=(body_feature_idx[i]);
				fixed_feature_weights[beforecount + i] = (body_feature_weights[i]);
				fixed_feature_pos[beforecount + i] = (body_feature_pos[i]);
			}
			body_feature_idx = idxset;
			body_feature_weights = weightset;
			body_feature_pos = posset;
		}
		decimal GetStrTerm(integer, integer, integer);
		decimal GetStrTerm1(integer from_id, integer to_id, integer e_id);
		decimal Es();

		decimal CalDiheral(const MyMesh::Point&, const MyMesh::Point&, const MyMesh::Point&, const MyMesh::Point&);
		decimal CalWLZHDiheral(const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&);
		decimal GetBenTerm(integer, integer, integer, integer, integer);
		
		decimal Eb();

		std::vector<decimal> Esb();
		decimal GetBenTermNew1(integer, integer, integer, integer, integer, const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&);
		decimal GetBenTermNew(integer, integer, integer, integer, integer, const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&);
		decimal EbNew();

		std::vector<decimal> EsbNew();
		std::vector<decimal> EsbNew1();
		decimal Ev();

		decimal GetWWTerm(integer);
		decimal Ew();

		void EgInit();

		//void CreateBasis(integer);
		//void BasisPCA( const Eigen::MatrixXd& );
		void BuildWntByBasis();
		void CalculateWeightForBasis(std::string);
		void UpdateWntVecWithBasis();
		void UpdateWntVecWithBasis2();
		void UpdateWntVecWithBasis3();
		void UpdateWntVecWithBasis3(Eigen::VectorXd cm);
		void UpdateWntVecWithBasis7();
		void BuildWntVecByBasis();
		void SetRBFlag(bool);
		void SeparateVerts(const Eigen::MatrixXd&);
		void BuildWntVecByBasis2();

		//void BasisSVD( const Eigen::MatrixXd& );
		////void createSparseBasis(Eigen::MatrixXd&, Eigen::MatrixXd&, Eigen::MatrixXd&, Eigen::MatrixXd&, std::vector<decimal>&, std::vector<decimal>&);

		void UpdateWntVec();
		void Transform2Example();

		void PrintLengthDihedral();
		void CreateDiLenMat();
		void CreateDiLenMat2();

		void ColorWeightedEdge(Eigen::MatrixXd&, decimal);
		void ColorWeightedEdge1(Eigen::MatrixXd&, decimal);
		void ColorWeightedEdge2(Eigen::MatrixXd&, decimal);
		void ColorWeightedEdge3(Eigen::MatrixXd&, decimal);
		void ColorWeightedEdge4(Eigen::MatrixXd&, decimal);
		void ColorComponent(Eigen::MatrixXd&, decimal);
		void ColorSpecificComponent(Eigen::MatrixXd&, integer, decimal);
		void ColorSpecificComponent2(Eigen::MatrixXd&, integer, decimal);
		bool ColorBwDecision();

		void BuildJmat(SMatXd&, SMatXd&, bool);

		void Deform(std::vector<integer>&, Eigen::MatrixXd&, Eigen::MatrixXd&);
		void Deform1(std::vector<integer>&, Eigen::MatrixXd&, Eigen::MatrixXd&);
		void DeformLD(std::vector<integer>&, Eigen::MatrixXd&, Eigen::MatrixXd&);
		void Deform_new(std::vector<integer>&, Eigen::MatrixXd&, Eigen::MatrixXd&);
		void DeformInitLD(std::vector<integer>&, Eigen::MatrixXd&, Eigen::MatrixXd&);
		void DeformInit(std::vector<integer>&, Eigen::MatrixXd&, Eigen::MatrixXd&);
		void DeformInit5(std::vector<integer>&, Eigen::MatrixXd&, Eigen::MatrixXd&);
		void DeformInit1(std::vector<integer>&, Eigen::MatrixXd&, Eigen::MatrixXd&);
		void ModMesh();
		void ModOriginalMesh(MyMesh&);
		bool SaveDeformedMesh();
		void ModVMat(Eigen::MatrixXd&);
		void VmatModAnotherMesh(const Eigen::MatrixXd&);
		void ClearAfterDeform();

		////cmsparse
		void Subdivision();
		void LoadData();
		void ReceiveBasis();
		void ReceiveBasis2();
		void ReceiveBasis2(string p_path);
		void ReceiveBasis3();
		void ReceiveBasis4();
		void ReceiveBasisWithName(std::string, integer);
		integer GetBasisCols();
		bool GetRBF();

		////subspace design
		void SetSubspaceW(const Eigen::MatrixXd&, decimal);
		bool GetUseSim();

		////
		void SetModelMod(integer);

		////weight cal
		void CalWeightByBasis();

		////calculate blend weight
		void CalWunknownByBasis();
		void CalWunknownByBasis2();

		////alpha change
		void Alpha_C_Change(decimal);
		void Alpha_S_Change(decimal);
		void Alpha_B_Change(decimal);
		void Alpha_V_Change(decimal);
		void Alpha_W_Change(decimal);
		void Clear4ReCal();
		bool Is_Alpha_Change();

		////error test
		//void ErrorEstimation();
		//decimal SpatialErrorEstimation();
		//decimal TemporalErrorEstimation(integer, decimal);
		//void STEDEstimation();

		////calculate optimal step length
		decimal CalStepLength(decimal, SMatXd&, Eigen::VectorXd&, Eigen::VectorXd&, Eigen::VectorXd&);

		decimal CSL_GetConTerm(integer, Eigen::VectorXd&);
		decimal CSL_Ec(Eigen::VectorXd&);

		decimal CSL_GetStrTerm(integer, integer, integer, Eigen::VectorXd&);

		decimal CSL_GetBenTerm(integer, integer, integer, integer, integer, Eigen::VectorXd&);

		std::vector<decimal> CSL_Esb(Eigen::VectorXd&);

		decimal CSL_Ev(Eigen::VectorXd&);

		decimal CSL_GetWWTerm(integer, Eigen::VectorXd&);
		decimal CSL_Ew(Eigen::VectorXd&);


		//Pose weight estimation
		void MeshWeightEst();
		void PoseWeightEst();
		void PoseAndWeightCompare(std::string, std::string, integer, integer, std::string);
	public:
		Eigen::VectorXd cm1;
		double econ_threshold;
		MyMesh mesh;
		MyMesh ori_mesh;
		MyMesh sim_mesh;
		MyMesh ex_mesh;
		MyMesh* eg_mesh;
		integer iter;
		integer max_iter;
		decimal step;
		decimal ori_step;

		bool use_sim_or_not;

		Eigen::VectorXd unknown;
		Eigen::VectorXd l_bw_unknown;
		Eigen::VectorXd d_bw_unknown;
		Eigen::VectorXd pre_unknown;

		decimal epsilon1;
		decimal epsilon2;

		int_vector idx_move;

		//importance of each term
		decimal alpha_c;
		decimal alpha_s;
		decimal alpha_b;
		decimal alpha_v;
		decimal alpha_w;

		//lamda_weight (importance for each term)
		decimal swc;
		decimal swv;
		dec_vector sws_set;
		dec_vector swb_set;
		decimal sww;

		//lower case f
		Eigen::VectorXd f;
		Eigen::VectorXd fc;
		Eigen::VectorXd fs;
		Eigen::VectorXd fb;
		Eigen::VectorXd fw;
		decimal fv;

		Eigen::VectorXd fsb;

		//constraint
		integer con_ele_cnt; // total constraint-vertex number
		MyPointSet cps;

		int_vector con_ori_map;
		integer move_pt_num;
		std::vector<int> pt_move_or_not;

		//stretching
		integer str_ele_cnt; // total edge-vertex number
		dec_vector ori_length_set;

		//bending
		integer ben_ele_cnt; // total edge-vertex number
		dec_vector ori_wlz_dihedral_set;
		dec_vector ori_dihedral_set;
		dec_vector ori_area_set;

		//s&b
		integer sb_ele_cnt; //2 times total edge-vertex number

							//volume
		Eigen::VectorXd u_vec;
		decimal ori_volume;

		//wanted paras
		dec_vector wnt_length_set;
		dec_vector wnt_wlz_dihedral_set;
		dec_vector wnt_dihedral_set;
		decimal wnt_volume;

		//example
		std::string example_str;
		bool use_eg_or_not;
		integer eg_num;
		std::vector<dec_vector> eg_length_vec;
		std::vector<dec_vector> eg_dihedral_vec;
		std::vector<decimal> eg_volume_vec;

		//basis
		bool use_basis_or_not;
		integer basis_col_num;
		integer create_basis_mesh_num;
		decimal contribution;
		Eigen::VectorXd data_avg;
		Eigen::MatrixXd basis;
		Eigen::VectorXd ori_len_di_mat;
		Eigen::MatrixXd sparse_basis;

		Eigen::MatrixXd local_basis;

		std::vector<dec_vector> length_basis_vec;
		std::vector<dec_vector> dihedral_basis_vec;

		bool rb_flag;// receive_basis_flag, 1 for have been set, 0 for not set

					 //use weight constraint for basis
		integer ww_ele_cnt;
		bool use_ww_or_not;

		//geodesic distance
		decimal** edge_geodis_mat;
		decimal** normalized_edge_geodis_mat;
		decimal largest_edge_geodis;

		//sparse container
		std::vector<T> triple;
		std::vector<T> j_triple;

		integer row_cnt; // total row

		decimal pre_energy;

		//cmsparse
		MyMesh subdivMesh;

		OpenMesh::VPropHandleT< MyMesh::VertexHandle > vhandle;
		OpenMesh::EPropHandleT< MyMesh::VertexHandle > ehandle;
		OpenMesh::EPropHandleT<integer> edgeIndex;
	public:
		string current_str;
		//subspace design
		Eigen::MatrixXd W_mat;
		SMatXd huge_W_mat;
		decimal limit;
		Eigen::MatrixXd high_v_mat;
		void outMesh1(Eigen::MatrixXd cordComponent, char* path);
		std::vector<std::vector<integer>> W_idx_set;
		dec_vector sign_flag;

		integer model_mod;
		int fitmode;
		//alpha change
		bool alpha_c_ch;
		bool alpha_s_ch;
		bool alpha_b_ch;
		bool alpha_v_ch;
		bool alpha_w_ch;
		bool need_ch[5];
		const static integer NEED_LEN = 5;
		decimal need_fac;
		integer ch_con;
		bool out_flag;
		//bool alpha_ch;

		//calculate blend weight
		bool cal_bw_or_not;
		std::vector<std::vector<int>> body_feature_idx;
		std::vector<std::vector<double>> body_feature_weights;
		std::vector<Eigen::Vector3d> body_feature_pos;
		std::vector<std::vector<int>> fixed_feature_idx;
		std::vector<std::vector<double>> fixed_feature_weights;
		std::vector<Eigen::Vector3d> fixed_feature_pos;
		std::vector<int> selectbodyfeature;
	};

	class MyMeshHandler
	{
	public:
		int fitmode;
		MyDeformationHandler* p_DeformHandler;

		MyMeshHandler(std::string, int, char**);
		MyMeshHandler(std::string, std::string, int, char**);
		MyMeshHandler(std::string, std::vector<std::vector<int>> , std::vector<std::vector<double>> );
		~MyMeshHandler();

		int argc;
		char** argv;

		bool SaveMesh();

		void ShowVertex() const;
		void ShowOneVertex(const int) const;
		void ShowEdge() const;
		void ShowHalfEdge() const;
		void ShowOneHalfEdge(const int) const;
		void ShowVertexNormal() const;
		void ShowFaceNormal() const;
		void ShowWedge();
		void ShowFaceHalfEdge();
		void ShowOneRingVertex(const int);

		void OutputInfo();

		void SetCollapse(const int, const int);
		void RemoveOutliers();
		void Bilateral(decimal, decimal, integer);

		void Deformation(std::vector<integer>&, Eigen::MatrixXd&, Eigen::MatrixXd&);
		void DeformationLD(std::vector<integer>&, Eigen::MatrixXd&, Eigen::MatrixXd&);

		MyMesh& GetMesh();
		MyPoint GetPoint(const int);
		MyPointSet GetPS();

		MyEdge GetEdge(const int);
		MyEdgeSet GetES();

		MyHalfEdge GetHalfEdge(const int);
		MyHalfEdgeSet GetHES();

		MyFace GetFace(const int);
		MyFaceSet GetFS();

		MyMesh& GetSimMesh();

		void HandlerInit(std::string, int, char**);
		void HandlerInit(std::string, std::string, int, char**);
		void HandlerInit(std::string, std::vector<std::vector<int>>, std::vector<std::vector<double>>);
		void HandlerInit(std::string);
	private:
		MyMesh mesh;
		MyMesh sim_mesh;
		MyMesh sam_mesh;
		std::string mesh_nm;
	};

}