#include "SPLOC.h"

using namespace std;
using namespace Eigen;

DeformationGraph::DeformationGraph(PGMesh *_mesh, std::vector<PGMesh::VertexHandle> &_vertexHandle, std::vector<int> &_vertexIndex, int seqSize_)
{
	mesh_ = _mesh;
	vertexHandle_ = _vertexHandle;
	vertexIndex_ = _vertexIndex;
	/*if(! mesh_->get_property_handle(rotation, "The rotation matrix of each node"))
	{
	mesh_->add_property(rotation, "The rotation matrix of each node");
	}
	if(! mesh_->get_property_handle(translation, "The translation vector of each node"))
	{
	mesh_->add_property(translation, "The translation vector of each node");
	}
	if(! mesh_->get_property_handle(nodePosition, "The nodePosition of each node"))
	{
	mesh_->add_property(nodePosition, "The nodePosition of each node");
	}*/
	if (!mesh_->get_property_handle(vh, "vh"))
	{
		mesh_->add_property(vh, "vh");
	}
	if (!mesh_->get_property_handle(vnum, "vnum"))
	{
		mesh_->add_property(vnum, "vnum");
	}
	sequenceSize = seqSize_;
	translation.resize(mesh_->n_vertices());
	nodePosition.resize(mesh_->n_vertices());
	rotation.resize(mesh_->n_vertices());
	for (int i = 0; i < mesh_->n_vertices(); i++)
	{
		translation[i].resize(sequenceSize);
		nodePosition[i].resize(sequenceSize);
		rotation[i].resize(sequenceSize);
		mesh_->property(vh, mesh_->vertex_handle(i)) = mesh_->vertex_handle(i);
		mesh_->property(vnum, mesh_->vertex_handle(i)) = i;
	}
}
/*排序函数*/
bool SparseLocalDecomposition::sortByDist(const Node  *v1, const Node  *v2)
{
	return v1->geodist < v2->geodist;//升序排列 
}

SparseLocalDecomposition::SparseLocalDecomposition(std::vector<PGMesh*> &_meshSeq)
{
	meshSeq_ = _meshSeq;
	sequenceSize = meshSeq_.size();
	ref_ = meshSeq_[0];
	/*Simplify the reference mesh*/
	//refMesh = *ref_;
	//decimate = new OpenMesh::Decimater::DecimaterT<PGMesh>(refMesh);

	//typedef OpenMesh::Decimater::ModQuadricT< OpenMesh::Decimater::DecimaterT<PGMesh> >::Handle HModQuadric;			
	//HModQuadric hModQuadric; 			 // use a quadric module			 
	//decimate->add( hModQuadric );   // register module at the decimater

	//decimate->t = geoDistToNodes[vit.handle().idx() * vertexIndexRemaining.size() + i];
	/** */
	//preprocessing2();
	preprocessing2();

}

void SparseLocalDecomposition::preprocessing2()
{
	std::ifstream fpp("template_feature_lnd.txt");

	int l;
	//double a, b, c;
	std::vector<int> init_index(800);
	for (int i = 0; i < 800; i++)
	{
		fpp >> l;
		init_index[i]=(l);
		//fppl >> a;
		//fppl >> b;
		//fppl >> c;
		//std::vector<double> ll;
		//ll.push_back(a);
		//ll.push_back(b);
		//ll.push_back(c);
		//feature_location.push_back(ll);
	}
	OpenMesh::IO::read_mesh(refMesh, "body_graph.obj");
	std::vector<int> vertexIndexRemaining;
	//std::vector<int> vertexIndexRemaining = QEMSimp->simplificating2(init_index);


	std::ifstream fff("body_node_idx.txt", ios::out);
	for (int i = 0; i < refMesh.n_vertices(); i++)
	{
		double l1;
		fff >> l1;
		vertexIndexRemaining.push_back(l1);
	}
	fff.close();
	int lllll = 0;
	for (PGMesh::VertexIter vit = refMesh.vertices_begin(); vit != refMesh.vertices_end(); ++vit)
	{
		PGMesh::Point  &p = (refMesh.point(*vit));
		PGMesh::VertexHandle vit1 = ref_->vertex_handle(lllll);
		PGMesh::Point  &p1 = (refMesh.point(vit1));
		p[0] = p1[0];
		p[1] = p1[1];
		p[2] = p1[2];
		lllll++;
	}
	std::vector<PGMesh::VertexHandle> vertexHandleRemaining(vertexIndexRemaining.size());
	for (int i = 0; i < init_index.size(); i++)
	{

		feature_index.push_back(init_index[i]);

	}


	/*建立deformation graph*/
	defgraph_ = new DeformationGraph(&refMesh, vertexHandleRemaining, vertexIndexRemaining, sequenceSize);
	/*Update the weights of each vertex on the reference mesh*/
	if (!ref_->get_property_handle(knNodes, "The k-nearnest nodes"))
		ref_->add_property(knNodes, "The k-nearnest nodes");
	if (!ref_->get_property_handle(knWeights, "The weights of k-nearnest nodes"))
		ref_->add_property(knWeights, "The weights of k-nearnest nodes");
	/*找出模型个点到各个结点的测地距离*/
	std::vector <double> geoToNode;
	std::vector <int> geoToNode_idx;
	ifstream geoIf("geoDis_body.txt");
	ifstream geoIf_idx("geoDis_body_idx.txt");
	for (int i = 0; i < ref_->n_vertices() * 5; i++)
	{
		double llll;
		geoIf >> llll;
		geoToNode.push_back(llll);
		double  nnnn;
		geoIf_idx >> nnnn;
		geoToNode_idx.push_back(int(nnnn));
	}

	for (PGMesh::VertexIter vit = ref_->vertices_begin(); vit != ref_->vertices_end(); ++vit)
	{
		std::vector<Node> vToNodes(5);
		for (int i = 0; i < 5; i++)
		{
			Node v;
			v.nodeNo = i;
			v.geodist = geoToNode[i  + 5 * vit->idx()];
			v.vh = refMesh.vertex_handle(geoToNode_idx[i + 5 * vit->idx()]);
			vToNodes[i]= (v);
		}
		/*对距离进行排序*/

		//std::sort(vToNodes.begin(), vToNodes.end());
		(ref_->property(knNodes, *vit)).push_back(vToNodes[0].vh);
		(ref_->property(knNodes, *vit)).push_back(vToNodes[1].vh);
		(ref_->property(knNodes, *vit)).push_back(vToNodes[2].vh);
		(ref_->property(knNodes, *vit)).push_back(vToNodes[3].vh);
		ref_->property(knWeights, *vit).push_back((1 - vToNodes[0].geodist / vToNodes[4].geodist) * (1 - vToNodes[0].geodist / vToNodes[4].geodist));
		ref_->property(knWeights, *vit).push_back((1 - vToNodes[1].geodist / vToNodes[4].geodist) * (1 - vToNodes[1].geodist / vToNodes[4].geodist));
		ref_->property(knWeights, *vit).push_back((1 - vToNodes[2].geodist / vToNodes[4].geodist) * (1 - vToNodes[2].geodist / vToNodes[4].geodist));
		ref_->property(knWeights, *vit).push_back((1 - vToNodes[3].geodist / vToNodes[4].geodist) * (1 - vToNodes[3].geodist / vToNodes[4].geodist));
	}

	for (int i = 0; i < sequenceSize; i++)
	{
		if (!meshSeq_[i]->get_property_handle(Yaxis, "the first step pose align axis"))
			meshSeq_[i]->add_property(Yaxis, "the first step pose align axis");
	}
	int startNo = 4000;
	int endNo =4000;

//#pragma omp parallel for
	for (int j = startNo; j <= endNo; j++)
	{
		int id_idx = j / 21;
		int pos_idx = j % 21;
		char *tmp1 = new char[50];
		char *tmp2 = new char[50];
		char *tmp3 = new char[50];
		char *tmp4 = new char[50];
		//sprintf_s(tmp1,50, "body_scan_lm_pos\\%d_%d_pose.txt", id_idx, pos_idx);
		//sprintf_s(tmp2, 50, "body_scan_lm_obj\\%d_%d_pose.obj", id_idx, pos_idx);
		//sprintf_s(tmp3, 50, "body_target_obj\\%d_%d_pose.obj", id_idx, pos_idx);
		//sprintf_s(tmp4, 50, "body_final_result\\%d_%d_pose.obj", id_idx, pos_idx);
		sprintf_s(tmp1, 50, "body_scan_lm_pos\\csr%04da.txt", j);
		sprintf_s(tmp2, 50, "D:\\result\\deformationgraph\\csr%04da.obj", j);
		sprintf_s(tmp3, 50, "body_target_obj\\csr%04da.obj", j);
		sprintf_s(tmp4, 50, "body_final_result\\csr%04da.obj", j);
		string str1 = tmp1;
		string str2 = tmp2;
		string str3 = tmp3;
		string str4 = tmp4;
		std::ifstream fppl(str1);
		if (!fppl)
		{
			continue;
		}
		double a, b, c;
		a = 0;
		b = 0;
		c = 0;
		std::vector<std::vector<double>> feature_location1;
		PGMesh *targetMesh = new PGMesh();
		//OpenMesh::IO::read_mesh(*targetMesh, str3);
		//step0: build the kd-tree of the target mesh
		//Eigen::Matrix<double, Dynamic, Dynamic>  mat(targetMesh->n_vertices(), 3);
		//int xxxx = 0;
		//for (PGMesh::VertexIter vit = targetMesh->vertices_begin(); xxxx < targetMesh->n_vertices(); ++vit)
		//{
		//	PointCloud<double>::Point pp;
		//	PGMesh::Point &pt = (targetMesh)->point(*vit);
		//	mat(xxxx, 0) = pt[0];
		//	mat(xxxx, 1) = pt[1];
		//	mat(xxxx, 2) = pt[2];
		//	xxxx++;
		//}
		//const int max_leaf = 10;
		//my_kd_tree_t mat_index(3, mat, max_leaf);
		//mat_index.index->buildIndex();
		
		//step1: sparse feature point alignment
		feature_location1.resize(800);
		string strTmp;
		std::vector<int> validIndex(0);
		for (int i = 0; i < 800; i++)
		{
			bool h;
			int calculate = 0;
			std::vector<double> ll;
			fppl >> strTmp;
			while (calculate < 3)
			{
				if (JudgeNum(strTmp, a) == 0)
					continue;
				if (JudgeNum(strTmp, a) == 1)
				{
					h = false;
					calculate++;
				}
				if (JudgeNum(strTmp, a) == 2)
				{
					h = true;
					ll.push_back(a);
					calculate++;
				}
				if (calculate < 3)
				{
					fppl >> strTmp;
				}
			}
			if (h)
				validIndex.push_back(i);
			feature_location1[i] = (ll);
		}
		std::vector<Eigen::VectorXd> deformed = deformWithFeaturePoints(init_index, feature_location1, str2, validIndex);
		//step2: Dence correspondence alignment
		//std::vector<int> feature_index1;
		//feature_index1.resize(0);
		//std::vector<std::vector<double>> feature_location2;
		//feature_location2.resize(0);
		//vector<int>::iterator it;
		//for (int i = 0; i < deformed.size(); i++)
		//{
		//	it = find(feature_index.begin(), feature_index.end(), i);
		//	if (it != feature_index.end())
		//	{
		//		feature_index1.push_back(i);
		//		feature_location2.push_back(feature_location1[int(it - feature_index.begin())]);
		//	}
		//	else
		//	{
		//		//找到最近点的坐标
		//		
		//		 double * dis = new double(1);
		//		 size_t * result = new size_t(1);
		//		 double *qry = new double(3);
		//		 qry[0] = deformed[i](0);
		//		 qry[1] = deformed[i](1);
		//		 qry[2] = deformed[i](2);
		//		 mat_index.query(qry, 1, result, dis);
		//		 feature_index1.push_back(i);
		//		 std::vector<double> po(3);
		//		 PGMesh::VertexHandle c = targetMesh->vertex_handle(int(result[0]));
		//		 PGMesh::Point &p = targetMesh->point(c);
		//		 po[0] = p[0];
		//		 po[1] = p[1];
		//		 po[2] = p[2];
		//		 feature_location2.push_back(po);
		//	}
		//}
		//std::vector<Eigen::VectorXd> deformed1 = deformWithFeaturePoints(feature_index1, feature_location2, str4);
	}
}

int SparseLocalDecomposition::JudgeNum(string str, double& iTmp)
{
	int bNum = 2;
	string::size_type szSize = str.size();
	if (szSize == 0)
		bNum = 0;
	for (int i = 0; i<szSize; ++i)
	{
		char ch = str.at(i);
		if (ch == ' ')
		{
			bNum = 0;
			break;
		}
		else if (((ch < '0') || (ch > '9')) && (ch != 'e') && (ch != '+') && (ch != '.') && (ch != '-'))
		{
			
			bNum = 1;
			break;
		}
	}
	if (bNum >1)
	{
		istringstream iss(str);
		iss >> iTmp;
	}
	return bNum;
}

std::vector<Eigen::VectorXd> SparseLocalDecomposition::deformWithFeaturePoints(std::vector<int> feature_index, std::vector<std::vector<double>> feature_location, string path, std::vector<int> v_index)
{
	PGMesh::VertexIter vit = defgraph_->mesh_->vertices_begin();
	for (; vit != defgraph_->mesh_->vertices_end(); ++vit)
	{
		int a = vit->idx();
		//outp<<"a: "<<a<<std::endl;
		int b = defgraph_->vertexIndex_[a];
		//outp<<"b: "<<b<<std::endl;
		PGMesh::VertexHandle c = meshSeq_[0]->vertex_handle(b);

		PGMesh::Point &p = meshSeq_[0]->point(c);
		Eigen::VectorXd position(3);
		position(0) = p[0];
		position(1) = p[1];
		position(2) = p[2];
		(defgraph_->nodePosition[vit->idx()])[0] = position;
		//outp<<"position: "<<std::endl<<position<<std::endl;
		/*Eigen::VectorXd position(3);
		PGMesh::Point &p = defgraph_->mesh_->point(vit.handle());
		position(0) = p[0];
		position(1) = p[1];
		position(2) = p[2];
		(defgraph_->nodePosition[vit.handle().idx()])[seqNo] = position;*/
	}
	int n_feature_Points = v_index.size();
	/*Initialization of the solver*/
	Eigen::SparseMatrix<double> A_rot(6 * defgraph_->mesh_->n_vertices(), 12 * defgraph_->mesh_->n_vertices());
	Eigen::SparseMatrix<double> A_reg(3 * defgraph_->mesh_->n_halfedges(), 12 * defgraph_->mesh_->n_vertices());
	//Eigen::SparseMatrix<double> A_reg(9 * defgraph_->mesh_->n_vertices(),12*defgraph_->mesh_->n_vertices());
	Eigen::SparseMatrix<double> A_con(3 * n_feature_Points, 12 * defgraph_->mesh_->n_vertices());
	Eigen::VectorXd b_reg(12 * defgraph_->mesh_->n_vertices());
	Eigen::VectorXd f_rot(6 * defgraph_->mesh_->n_vertices());
	Eigen::VectorXd f_reg(3 * defgraph_->mesh_->n_halfedges());
	//Eigen::VectorXd f_reg(12 * defgraph_->mesh_->n_vertices());
	Eigen::VectorXd f_con(3 * n_feature_Points);

	int currentIter = 0;
	int totalIter = 5;
	//double alpha_con = 0.003;
	for (int i = 0; i< defgraph_->mesh_->n_vertices(); i++)
	{
		b_reg(12 * i + 0) = 1;
		b_reg(12 * i + 1) = 0;
		b_reg(12 * i + 2) = 0;
		b_reg(12 * i + 3) = 0;
		b_reg(12 * i + 4) = 1;
		b_reg(12 * i + 5) = 0;
		b_reg(12 * i + 6) = 0;
		b_reg(12 * i + 7) = 0;
		b_reg(12 * i + 8) = 1;
		b_reg(12 * i + 9) = 0;
		b_reg(12 * i + 10) = 0;
		b_reg(12 * i + 11) = 0;
	}
	/*Do the iteration*/
	while (1)
	{
		if (currentIter >= totalIter)
			break;
		currentIter++;
		/*Set the three sparse matrices zero*/
		A_rot.setZero();
		A_reg.setZero();
		A_con.setZero();
		/*calculate the value of f_rot*/
		for (int i = 0; i < defgraph_->mesh_->n_vertices(); i++)
		{
			f_rot(i * 6 + 0) = b_reg(i * 12 + 0) * b_reg(i * 12 + 1) + b_reg(i * 12 + 3) * b_reg(i * 12 + 4) + b_reg(i * 12 + 6) * b_reg(i * 12 + 7);//
			f_rot(i * 6 + 1) = b_reg(i * 12 + 0) * b_reg(i * 12 + 2) + b_reg(i * 12 + 3) * b_reg(i * 12 + 5) + b_reg(i * 12 + 6) * b_reg(i * 12 + 8);//
			f_rot(i * 6 + 2) = b_reg(i * 12 + 1) * b_reg(i * 12 + 2) + b_reg(i * 12 + 4) * b_reg(i * 12 + 5) + b_reg(i * 12 + 7) * b_reg(i * 12 + 8);//
			f_rot(i * 6 + 3) = b_reg(i * 12 + 0) * b_reg(i * 12 + 0) + b_reg(i * 12 + 3) * b_reg(i * 12 + 3) + b_reg(i * 12 + 6) * b_reg(i * 12 + 6) - 1;
			f_rot(i * 6 + 4) = b_reg(i * 12 + 1) * b_reg(i * 12 + 1) + b_reg(i * 12 + 4) * b_reg(i * 12 + 4) + b_reg(i * 12 + 7) * b_reg(i * 12 + 7) - 1;
			f_rot(i * 6 + 5) = b_reg(i * 12 + 2) * b_reg(i * 12 + 2) + b_reg(i * 12 + 5) * b_reg(i * 12 + 5) + b_reg(i * 12 + 8) * b_reg(i * 12 + 8) - 1;
		}
		/*calculate the Jacobian matrtix of f_rot*/
		for (int i = 0; i < defgraph_->mesh_->n_vertices(); i++)
		{
			A_rot.insert(i * 6 + 0, i * 12 + 0) = b_reg(i * 12 + 1);//
			A_rot.insert(i * 6 + 0, i * 12 + 1) = b_reg(i * 12 + 0);//
			A_rot.insert(i * 6 + 0, i * 12 + 3) = b_reg(i * 12 + 4);//
			A_rot.insert(i * 6 + 0, i * 12 + 4) = b_reg(i * 12 + 3);//
			A_rot.insert(i * 6 + 0, i * 12 + 6) = b_reg(i * 12 + 7);//
			A_rot.insert(i * 6 + 0, i * 12 + 7) = b_reg(i * 12 + 6);//
			A_rot.insert(i * 6 + 1, i * 12 + 0) = b_reg(i * 12 + 2);//
			A_rot.insert(i * 6 + 1, i * 12 + 2) = b_reg(i * 12 + 0);//
			A_rot.insert(i * 6 + 1, i * 12 + 3) = b_reg(i * 12 + 5);//
			A_rot.insert(i * 6 + 1, i * 12 + 5) = b_reg(i * 12 + 3);//
			A_rot.insert(i * 6 + 1, i * 12 + 6) = b_reg(i * 12 + 8);//
			A_rot.insert(i * 6 + 1, i * 12 + 8) = b_reg(i * 12 + 6);//
			A_rot.insert(i * 6 + 2, i * 12 + 1) = b_reg(i * 12 + 2);//
			A_rot.insert(i * 6 + 2, i * 12 + 2) = b_reg(i * 12 + 1);//
			A_rot.insert(i * 6 + 2, i * 12 + 4) = b_reg(i * 12 + 5);//
			A_rot.insert(i * 6 + 2, i * 12 + 5) = b_reg(i * 12 + 4);//
			A_rot.insert(i * 6 + 2, i * 12 + 7) = b_reg(i * 12 + 8);//
			A_rot.insert(i * 6 + 2, i * 12 + 8) = b_reg(i * 12 + 7); //
			A_rot.insert(i * 6 + 3, i * 12 + 0) = 2 * b_reg(i * 12 + 0); //
			A_rot.insert(i * 6 + 3, i * 12 + 3) = 2 * b_reg(i * 12 + 3); //
			A_rot.insert(i * 6 + 3, i * 12 + 6) = 2 * b_reg(i * 12 + 6); //
			A_rot.insert(i * 6 + 4, i * 12 + 1) = 2 * b_reg(i * 12 + 1); //
			A_rot.insert(i * 6 + 4, i * 12 + 4) = 2 * b_reg(i * 12 + 4);//
			A_rot.insert(i * 6 + 4, i * 12 + 7) = 2 * b_reg(i * 12 + 7);//
			A_rot.insert(i * 6 + 5, i * 12 + 2) = 2 * b_reg(i * 12 + 2);//
			A_rot.insert(i * 6 + 5, i * 12 + 5) = 2 * b_reg(i * 12 + 5);//
			A_rot.insert(i * 6 + 5, i * 12 + 8) = 2 * b_reg(i * 12 + 8);//
		}

		//std::fstream op5("/Users/wangyupan/human model dataset/mesh/f_reg_vkvj.txt", ios::out);

		for (PGMesh::HalfedgeIter hit = defgraph_->mesh_->halfedges_begin(); hit != defgraph_->mesh_->halfedges_end(); ++hit)
		{
			/*calculate the value of f_reg*/
			PGMesh::VertexHandle vk = defgraph_->mesh_->from_vertex_handle(*hit);
			PGMesh::VertexHandle vj = defgraph_->mesh_->to_vertex_handle(*hit);
			Eigen::VectorXd lk = (defgraph_->nodePosition[vk.idx()])[0];
			Eigen::VectorXd lj = (defgraph_->nodePosition[vj.idx()])[0];
			//op5 << "vkvjlklj: " << std::endl << vk << "," << vj << "," << std::endl << lk << std::endl << lj << std::endl;
			f_reg(3 * hit->idx() + 0) = b_reg(12 * vj.idx() + 0) * (lk(0) - lj(0)) + b_reg(12 * vj.idx() + 1) * (lk(1) - lj(1))
				+ b_reg(12 * vj.idx() + 2) * (lk(2) - lj(2)) + lj(0) + b_reg(12 * vj.idx() + 9) - lk(0) - b_reg(12 * vk.idx() + 9);
			f_reg(3 * hit->idx() + 1) = b_reg(12 * vj.idx() + 3) * (lk(0) - lj(0)) + b_reg(12 * vj.idx() + 4) * (lk(1) - lj(1))
				+ b_reg(12 * vj.idx() + 5) * (lk(2) - lj(2)) + lj(1) + b_reg(12 * vj.idx() + 10) - lk(1) - b_reg(12 * vk.idx() + 10);
			f_reg(3 * hit->idx() + 2) = b_reg(12 * vj.idx() + 6) * (lk(0) - lj(0)) + b_reg(12 * vj.idx() + 7) * (lk(1) - lj(1))
				+ b_reg(12 * vj.idx() + 8) * (lk(2) - lj(2)) + lj(2) + b_reg(12 * vj.idx() + 11) - lk(2) - b_reg(12 * vk.idx() + 11);
			/*calculate the Jacobian matrix of f_reg*/
			A_reg.insert(3 * hit->idx() + 0, 12 * vj.idx() + 0) = (lk(0) - lj(0));
			A_reg.insert(3 * hit->idx() + 0, 12 * vj.idx() + 1) = (lk(1) - lj(1));
			A_reg.insert(3 * hit->idx() + 0, 12 * vj.idx() + 2) = (lk(2) - lj(2));
			A_reg.insert(3 * hit->idx() + 0, 12 * vj.idx() + 9) = 1;
			A_reg.insert(3 * hit->idx() + 0, 12 * vk.idx() + 9) = -1;
			A_reg.insert(3 * hit->idx() + 1, 12 * vj.idx() + 3) = (lk(0) - lj(0));
			A_reg.insert(3 * hit->idx() + 1, 12 * vj.idx() + 4) = (lk(1) - lj(1));
			A_reg.insert(3 * hit->idx() + 1, 12 * vj.idx() + 5) = (lk(2) - lj(2));
			A_reg.insert(3 * hit->idx() + 1, 12 * vj.idx() + 10) = 1;
			A_reg.insert(3 * hit->idx() + 1, 12 * vk.idx() + 10) = -1;
			A_reg.insert(3 * hit->idx() + 2, 12 * vj.idx() + 6) = (lk(0) - lj(0));
			A_reg.insert(3 * hit->idx() + 2, 12 * vj.idx() + 7) = (lk(1) - lj(1));
			A_reg.insert(3 * hit->idx() + 2, 12 * vj.idx() + 8) = (lk(2) - lj(2));
			A_reg.insert(3 * hit->idx() + 2, 12 * vj.idx() + 11) = 1;
			A_reg.insert(3 * hit->idx() + 2, 12 * vk.idx() + 11) = -1;
		}

		/*calculate the value of f_con*/
		std::fstream op("E://opp.txt", ios::out);
		for (int i = 0; i < n_feature_Points; i++)

			//for(int j = 0; j <  int (defgraph_->mesh_->n_vertices() / 2); j++)
		{
			//int i = 2 * j;
			Eigen::VectorXd tildeV(3);
			int conIdx = feature_index[v_index[i]];
			op << "conIdx: " << conIdx << std::endl;
			PGMesh::VertexHandle v0 = ref_->property(knNodes, ref_->vertex_handle(conIdx))[0];
			PGMesh::VertexHandle v1 = ref_->property(knNodes, ref_->vertex_handle(conIdx))[1];
			PGMesh::VertexHandle v2 = ref_->property(knNodes, ref_->vertex_handle(conIdx))[2];
			PGMesh::VertexHandle v3 = ref_->property(knNodes, ref_->vertex_handle(conIdx))[3];
			//op << "v0,v1,v2,v3: " << v0 << " " << v1 << " " << v2 << " " << v3 << std::endl;
			//PGMesh::VertexHandle v4 = ref_->property(knNodes, ref_->vertex_handle(conIdx))[4];
			PGMesh::Point pp = meshSeq_[0]->point(meshSeq_[0]->vertex_handle(conIdx));
			PGMesh::Point pobj = PGMesh::Point(feature_location[v_index[i]][0], feature_location[v_index[i]][1], feature_location[v_index[i]][2]);
			Eigen::VectorXd p0 = (defgraph_->nodePosition[v0.idx()])[0];
			Eigen::VectorXd p1 = (defgraph_->nodePosition[v1.idx()])[0];
			Eigen::VectorXd p2 = (defgraph_->nodePosition[v2.idx()])[0];
			Eigen::VectorXd p3 = (defgraph_->nodePosition[v3.idx()])[0];
			op << "pobj: " << pobj[0] << " " << pobj[1] << " " << pobj[2] << " ;" << std::endl;

			//op << "p0,p1,p2,p3: " << p0 << "; " << p1 << " ;" << p2 << "; " << p3 << std::endl;

			//Eigen::VectorXd p4 = (defgraph_->nodePosition[v4.idx()])[seqNo];
			Eigen::VectorXd p(3);
			p << pp[0], pp[1], pp[2];
			Eigen::VectorXd ppp1(3);
			ppp1 << pobj[0], pobj[1], pobj[2];
			Eigen::MatrixXd R0(3, 3);
			Eigen::MatrixXd R1(3, 3);
			Eigen::MatrixXd R2(3, 3);
			Eigen::MatrixXd R3(3, 3);
			//Eigen::MatrixXd R4(3,3);
			Eigen::VectorXd t0(3);
			Eigen::VectorXd t1(3);
			Eigen::VectorXd t2(3);
			Eigen::VectorXd t3(3);
			//Eigen::VectorXd t4(3);
			R0 << b_reg(12 * v0.idx() + 0), b_reg(12 * v0.idx() + 1), b_reg(12 * v0.idx() + 2),
				b_reg(12 * v0.idx() + 3), b_reg(12 * v0.idx() + 4), b_reg(12 * v0.idx() + 5),
				b_reg(12 * v0.idx() + 6), b_reg(12 * v0.idx() + 7), b_reg(12 * v0.idx() + 8);
			R1 << b_reg(12 * v1.idx() + 0), b_reg(12 * v1.idx() + 1), b_reg(12 * v1.idx() + 2),
				b_reg(12 * v1.idx() + 3), b_reg(12 * v1.idx() + 4), b_reg(12 * v1.idx() + 5),
				b_reg(12 * v1.idx() + 6), b_reg(12 * v1.idx() + 7), b_reg(12 * v1.idx() + 8);
			R2 << b_reg(12 * v2.idx() + 0), b_reg(12 * v2.idx() + 1), b_reg(12 * v2.idx() + 2),
				b_reg(12 * v2.idx() + 3), b_reg(12 * v2.idx() + 4), b_reg(12 * v2.idx() + 5),
				b_reg(12 * v2.idx() + 6), b_reg(12 * v2.idx() + 7), b_reg(12 * v2.idx() + 8);
			R3 << b_reg(12 * v3.idx() + 0), b_reg(12 * v3.idx() + 1), b_reg(12 * v3.idx() + 2),
				b_reg(12 * v3.idx() + 3), b_reg(12 * v3.idx() + 4), b_reg(12 * v3.idx() + 5),
				b_reg(12 * v3.idx() + 6), b_reg(12 * v3.idx() + 7), b_reg(12 * v3.idx() + 8);
			//R4 << b_reg(12 * v4.idx() + 0), b_reg(12 * v4.idx() + 1), b_reg(12 * v4.idx() + 2),
			//	b_reg(12 * v4.idx() + 3), b_reg(12 * v4.idx() + 4), b_reg(12 * v4.idx() + 5),
			//	b_reg(12 * v4.idx() + 6), b_reg(12 * v4.idx() + 7), b_reg(12 * v4.idx() + 8);
			t0 << b_reg(12 * v0.idx() + 9), b_reg(12 * v0.idx() + 10), b_reg(12 * v0.idx() + 11);
			t1 << b_reg(12 * v1.idx() + 9), b_reg(12 * v1.idx() + 10), b_reg(12 * v1.idx() + 11);
			t2 << b_reg(12 * v2.idx() + 9), b_reg(12 * v2.idx() + 10), b_reg(12 * v2.idx() + 11);
			t3 << b_reg(12 * v3.idx() + 9), b_reg(12 * v3.idx() + 10), b_reg(12 * v3.idx() + 11);
			//t4 << b_reg(12 * v4.idx() + 9), b_reg(12 * v4.idx() + 10), b_reg(12 * v4.idx() + 11);
			double w0 = ref_->property(knWeights, ref_->vertex_handle(conIdx))[0];
			double w1 = ref_->property(knWeights, ref_->vertex_handle(conIdx))[1];
			double w2 = ref_->property(knWeights, ref_->vertex_handle(conIdx))[2];
			double w3 = ref_->property(knWeights, ref_->vertex_handle(conIdx))[3];
			//double w4 = ref_->property(knWeights, ref_->vertex_handle(conIdx))[4];
			tildeV = (w0 * (R0 * (p - p0) + p0 + t0) + w1 * (R1 * (p - p1) + p1 + t1) + w2 * (R2 * (p - p2) + p2 + t2) + w3 * (R3 * (p - p3) + p3 + t3)) / (w0 + w1 + w2 + w3);
			Eigen::VectorXd dif = tildeV - ppp1;
			f_con(3 * i + 0) = dif(0);
			f_con(3 * i + 1) = dif(1);
			f_con(3 * i + 2) = dif(2);
			/*calculate the Jacobian matrix of f_con*/
			A_con.insert(3 * i + 0, 12 * v0.idx() + 0) = (w0 * (p - p0)(0)) / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 0, 12 * v0.idx() + 1) = (w0 * (p - p0)(1)) / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 0, 12 * v0.idx() + 2) = (w0 * (p - p0)(2)) / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 0, 12 * v0.idx() + 9) = w0 / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 0, 12 * v1.idx() + 0) = (w1 * (p - p0)(0)) / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 0, 12 * v1.idx() + 1) = w1 * (p - p0)(1) / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 0, 12 * v1.idx() + 2) = w1 * (p - p0)(2) / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 0, 12 * v1.idx() + 9) = w1 / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 0, 12 * v2.idx() + 0) = w2 * (p - p0)(0) / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 0, 12 * v2.idx() + 1) = w2 * (p - p0)(1) / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 0, 12 * v2.idx() + 2) = w2 * (p - p0)(2) / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 0, 12 * v2.idx() + 9) = w2 / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 0, 12 * v3.idx() + 0) = w3 * (p - p0)(0) / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 0, 12 * v3.idx() + 1) = w3 * (p - p0)(1) / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 0, 12 * v3.idx() + 2) = w3 * (p - p0)(2) / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 0, 12 * v3.idx() + 9) = w3 / (w0 + w1 + w2 + w3);
			//A_con.insert(3 * i + 0, 12 * v4.idx() + 0) = w4 * (p - p0)(0) / (w0 + w1 + w2 + w3 + w4);
			//A_con.insert(3 * i + 0, 12 * v4.idx() + 1) = w4 * (p - p0)(1) / (w0 + w1 + w2 + w3 + w4);
			//A_con.insert(3 * i + 0, 12 * v4.idx() + 2) = w4 * (p - p0)(2) / (w0 + w1 + w2 + w3 + w4);
			//A_con.insert(3 * i + 0, 12 * v4.idx() + 9) = w4 / (w0 + w1 + w2 + w3 + w4);

			A_con.insert(3 * i + 1, 12 * v0.idx() + 3) = w0 * (p - p0)(0) / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 1, 12 * v0.idx() + 4) = w0 * (p - p0)(1) / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 1, 12 * v0.idx() + 5) = w0 * (p - p0)(2) / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 1, 12 * v0.idx() + 10) = w0 / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 1, 12 * v1.idx() + 3) = w1 * (p - p0)(0) / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 1, 12 * v1.idx() + 4) = w1 * (p - p0)(1) / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 1, 12 * v1.idx() + 5) = w1 * (p - p0)(2) / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 1, 12 * v1.idx() + 10) = w1 / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 1, 12 * v2.idx() + 3) = w2 * (p - p0)(0) / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 1, 12 * v2.idx() + 4) = w2 * (p - p0)(1) / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 1, 12 * v2.idx() + 5) = w2 * (p - p0)(2) / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 1, 12 * v2.idx() + 10) = w2 / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 1, 12 * v3.idx() + 3) = w3 * (p - p0)(0) / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 1, 12 * v3.idx() + 4) = w3 * (p - p0)(1) / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 1, 12 * v3.idx() + 5) = w3 * (p - p0)(2) / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 1, 12 * v3.idx() + 10) = w3 / (w0 + w1 + w2 + w3);
			/*A_con.insert(3 * i + 1, 12 * v4.idx() + 3) = w4 * (p - p0)(0) / (w0 + w1 + w2+ w3 + w4);
			A_con.insert(3 * i + 1, 12 * v4.idx() + 4) = w4 * (p - p0)(1) / (w0 + w1 + w2+ w3 + w4);
			A_con.insert(3 * i + 1, 12 * v4.idx() + 5) = w4 * (p - p0)(2) / (w0 + w1 + w2+ w3 + w4);
			A_con.insert(3 * i + 1, 12 * v4.idx() + 10) = w4 / (w0 + w1 + w2+ w3 + w4);*/

			A_con.insert(3 * i + 2, 12 * v0.idx() + 6) = w0 * (p - p0)(0) / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 2, 12 * v0.idx() + 7) = w0 * (p - p0)(1) / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 2, 12 * v0.idx() + 8) = w0 * (p - p0)(2) / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 2, 12 * v0.idx() + 11) = w0 / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 2, 12 * v1.idx() + 6) = w1 * (p - p0)(0) / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 2, 12 * v1.idx() + 7) = w1 * (p - p0)(1) / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 2, 12 * v1.idx() + 8) = w1 * (p - p0)(2) / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 2, 12 * v1.idx() + 11) = w1 / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 2, 12 * v2.idx() + 6) = w2 * (p - p0)(0) / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 2, 12 * v2.idx() + 7) = w2 * (p - p0)(1) / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 2, 12 * v2.idx() + 8) = w2 * (p - p0)(2) / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 2, 12 * v2.idx() + 11) = w2 / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 2, 12 * v3.idx() + 6) = w3 * (p - p0)(0) / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 2, 12 * v3.idx() + 7) = w3 * (p - p0)(1) / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 2, 12 * v3.idx() + 8) = w3 * (p - p0)(2) / (w0 + w1 + w2 + w3);
			A_con.insert(3 * i + 2, 12 * v3.idx() + 11) = w3 / (w0 + w1 + w2 + w3);
			/*A_con.insert(3 * i + 2, 12 * v4.idx() + 6) = w4 * (p - p0)(0) / (w0 + w1 + w2 + w3 + w4);
			A_con.insert(3 * i + 2, 12 * v4.idx() + 7) = w4 * (p - p0)(1) / (w0 + w1 + w2 + w3 + w4);
			A_con.insert(3 * i + 2, 12 * v4.idx() + 8) = w4 * (p - p0)(2) / (w0 + w1 + w2 + w3 + w4);
			A_con.insert(3 * i + 2, 12 * v4.idx() + 11) = w4 / (w0 + w1 + w2 + w3 + w4);*/
		}

		double w_rot = 5.0;
		double w_reg = 25.0;
		double w_con = 5.0;
		//Eigen::CholmodSimplicialLDLT< Eigen::SparseMatrix<double> >solver;
		//SparseQR < SparseMatrix < double > , AMDOrdering < int > > solver ;
		SimplicialCholesky<SparseMatrix<double>> solver;
		solver.compute(w_rot * A_rot.transpose() * A_rot + w_reg * A_reg.transpose() * A_reg + w_con * A_con.transpose() * A_con);
		//solver.compute(w_rot * A_rot.transpose() * A_rot + w_con * A_con.transpose() * A_con);
		if (solver.info() != Success)
		{
			/// decomposit ion failed
			std::cout << "Decomposition failed" << std::endl;
			break;
		}

		VectorXd x_update;
		VectorXd Axb;


		Axb = -1.0f * (w_rot  * A_rot.transpose() * f_rot + w_reg * A_reg.transpose() * f_reg + w_con * A_con.transpose() * f_con);
		//Axb = -1.0f * (w_rot  * A_rot.transpose() * f_rot + w_con * A_con.transpose() * f_con);

		x_update = solver.solve(Axb);
		b_reg = b_reg + x_update;
		f_error = f_rot.squaredNorm() + f_reg.squaredNorm() + f_con.squaredNorm();

	}
	std::vector<Eigen::MatrixXd> rotation1(defgraph_->mesh_->n_vertices());
	std::vector<Eigen::VectorXd> translation1(defgraph_->mesh_->n_vertices());
	for (int i = 0; i < defgraph_->mesh_->n_vertices(); i++)
	{
		Eigen::MatrixXd rot(3, 3);
		rot << b_reg(12 * i + 0), b_reg(12 * i + 1), b_reg(12 * i + 2),
			b_reg(12 * i + 3), b_reg(12 * i + 4), b_reg(12 * i + 5),
			b_reg(12 * i + 6), b_reg(12 * i + 7), b_reg(12 * i + 8);
		/*rot << 1, 0, 0,
		0, 1, 0,
		0, 0, 1;*/
		rotation1[i] = rot;
		Eigen::VectorXd trans(3);
		trans << b_reg(12 * i + 9), b_reg(12 * i + 10), b_reg(12 * i + 11);
		//trans << 0, 0, 0;
		translation1[i] = trans;
	}

	/*reconstruction of the seqNo-th frame*/
	//std::fstream outf("F:\\reconstruction.txt",ios::out);
	std::vector<Eigen::VectorXd> Yaxis1(meshSeq_[0]->n_vertices());
	for (int i = 0; i< meshSeq_[0]->n_vertices(); i++)
	{

		PGMesh::VertexHandle v0 = ref_->property(knNodes, ref_->vertex_handle(i))[0];
		PGMesh::VertexHandle v1 = ref_->property(knNodes, ref_->vertex_handle(i))[1];
		PGMesh::VertexHandle v2 = ref_->property(knNodes, ref_->vertex_handle(i))[2];
		PGMesh::VertexHandle v3 = ref_->property(knNodes, ref_->vertex_handle(i))[3];
		//outf<<"v0,...,v3"<<v0<<" "<<v1<<" "<<v2<<" "<<v3<<std::endl;
		//PGMesh::VertexHandle v4 = ref_->property(knNodes, ref_->vertex_handle(i))[4];
		PGMesh::Point pp = meshSeq_[0]->point(meshSeq_[0]->vertex_handle(i));
		Eigen::VectorXd p0 = (defgraph_->nodePosition[v0.idx()])[0];
		Eigen::VectorXd p1 = (defgraph_->nodePosition[v1.idx()])[0];
		Eigen::VectorXd p2 = (defgraph_->nodePosition[v2.idx()])[0];
		Eigen::VectorXd p3 = (defgraph_->nodePosition[v3.idx()])[0];


		//Eigen::VectorXd p4 = (defgraph_->nodePosition[v4.idx()])[seqNo];
		Eigen::VectorXd p(3);
		p << pp[0], pp[1], pp[2];

		//outf<<"p,p0,...,p3"<<p0<<std::endl<<p<<std::endl<<p1<<std::endl<<p2<<std::endl<<p3<<std::endl;
		Eigen::MatrixXd R0 = rotation1[v0.idx()];
		Eigen::MatrixXd R1 = rotation1[v1.idx()];
		Eigen::MatrixXd R2 = rotation1[v2.idx()];
		Eigen::MatrixXd R3 = rotation1[v3.idx()];
		//outf<<"R0,...,R3"<<R0<<std::endl<<R1<<std::endl<<R2<<std::endl<<R3<<std::endl;
		//Eigen::MatrixXd R4 = defgraph_->rotation[v4.idx()][seqNo];
		Eigen::VectorXd t0 = translation1[v0.idx()];
		Eigen::VectorXd t1 = translation1[v1.idx()];
		Eigen::VectorXd t2 = translation1[v2.idx()];
		Eigen::VectorXd t3 = translation1[v3.idx()];
		//outf<<"t0,...,t3"<<t0<<std::endl<<t1<<std::endl<<t2<<std::endl<<t3<<std::endl;
		//Eigen::VectorXd t4 = defgraph_->translation[v4.idx()][seqNo];
		double w0 = ref_->property(knWeights, ref_->vertex_handle(i))[0];
		double w1 = ref_->property(knWeights, ref_->vertex_handle(i))[1];
		double w2 = ref_->property(knWeights, ref_->vertex_handle(i))[2];
		double w3 = ref_->property(knWeights, ref_->vertex_handle(i))[3];
		//double w4 = ref_->property(knWeights, ref_->vertex_handle(i))[4];
		Eigen::VectorXd llll(3);
		llll = ((w0 * (R0 * (p - p0) + p0 + t0) + w1 * (R1 * (p - p1) + p1 + t1) + w2 * (R2 * (p - p2) + p2 + t2) + w3 * (R3 * (p - p3) + p3 + t3)) / (w0 + w1 + w2 + w3));
		Yaxis1[i] = llll;
	}
	//outf.close();
	/*output the aligned mesh to a file: seqNo.obj*/
	//char *pPath = new char[16];
	//sprintf(pPath,"deformed%d.obj",0);
	//string offname = pPath;
	ofstream offfile(path);
	//delete []pPath;
	for (int i = 0; i< meshSeq_[0]->n_vertices(); i++)
	{
		Eigen::VectorXd v = Yaxis1[i];
		offfile << "v" << " " << v(0) << " " << v(1) << " " << v(2) << endl;
	}
	for (PGMesh::FaceIter fit = meshSeq_[0]->faces_begin(); fit != meshSeq_[0]->faces_end(); ++fit)
	{
		offfile << "f" << " ";
		for (PGMesh::FVIter fvit = meshSeq_[0]->fv_begin(*fit); fvit != meshSeq_[0]->fv_end(*fit); ++fvit)
		{
			offfile << fvit->idx() + 1 << " ";
		}
		offfile << endl;
	}
	return Yaxis1;
}

void SparseLocalDecomposition::rigidAlignment(Eigen::MatrixXd from, Eigen::MatrixXd to, Eigen::MatrixXd * R, double * s, Eigen::VectorXd *t)
/*
        from ---------------  n * 3 matrix
		to   ---------------  n * 3 matrix
		R    ---------------  3 * 3 rotational matrix
		s    ---------------  double  scale factor
		t    ---------------  3 * 1 vector  translation vector
*/
{
	//compute the transformation from 'from' to 'to'
	Eigen::MatrixXd trans(3, 3);

}