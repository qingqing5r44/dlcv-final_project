#include "edit.h"

Eigen::Matrix3d CalTensor(Eigen::Vector3d w)
{
	Eigen::Matrix3d t = Eigen::Matrix3d::Zero();
	t(0, 1) = -w(2);
	t(1, 0) = w(2);
	t(0, 2) = w(1);
	t(2, 0) = -w(1);
	t(1, 2) = -w(0);
	t(2, 1) = w(0);

	return t;
}

std::vector<Eigen::Matrix3d> CalTensor(Eigen::Matrix3d w)
{
	Eigen::Vector3d v1 = w.block(0, 0, 3, 1);
	Eigen::Vector3d v2 = w.block(0, 1, 3, 1);
	Eigen::Vector3d v3 = w.block(0, 2, 3, 1);

	std::vector<Eigen::Matrix3d> t;
	t.push_back(CalTensor(v1));
	t.push_back(CalTensor(v2));
	t.push_back(CalTensor(v3));

	return t;
}

decimal CalTriangleArea(Eigen::Vector3d p1, Eigen::Vector3d p2, Eigen::Vector3d p3)
{
	decimal t1_side[3];
	t1_side[0] = (p1 - p2).norm();
	t1_side[1] = (p3 - p2).norm();
	t1_side[2] = (p1 - p3).norm();

	decimal t1 = (t1_side[0] + t1_side[1] + t1_side[2]) / 2; //半周长;  
	decimal area = std::sqrt(t1*(t1 - t1_side[0])*(t1 - t1_side[1])*(t1 - t1_side[2]));
	return area;
}

mns::MyMeshHandler::MyMeshHandler(std::string nm, int argc, char** argv)
{
	HandlerInit(nm, argc, argv);
}

mns::MyMeshHandler::MyMeshHandler(std::string nm, std::string parmpath, int argc, char** argv)
{
	HandlerInit(nm, parmpath, argc, argv);
}


mns::MyMeshHandler::MyMeshHandler(std::string nm, std::vector<std::vector<int>> body_feature_idx, std::vector<std::vector<double>> body_feature_weights)
{
	HandlerInit(nm, body_feature_idx, body_feature_weights);
}

void mns::MyMeshHandler::HandlerInit(std::string nm)
{
	this->mesh_nm = nm;
	bool res = OpenMesh::IO::read_mesh(this->mesh, nm);
	if (res == false)
	{
		std::cout << "read error\n";
	}
	this->mesh.request_face_normals();
	this->mesh.request_vertex_normals();
	this->mesh.update_normals();

	integer max_it = 100;
	decimal ep1 = 0.01;
	decimal ep2 = 0.1;
	decimal stp = 0.25;
	bool bw_or_not = 0;
	bool uw_or_not = 0;
	bool sim_or_not = 0;
	bool eg_or_not = 0;
	bool ba_or_not = 0;
	std::string e_str = "";
	integer n_eg = 0;
	integer b_m_num = 0;
	integer b_col = 0;
	decimal contri = 0;
	decimal a_c = 800000;
	decimal a_s = 1;
	decimal a_b = 1.2;
	decimal a_v = 0;
	decimal a_w = 0;

	this->p_DeformHandler = new MyDeformationHandler(this->mesh, this->GetSimMesh(), uw_or_not, sim_or_not, max_it, stp, bw_or_not, eg_or_not, ba_or_not, e_str, n_eg, ep1, ep2, b_m_num, b_col, contri, a_c, a_s, a_b, a_v, a_w);

	integer mod = -1;
	for (integer i = 0; i < nm.size(); i++)
	{
		if (nm[i] == '.')
		{
			if ((nm[i + 1] == 'p' && nm[i + 2] == 'l' && nm[i + 3] == 'y') || (nm[i + 1] == 'P' && nm[i + 2] == 'L' && nm[i + 3] == 'Y'))
			{
				mod = 1;
			}
			if ((nm[i + 1] == 'o' && nm[i + 2] == 'f' && nm[i + 3] == 'f') || (nm[i + 1] == 'O' && nm[i + 2] == 'F' && nm[i + 3] == 'F'))
			{
				mod = 2;
			}
			if ((nm[i + 1] == 'o' && nm[i + 2] == 'b' && nm[i + 3] == 'j') || (nm[i + 1] == 'O' && nm[i + 2] == 'B' && nm[i + 3] == 'J'))
			{
				mod = 3;
			}
		}
	}
	this->p_DeformHandler->SetModelMod(mod);
}

void mns::MyMeshHandler::HandlerInit(std::string nm, std::vector<std::vector<int>> body_feature_idx, std::vector<std::vector<double>> body_feature_weights)
{
	this->mesh_nm = nm;
	bool res = OpenMesh::IO::read_mesh(this->mesh, nm);
	if (res == false)
	{
		std::cout << "read error\n";
	}
	this->mesh.request_face_normals();
	this->mesh.request_vertex_normals();
	this->mesh.update_normals();

	integer max_it = 100;
	decimal ep1 = 0.01;
	decimal ep2 = 0.1;
	decimal stp = 0.5;
	bool bw_or_not = 0;
	bool uw_or_not = 0;
	bool sim_or_not = 0;
	bool eg_or_not = 0;
	bool ba_or_not = 0;
	std::string e_str = "";
	integer n_eg = 0;
	integer b_m_num = 0;
	integer b_col = 0;
	decimal contri = 0;
	decimal a_c;
	a_c = 500000;
	decimal a_s = 1;
	decimal a_b = 1.2;
	decimal a_v = 0;
	decimal a_w = 0;

	this->p_DeformHandler = new MyDeformationHandler(body_feature_idx, body_feature_weights, this->mesh, this->GetSimMesh(), uw_or_not, sim_or_not, max_it, stp, bw_or_not, eg_or_not, ba_or_not, e_str, n_eg, ep1, ep2, b_m_num, b_col, contri, a_c, a_s, a_b, a_v, a_w);

	integer mod = -1;
	for (integer i = 0; i < nm.size(); i++)
	{
		if (nm[i] == '.')
		{
			if ((nm[i + 1] == 'p' && nm[i + 2] == 'l' && nm[i + 3] == 'y') || (nm[i + 1] == 'P' && nm[i + 2] == 'L' && nm[i + 3] == 'Y'))
			{
				mod = 1;
			}
			if ((nm[i + 1] == 'o' && nm[i + 2] == 'f' && nm[i + 3] == 'f') || (nm[i + 1] == 'O' && nm[i + 2] == 'F' && nm[i + 3] == 'F'))
			{
				mod = 2;
			}
			if ((nm[i + 1] == 'o' && nm[i + 2] == 'b' && nm[i + 3] == 'j') || (nm[i + 1] == 'O' && nm[i + 2] == 'B' && nm[i + 3] == 'J'))
			{
				mod = 3;
			}
		}
	}
	this->p_DeformHandler->SetModelMod(mod);
}

void mns::MyMeshHandler::HandlerInit(std::string nm, std::string parmpath, int argc, char** argv)
{
	this->argc = argc;
	this->argv = argv;
	this->mesh_nm = nm;
	bool res = OpenMesh::IO::read_mesh(this->mesh, nm);
	if (res == false)
	{
		std::cout << "read error\n";
	}
	this->mesh.request_face_normals();
	this->mesh.request_vertex_normals();
	this->mesh.update_normals();

	std::ifstream ifile(parmpath, std::ios::in);
	integer max_it;
	ifile >> max_it;
	decimal ep1;
	ifile >> ep1;
	decimal ep2;
	ifile >> ep2;
	decimal stp;
	ifile >> stp;
	bool bw_or_not = 0;
	ifile >> bw_or_not;
	bool uw_or_not = 0;
	ifile >> uw_or_not;
	bool sim_or_not = 0;
	ifile >> sim_or_not;
	bool eg_or_not = 0;
	ifile >> eg_or_not;
	bool ba_or_not = 0;
	ifile >> ba_or_not;
	std::string e_str;
	ifile >> e_str;
	integer n_eg;
	ifile >> n_eg;
	integer b_m_num;
	ifile >> b_m_num;
	integer b_col;
	ifile >> b_col;
	decimal contri;
	ifile >> contri;
	decimal a_c;
	ifile >> a_c;
	decimal a_s;
	ifile >> a_s;
	decimal a_b;
	ifile >> a_b;
	decimal a_v;
	ifile >> a_v;
	decimal a_w;
	ifile >> a_w;
	ifile.close();

	this->p_DeformHandler = new MyDeformationHandler(this->mesh, this->GetSimMesh(), uw_or_not, sim_or_not, max_it, stp, bw_or_not, eg_or_not, ba_or_not, e_str, n_eg, ep1, ep2, b_m_num, b_col, contri, a_c, a_s, a_b, a_v, a_w);

	integer mod = -1;
	for (integer i = 0; i < nm.size(); i++)
	{
		if (nm[i] == '.')
		{
			if ((nm[i + 1] == 'p' && nm[i + 2] == 'l' && nm[i + 3] == 'y') || (nm[i + 1] == 'P' && nm[i + 2] == 'L' && nm[i + 3] == 'Y'))
			{
				mod = 1;
			}
			if ((nm[i + 1] == 'o' && nm[i + 2] == 'f' && nm[i + 3] == 'f') || (nm[i + 1] == 'O' && nm[i + 2] == 'F' && nm[i + 3] == 'F'))
			{
				mod = 2;
			}
			if ((nm[i + 1] == 'o' && nm[i + 2] == 'b' && nm[i + 3] == 'j') || (nm[i + 1] == 'O' && nm[i + 2] == 'B' && nm[i + 3] == 'J'))
			{
				mod = 3;
			}
		}
	}
	this->p_DeformHandler->SetModelMod(mod);
}


void mns::MyMeshHandler::HandlerInit(std::string nm, int argc, char** argv)
{
	this->argc = argc;
	this->argv = argv;
	this->mesh_nm = nm;
	bool res = OpenMesh::IO::read_mesh(this->mesh, nm);
	if (res == false)
	{
		std::cout << "read error\n";
	}
	this->mesh.request_face_normals();
	this->mesh.request_vertex_normals();
	this->mesh.update_normals();

	std::ifstream ifile("parameters.txt", std::ios::in);
	integer max_it;
	ifile >> max_it;
	decimal ep1;
	ifile >> ep1;
	decimal ep2;
	ifile >> ep2;
	decimal stp;
	ifile >> stp;
	bool bw_or_not = 0;
	ifile >> bw_or_not;
	bool uw_or_not = 0;
	ifile >> uw_or_not;        
	bool sim_or_not = 0;
	ifile >> sim_or_not;
	bool eg_or_not = 0;
	ifile >> eg_or_not;
	bool ba_or_not = 0;
	ifile >> ba_or_not;
	std::string e_str;
	ifile >> e_str;
	integer n_eg;
	ifile >> n_eg;
	integer b_m_num;
	ifile >> b_m_num;
	integer b_col;
	ifile >> b_col;
	decimal contri;
	ifile >> contri;
	decimal a_c;
	ifile >> a_c;
	decimal a_s;
	ifile >> a_s;
	decimal a_b;
	ifile >> a_b;
	decimal a_v;
	ifile >> a_v;
	decimal a_w;
	ifile >> a_w;
	ifile.close();

	this->p_DeformHandler = new MyDeformationHandler(this->mesh, this->GetSimMesh(), uw_or_not, sim_or_not, max_it, stp, bw_or_not, eg_or_not, ba_or_not, e_str, n_eg, ep1, ep2, b_m_num, b_col, contri, a_c, a_s, a_b, a_v, a_w);

	integer mod = -1;
	for (integer i = 0; i < nm.size(); i++)
	{
		if (nm[i] == '.')
		{
			if ((nm[i + 1] == 'p' && nm[i + 2] == 'l' && nm[i + 3] == 'y') || (nm[i + 1] == 'P' && nm[i + 2] == 'L' && nm[i + 3] == 'Y'))
			{
				mod = 1;
			}
			if ((nm[i + 1] == 'o' && nm[i + 2] == 'f' && nm[i + 3] == 'f') || (nm[i + 1] == 'O' && nm[i + 2] == 'F' && nm[i + 3] == 'F'))
			{
				mod = 2;
			}
			if ((nm[i + 1] == 'o' && nm[i + 2] == 'b' && nm[i + 3] == 'j') || (nm[i + 1] == 'O' && nm[i + 2] == 'B' && nm[i + 3] == 'J'))
			{
				mod = 3;
			}
		}
	}
	this->p_DeformHandler->SetModelMod(mod);
}

mns::MyMeshHandler::~MyMeshHandler()
{
	//delete p_DeformHandler;
	//delete p_ViewerHandler;
	std::cout << "destruct\n";
}

bool mns::MyMeshHandler::SaveMesh()
{
	bool res = OpenMesh::IO::write_mesh(this->mesh, "res.ply");
	return res;
}

void mns::MyMeshHandler::ShowVertex() const
{
	for (auto it = this->mesh.vertices_begin(); it != this->mesh.vertices_end(); ++it)
	{
		auto point = this->mesh.point(*it);
		std::cout << "x:" << point.data()[0] << " y:" << point.data()[1] << " z:" << point.data()[2] << std::endl;
	}
}

void mns::MyMeshHandler::ShowOneVertex(const int id) const
{
	auto it = this->mesh.vertices_begin();
	for (int k = 0; k < id; k++)
	{
		it++;
	}

	auto point = this->mesh.point(*it);
	std::cout << "x:" << point.data()[0] << " y:" << point.data()[1] << " z:" << point.data()[2] << std::endl;

}

void mns::MyMeshHandler::ShowWedge()
{
	auto it = this->mesh.halfedges_begin();
	auto f1 = this->mesh.face_handle(*it);
	auto op_he = this->mesh.opposite_halfedge_handle(*it);
	auto f2 = this->mesh.face_handle(op_he);
	for (auto fv_it1 = this->mesh.fv_begin(f1); fv_it1 != this->mesh.fv_end(f1); ++fv_it1)
	{
		auto v = *fv_it1;
		std::cout << "f1: " << v.idx() << " ";
	}
	std::cout << std::endl;
	for (auto fv_it2 = this->mesh.fv_begin(f2); fv_it2 != this->mesh.fv_end(f2); ++fv_it2)
	{
		auto v = *fv_it2;
		std::cout << "f2: " << v.idx() << " ";
	}
}

void mns::MyMeshHandler::ShowEdge() const
{
	unsigned int count = 0;
	for (auto it = this->mesh.edges_begin(); it != this->mesh.edges_end(); ++it)
	{
		auto edge = this->mesh.edge_handle(count++);
		std::cout << edge.idx();
	}
}

void mns::MyMeshHandler::ShowHalfEdge() const
{
	for (auto it = this->mesh.halfedges_begin(); it != this->mesh.halfedges_end(); ++it)
	{
		auto fromVertex = this->mesh.from_vertex_handle(*it);
		auto toVertex = this->mesh.to_vertex_handle(*it);
		std::cout << "half edge:" << it->idx() << " from:" << fromVertex.idx() << " to:" << toVertex.idx() << std::endl;
	}
}

void mns::MyMeshHandler::ShowOneHalfEdge(const int id) const
{
	auto it = this->mesh.halfedges_begin();
	for (int k = 0; k < id; k++)
	{
		it++;
	}
	auto fromVertex = this->mesh.from_vertex_handle(*it);
	auto toVertex = this->mesh.to_vertex_handle(*it);
	std::cout << "half edge:" << it->idx() << " from:" << fromVertex.idx() << " to:" << toVertex.idx() << std::endl;
}

void mns::MyMeshHandler::ShowVertexNormal() const
{
	for (auto it = this->mesh.vertices_begin(); it != this->mesh.vertices_end(); ++it)
	{
		auto vertex = *it;
		MyMesh::Normal normal = mesh.normal(vertex);
		double x = normal.data()[0];
		double y = normal.data()[1];
		double z = normal.data()[2];
		std::cout << "Vertex:" << it->idx() << " x:" << x << " y:" << y << " z:" << z << std::endl;
	}
}

void mns::MyMeshHandler::ShowFaceNormal() const
{
	for (auto it = this->mesh.faces_begin(); it != this->mesh.faces_end(); ++it)
	{
		auto face = *it;
		MyMesh::Normal normal = mesh.normal(face);
		decimal x = normal.data()[0];
		decimal y = normal.data()[1];
		decimal z = normal.data()[2];
		std::cout << "face:" << it->idx() << " x:" << x << " y:" << y << " z:" << z << std::endl;
	}
}

void mns::MyMeshHandler::ShowFaceHalfEdge()
{
	for (auto it = this->mesh.faces_begin(); it != this->mesh.faces_end(); ++it)
	{
		auto face = *it;
		MyMesh::FaceHalfedgeIter fh_it = mesh.fh_begin(face);

		for (; fh_it.is_valid(); ++fh_it) {
			std::cout << "Halfedge has handle " << *fh_it << std::endl;
		}
	}
}

void mns::MyMeshHandler::ShowOneRingVertex(const int id)
{
	auto it = this->mesh.vertices_begin();
	for (int k = 0; k < id; k++)
	{
		it++;
	}

	auto vertex = *it;
	for (auto v_it = this->mesh.vv_begin(vertex); v_it != this->mesh.vv_end(vertex); v_it++)
	{
		auto n_vertex = *v_it;
		std::cout << "1-ring vertex: " << n_vertex.idx() << std::endl;
	}
}

void mns::MyMeshHandler::OutputInfo()
{
	std::ofstream info_file("mesh_info.txt");
	info_file << "name: " << this->mesh_nm << std::endl;
	info_file << "#vertices: " << this->mesh.n_vertices() << std::endl;
	info_file << "#edges: " << this->mesh.n_edges() << std::endl;
	info_file << "#faces: " << this->mesh.n_faces() << std::endl;
	info_file.close();
}

void mns::MyMeshHandler::SetCollapse(const int from_id, const int to_id)
{
	auto f_it = this->mesh.vertices_begin();
	for (int k = 0; k < from_id; k++)
	{
		f_it++;
	}
	auto from_v = *f_it;

	auto t_it = this->mesh.vertices_begin();
	for (int k = 0; k < to_id; k++)
	{
		t_it++;
	}
	auto to_v = *t_it;

	if (!this->mesh.has_vertex_status())
		this->mesh.request_vertex_status();
	if (!this->mesh.has_face_status())
		this->mesh.request_face_status();
	if (!this->mesh.has_edge_status())
		this->mesh.request_edge_status();

	for (auto it = mesh.halfedges_begin(); it != mesh.halfedges_end(); ++it)
	{
		if (this->mesh.OpenMesh::PolyConnectivity::to_vertex_handle(*it) == to_v &&
			this->mesh.OpenMesh::PolyConnectivity::from_vertex_handle(*it) == from_v)
		{

			if (!(this->mesh.OpenMesh::PolyConnectivity::is_collapse_ok(*it)))
			{
				std::cout << "cannot collapse";
				return;
			}
			// Collapse edge
			mesh.collapse(*it);
			break;
		}
	}
	if (this->mesh.has_vertex_status())
		this->mesh.release_vertex_status();/**/
	if (this->mesh.has_face_status())
		this->mesh.release_face_status();
	if (this->mesh.has_edge_status())
		this->mesh.release_edge_status();
}

void mns::MyMeshHandler::RemoveOutliers()
{
	if (!this->mesh.has_vertex_status())
		this->mesh.request_vertex_status();
	if (!this->mesh.has_face_status())
		this->mesh.request_face_status();
	if (!this->mesh.has_edge_status())
		this->mesh.request_edge_status();
	if (!this->mesh.has_halfedge_status())
		this->mesh.has_halfedge_status();
	for (auto it = this->mesh.faces_begin(); it != this->mesh.faces_end(); ++it)
	{
		auto facet = this->mesh.face(*it);
		auto face = *it;
		for (auto fh_it = this->mesh.fh_begin(face); fh_it != this->mesh.fh_end(face); ++fh_it)
		{
			std::cout << fh_it->idx();
			auto he = this->mesh.halfedge(*fh_it);
		}
	}
	for (auto it = this->mesh.vertices_begin(); it != this->mesh.vertices_end(); ++it)
	{
		auto vertex = *it;
		int count = 0;
		for (auto v_it = this->mesh.vv_begin(vertex); v_it != this->mesh.vv_end(vertex); v_it++)
		{
			count++;
		}
		if (count == 0)
		{
			this->mesh.delete_vertex(vertex, false);
		}
	}

	this->mesh.garbage_collection();

	if (this->mesh.has_vertex_status())
		this->mesh.release_vertex_status();/**/
	if (this->mesh.has_face_status())
		this->mesh.release_face_status();
	if (this->mesh.has_edge_status())
		this->mesh.release_edge_status();
	if (this->mesh.has_halfedge_status())
		this->mesh.release_halfedge_status();
}

void mns::MyMeshHandler::Bilateral(decimal sigma_s, decimal sigma_c, integer iter)
{
	decimal s_coeff = -0.5 / (sigma_s*sigma_s);
	decimal c_coeff = -0.5 / (sigma_c*sigma_c);
	for (int it = 0; it < iter; it++)
	{
		for (auto v_it = this->mesh.vertices_begin(); v_it != this->mesh.vertices_end(); ++v_it)
		{
			auto vertex = *v_it;
			MyMesh::Point v = this->mesh.point(vertex);
			MyMesh::Normal n = this->mesh.normal(vertex);

			decimal sum = 0;
			decimal weight_sum = 0;

			for (auto vv_it = this->mesh.vv_begin(vertex); vv_it != this->mesh.vv_end(vertex); ++vv_it)
			{
				auto neighbor_vertex = *vv_it;
				MyMesh::Point q = this->mesh.point(neighbor_vertex);

				MyMesh::Point diff = v - q;
				double sq_t = diff.length()*diff.length();
				decimal wc = std::exp(sq_t*s_coeff);

				decimal h = n | diff;
				decimal sq_h = h*h;
				decimal ws = std::exp(sq_h*c_coeff);

				decimal weight = wc*ws;
				sum += weight*h;
				weight_sum += weight;
			}
			decimal movement = sum / weight_sum;
			MyMesh::Point new_v = v + n*movement;
			mesh.set_point(vertex, new_v);
		}
		this->mesh.update_normals();
	}
}
mns::MyPoint mns::MyMeshHandler::GetPoint(const int id)
{
	mns::MyPoint temp;

	auto it = this->mesh.vertices_begin();
	for (int k = 0; k < id; k++)
	{
		it++;
	}
	auto vertex = *it;
	auto pt = this->mesh.point(*it);

	temp.idx = id;
	temp.x = pt.data()[0];
	temp.y = pt.data()[1];
	temp.z = pt.data()[2];

	MyMesh::Normal normal = mesh.normal(vertex);
	temp.nx = normal.data()[0];
	temp.ny = normal.data()[1];
	temp.nz = normal.data()[2];

	for (auto v_it = this->mesh.vv_begin(vertex); v_it != this->mesh.vv_end(vertex); v_it++)
	{
		auto n_vertex = *v_it;
		temp.neighbor.push_back((integer)n_vertex.idx());
	}
	return temp;
}

mns::MyPointSet mns::MyMeshHandler::GetPS()
{
	mns::MyPointSet tempPS;

	for (auto it = this->mesh.vertices_begin(); it != this->mesh.vertices_end(); ++it)
	{
		mns::MyPoint temp;
		auto vertex = *it;
		auto pt = this->mesh.point(*it);

		temp.idx = vertex.idx();
		temp.x = pt.data()[0];
		temp.y = pt.data()[1];
		temp.z = pt.data()[2];

		MyMesh::Normal normal = mesh.normal(vertex);
		temp.nx = normal.data()[0];
		temp.ny = normal.data()[1];
		temp.nz = normal.data()[2];

		for (auto v_it = this->mesh.vv_begin(vertex); v_it != this->mesh.vv_end(vertex); v_it++)
		{
			auto n_vertex = *v_it;
			temp.neighbor.push_back((integer)n_vertex.idx());
		}

		tempPS.set.push_back(temp);
	}
	return tempPS;
}

mns::MyEdge mns::MyMeshHandler::GetEdge(const int id)
{
	mns::MyEdge temp;

	auto it = this->mesh.halfedges_begin();
	for (int k = 0; k < 2 * id; ++k)
	{
		it++;
	}
	auto fromVertex = this->mesh.from_vertex_handle(*it);
	auto toVertex = this->mesh.to_vertex_handle(*it);

	temp.length = this->mesh.calc_edge_length(*it);

	if (fromVertex.idx() < toVertex.idx())
	{
		temp.points.push_back((integer)fromVertex.idx());
		temp.points.push_back((integer)toVertex.idx());
	}
	else
	{
		temp.points.push_back((integer)toVertex.idx());
		temp.points.push_back((integer)fromVertex.idx());
	}

	return temp;
}

mns::MyEdgeSet mns::MyMeshHandler::GetES()
{
	mns::MyEdgeSet tempES;
	for (auto it = this->mesh.halfedges_begin(); it != this->mesh.halfedges_end(); )
	{
		MyEdge temp;

		auto fromVertex = this->mesh.from_vertex_handle(*it);
		auto toVertex = this->mesh.to_vertex_handle(*it);

		temp.length = this->mesh.calc_edge_length(*it);

		if (fromVertex.idx() < toVertex.idx())
		{
			temp.points.push_back((integer)fromVertex.idx());
			temp.points.push_back((integer)toVertex.idx());
		}
		else
		{
			temp.points.push_back((integer)toVertex.idx());
			temp.points.push_back((integer)fromVertex.idx());
		}

		tempES.set.push_back(temp);
		++it;
		++it;
	}
	return tempES;
}

mns::MyHalfEdge mns::MyMeshHandler::GetHalfEdge(const int id)
{
	mns::MyHalfEdge temp;

	auto it = this->mesh.halfedges_begin();
	for (int k = 0; k < id; ++k)
	{
		it++;
	}
	auto fromVertex = this->mesh.from_vertex_handle(*it);
	auto toVertex = this->mesh.to_vertex_handle(*it);

	temp.length = this->mesh.calc_edge_length(*it);

	temp.from_vertex = (integer)fromVertex.idx();
	temp.to_vertex = (integer)toVertex.idx();

	return temp;
}

mns::MyHalfEdgeSet mns::MyMeshHandler::GetHES()
{
	mns::MyHalfEdgeSet tempHES;
	for (auto it = this->mesh.halfedges_begin(); it != this->mesh.halfedges_end(); ++it)
	{
		MyHalfEdge temp;

		auto fromVertex = this->mesh.from_vertex_handle(*it);
		auto toVertex = this->mesh.to_vertex_handle(*it);

		temp.length = this->mesh.calc_edge_length(*it);

		temp.from_vertex = (integer)fromVertex.idx();
		temp.to_vertex = (integer)toVertex.idx();

		tempHES.set.push_back(temp);
	}
	return tempHES;
}

mns::MyFace mns::MyMeshHandler::GetFace(const int id)
{
	mns::MyFace temp;
	auto it = this->mesh.faces_begin();
	for (int k = 0; k < id; k++)
	{
		it++;
	}
	auto face = *it;

	int_vector idx_vec;
	pt_vector pts_vec;
	for (auto fv_it = this->mesh.fv_begin(face); fv_it != this->mesh.fv_end(face); ++fv_it)
	{
		auto vertex = *fv_it;
		auto pt = this->mesh.point(vertex);
		idx_vec.push_back((integer)(vertex.idx()));
		pts_vec.push_back(pt);
	}

	std::sort(idx_vec.begin(), idx_vec.end());
	for (int i = 0; i < idx_vec.size(); i++)
	{
		temp.points.push_back(idx_vec[i]);
	}

	MyMesh::Normal normal = this->mesh.normal(face);
	temp.nx = normal.data()[0];
	temp.ny = normal.data()[1];
	temp.nz = normal.data()[2];

	decimal side[3];
	side[0] = (pts_vec[0] - pts_vec[1]).length();
	side[1] = (pts_vec[2] - pts_vec[1]).length();
	side[2] = (pts_vec[0] - pts_vec[2]).length();

	decimal p = (side[0] + side[1] + side[2]) / 2; //半周长;  
	temp.area = std::sqrt(p*(p - side[0])*(p - side[1])*(p - side[2]));

	return temp;
}

mns::MyFaceSet mns::MyMeshHandler::GetFS()
{
	mns::MyFaceSet tempFS;
	for (auto it = this->mesh.faces_begin(); it != this->mesh.faces_end(); ++it)
	{
		mns::MyFace temp;
		auto face = *it;

		int_vector idx_vec;
		pt_vector pts_vec;
		for (auto fv_it = this->mesh.fv_begin(face); fv_it != this->mesh.fv_end(face); ++fv_it)
		{
			auto vertex = *fv_it;
			auto pt = this->mesh.point(vertex);
			idx_vec.push_back((integer)(vertex.idx()));
			pts_vec.push_back(pt);
		}

		std::sort(idx_vec.begin(), idx_vec.end());
		for (int i = 0; i < idx_vec.size(); i++)
		{
			temp.points.push_back(idx_vec[i]);
		}

		MyMesh::Normal normal = this->mesh.normal(face);
		temp.nx = normal.data()[0];
		temp.ny = normal.data()[1];
		temp.nz = normal.data()[2];

		decimal side[3];
		side[0] = (pts_vec[0] - pts_vec[1]).length();
		side[1] = (pts_vec[2] - pts_vec[1]).length();
		side[2] = (pts_vec[0] - pts_vec[2]).length();

		decimal p = (side[0] + side[1] + side[2]) / 2; //半周长;  
		temp.area = std::sqrt(p*(p - side[0])*(p - side[1])*(p - side[2]));

		tempFS.set.push_back(temp);
	}

	return tempFS;
}

MyMesh& mns::MyMeshHandler::GetSimMesh()
{
	std::ifstream ifile("nm.txt", std::ios::in);
	std::string nm;
	ifile >> nm;
	std::string test_nm;
	ifile >> test_nm;
	std::string sim_nm;
	ifile >> sim_nm;
	ifile.close();

	bool res = OpenMesh::IO::read_mesh(this->sim_mesh, sim_nm);
	if (res == false)
	{
		std::cout << "read error\n";
	}
	this->sim_mesh.request_face_normals();
	this->sim_mesh.request_vertex_normals();
	this->sim_mesh.update_normals();

	return this->sim_mesh;
}

MyMesh& mns::MyMeshHandler::GetMesh()
{
	return mesh;
}

void mns::MyMeshHandler::Deformation(std::vector<integer>& s_vid_vec, Eigen::MatrixXd& v_mat, Eigen::MatrixXd& con_v_mat)
{
	this->p_DeformHandler->current_str = "conset0.txt";
	this->p_DeformHandler->DeformInit(s_vid_vec, v_mat, con_v_mat);
	//this->p_DeformHandler->Deform(s_vid_vec, v_mat, con_v_mat);
	/*this->p_DeformHandler->current_str = "conset1.txt";
	this->p_DeformHandler->DeformInit1(s_vid_vec, v_mat, con_v_mat);*/
	//this->p_DeformHandler->Deform(s_vid_vec, v_mat, con_v_mat);
	this->p_DeformHandler->current_str = "conset.txt";
	this->p_DeformHandler->DeformInit1(s_vid_vec, v_mat, con_v_mat);
	this->p_DeformHandler->Deform(s_vid_vec, v_mat, con_v_mat);
	/*this->p_DeformHandler->current_str = "conset1.txt";
	this->p_DeformHandler->DeformInit1(s_vid_vec, v_mat, con_v_mat);
	this->p_DeformHandler->Deform(s_vid_vec, v_mat, con_v_mat);
	this->p_DeformHandler->current_str = "conset2.txt";
	this->p_DeformHandler->DeformInit1(s_vid_vec, v_mat, con_v_mat);
	this->p_DeformHandler->Deform(s_vid_vec, v_mat, con_v_mat);*/
	this->p_DeformHandler->ModOriginalMesh(this->mesh);
	this->p_DeformHandler->ModMesh();
	
	this->p_DeformHandler->ModVMat(v_mat);
	this->p_DeformHandler->ClearAfterDeform();
}

void mns::MyMeshHandler::DeformationLD(std::vector<integer>& s_vid_vec, Eigen::MatrixXd& v_mat, Eigen::MatrixXd& con_v_mat)
{
	this->p_DeformHandler->DeformLD(s_vid_vec, v_mat, con_v_mat);
	//this->p_DeformHandler->ModOriginalMesh(this->mesh);
	this->p_DeformHandler->ModMesh();
	this->p_DeformHandler->ModVMat(v_mat);
	this->p_DeformHandler->ClearAfterDeform();
}

mns::MyDeformationHandler::MyDeformationHandler(std::vector<std::vector<int>> _body_feature_idx, std::vector<std::vector<double>> _body_feature_weights, MyMesh& mesh, MyMesh& s_mesh, const bool use_ww, const bool use_sim, const integer max_it, const decimal stp, const bool bw_or_not, const bool eg_or_not, const bool ba_or_not, const std::string e_str, const integer n_eg, const decimal ep1, const decimal ep2, const integer b_m_num, const integer b_col, const decimal contri, const decimal a_c, const decimal a_s, const decimal a_b, const decimal a_v, const decimal a_w)
	:epsilon1(ep1), epsilon2(ep2), step(stp), use_eg_or_not(eg_or_not), use_basis_or_not(ba_or_not), example_str(e_str), eg_num(n_eg), max_iter(max_it), create_basis_mesh_num(b_m_num), basis_col_num(b_col), contribution(contri), alpha_c(a_c), alpha_s(a_s), alpha_b(a_b), alpha_v(a_v), alpha_w(a_w), mesh(mesh), ori_mesh(mesh), sim_mesh(s_mesh), use_sim_or_not(use_sim), use_ww_or_not(use_ww), cal_bw_or_not(bw_or_not)
{
	/*this->con_ele_cnt = 0;
	this->str_ele_cnt = 0;
	this->ben_ele_cnt = 0;
	this->row_cnt = 0;
	this->pre_energy = 1000000;
	this->iter = 0;*/
	this->body_feature_idx = _body_feature_idx;
	this->body_feature_weights = _body_feature_weights;
	this->fixed_feature_idx.resize(0);
	this->fixed_feature_pos.resize(0);
	this->fixed_feature_weights.resize(0);
	this->rb_flag = false;
	this->alpha_c_ch = false;
	this->alpha_s_ch = false;
	this->alpha_b_ch = false;
	this->alpha_v_ch = false;
	this->alpha_w_ch = false;
	this->ch_con = 0;
	this->out_flag = false;
	this->ori_step = step;
	std::ifstream a_ifile("alpha_change.txt", std::ios::in);
	for (integer i = 0; i < NEED_LEN; i++)
	{
		bool tmp = false;
		a_ifile >> tmp;
		need_ch[i] = tmp;
	}
	a_ifile >> need_fac;
	a_ifile.close();
}


mns::MyDeformationHandler::MyDeformationHandler(MyMesh& mesh, MyMesh& s_mesh, const bool use_ww, const bool use_sim, const integer max_it, const decimal stp, const bool bw_or_not, const bool eg_or_not, const bool ba_or_not, const std::string e_str, const integer n_eg, const decimal ep1, const decimal ep2, const integer b_m_num, const integer b_col, const decimal contri, const decimal a_c, const decimal a_s, const decimal a_b, const decimal a_v, const decimal a_w)
	:epsilon1(ep1), epsilon2(ep2), step(stp), use_eg_or_not(eg_or_not), use_basis_or_not(ba_or_not), example_str(e_str), eg_num(n_eg), max_iter(max_it), create_basis_mesh_num(b_m_num), basis_col_num(b_col), contribution(contri), alpha_c(a_c), alpha_s(a_s), alpha_b(a_b), alpha_v(a_v), alpha_w(a_w), mesh(mesh), ori_mesh(mesh), sim_mesh(s_mesh), use_sim_or_not(use_sim), use_ww_or_not(use_ww), cal_bw_or_not(bw_or_not)
{
	/*this->con_ele_cnt = 0;
	this->str_ele_cnt = 0;
	this->ben_ele_cnt = 0;
	this->row_cnt = 0;
	this->pre_energy = 1000000;
	this->iter = 0;*/
	this->rb_flag = false;
	this->alpha_c_ch = false;
	this->alpha_s_ch = false;
	this->alpha_b_ch = false;
	this->alpha_v_ch = false;
	this->alpha_w_ch = false;
	this->ch_con = 0;
	this->out_flag = false;
	this->ori_step = step;
	std::ifstream a_ifile("alpha_change.txt", std::ios::in);
	for (integer i = 0; i < NEED_LEN; i++)
	{
		bool tmp = false;
		a_ifile >> tmp;
		need_ch[i] = tmp;
	}
	a_ifile >> need_fac;
	a_ifile.close();
}

mns::MyDeformationHandler::~MyDeformationHandler()
{
	ori_length_set.clear();
	ori_area_set.clear();
	ori_dihedral_set.clear();
	wnt_dihedral_set.clear();
	wnt_length_set.clear();
	idx_move.clear();
	con_ori_map.clear();
	/*for (integer k = 0; k < this->mesh.n_edges(); k++)
	{
	delete[]edge_geodis_mat[k];
	delete[]normalized_edge_geodis_mat[k];
	}
	delete[]edge_geodis_mat;
	delete[]normalized_edge_geodis_mat;*/
}

void mns::MyDeformationHandler::DeformInit1(std::vector<integer>&, Eigen::MatrixXd&, Eigen::MatrixXd&)
{
	this->pre_energy = 100000000000000;
	//constraint initialization
	//CreateConstraintPointSet(GetConSetByV(s_vid_vec, con_v_mat));
	triple.clear();
	con_ele_cnt = 0;
	sb_ele_cnt = 0;
}

void mns::MyDeformationHandler::DeformInit5(std::vector<integer>& s_vid_vec, Eigen::MatrixXd& v_mat, Eigen::MatrixXd& con_v_mat)
{
	this->con_ele_cnt = 0;
	this->str_ele_cnt = 0;
	this->ben_ele_cnt = 0;
	this->sb_ele_cnt = 0;
	this->row_cnt = 0;
	this->pre_energy = 100000000000000;
	this->iter = 0;

	//std::ofstream mydi_ofile("mydihedral.txt");
	//std::ofstream wdi_ofile("wdihedral.txt");

	for (auto he_it = this->mesh.halfedges_begin(); he_it != this->mesh.halfedges_end(); ++he_it)
	{
		auto e = *he_it;
		auto fromVertex = this->mesh.from_vertex_handle(e);
		auto toVertex = this->mesh.to_vertex_handle(e);
		auto fromPt = this->mesh.point(fromVertex);
		auto toPt = this->mesh.point(toVertex);
		auto edgee = this->mesh.edge_handle(e);
		int e_idx = edgee.idx();
		//stretching initialization		
		decimal l = this->mesh.calc_edge_length(e);
		this->ori_length_set.push_back(l);
		this->wnt_length_set.push_back(l * (1 + cm1(2 * e_idx + 1)));

		//bend initialization
		auto f1 = this->mesh.face_handle(e);
		auto op_he = this->mesh.opposite_halfedge_handle(e);
		auto f2 = this->mesh.face_handle(op_he);
		pt_vector pt_set;
		pt_vector t1_set;
		pt_vector t2_set;

		if (f1.idx() == -1 || f2.idx() == -1)
		{
			this->ori_dihedral_set.push_back(0);
			this->wnt_dihedral_set.push_back(0);
			ori_area_set.push_back(0);

			this->ori_wlz_dihedral_set.push_back(0);
			this->wnt_wlz_dihedral_set.push_back(0);
			sign_flag.push_back(1);

			//edge, not half edge
			++he_it;
			continue;
		}

		pt_set.push_back(fromPt);
		t1_set.push_back(fromPt);
		t2_set.push_back(fromPt);
		pt_set.push_back(toPt);
		t1_set.push_back(toPt);
		t2_set.push_back(toPt);

		for (auto fv_it1 = this->mesh.fv_begin(f1); fv_it1 != this->mesh.fv_end(f1); ++fv_it1)
		{
			auto v = *fv_it1;
			auto pt = this->mesh.point(v);
			if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
			{
				pt_set.push_back(pt);
				t1_set.push_back(pt);
				break;
			}
		}
		integer temp_cnt = 0;
		for (auto fv_it2 = this->mesh.fv_begin(f2); fv_it2 != this->mesh.fv_end(f2); ++fv_it2)
		{
			auto v = *fv_it2;
			auto pt = this->mesh.point(v);
			if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
			{
				t2_set.push_back(pt);
				pt_set.push_back(pt);
				break;
			}
		}
		decimal t1_side[3];
		t1_side[0] = (t1_set[0] - t1_set[1]).length();
		t1_side[1] = (t1_set[2] - t1_set[1]).length();
		t1_side[2] = (t1_set[0] - t1_set[2]).length();

		decimal p1 = (t1_side[0] + t1_side[1] + t1_side[2]) / 2; //半周长;  
		decimal area1 = std::sqrt(p1*(p1 - t1_side[0])*(p1 - t1_side[1])*(p1 - t1_side[2]));

		decimal t2_side[3];
		t2_side[0] = (t2_set[0] - t2_set[1]).length();
		t2_side[1] = (t2_set[2] - t2_set[1]).length();
		t2_side[2] = (t2_set[0] - t2_set[2]).length();

		decimal p2 = (t2_side[0] + t2_side[1] + t2_side[2]) / 2; //半周长;  
		decimal area2 = std::sqrt(p2*(p2 - t2_side[0])*(p2 - t2_side[1])*(p2 - t2_side[2]));
		ori_area_set.push_back(area1 + area2);
		Eigen::Vector3d pt0(pt_set[0][0], pt_set[0][1], pt_set[0][2]);
		Eigen::Vector3d pt1(pt_set[1][0], pt_set[1][1], pt_set[1][2]);
		Eigen::Vector3d pt2(pt_set[2][0], pt_set[2][1], pt_set[2][2]);
		Eigen::Vector3d pt3(pt_set[3][0], pt_set[3][1], pt_set[3][2]);
		//decimal angle = 3.1415926 - CalWLZHDiheral(pt0, pt1, pt2, pt3);
		//decimal angle = CalDiheral(pt_set[0], pt_set[1], pt_set[2], pt_set[3]);
		/*this->ori_dihedral_set.push_back(angle);
		this->wnt_dihedral_set.push_back(angle);*/

		//mydi_ofile << angle << std::endl;

		//sign_flag computation
		Eigen::Vector3d normal1;
		Eigen::Vector3d normal2;
		{
			auto fhe_it1 = this->mesh.fh_begin(f1);
			auto he1 = *fhe_it1;
			auto he1idx = he1.idx();
			auto fromVertex1 = this->mesh.from_vertex_handle(he1);
			auto toVertex1 = this->mesh.to_vertex_handle(he1);
			auto fromPt1 = this->mesh.point(fromVertex1);
			auto toPt1 = this->mesh.point(toVertex1);
			auto e1 = toPt1 - fromPt1;
			Eigen::Vector3d ev1;
			ev1 << e1.data()[0], e1.data()[1], e1.data()[2];
			++fhe_it1;
			auto he2 = *fhe_it1;
			auto he2idx = he2.idx();
			auto fromVertex2 = this->mesh.from_vertex_handle(he2);
			auto toVertex2 = this->mesh.to_vertex_handle(he2);
			auto fromPt2 = this->mesh.point(fromVertex2);
			auto toPt2 = this->mesh.point(toVertex2);
			auto e2 = toPt2 - fromPt2;
			Eigen::Vector3d ev2;
			ev2 << e2.data()[0], e2.data()[1], e2.data()[2];
			normal1 = ev1.cross(ev2);
			normal1 = normal1 / normal1.norm();
		}
		{
			auto fhe_it2 = this->mesh.fh_begin(f2);
			auto he1 = *fhe_it2;
			auto he1idx = he1.idx();
			auto fromVertex1 = this->mesh.from_vertex_handle(he1);
			auto toVertex1 = this->mesh.to_vertex_handle(he1);
			auto fromPt1 = this->mesh.point(fromVertex1);
			auto toPt1 = this->mesh.point(toVertex1);
			auto e1 = toPt1 - fromPt1;
			Eigen::Vector3d ev1;
			ev1 << e1.data()[0], e1.data()[1], e1.data()[2];
			++fhe_it2;
			auto he2 = *fhe_it2;
			auto he2idx = he2.idx();
			auto fromVertex2 = this->mesh.from_vertex_handle(he2);
			auto toVertex2 = this->mesh.to_vertex_handle(he2);
			auto fromPt2 = this->mesh.point(fromVertex2);
			auto toPt2 = this->mesh.point(toVertex2);
			auto e2 = toPt2 - fromPt2;
			Eigen::Vector3d ev2;
			ev2 << e2.data()[0], e2.data()[1], e2.data()[2];
			normal2 = ev1.cross(ev2);
			normal2 = normal2 / normal2.norm();
		}

		auto e_ij_tmp = toPt - fromPt;
		Eigen::Vector3d e_ij;
		e_ij << e_ij_tmp.data()[0], e_ij_tmp.data()[1], e_ij_tmp.data()[2];
		e_ij = e_ij / (e_ij.norm());

		decimal phi = CalWLZHDiheral(normal1, normal2, e_ij);
		this->ori_wlz_dihedral_set.push_back(phi);
		this->wnt_wlz_dihedral_set.push_back(phi);
		this->ori_dihedral_set.push_back(phi);
		//cout << ori_dihedral_set.size();
		this->wnt_dihedral_set.push_back(phi + cm1(2 * e_idx + 0));
		//cout << wnt_dihedral_set.size();

		if (phi > 0)
		{
			sign_flag.push_back(1);
		}
		else
		{
			sign_flag.push_back(-1);
		}

		//wdi_ofile << phi << std::endl;

		//edge, not half edge
		++he_it;
	}
	//mydi_ofile.close();
	//wdi_ofile.close();

	for (integer i = 0; i < this->ori_dihedral_set.size(); i++)
	{
		decimal ws = alpha_s / (ori_length_set[i] * ori_length_set[i]);
		decimal sqr_ws = std::sqrt(ws);
		this->sws_set.push_back(sqr_ws);
		decimal wb = (alpha_b * ori_length_set[i] * ori_length_set[i]) / ori_area_set[i];
		decimal sqr_wb = std::sqrt(wb);
		this->swb_set.push_back(sqr_wb);
		if (_isnan(wb))
		{
			cout << sqr_wb << endl;
			cout << ori_area_set[i] << endl;
			cout << sqr_ws << endl;
			cout << ori_length_set[i] << endl;
			cout << alpha_b << endl;
			cout << ori_area_set[i] << endl;
		}
	}

	//volume initialization
	ori_volume = 0.0;
	if (alpha_v != 0)
	{
		for (auto f_it = this->mesh.faces_begin(); f_it != this->mesh.faces_end(); ++f_it)
		{
			auto face = *f_it;
			pt_vector pt_vec;
			for (auto fv_it = this->mesh.fv_begin(face); fv_it != this->mesh.fv_end(face); ++fv_it)
			{
				auto v = *fv_it;
				auto pt = this->mesh.point(v);
				pt_vec.push_back(pt);
			}
			decimal temp_vol = (1.0 / 6.0) * ((pt_vec[0] % pt_vec[1]) | pt_vec[2]);
			ori_volume += temp_vol;
		}
		wnt_volume = ori_volume;
		decimal wv = alpha_v / (ori_volume * ori_volume);
		swv = std::sqrt(wv);
	}

	//average
	Eigen::VectorXd avg = Eigen::VectorXd::Zero(2 * this->mesh.n_edges());
#pragma omp parallel
	{
#pragma omp for
		for (integer i = 0; i < this->mesh.n_edges(); i++)
		{
			avg(2 * i + 1) = ori_dihedral_set[i];
			avg(2 * i + 0) = ori_length_set[i];
		}
	}

	this->data_avg = avg;
	Eigen::VectorXd tmp_avg = Eigen::VectorXd::Zero(2 * this->mesh.n_edges() + 1);
	if (alpha_v != 0)
	{
		tmp_avg.block(0, 0, avg.rows(), 1) = avg;
		tmp_avg(avg.rows()) = ori_volume;
		this->data_avg.resize(2 * this->mesh.n_edges() + 1);
		this->data_avg = tmp_avg;
	}


	//constraint initialization
	//CreateConstraintPointSet(GetConSetByV(s_vid_vec, con_v_mat));
	std::ifstream input_ifile("conset_nm.txt", std::ios::in);
	integer choice;
	input_ifile >> choice;
	input_ifile.close();
	if (choice == 1)
	{
		CreateConstraintPointSet(GetConSetByMesh());
	}
	else if (choice == 2)
	{
		if (cps.set.size() == 0)
		{
			std::cout << "no point constrins\n";
		}
	}
	else
	{
		mns::MyPointSet p1 = GetConSetByInput(current_str);
		CreateConstraintPointSet1(p1);
	}
	//PrintPtControl();
	decimal con_set_num = ((decimal)cps.set.size());
	decimal wc = alpha_c / con_set_num;
	swc = std::sqrt(wc);
	integer cnt = 0;
	bool cnt_flag = true;
	//this->move_pt_num = 0;

	//build unknown	
	if (this->use_eg_or_not == false && this->use_basis_or_not == false)
	{
		if (this->use_sim_or_not == false)
		{
			//soft constraint(soft con): unknown includes constraint point, constraint term must not be 0
			unknown = Eigen::VectorXd::Zero(3 * this->mesh.n_vertices());

			//hard const(hard con): unknown excludes constraint point, constraint term must be 0
			//unknown = Eigen::VectorXd::Zero(3 * (this->mesh.n_vertices() - this->cps.set.size()));

			int i = 0;
			for (auto it = this->mesh.vertices_begin(); it != this->mesh.vertices_end(); ++it, ++i)
			{
				auto point = this->mesh.point(*it);
				//cal moved pt
				if (i == cps.set[cnt].idx && cnt_flag)
				{
					int check_res = (CheckMovedPt(point, Eigen::RowVector3d(cps.set[cnt].x, cps.set[cnt].y, cps.set[cnt].z)));
					this->move_pt_num += check_res;
					this->pt_move_or_not.push_back(check_res);
					++cnt;
					if (cnt == cps.set.size())
					{
						cnt--;
						cnt_flag = false;
					}
				}

				//soft con
				unknown(3 * i + 0) = point.data()[0];
				unknown(3 * i + 1) = point.data()[1];
				unknown(3 * i + 2) = point.data()[2];

				//hard con
				/*if (idx_move[i] != -1)
				{
				unknown(3 * (i - idx_move[i]) + 0) = point.data()[0];
				unknown(3 * (i - idx_move[i]) + 1) = point.data()[1];
				unknown(3 * (i - idx_move[i]) + 2) = point.data()[2];
				}
				else
				{
				continue;
				}*/
			}
		}
		if (this->use_sim_or_not == true)
		{
			//soft constraint(soft con): unknown includes constraint point, constraint term must not be 0
			unknown = Eigen::VectorXd::Zero(3 * this->sim_mesh.n_vertices());
			high_v_mat = Eigen::MatrixXd::Zero(mesh.n_vertices(), 3);

			//hard const(hard con): unknown excludes constraint point, constraint term must be 0(not ready);

			int i = 0;
			for (auto it = this->sim_mesh.vertices_begin(); it != this->sim_mesh.vertices_end(); ++it, ++i)
			{
				auto point = this->sim_mesh.point(*it);
				//soft con
				unknown(3 * i + 0) = point.data()[0];
				unknown(3 * i + 1) = point.data()[1];
				unknown(3 * i + 2) = point.data()[2];

				//hard con(not ready)
			}
			i = 0;
			for (auto it = this->mesh.vertices_begin(); it != this->mesh.vertices_end(); ++it, ++i)
			{
				auto point = this->mesh.point(*it);
				//cal moved pt
				if (i == cps.set[cnt].idx && cnt_flag)
				{
					int check_res = (CheckMovedPt(point, Eigen::RowVector3d(cps.set[cnt].x, cps.set[cnt].y, cps.set[cnt].z)));
					this->move_pt_num += check_res;
					this->pt_move_or_not.push_back(check_res);
					++cnt;
					if (cnt == cps.set.size())
					{
						cnt--;
						cnt_flag = false;
					}
				}

				//soft con
				high_v_mat(i, 0) = point.data()[0];
				high_v_mat(i, 1) = point.data()[1];
				high_v_mat(i, 2) = point.data()[2];

				//hard con
				/*if (idx_move[i] != -1)
				{
				unknown(3 * (i - idx_move[i]) + 0) = point.data()[0];
				unknown(3 * (i - idx_move[i]) + 1) = point.data()[1];
				unknown(3 * (i - idx_move[i]) + 2) = point.data()[2];
				}
				else
				{
				continue;
				}*/
			}
		}
	}

	//example-driven initialization inside build unknown
	if (this->use_eg_or_not == true && this->use_basis_or_not == false)
	{
		EgInit();
		if (this->use_sim_or_not == false)
		{
			//soft con
			unknown = Eigen::VectorXd::Zero(3 * mesh.n_vertices() + this->eg_num);

			//hard con
			//unknown = Eigen::VectorXd::Zero(3 * (this->mesh.n_vertices() - this->cps.set.size()) + this->eg_num);

			int i = 0;
			for (auto it = this->mesh.vertices_begin(); it != this->mesh.vertices_end(); ++it, ++i)
			{
				auto point = this->mesh.point(*it);
				//cal moved pt
				if (i == cps.set[cnt].idx && cnt_flag)
				{
					int check_res = (CheckMovedPt(point, Eigen::RowVector3d(cps.set[cnt].x, cps.set[cnt].y, cps.set[cnt].z)));
					this->move_pt_num += check_res;
					this->pt_move_or_not.push_back(check_res);
					++cnt;
					if (cnt == cps.set.size())
					{
						cnt--;
						cnt_flag = false;
					}
				}

				//soft con
				unknown(3 * i + 0) = point.data()[0];
				unknown(3 * i + 1) = point.data()[1];
				unknown(3 * i + 2) = point.data()[2];

				//hard con
				/*if (idx_move[i] != -1)
				{
				unknown(3 * (i - idx_move[i]) + 0) = point.data()[0];
				unknown(3 * (i - idx_move[i]) + 1) = point.data()[1];
				unknown(3 * (i - idx_move[i]) + 2) = point.data()[2];
				}
				else
				{
				continue;
				}*/
			}

			//soft con
			for (integer k = 3 * i; k < 3 * i + this->eg_num; k++)
				//hard con
				//for (integer k = 3 * (i - this->cps.set.size()); k < 3 * (i - this->cps.set.size()) + this->eg_num; k++)
			{
				decimal tmp = 1.0 / ((decimal)eg_num);
				//decimal tmp = 0.0;
				unknown(k) = tmp;
			}
		}
		if (this->use_sim_or_not == true)
		{
			//soft con
			unknown = Eigen::VectorXd::Zero(3 * mesh.n_vertices() + this->eg_num);
			high_v_mat = Eigen::MatrixXd::Zero(mesh.n_vertices(), 3);

			//hard con(not ready)

			int i = 0;
			for (auto it = this->sim_mesh.vertices_begin(); it != this->sim_mesh.vertices_end(); ++it, ++i)
			{
				auto point = this->sim_mesh.point(*it);
				//soft con
				unknown(3 * i + 0) = point.data()[0];
				unknown(3 * i + 1) = point.data()[1];
				unknown(3 * i + 2) = point.data()[2];

				//hard con(not ready)
			}
			int j = 0;
			for (auto it = this->mesh.vertices_begin(); it != this->mesh.vertices_end(); ++it, ++j)
			{
				auto point = this->mesh.point(*it);
				//cal moved pt
				if (i == cps.set[cnt].idx && cnt_flag)
				{
					int check_res = (CheckMovedPt(point, Eigen::RowVector3d(cps.set[cnt].x, cps.set[cnt].y, cps.set[cnt].z)));
					this->move_pt_num += check_res;
					this->pt_move_or_not.push_back(check_res);
					++cnt;
					if (cnt == cps.set.size())
					{
						cnt--;
						cnt_flag = false;
					}
				}

				//soft con
				high_v_mat(j, 0) = point.data()[0];
				high_v_mat(j, 1) = point.data()[1];
				high_v_mat(j, 2) = point.data()[2];

				//hard con
				/*if (idx_move[i] != -1)
				{
				unknown(3 * (i - idx_move[i]) + 0) = point.data()[0];
				unknown(3 * (i - idx_move[i]) + 1) = point.data()[1];
				unknown(3 * (i - idx_move[i]) + 2) = point.data()[2];
				}
				else
				{
				continue;
				}*/
			}
			//soft con
			for (integer k = 3 * i; k < 3 * i + this->eg_num; k++)
				//hard con
				//for (integer k = 3 * (i - this->cps.set.size()); k < 3 * (i - this->cps.set.size()) + this->eg_num; k++)
			{
				decimal tmp = 1.0 / ((decimal)eg_num);
				//decimal tmp = 0.0;
				unknown(k) = tmp;
			}
		}
	}
	//basis initialization inside build unknown
	if (this->use_eg_or_not == false && this->use_basis_or_not == true)
	{
		//CreateBasis();
		//receive basis
		if (this->rb_flag == false)
		{
			ReceiveBasis();
		}
		if (this->use_sim_or_not == false && this->use_ww_or_not == false)//1
		{
			//soft con
			//unknown = Eigen::VectorXd::Zero(3 * mesh.n_vertices() + this->basis.cols());
			unknown = Eigen::VectorXd::Zero(3 * mesh.n_vertices() + 2 * this->basis.cols());
			MyMesh *mesh2 = new MyMesh();
			OpenMesh::IO::read_mesh(*mesh2, "refMesh2.off");
			//hard con
			//unknown = Eigen::VectorXd::Zero(3 * (this->mesh.n_vertices() - this->cps.set.size()) + this->basis.cols());

			int i = 0;
			for (auto it = mesh2->vertices_begin(); it != mesh2->vertices_end(); ++it, ++i)
			{
				auto point = mesh2->point(*it);

				//cal moved pt
				/*if (i == cps.set[cnt].idx && cnt_flag)
				{
				int check_res = (CheckMovedPt(point, Eigen::RowVector3d(cps.set[cnt].x, cps.set[cnt].y, cps.set[cnt].z)));
				this->move_pt_num += check_res;
				this->pt_move_or_not.push_back(check_res);
				++cnt;
				if (cnt == cps.set.size())
				{
				cnt--;
				cnt_flag = false;
				}
				}*/

				//soft con
				unknown(3 * i + 0) = point.data()[0];
				unknown(3 * i + 1) = point.data()[1];
				unknown(3 * i + 2) = point.data()[2];
				//cout << point.data()[0] << " " << point.data()[1] << " " << point.data()[2] << endl;

				//hard con
				/*if (idx_move[i] != -1)
				{
				unknown(3 * (i - idx_move[i]) + 0) = point.data()[0];
				unknown(3 * (i - idx_move[i]) + 1) = point.data()[1];
				unknown(3 * (i - idx_move[i]) + 2) = point.data()[2];
				}
				else
				{
				continue;
				}*/
			}

			//soft con
			//for (integer k = 3 * i; k < 3 * i + 2 * this->basis.cols(); k++)
			//	//hard con
			//	//for (integer k = 3 * (i - this->cps.set.size()); k < 3 * (i - this->cps.set.size()) + this->basis.cols(); k++)
			//{
			//	decimal tmp = 1.0 / ((decimal)this->basis.cols());
			//	//decimal tmp = 0;
			//	unknown(k) = tmp;
			//	cout << unknown(k) << endl;
			//}
			ifstream readtxt("initialize.txt");
			for (integer k = 3 * i; k < 3 * i + 2 * this->basis.cols(); k++)
			{
				readtxt >> unknown(k);
				//cout << unknown(k) << endl;
			}
		}
		if (this->use_sim_or_not == true && this->use_ww_or_not == false)//2
		{
			//soft con
			unknown = Eigen::VectorXd::Zero(3 * sim_mesh.n_vertices() + 2 * this->basis.cols());
			high_v_mat = Eigen::MatrixXd::Zero(mesh.n_vertices(), 3);

			//hard con(not ready)

			int i = 0;
			for (auto it = this->sim_mesh.vertices_begin(); it != this->sim_mesh.vertices_end(); ++it, ++i)
			{
				auto point = this->sim_mesh.point(*it);
				//soft con
				unknown(3 * i + 0) = point.data()[0];
				unknown(3 * i + 1) = point.data()[1];
				unknown(3 * i + 2) = point.data()[2];

				//hard con(not ready)
			}
			int j = 0;
			for (auto it = this->mesh.vertices_begin(); it != this->mesh.vertices_end(); ++it, ++j)
			{
				auto point = this->mesh.point(*it);
				//cal moved pt
				if (i == cps.set[cnt].idx && cnt_flag)
				{
					int check_res = (CheckMovedPt(point, Eigen::RowVector3d(cps.set[cnt].x, cps.set[cnt].y, cps.set[cnt].z)));
					this->move_pt_num += check_res;
					this->pt_move_or_not.push_back(check_res);
					++cnt;
					if (cnt == cps.set.size())
					{
						cnt--;
						cnt_flag = false;
					}
				}

				//soft con
				high_v_mat(j, 0) = point.data()[0];
				high_v_mat(j, 1) = point.data()[1];
				high_v_mat(j, 2) = point.data()[2];

				//hard con
				/*if (idx_move[i] != -1)
				{
				unknown(3 * (i - idx_move[i]) + 0) = point.data()[0];
				unknown(3 * (i - idx_move[i]) + 1) = point.data()[1];
				unknown(3 * (i - idx_move[i]) + 2) = point.data()[2];
				}
				else
				{
				continue;
				}*/
			}
			//soft con
			for (integer k = 3 * i; k < 3 * i + 2 * this->basis.cols(); k++)
				//hard con(not ready)
			{
				decimal tmp = 1.0 / (2 * (decimal)this->basis.cols());
				//decimal tmp = 0.0;
				unknown(k) = tmp;
			}
		}
		if (this->use_sim_or_not == false && this->use_ww_or_not == true)//3
		{
			MyMesh *mesh2 = new MyMesh();
			OpenMesh::IO::read_mesh(*mesh2, "out_ini.obj");
			//soft con
			//unknown = Eigen::VectorXd::Zero(3 * mesh.n_vertices() + this->basis.cols());
			unknown = Eigen::VectorXd::Zero(3 * mesh.n_vertices() + 4 * this->basis.cols());

			//hard con
			//unknown = Eigen::VectorXd::Zero(3 * (this->mesh.n_vertices() - this->cps.set.size()) + this->basis.cols());

			int i = 0;
			for (auto it = mesh2->vertices_begin(); it != mesh2->vertices_end(); ++it, ++i)
			{
				auto point = mesh2->point(*it);
				//cal moved pt
				if (i == cps.set[cnt].idx && cnt_flag)
				{
					int check_res = (CheckMovedPt(point, Eigen::RowVector3d(cps.set[cnt].x, cps.set[cnt].y, cps.set[cnt].z)));
					this->move_pt_num += check_res;
					this->pt_move_or_not.push_back(check_res);
					++cnt;
					if (cnt == cps.set.size())
					{
						cnt--;
						cnt_flag = false;
					}
				}

				//soft con
				unknown(3 * i + 0) = point.data()[0];
				unknown(3 * i + 1) = point.data()[1];
				unknown(3 * i + 2) = point.data()[2];

				//hard con
				/*if (idx_move[i] != -1)
				{
				unknown(3 * (i - idx_move[i]) + 0) = point.data()[0];
				unknown(3 * (i - idx_move[i]) + 1) = point.data()[1];
				unknown(3 * (i - idx_move[i]) + 2) = point.data()[2];
				}
				else
				{
				continue;
				}*/
			}

			//soft con
			for (integer k = 3 * i; k < 3 * i + 2 * this->basis.cols(); k++)
				//hard con
				//for (integer k = 3 * (i - this->cps.set.size()); k < 3 * (i - this->cps.set.size()) + this->basis.cols(); k++)
			{
				decimal tmp = 1.0 / ((decimal)this->basis.cols());
				//decimal tmp = 0;
				unknown(k) = tmp;
			}
			for (integer k = 3 * i + 2 * this->basis.cols(); k < 3 * i + 4 * this->basis.cols(); k++)
				//hard con(not ready)
			{
				//decimal tmp = 1.0 / (4.0 * (decimal)this->basis.cols());
				decimal tmp = 5.0 / ((decimal)this->basis.cols());
				unknown(k) = tmp;
			}
		}
		if (this->use_sim_or_not == true && this->use_ww_or_not == true)//4
		{
			//soft con
			unknown = Eigen::VectorXd::Zero(3 * sim_mesh.n_vertices() + 4 * this->basis.cols());
			high_v_mat = Eigen::MatrixXd::Zero(mesh.n_vertices(), 3);

			//hard con(not ready)

			int i = 0;
			for (auto it = this->sim_mesh.vertices_begin(); it != this->sim_mesh.vertices_end(); ++it, ++i)
			{
				auto point = this->sim_mesh.point(*it);
				//soft con
				unknown(3 * i + 0) = point.data()[0];
				unknown(3 * i + 1) = point.data()[1];
				unknown(3 * i + 2) = point.data()[2];

				//hard con(not ready)
			}
			int j = 0;
			for (auto it = this->mesh.vertices_begin(); it != this->mesh.vertices_end(); ++it, ++j)
			{
				auto point = this->mesh.point(*it);
				//cal moved pt
				if (i == cps.set[cnt].idx && cnt_flag)
				{
					int check_res = (CheckMovedPt(point, Eigen::RowVector3d(cps.set[cnt].x, cps.set[cnt].y, cps.set[cnt].z)));
					this->move_pt_num += check_res;
					this->pt_move_or_not.push_back(check_res);
					++cnt;
					if (cnt == cps.set.size())
					{
						cnt--;
						cnt_flag = false;
					}
				}

				//soft con
				high_v_mat(j, 0) = point.data()[0];
				high_v_mat(j, 1) = point.data()[1];
				high_v_mat(j, 2) = point.data()[2];

				//hard con
				/*if (idx_move[i] != -1)
				{
				unknown(3 * (i - idx_move[i]) + 0) = point.data()[0];
				unknown(3 * (i - idx_move[i]) + 1) = point.data()[1];
				unknown(3 * (i - idx_move[i]) + 2) = point.data()[2];
				}
				else
				{
				continue;
				}*/
			}
			//soft con
			for (integer k = 3 * i; k < 3 * i + 2 * this->basis.cols(); k++)
				//hard con(not ready)
			{
				//decimal tmp = 1.0 / (4.0 * (decimal)this->basis.cols());
				decimal tmp = 0.0;
				unknown(k) = tmp;
			}
			for (integer k = 3 * i + 2 * this->basis.cols(); k < 3 * i + 4 * this->basis.cols(); k++)
				//hard con(not ready)
			{
				//decimal tmp = 1.0 / (4.0 * (decimal)this->basis.cols());
				decimal tmp = 0.1;
				unknown(k) = tmp;
			}
		}
	}

	if (this->use_eg_or_not == true && this->use_basis_or_not == true)
	{
		if (this->rb_flag == false)
		{
			ReceiveBasis();
		}
		if (this->use_sim_or_not == false && this->use_ww_or_not == false)//1
		{
			unknown = Eigen::VectorXd::Zero(3 * mesh.n_vertices());
			l_bw_unknown = Eigen::VectorXd::Zero(this->basis.cols());
			d_bw_unknown = Eigen::VectorXd::Zero(this->basis.cols());

			int i = 0;
			for (auto it = this->mesh.vertices_begin(); it != this->mesh.vertices_end(); ++it, ++i)
			{
				auto point = this->mesh.point(*it);
				//cal moved pt
				if (i == cps.set[cnt].idx && cnt_flag)
				{
					int check_res = (CheckMovedPt(point, Eigen::RowVector3d(cps.set[cnt].x, cps.set[cnt].y, cps.set[cnt].z)));
					this->move_pt_num += check_res;
					this->pt_move_or_not.push_back(check_res);
					++cnt;
					if (cnt == cps.set.size())
					{
						cnt--;
						cnt_flag = false;
					}
				}

				//soft con
				unknown(3 * i + 0) = point.data()[0];
				unknown(3 * i + 1) = point.data()[1];
				unknown(3 * i + 2) = point.data()[2];

				//hard con
				/*if (idx_move[i] != -1)
				{
				unknown(3 * (i - idx_move[i]) + 0) = point.data()[0];
				unknown(3 * (i - idx_move[i]) + 1) = point.data()[1];
				unknown(3 * (i - idx_move[i]) + 2) = point.data()[2];
				}
				else
				{
				continue;
				}*/
			}

			//soft con
			for (integer k = 0; k < this->basis.cols(); k++)
				//hard con
				//for (integer k = 3 * (i - this->cps.set.size()); k < 3 * (i - this->cps.set.size()) + this->basis.cols(); k++)
			{
				decimal tmp = 1.0 / ((decimal)this->basis.cols());
				//decimal tmp = 0;
				l_bw_unknown(k) = tmp;
				d_bw_unknown(k) = tmp;
			}
		}
		if (this->use_sim_or_not == true && this->use_ww_or_not == false)//2
		{
			//soft con
			unknown = Eigen::VectorXd::Zero(3 * sim_mesh.n_vertices());
			l_bw_unknown = Eigen::VectorXd::Zero(this->basis.cols());
			d_bw_unknown = Eigen::VectorXd::Zero(this->basis.cols());
			high_v_mat = Eigen::MatrixXd::Zero(mesh.n_vertices(), 3);

			//hard con(not ready)

			int i = 0;
			for (auto it = this->sim_mesh.vertices_begin(); it != this->sim_mesh.vertices_end(); ++it, ++i)
			{
				auto point = this->sim_mesh.point(*it);
				//soft con
				unknown(3 * i + 0) = point.data()[0];
				unknown(3 * i + 1) = point.data()[1];
				unknown(3 * i + 2) = point.data()[2];

				//hard con(not ready)
			}
			int j = 0;
			for (auto it = this->mesh.vertices_begin(); it != this->mesh.vertices_end(); ++it, ++j)
			{
				auto point = this->mesh.point(*it);
				//cal moved pt
				if (i == cps.set[cnt].idx && cnt_flag)
				{
					int check_res = (CheckMovedPt(point, Eigen::RowVector3d(cps.set[cnt].x, cps.set[cnt].y, cps.set[cnt].z)));
					this->move_pt_num += check_res;
					this->pt_move_or_not.push_back(check_res);
					++cnt;
					if (cnt == cps.set.size())
					{
						cnt--;
						cnt_flag = false;
					}
				}

				//soft con
				high_v_mat(j, 0) = point.data()[0];
				high_v_mat(j, 1) = point.data()[1];
				high_v_mat(j, 2) = point.data()[2];

				//hard con
				/*if (idx_move[i] != -1)
				{
				unknown(3 * (i - idx_move[i]) + 0) = point.data()[0];
				unknown(3 * (i - idx_move[i]) + 1) = point.data()[1];
				unknown(3 * (i - idx_move[i]) + 2) = point.data()[2];
				}
				else
				{
				continue;
				}*/
			}
			//soft con
			for (integer k = 0; k < this->basis.cols(); k++)
				//hard con(not ready)
			{
				decimal tmp = 1.0 / (2 * (decimal)this->basis.cols());
				//decimal tmp = 0.0;
				l_bw_unknown(k) = tmp;
				d_bw_unknown(k) = tmp;
			}
		}
	}
	//cout << unknown << endl;
	pre_unknown = unknown;
}


void mns::MyDeformationHandler::DeformInit(std::vector<integer>& s_vid_vec, Eigen::MatrixXd& v_mat, Eigen::MatrixXd& con_v_mat)
{
	this->con_ele_cnt = 0;
	this->str_ele_cnt = 0;
	this->ben_ele_cnt = 0;
	this->sb_ele_cnt = 0;
	this->row_cnt = 0;
	this->pre_energy = 100000000000000;
	this->iter = 0;

	//std::ofstream mydi_ofile("mydihedral.txt");
	//std::ofstream wdi_ofile("wdihedral.txt");

	for (auto he_it = this->mesh.halfedges_begin(); he_it != this->mesh.halfedges_end(); ++he_it)
	{
		auto e = *he_it;
		auto fromVertex = this->mesh.from_vertex_handle(e);
		auto toVertex = this->mesh.to_vertex_handle(e);
		auto fromPt = this->mesh.point(fromVertex);
		auto toPt = this->mesh.point(toVertex);

		//stretching initialization		
		decimal l = this->mesh.calc_edge_length(e);
		this->ori_length_set.push_back(l);
		this->wnt_length_set.push_back(l);

		//bend initialization
		auto f1 = this->mesh.face_handle(e);
		auto op_he = this->mesh.opposite_halfedge_handle(e);
		auto f2 = this->mesh.face_handle(op_he);
		pt_vector pt_set;
		pt_vector t1_set;
		pt_vector t2_set;

		if (f1.idx() == -1 || f2.idx() == -1)
		{
			this->ori_dihedral_set.push_back(0);
			this->wnt_dihedral_set.push_back(0);
			ori_area_set.push_back(0);

			this->ori_wlz_dihedral_set.push_back(0);
			this->wnt_wlz_dihedral_set.push_back(0);
			sign_flag.push_back(1);

			//edge, not half edge
			++he_it;
			continue;
		}

		pt_set.push_back(fromPt);
		t1_set.push_back(fromPt);
		t2_set.push_back(fromPt);
		pt_set.push_back(toPt);
		t1_set.push_back(toPt);
		t2_set.push_back(toPt);

		for (auto fv_it1 = this->mesh.fv_begin(f1); fv_it1 != this->mesh.fv_end(f1); ++fv_it1)
		{
			auto v = *fv_it1;
			auto pt = this->mesh.point(v);
			if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
			{
				pt_set.push_back(pt);
				t1_set.push_back(pt);
				break;
			}
		}
		integer temp_cnt = 0;
		for (auto fv_it2 = this->mesh.fv_begin(f2); fv_it2 != this->mesh.fv_end(f2); ++fv_it2)
		{
			auto v = *fv_it2;
			auto pt = this->mesh.point(v);
			if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
			{
				t2_set.push_back(pt);
				pt_set.push_back(pt);
				break;
			}
		}
		decimal t1_side[3];
		t1_side[0] = (t1_set[0] - t1_set[1]).length();
		t1_side[1] = (t1_set[2] - t1_set[1]).length();
		t1_side[2] = (t1_set[0] - t1_set[2]).length();

		decimal p1 = (t1_side[0] + t1_side[1] + t1_side[2]) / 2; //半周长;  
		decimal area1 = std::sqrt(p1*(p1 - t1_side[0])*(p1 - t1_side[1])*(p1 - t1_side[2]));

		decimal t2_side[3];
		t2_side[0] = (t2_set[0] - t2_set[1]).length();
		t2_side[1] = (t2_set[2] - t2_set[1]).length();
		t2_side[2] = (t2_set[0] - t2_set[2]).length();

		decimal p2 = (t2_side[0] + t2_side[1] + t2_side[2]) / 2; //半周长;  
		decimal area2 = std::sqrt(p2*(p2 - t2_side[0])*(p2 - t2_side[1])*(p2 - t2_side[2]));
		ori_area_set.push_back(area1 + area2);
		Eigen::Vector3d pt0(pt_set[0][0], pt_set[0][1], pt_set[0][2]);
		Eigen::Vector3d pt1(pt_set[1][0], pt_set[1][1], pt_set[1][2]);
		Eigen::Vector3d pt2(pt_set[2][0], pt_set[2][1], pt_set[2][2]);
		Eigen::Vector3d pt3(pt_set[3][0], pt_set[3][1], pt_set[3][2]);
		//decimal angle = 3.1415926 - CalWLZHDiheral(pt0, pt1, pt2, pt3);
		//decimal angle = CalDiheral(pt_set[0], pt_set[1], pt_set[2], pt_set[3]);
		/*this->ori_dihedral_set.push_back(angle);
		this->wnt_dihedral_set.push_back(angle);*/

		//mydi_ofile << angle << std::endl;

		//sign_flag computation
		Eigen::Vector3d normal1;
		Eigen::Vector3d normal2;
		{
			auto fhe_it1 = this->mesh.fh_begin(f1);
			auto he1 = *fhe_it1;
			auto he1idx = he1.idx();
			auto fromVertex1 = this->mesh.from_vertex_handle(he1);
			auto toVertex1 = this->mesh.to_vertex_handle(he1);
			auto fromPt1 = this->mesh.point(fromVertex1);
			auto toPt1 = this->mesh.point(toVertex1);
			auto e1 = toPt1 - fromPt1;
			Eigen::Vector3d ev1;
			ev1 << e1.data()[0], e1.data()[1], e1.data()[2];
			++fhe_it1;
			auto he2 = *fhe_it1;
			auto he2idx = he2.idx();
			auto fromVertex2 = this->mesh.from_vertex_handle(he2);
			auto toVertex2 = this->mesh.to_vertex_handle(he2);
			auto fromPt2 = this->mesh.point(fromVertex2);
			auto toPt2 = this->mesh.point(toVertex2);
			auto e2 = toPt2 - fromPt2;
			Eigen::Vector3d ev2;
			ev2 << e2.data()[0], e2.data()[1], e2.data()[2];
			normal1 = ev1.cross(ev2);
			normal1 = normal1 / normal1.norm();
		}
		{
			auto fhe_it2 = this->mesh.fh_begin(f2);
			auto he1 = *fhe_it2;
			auto he1idx = he1.idx();
			auto fromVertex1 = this->mesh.from_vertex_handle(he1);
			auto toVertex1 = this->mesh.to_vertex_handle(he1);
			auto fromPt1 = this->mesh.point(fromVertex1);
			auto toPt1 = this->mesh.point(toVertex1);
			auto e1 = toPt1 - fromPt1;
			Eigen::Vector3d ev1;
			ev1 << e1.data()[0], e1.data()[1], e1.data()[2];
			++fhe_it2;
			auto he2 = *fhe_it2;
			auto he2idx = he2.idx();
			auto fromVertex2 = this->mesh.from_vertex_handle(he2);
			auto toVertex2 = this->mesh.to_vertex_handle(he2);
			auto fromPt2 = this->mesh.point(fromVertex2);
			auto toPt2 = this->mesh.point(toVertex2);
			auto e2 = toPt2 - fromPt2;
			Eigen::Vector3d ev2;
			ev2 << e2.data()[0], e2.data()[1], e2.data()[2];
			normal2 = ev1.cross(ev2);
			normal2 = normal2 / normal2.norm();
		}

		auto e_ij_tmp = toPt - fromPt;
		Eigen::Vector3d e_ij;
		e_ij << e_ij_tmp.data()[0], e_ij_tmp.data()[1], e_ij_tmp.data()[2];
		e_ij = e_ij / (e_ij.norm());

		decimal phi = CalWLZHDiheral(normal1, normal2, e_ij);
		this->ori_wlz_dihedral_set.push_back( phi);
		this->wnt_wlz_dihedral_set.push_back( phi);
		this->ori_dihedral_set.push_back( phi);
		//cout << ori_dihedral_set.size();
		this->wnt_dihedral_set.push_back(phi);
		//cout << wnt_dihedral_set.size();

		if (phi > 0)
		{
			sign_flag.push_back(1);
		}
		else
		{
			sign_flag.push_back(-1);
		}

		//wdi_ofile << phi << std::endl;

		//edge, not half edge
		++he_it;
	}
	//mydi_ofile.close();
	//wdi_ofile.close();

	for (integer i = 0; i < this->ori_dihedral_set.size(); i++)
	{
		decimal ws = alpha_s / (ori_length_set[i] * ori_length_set[i]);
		decimal sqr_ws = std::sqrt(ws);
		this->sws_set.push_back(sqr_ws);
		decimal wb = (alpha_b * ori_length_set[i] * ori_length_set[i]) / ori_area_set[i];
		decimal sqr_wb = std::sqrt(wb);
		this->swb_set.push_back(sqr_wb);
		if (_isnan(wb))
		{
			cout << sqr_wb << endl;
			cout << ori_area_set[i] << endl;
			cout << sqr_ws<<endl;
			cout << ori_length_set[i]<<endl;
			cout << alpha_b << endl;
			cout << ori_area_set[i] << endl;
		}
	}

	//volume initialization
	ori_volume = 0.0;
	if (alpha_v != 0)
	{
		for (auto f_it = this->mesh.faces_begin(); f_it != this->mesh.faces_end(); ++f_it)
		{
			auto face = *f_it;
			pt_vector pt_vec;
			for (auto fv_it = this->mesh.fv_begin(face); fv_it != this->mesh.fv_end(face); ++fv_it)
			{
				auto v = *fv_it;
				auto pt = this->mesh.point(v);
				pt_vec.push_back(pt);
			}
			decimal temp_vol = (1.0 / 6.0) * ((pt_vec[0] % pt_vec[1]) | pt_vec[2]);
			ori_volume += temp_vol;
		}
		wnt_volume = ori_volume;
		decimal wv = alpha_v / (ori_volume * ori_volume);
		swv = std::sqrt(wv);
	}

	//average
	Eigen::VectorXd avg = Eigen::VectorXd::Zero(2 * this->mesh.n_edges());
#pragma omp parallel
	{
#pragma omp for
		for (integer i = 0; i < this->mesh.n_edges(); i++)
		{
			avg(2 * i + 1) = ori_dihedral_set[i];
			avg(2 * i + 0) = ori_length_set[i];
		}
	}

	this->data_avg = avg;
	Eigen::VectorXd tmp_avg = Eigen::VectorXd::Zero(2 * this->mesh.n_edges() + 1);
	if (alpha_v != 0)
	{
		tmp_avg.block(0, 0, avg.rows(), 1) = avg;
		tmp_avg(avg.rows()) = ori_volume;
		this->data_avg.resize(2 * this->mesh.n_edges() + 1);
		this->data_avg = tmp_avg;
	}


	//constraint initialization
	//CreateConstraintPointSet(GetConSetByV(s_vid_vec, con_v_mat));
	std::ifstream input_ifile("conset_nm_test.txt", std::ios::in);//conset_nm.txt改为conset_nm_test.txt
	integer choice;
	input_ifile >> choice;
	input_ifile.close();
	if (choice == 1)
	{
		CreateConstraintPointSet(GetConSetByMesh());
	}
	else if (choice == 2)
	{
		if (cps.set.size() == 0)
		{
			std::cout << "no point constrins\n";
		}
	}
	else
	{
		mns::MyPointSet p1 = GetConSetByInput(current_str);
		CreateConstraintPointSet1(p1);
	}
	//PrintPtControl();
	decimal con_set_num = ((decimal)cps.set.size());
	decimal wc = alpha_c / con_set_num;
	swc = std::sqrt(wc);
	integer cnt = 0;
	bool cnt_flag = true;
	//this->move_pt_num = 0;

	//build unknown	
	if (this->use_eg_or_not == false && this->use_basis_or_not == false)
	{
		if (this->use_sim_or_not == false)
		{
			//soft constraint(soft con): unknown includes constraint point, constraint term must not be 0
			unknown = Eigen::VectorXd::Zero(3 * this->mesh.n_vertices());

			//hard const(hard con): unknown excludes constraint point, constraint term must be 0
			//unknown = Eigen::VectorXd::Zero(3 * (this->mesh.n_vertices() - this->cps.set.size()));

			int i = 0;
			for (auto it = this->mesh.vertices_begin(); it != this->mesh.vertices_end(); ++it, ++i)
			{
				auto point = this->mesh.point(*it);
				//cal moved pt
				if (i == cps.set[cnt].idx && cnt_flag)
				{
					int check_res = (CheckMovedPt(point, Eigen::RowVector3d(cps.set[cnt].x, cps.set[cnt].y, cps.set[cnt].z)));
					this->move_pt_num += check_res;
					this->pt_move_or_not.push_back(check_res);
					++cnt;
					if (cnt == cps.set.size())
					{
						cnt--;
						cnt_flag = false;
					}
				}

				//soft con
				unknown(3 * i + 0) = point.data()[0];
				unknown(3 * i + 1) = point.data()[1];
				unknown(3 * i + 2) = point.data()[2];

				//hard con
				/*if (idx_move[i] != -1)
				{
				unknown(3 * (i - idx_move[i]) + 0) = point.data()[0];
				unknown(3 * (i - idx_move[i]) + 1) = point.data()[1];
				unknown(3 * (i - idx_move[i]) + 2) = point.data()[2];
				}
				else
				{
				continue;
				}*/
			}
		}
		if (this->use_sim_or_not == true)
		{
			//soft constraint(soft con): unknown includes constraint point, constraint term must not be 0
			unknown = Eigen::VectorXd::Zero(3 * this->sim_mesh.n_vertices());
			high_v_mat = Eigen::MatrixXd::Zero(mesh.n_vertices(), 3);

			//hard const(hard con): unknown excludes constraint point, constraint term must be 0(not ready);

			int i = 0;
			for (auto it = this->sim_mesh.vertices_begin(); it != this->sim_mesh.vertices_end(); ++it, ++i)
			{
				auto point = this->sim_mesh.point(*it);
				//soft con
				unknown(3 * i + 0) = point.data()[0];
				unknown(3 * i + 1) = point.data()[1];
				unknown(3 * i + 2) = point.data()[2];

				//hard con(not ready)
			}
			i = 0;
			for (auto it = this->mesh.vertices_begin(); it != this->mesh.vertices_end(); ++it, ++i)
			{
				auto point = this->mesh.point(*it);
				//cal moved pt
				if (i == cps.set[cnt].idx && cnt_flag)
				{
					int check_res = (CheckMovedPt(point, Eigen::RowVector3d(cps.set[cnt].x, cps.set[cnt].y, cps.set[cnt].z)));
					this->move_pt_num += check_res;
					this->pt_move_or_not.push_back(check_res);
					++cnt;
					if (cnt == cps.set.size())
					{
						cnt--;
						cnt_flag = false;
					}
				}

				//soft con
				high_v_mat(i, 0) = point.data()[0];
				high_v_mat(i, 1) = point.data()[1];
				high_v_mat(i, 2) = point.data()[2];

				//hard con
				/*if (idx_move[i] != -1)
				{
				unknown(3 * (i - idx_move[i]) + 0) = point.data()[0];
				unknown(3 * (i - idx_move[i]) + 1) = point.data()[1];
				unknown(3 * (i - idx_move[i]) + 2) = point.data()[2];
				}
				else
				{
				continue;
				}*/
			}
		}
	}

	//example-driven initialization inside build unknown
	if (this->use_eg_or_not == true && this->use_basis_or_not == false)
	{
		EgInit();
		if (this->use_sim_or_not == false)
		{
			//soft con
			unknown = Eigen::VectorXd::Zero(3 * mesh.n_vertices() + this->eg_num);

			//hard con
			//unknown = Eigen::VectorXd::Zero(3 * (this->mesh.n_vertices() - this->cps.set.size()) + this->eg_num);

			int i = 0;
			for (auto it = this->mesh.vertices_begin(); it != this->mesh.vertices_end(); ++it, ++i)
			{
				auto point = this->mesh.point(*it);
				//cal moved pt
				if (i == cps.set[cnt].idx && cnt_flag)
				{
					int check_res = (CheckMovedPt(point, Eigen::RowVector3d(cps.set[cnt].x, cps.set[cnt].y, cps.set[cnt].z)));
					this->move_pt_num += check_res;
					this->pt_move_or_not.push_back(check_res);
					++cnt;
					if (cnt == cps.set.size())
					{
						cnt--;
						cnt_flag = false;
					}
				}

				//soft con
				unknown(3 * i + 0) = point.data()[0];
				unknown(3 * i + 1) = point.data()[1];
				unknown(3 * i + 2) = point.data()[2];

				//hard con
				/*if (idx_move[i] != -1)
				{
				unknown(3 * (i - idx_move[i]) + 0) = point.data()[0];
				unknown(3 * (i - idx_move[i]) + 1) = point.data()[1];
				unknown(3 * (i - idx_move[i]) + 2) = point.data()[2];
				}
				else
				{
				continue;
				}*/
			}

			//soft con
			for (integer k = 3 * i; k < 3 * i + this->eg_num; k++)
				//hard con
				//for (integer k = 3 * (i - this->cps.set.size()); k < 3 * (i - this->cps.set.size()) + this->eg_num; k++)
			{
				decimal tmp = 1.0 / ((decimal)eg_num);
				//decimal tmp = 0.0;
				unknown(k) = tmp;
			}
		}
		if (this->use_sim_or_not == true)
		{
			//soft con
			unknown = Eigen::VectorXd::Zero(3 * mesh.n_vertices() + this->eg_num);
			high_v_mat = Eigen::MatrixXd::Zero(mesh.n_vertices(), 3);

			//hard con(not ready)

			int i = 0;
			for (auto it = this->sim_mesh.vertices_begin(); it != this->sim_mesh.vertices_end(); ++it, ++i)
			{
				auto point = this->sim_mesh.point(*it);
				//soft con
				unknown(3 * i + 0) = point.data()[0];
				unknown(3 * i + 1) = point.data()[1];
				unknown(3 * i + 2) = point.data()[2];

				//hard con(not ready)
			}
			int j = 0;
			for (auto it = this->mesh.vertices_begin(); it != this->mesh.vertices_end(); ++it, ++j)
			{
				auto point = this->mesh.point(*it);
				//cal moved pt
				if (i == cps.set[cnt].idx && cnt_flag)
				{
					int check_res = (CheckMovedPt(point, Eigen::RowVector3d(cps.set[cnt].x, cps.set[cnt].y, cps.set[cnt].z)));
					this->move_pt_num += check_res;
					this->pt_move_or_not.push_back(check_res);
					++cnt;
					if (cnt == cps.set.size())
					{
						cnt--;
						cnt_flag = false;
					}
				}

				//soft con
				high_v_mat(j, 0) = point.data()[0];
				high_v_mat(j, 1) = point.data()[1];
				high_v_mat(j, 2) = point.data()[2];

				//hard con
				/*if (idx_move[i] != -1)
				{
				unknown(3 * (i - idx_move[i]) + 0) = point.data()[0];
				unknown(3 * (i - idx_move[i]) + 1) = point.data()[1];
				unknown(3 * (i - idx_move[i]) + 2) = point.data()[2];
				}
				else
				{
				continue;
				}*/
			}
			//soft con
			for (integer k = 3 * i; k < 3 * i + this->eg_num; k++)
				//hard con
				//for (integer k = 3 * (i - this->cps.set.size()); k < 3 * (i - this->cps.set.size()) + this->eg_num; k++)
			{
				decimal tmp = 1.0 / ((decimal)eg_num);
				//decimal tmp = 0.0;
				unknown(k) = tmp;
			}
		}
	}
	//basis initialization inside build unknown
	if (this->use_eg_or_not == false && this->use_basis_or_not == true)
	{
		//CreateBasis();
		//receive basis
		if (this->rb_flag == false)
		{
			ReceiveBasis();
		}
		if (this->use_sim_or_not == false && this->use_ww_or_not == false)//1
		{
			//soft con
			//unknown = Eigen::VectorXd::Zero(3 * mesh.n_vertices() + this->basis.cols());
			unknown = Eigen::VectorXd::Zero(3 * mesh.n_vertices() + 2 * this->basis.cols());
			MyMesh *mesh2 = new MyMesh();
			OpenMesh::IO::read_mesh(*mesh2, "refMesh2.off");
			//hard con
			//unknown = Eigen::VectorXd::Zero(3 * (this->mesh.n_vertices() - this->cps.set.size()) + this->basis.cols());

			int i = 0;
			for (auto it = mesh2->vertices_begin(); it != mesh2->vertices_end(); ++it, ++i)
			{
				auto point = mesh2->point(*it);

				//cal moved pt
				/*if (i == cps.set[cnt].idx && cnt_flag)
				{
					int check_res = (CheckMovedPt(point, Eigen::RowVector3d(cps.set[cnt].x, cps.set[cnt].y, cps.set[cnt].z)));
					this->move_pt_num += check_res;
					this->pt_move_or_not.push_back(check_res);
					++cnt;
					if (cnt == cps.set.size())
					{
						cnt--;
						cnt_flag = false;
					}
				}*/

				//soft con
				unknown(3 * i + 0) = point.data()[0];
				unknown(3 * i + 1) = point.data()[1];
				unknown(3 * i + 2) = point.data()[2];
				//cout << point.data()[0] << " " << point.data()[1] << " " << point.data()[2] << endl;

				//hard con
				/*if (idx_move[i] != -1)
				{
				unknown(3 * (i - idx_move[i]) + 0) = point.data()[0];
				unknown(3 * (i - idx_move[i]) + 1) = point.data()[1];
				unknown(3 * (i - idx_move[i]) + 2) = point.data()[2];
				}
				else
				{
				continue;
				}*/
			}

			//soft con
			//for (integer k = 3 * i; k < 3 * i + 2 * this->basis.cols(); k++)
			//	//hard con
			//	//for (integer k = 3 * (i - this->cps.set.size()); k < 3 * (i - this->cps.set.size()) + this->basis.cols(); k++)
			//{
			//	decimal tmp = 1.0 / ((decimal)this->basis.cols());
			//	//decimal tmp = 0;
			//	unknown(k) = tmp;
			//	cout << unknown(k) << endl;
			//}
			ifstream readtxt("initialize.txt");
			for (integer k = 3 * i; k < 3 * i + 2 * this->basis.cols(); k++)
			{
				readtxt >> unknown(k);
				//cout << unknown(k) << endl;
			}
		}
		if (this->use_sim_or_not == true && this->use_ww_or_not == false)//2
		{
			//soft con
			unknown = Eigen::VectorXd::Zero(3 * sim_mesh.n_vertices() + 2 * this->basis.cols());
			high_v_mat = Eigen::MatrixXd::Zero(mesh.n_vertices(), 3);

			//hard con(not ready)

			int i = 0;
			for (auto it = this->sim_mesh.vertices_begin(); it != this->sim_mesh.vertices_end(); ++it, ++i)
			{
				auto point = this->sim_mesh.point(*it);
				//soft con
				unknown(3 * i + 0) = point.data()[0];
				unknown(3 * i + 1) = point.data()[1];
				unknown(3 * i + 2) = point.data()[2];

				//hard con(not ready)
			}
			int j = 0;
			for (auto it = this->mesh.vertices_begin(); it != this->mesh.vertices_end(); ++it, ++j)
			{
				auto point = this->mesh.point(*it);
				//cal moved pt
				if (i == cps.set[cnt].idx && cnt_flag)
				{
					int check_res = (CheckMovedPt(point, Eigen::RowVector3d(cps.set[cnt].x, cps.set[cnt].y, cps.set[cnt].z)));
					this->move_pt_num += check_res;
					this->pt_move_or_not.push_back(check_res);
					++cnt;
					if (cnt == cps.set.size())
					{
						cnt--;
						cnt_flag = false;
					}
				}

				//soft con
				high_v_mat(j, 0) = point.data()[0];
				high_v_mat(j, 1) = point.data()[1];
				high_v_mat(j, 2) = point.data()[2];

				//hard con
				/*if (idx_move[i] != -1)
				{
				unknown(3 * (i - idx_move[i]) + 0) = point.data()[0];
				unknown(3 * (i - idx_move[i]) + 1) = point.data()[1];
				unknown(3 * (i - idx_move[i]) + 2) = point.data()[2];
				}
				else
				{
				continue;
				}*/
			}
			//soft con
			for (integer k = 3 * i; k < 3 * i + 2 * this->basis.cols(); k++)
				//hard con(not ready)
			{
				decimal tmp = 1.0 / (2 * (decimal)this->basis.cols());
				//decimal tmp = 0.0;
				unknown(k) = tmp;
			}
		}
		if (this->use_sim_or_not == false && this->use_ww_or_not == true)//3
		{
			MyMesh *mesh2 = new MyMesh();
			OpenMesh::IO::read_mesh(*mesh2, "out_ini.obj");
			//soft con
			//unknown = Eigen::VectorXd::Zero(3 * mesh.n_vertices() + this->basis.cols());
			unknown = Eigen::VectorXd::Zero(3 * mesh.n_vertices() + 4 * this->basis.cols());

			//hard con
			//unknown = Eigen::VectorXd::Zero(3 * (this->mesh.n_vertices() - this->cps.set.size()) + this->basis.cols());

			int i = 0;
			for (auto it = mesh2->vertices_begin(); it != mesh2->vertices_end(); ++it, ++i)
			{
				auto point = mesh2->point(*it);
				//cal moved pt
				if (i == cps.set[cnt].idx && cnt_flag)
				{
					int check_res = (CheckMovedPt(point, Eigen::RowVector3d(cps.set[cnt].x, cps.set[cnt].y, cps.set[cnt].z)));
					this->move_pt_num += check_res;
					this->pt_move_or_not.push_back(check_res);
					++cnt;
					if (cnt == cps.set.size())
					{
						cnt--;
						cnt_flag = false;
					}
				}

				//soft con
				unknown(3 * i + 0) = point.data()[0];
				unknown(3 * i + 1) = point.data()[1];
				unknown(3 * i + 2) = point.data()[2];

				//hard con
				/*if (idx_move[i] != -1)
				{
				unknown(3 * (i - idx_move[i]) + 0) = point.data()[0];
				unknown(3 * (i - idx_move[i]) + 1) = point.data()[1];
				unknown(3 * (i - idx_move[i]) + 2) = point.data()[2];
				}
				else
				{
				continue;
				}*/
			}

			//soft con
			for (integer k = 3 * i; k < 3 * i + 2 * this->basis.cols(); k++)
				//hard con
				//for (integer k = 3 * (i - this->cps.set.size()); k < 3 * (i - this->cps.set.size()) + this->basis.cols(); k++)
			{
				decimal tmp = 1.0 / ((decimal)this->basis.cols());
				//decimal tmp = 0;
				unknown(k) = tmp;
			}
			for (integer k = 3 * i + 2 * this->basis.cols(); k < 3 * i + 4 * this->basis.cols(); k++)
				//hard con(not ready)
			{
				//decimal tmp = 1.0 / (4.0 * (decimal)this->basis.cols());
				decimal tmp = 5.0 / ((decimal)this->basis.cols());
				unknown(k) = tmp;
			}
		}
		if (this->use_sim_or_not == true && this->use_ww_or_not == true)//4
		{
			//soft con
			unknown = Eigen::VectorXd::Zero(3 * sim_mesh.n_vertices() + 4 * this->basis.cols());
			high_v_mat = Eigen::MatrixXd::Zero(mesh.n_vertices(), 3);

			//hard con(not ready)

			int i = 0;
			for (auto it = this->sim_mesh.vertices_begin(); it != this->sim_mesh.vertices_end(); ++it, ++i)
			{
				auto point = this->sim_mesh.point(*it);
				//soft con
				unknown(3 * i + 0) = point.data()[0];
				unknown(3 * i + 1) = point.data()[1];
				unknown(3 * i + 2) = point.data()[2];

				//hard con(not ready)
			}
			int j = 0;
			for (auto it = this->mesh.vertices_begin(); it != this->mesh.vertices_end(); ++it, ++j)
			{
				auto point = this->mesh.point(*it);
				//cal moved pt
				if (i == cps.set[cnt].idx && cnt_flag)
				{
					int check_res = (CheckMovedPt(point, Eigen::RowVector3d(cps.set[cnt].x, cps.set[cnt].y, cps.set[cnt].z)));
					this->move_pt_num += check_res;
					this->pt_move_or_not.push_back(check_res);
					++cnt;
					if (cnt == cps.set.size())
					{
						cnt--;
						cnt_flag = false;
					}
				}

				//soft con
				high_v_mat(j, 0) = point.data()[0];
				high_v_mat(j, 1) = point.data()[1];
				high_v_mat(j, 2) = point.data()[2];

				//hard con
				/*if (idx_move[i] != -1)
				{
				unknown(3 * (i - idx_move[i]) + 0) = point.data()[0];
				unknown(3 * (i - idx_move[i]) + 1) = point.data()[1];
				unknown(3 * (i - idx_move[i]) + 2) = point.data()[2];
				}
				else
				{
				continue;
				}*/
			}
			//soft con
			for (integer k = 3 * i; k < 3 * i + 2 * this->basis.cols(); k++)
				//hard con(not ready)
			{
				//decimal tmp = 1.0 / (4.0 * (decimal)this->basis.cols());
				decimal tmp = 0.0;
				unknown(k) = tmp;
			}
			for (integer k = 3 * i + 2 * this->basis.cols(); k < 3 * i + 4 * this->basis.cols(); k++)
				//hard con(not ready)
			{
				//decimal tmp = 1.0 / (4.0 * (decimal)this->basis.cols());
				decimal tmp = 0.1;
				unknown(k) = tmp;
			}
		}
	}

	if (this->use_eg_or_not == true && this->use_basis_or_not == true)
	{
		if (this->rb_flag == false)
		{
			ReceiveBasis();
		}
		if (this->use_sim_or_not == false && this->use_ww_or_not == false)//1
		{
			unknown = Eigen::VectorXd::Zero(3 * mesh.n_vertices());
			l_bw_unknown = Eigen::VectorXd::Zero(this->basis.cols());
			d_bw_unknown = Eigen::VectorXd::Zero(this->basis.cols());

			int i = 0;
			for (auto it = this->mesh.vertices_begin(); it != this->mesh.vertices_end(); ++it, ++i)
			{
				auto point = this->mesh.point(*it);
				//cal moved pt
				if (i == cps.set[cnt].idx && cnt_flag)
				{
					int check_res = (CheckMovedPt(point, Eigen::RowVector3d(cps.set[cnt].x, cps.set[cnt].y, cps.set[cnt].z)));
					this->move_pt_num += check_res;
					this->pt_move_or_not.push_back(check_res);
					++cnt;
					if (cnt == cps.set.size())
					{
						cnt--;
						cnt_flag = false;
					}
				}

				//soft con
				unknown(3 * i + 0) = point.data()[0];
				unknown(3 * i + 1) = point.data()[1];
				unknown(3 * i + 2) = point.data()[2];

				//hard con
				/*if (idx_move[i] != -1)
				{
				unknown(3 * (i - idx_move[i]) + 0) = point.data()[0];
				unknown(3 * (i - idx_move[i]) + 1) = point.data()[1];
				unknown(3 * (i - idx_move[i]) + 2) = point.data()[2];
				}
				else
				{
				continue;
				}*/
			}

			//soft con
			for (integer k = 0; k < this->basis.cols(); k++)
				//hard con
				//for (integer k = 3 * (i - this->cps.set.size()); k < 3 * (i - this->cps.set.size()) + this->basis.cols(); k++)
			{
				decimal tmp = 1.0 / ((decimal)this->basis.cols());
				//decimal tmp = 0;
				l_bw_unknown(k) = tmp;
				d_bw_unknown(k) = tmp;
			}
		}
		if (this->use_sim_or_not == true && this->use_ww_or_not == false)//2
		{
			//soft con
			unknown = Eigen::VectorXd::Zero(3 * sim_mesh.n_vertices());
			l_bw_unknown = Eigen::VectorXd::Zero(this->basis.cols());
			d_bw_unknown = Eigen::VectorXd::Zero(this->basis.cols());
			high_v_mat = Eigen::MatrixXd::Zero(mesh.n_vertices(), 3);

			//hard con(not ready)

			int i = 0;
			for (auto it = this->sim_mesh.vertices_begin(); it != this->sim_mesh.vertices_end(); ++it, ++i)
			{
				auto point = this->sim_mesh.point(*it);
				//soft con
				unknown(3 * i + 0) = point.data()[0];
				unknown(3 * i + 1) = point.data()[1];
				unknown(3 * i + 2) = point.data()[2];

				//hard con(not ready)
			}
			int j = 0;
			for (auto it = this->mesh.vertices_begin(); it != this->mesh.vertices_end(); ++it, ++j)
			{
				auto point = this->mesh.point(*it);
				//cal moved pt
				if (i == cps.set[cnt].idx && cnt_flag)
				{
					int check_res = (CheckMovedPt(point, Eigen::RowVector3d(cps.set[cnt].x, cps.set[cnt].y, cps.set[cnt].z)));
					this->move_pt_num += check_res;
					this->pt_move_or_not.push_back(check_res);
					++cnt;
					if (cnt == cps.set.size())
					{
						cnt--;
						cnt_flag = false;
					}
				}

				//soft con
				high_v_mat(j, 0) = point.data()[0];
				high_v_mat(j, 1) = point.data()[1];
				high_v_mat(j, 2) = point.data()[2];

				//hard con
				/*if (idx_move[i] != -1)
				{
				unknown(3 * (i - idx_move[i]) + 0) = point.data()[0];
				unknown(3 * (i - idx_move[i]) + 1) = point.data()[1];
				unknown(3 * (i - idx_move[i]) + 2) = point.data()[2];
				}
				else
				{
				continue;
				}*/
			}
			//soft con
			for (integer k = 0; k < this->basis.cols(); k++)
				//hard con(not ready)
			{
				decimal tmp = 1.0 / (2 * (decimal)this->basis.cols());
				//decimal tmp = 0.0;
				l_bw_unknown(k) = tmp;
				d_bw_unknown(k) = tmp;
			}
		}
	}
	//cout << unknown << endl;
	pre_unknown = unknown;
}

void mns::MyDeformationHandler::ClearAfterDeform()
{
	cps.set.clear();
	triple.clear();
	ori_length_set.clear();
	ori_area_set.clear();
	ori_dihedral_set.clear();
	wnt_dihedral_set.clear();
	wnt_length_set.clear();
	idx_move.clear();
	con_ori_map.clear();
	sign_flag.clear();
	pt_move_or_not.clear();
}

void mns::MyDeformationHandler::EgInit()
{
	std::string eg_begin = this->example_str;
	std::string eg_end;
	std::string eg_idx0 = "1";
	std::string eg_idx1 = "0";
	std::string eg_idx2 = "0";
	std::string eg_idx3 = "0";
	if (model_mod == 1)
	{
		eg_end = ".ply";
	}
	if (model_mod == 2)
	{
		eg_end = ".off";
	}
	if (model_mod == 3)
	{
		eg_end = ".obj";
	}
	eg_mesh = new MyMesh[this->eg_num];
	for (integer k = 0; k < this->eg_num; k++)
	{
		std::string eg_idx;
		/*if (k<9)
		{
		eg_idx = eg_idx0;
		}*/
		if ((k + 1) % 10 == 0 && (k + 1) % 100 != 0)
		{
			eg_idx0 = "0";
			eg_idx1[0]++;
			//eg_idx = eg_idx1 + eg_idx0;
		}
		if ((k + 1) % 100 == 0 && (k + 1) % 1000 != 0)
		{
			eg_idx0 = "0";
			eg_idx1 = "0";
			eg_idx2[0]++;
			//eg_idx = eg_idx2 + eg_idx1 + eg_idx0;
		}
		eg_idx = eg_idx2 + eg_idx1 + eg_idx0;
		/*if ((k + 1) % 1000 == 0 && (k + 1) % 10000 != 0)
		{
		eg_idx0 = "0";
		eg_idx1 = "0";
		eg_idx2 = "0";
		eg_idx3[0]++;
		eg_idx = eg_idx3 + eg_idx2 + eg_idx1 + eg_idx0;
		}*/
		std::string eg_nm = eg_begin + eg_idx + eg_end;
		//MyMesh eg_mesh;
		bool res = OpenMesh::IO::read_mesh(eg_mesh[k], eg_nm);
		if (res == false)
		{
			std::cout << "read error\n";
		}

		dec_vector eg_length_set;
		dec_vector eg_dihedral_set;

		//length & dihedral init
		for (auto he_it = eg_mesh[k].halfedges_begin(); he_it != eg_mesh[k].halfedges_end(); ++he_it)
		{
			auto e = *he_it;
			auto fromVertex = eg_mesh[k].from_vertex_handle(e);
			auto toVertex = eg_mesh[k].to_vertex_handle(e);
			auto fromPt = eg_mesh[k].point(fromVertex);
			auto toPt = eg_mesh[k].point(toVertex);

			//stretching initialization		
			decimal l = eg_mesh[k].calc_edge_length(e);
			eg_length_set.push_back(l);

			//bend initialization
			auto f1 = eg_mesh[k].face_handle(e);
			auto op_he = eg_mesh[k].opposite_halfedge_handle(e);
			auto f2 = eg_mesh[k].face_handle(op_he);
			pt_vector pt_set;
			pt_vector t1_set;
			pt_vector t2_set;

			if (f1.idx() == -1 || f2.idx() == -1)
			{
				eg_dihedral_set.push_back(0);

				//edge, not half edge
				++he_it;
				continue;
			}

			pt_set.push_back(fromPt);
			t1_set.push_back(fromPt);
			t2_set.push_back(fromPt);
			pt_set.push_back(toPt);
			t1_set.push_back(toPt);
			t2_set.push_back(toPt);

			for (auto fv_it1 = eg_mesh[k].fv_begin(f1); fv_it1 != eg_mesh[k].fv_end(f1); ++fv_it1)
			{
				auto v = *fv_it1;
				auto pt = eg_mesh[k].point(v);
				if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
				{
					pt_set.push_back(pt);
					t1_set.push_back(pt);
					break;
				}
			}
			integer temp_cnt = 0;
			for (auto fv_it2 = eg_mesh[k].fv_begin(f2); fv_it2 != eg_mesh[k].fv_end(f2); ++fv_it2)
			{
				auto v = *fv_it2;
				auto pt = eg_mesh[k].point(v);
				if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
				{
					t2_set.push_back(pt);
					pt_set.push_back(pt);
					break;
				}
			}

			decimal angle = CalDiheral(pt_set[0], pt_set[1], pt_set[2], pt_set[3]);
			eg_dihedral_set.push_back(angle);

			//edge, not half edge
			++he_it;
		}

		this->eg_length_vec.push_back(eg_length_set);
		this->eg_dihedral_vec.push_back(eg_dihedral_set);

		decimal eg_volume = 0;
		if (alpha_v != 0)
		{
			for (auto f_it = eg_mesh[k].faces_begin(); f_it != eg_mesh[k].faces_end(); ++f_it)
			{
				auto face = *f_it;
				pt_vector pt_vec;
				for (auto fv_it = eg_mesh[k].fv_begin(face); fv_it != eg_mesh[k].fv_end(face); ++fv_it)
				{
					auto v = *fv_it;
					auto pt = eg_mesh[k].point(v);
					pt_vec.push_back(pt);
				}
				decimal temp_vol = (1.0 / 6.0) * ((pt_vec[0] % pt_vec[1]) | pt_vec[2]);
				eg_volume += temp_vol;
			}

			eg_volume_vec.push_back(eg_volume);
		}

		eg_idx0[0]++;
	}
	//delete[] eg_mesh;
}

void mns::MyDeformationHandler::BuildWntByBasis()
{
	Eigen::VectorXd weight = Eigen::VectorXd::Zero(this->basis.cols());

	std::ifstream w_ifile("weight.txt", std::ios::in);
	for (integer i = 0; i < this->basis.cols(); i++)
	{
		//weight(i) = 0.8;
		w_ifile >> weight(i);
	}
	w_ifile.close();
	Eigen::VectorXd x = basis*weight + data_avg;
	for (integer i = 0; i < this->mesh.n_edges(); i++)
	{
		wnt_length_set[i] = x(2 * i + 0);
		wnt_dihedral_set[i] = x(2 * i + 1);
	}
	if (this->alpha_v != 0)
	{
		wnt_volume = x(2 * this->mesh.n_edges());
	}
}

void mns::MyDeformationHandler::CalculateWeightForBasis(std::string nm)
{
	MyMesh input;
	bool res = OpenMesh::IO::read_mesh(input, nm);
	if (res == false)
	{
		std::cout << "read error\n";
	}
	dec_vector input_length_set;
	dec_vector input_dihedral_set;
	for (auto he_it = input.halfedges_begin(); he_it != input.halfedges_end(); ++he_it)
	{
		auto e = *he_it;
		auto fromVertex = input.from_vertex_handle(e);
		auto toVertex = input.to_vertex_handle(e);
		auto fromPt = input.point(fromVertex);
		auto toPt = input.point(toVertex);

		//stretching initialization		
		decimal l = input.calc_edge_length(e);
		input_length_set.push_back(l);

		//bend initialization
		auto f1 = input.face_handle(e);
		auto op_he = input.opposite_halfedge_handle(e);
		auto f2 = input.face_handle(op_he);
		pt_vector pt_set;
		pt_vector t1_set;
		pt_vector t2_set;

		if (f1.idx() == -1 || f2.idx() == -1)
		{
			input_dihedral_set.push_back(0);

			//edge, not half edge
			++he_it;
			continue;
		}

		pt_set.push_back(fromPt);
		t1_set.push_back(fromPt);
		t2_set.push_back(fromPt);
		pt_set.push_back(toPt);
		t1_set.push_back(toPt);
		t2_set.push_back(toPt);

		for (auto fv_it1 = input.fv_begin(f1); fv_it1 != input.fv_end(f1); ++fv_it1)
		{
			auto v = *fv_it1;
			auto pt = input.point(v);
			if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
			{
				pt_set.push_back(pt);
				t1_set.push_back(pt);
				break;
			}
		}
		integer temp_cnt = 0;
		for (auto fv_it2 = input.fv_begin(f2); fv_it2 != input.fv_end(f2); ++fv_it2)
		{
			auto v = *fv_it2;
			auto pt = input.point(v);
			if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
			{
				t2_set.push_back(pt);
				pt_set.push_back(pt);
				break;
			}
		}

		decimal angle = CalDiheral(pt_set[0], pt_set[1], pt_set[2], pt_set[3]);
		input_dihedral_set.push_back(angle);

		//edge, not half edge
		++he_it;
	}

	decimal input_volume = 0;
	if (this->alpha_v != 0)
	{
		for (auto f_it = input.faces_begin(); f_it != input.faces_end(); ++f_it)
		{
			auto face = *f_it;
			pt_vector pt_vec;
			for (auto fv_it = input.fv_begin(face); fv_it != input.fv_end(face); ++fv_it)
			{
				auto v = *fv_it;
				auto pt = input.point(v);
				pt_vec.push_back(pt);
			}
			decimal temp_vol = (1.0 / 6.0) * ((pt_vec[0] % pt_vec[1]) | pt_vec[2]);
			input_volume += temp_vol;
		}
	}

	Eigen::VectorXd len_di_vec;
	if (this->alpha_v == 0)
	{
		len_di_vec = Eigen::VectorXd::Zero(2 * input.n_edges());
	}
	else
	{
		len_di_vec = Eigen::VectorXd::Zero(2 * input.n_edges() + 1);
	}

	for (integer i = 0; i < input.n_edges(); i++)
	{
		len_di_vec(2 * i + 0) = input_length_set[i];
		len_di_vec(2 * i + 1) = input_dihedral_set[i];
	}
	if (this->alpha_v != 0)
	{
		len_di_vec(2 * input.n_edges()) = input_volume;
	}

	Eigen::VectorXd x = this->basis.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(len_di_vec);
	std::ofstream w_ofile("weight.txt", std::ios::out);
	w_ofile << x;
	w_ofile.close();
}

void mns::MyDeformationHandler::UpdateWntVec()
{
	if (use_eg_or_not == false)
	{
		return;
	}
	//length
	for (integer i = 0; i < wnt_length_set.size(); i++)
	{
		decimal tmp_length_diff = 0;
		if (this->use_eg_or_not == true)
		{
			for (integer k = 0; k < this->eg_num; k++)
			{
				//soft con
				tmp_length_diff += unknown(3 * this->mesh.n_vertices() + k) * (eg_length_vec[k][i] - ori_length_set[i]);

				//hard con
				//tmp_length_diff += unknown(3 * (this->mesh.n_vertices() - cps.set.size()) + k) * (eg_length_vec[k][i] - ori_length_set[i]);				
			}
		}

		wnt_length_set[i] = ori_length_set[i] + tmp_length_diff;
	}

	//dihedral
	for (integer i = 0; i < wnt_dihedral_set.size(); i++)
	{
		decimal tmp_dihedral_diff = 0;
		if (this->use_eg_or_not == true)
		{
			for (integer k = 0; k < this->eg_num; k++)
			{
				//soft con
				tmp_dihedral_diff += unknown(3 * this->mesh.n_vertices() + k) * (eg_dihedral_vec[k][i] - ori_dihedral_set[i]);

				//hard con
				//tmp_dihedral_diff += unknown(3 * (this->mesh.n_vertices() - cps.set.size()) + k) * (eg_dihedral_vec[k][i] - ori_dihedral_set[i]);
			}
		}

		wnt_dihedral_set[i] = ori_dihedral_set[i] + tmp_dihedral_diff;
	}

	//volume
	decimal tmp_volume_diff = 0;
	if (alpha_v != 0)
	{
		if (this->use_eg_or_not == true)
		{
			for (integer k = 0; k < this->eg_num; k++)
			{
				//soft con
				decimal alpha_k = unknown(3 * this->mesh.n_vertices() + k);

				//hard con
				//decimal alpha_k = unknown(3 * (this->mesh.n_vertices() - cps.set.size()) + k);

				tmp_volume_diff += alpha_k * (eg_volume_vec[k] - ori_volume);
			}
		}
		wnt_volume = ori_volume + tmp_volume_diff;
	}
}

void mns::MyDeformationHandler::UpdateWntVecWithBasis()
{
	if (use_basis_or_not == false)
	{
		return;
	}

	Eigen::VectorXd basis_weight = Eigen::VectorXd::Zero(this->basis.cols());

	for (integer i = 0; i < this->basis.cols(); i++)
	{
		//soft con
		basis_weight(i) = unknown(3 * this->mesh.n_vertices() + i);

		//hard con
		//basis_weight(i) = unknown(3 * (this->mesh.n_vertices() - cps.set.size()) + i);
	}
	//no origin
	//Eigen::VectorXd wanted_vec = this->basis*basis_weight + this->data_avg;

	//origin
	Eigen::VectorXd wanted_vec = this->basis*basis_weight + this->data_avg;

	for (integer i = 0; i < this->mesh.n_edges(); i++)
	{
		wnt_length_set[i] = wanted_vec(2 * i + 0);
		wnt_dihedral_set[i] = wanted_vec(2 * i + 1);
	}
	if (alpha_v != 0)
	{
		wnt_volume = wanted_vec(2 * this->mesh.n_edges());
	}
}
void mns::MyDeformationHandler::UpdateWntVecWithBasis7()
{
	if (use_basis_or_not == false)
	{
		return;
	}

	for (integer j = 0; j < this->mesh.n_edges(); j++)
	{
		decimal tmp_length_diff = 0;
		decimal tmp_dihedral_diff = 0;
		if (this->use_sim_or_not == false)
		{
#pragma omp parallel
			{
#pragma omp for reduction(+: tmp_length_diff, tmp_dihedral_diff)
				for (integer k = 0; k < basis.cols(); k++)
				{
					tmp_length_diff += unknown(3 * this->mesh.n_vertices() + 2 * k + 0) * length_basis_vec[k][j];
					tmp_dihedral_diff += unknown(3 * this->mesh.n_vertices() + 2 * k + 1) * dihedral_basis_vec[k][j];
				}
			}
		}
		else
		{
#pragma omp parallel
			{
#pragma omp for reduction(+: tmp_length_diff, tmp_dihedral_diff)
				for (integer k = 0; k < basis.cols(); k++)
				{
					tmp_length_diff += unknown(3 * this->sim_mesh.n_vertices() + 2 * k + 0) * length_basis_vec[k][j];
					tmp_dihedral_diff += unknown(3 * this->sim_mesh.n_vertices() + 2 * k + 1) * dihedral_basis_vec[k][j];
				}
			}
		}
		wnt_length_set[j] = (ori_length_set[j] + tmp_length_diff) * (1 + cm1(2 * j + 1));
		wnt_dihedral_set[j] =  ori_dihedral_set[j] + cm1( 2 * j + 0)+ tmp_dihedral_diff;
	}
	if (alpha_v != 0)
	{
		wnt_volume = ori_volume;
	}
}

void mns::MyDeformationHandler::UpdateWntVecWithBasis2()
{
	if (use_basis_or_not == false)
	{
		return;
	}

	for (integer j = 0; j < this->mesh.n_edges(); j++)
	{
		decimal tmp_length_diff = 0;
		decimal tmp_dihedral_diff = 0;
		if (this->use_sim_or_not == false)
		{
#pragma omp parallel
			{
#pragma omp for reduction(+: tmp_length_diff, tmp_dihedral_diff)
				for (integer k = 0; k < basis.cols(); k++)
				{
					tmp_length_diff += unknown(3 * this->mesh.n_vertices() + 2 * k + 0) * length_basis_vec[k][j];
					tmp_dihedral_diff += unknown(3 * this->mesh.n_vertices() + 2 * k + 1) * dihedral_basis_vec[k][j];
				}
			}
		}
		else
		{
#pragma omp parallel
			{
#pragma omp for reduction(+: tmp_length_diff, tmp_dihedral_diff)
				for (integer k = 0; k < basis.cols(); k++)
				{
					tmp_length_diff += unknown(3 * this->sim_mesh.n_vertices() + 2 * k + 0) * length_basis_vec[k][j];
					tmp_dihedral_diff += unknown(3 * this->sim_mesh.n_vertices() + 2 * k + 1) * dihedral_basis_vec[k][j];
				}
			}
		}
		wnt_length_set[j] = ori_length_set[j] * (1 + tmp_length_diff);
		//wnt_dihedral_set[j] = 3.1415926 - ( (3.1415926 -  ori_dihedral_set[j]) + tmp_dihedral_diff);
		wnt_dihedral_set[j] =ori_dihedral_set[j] + tmp_dihedral_diff;
	}
	if (alpha_v != 0)
	{
		wnt_volume = ori_volume;
	}
}

void mns::MyDeformationHandler::UpdateWntVecWithBasis3(Eigen::VectorXd cm)
{
	for (integer j = 0; j < this->mesh.n_edges(); j++)
	{
		wnt_length_set[j] = cm(2 * j + 1);
		wnt_dihedral_set[j] = cm(2 * j + 0);
	}
}

void mns::MyDeformationHandler::UpdateWntVecWithBasis3()
{
	if (use_basis_or_not == false)
	{
		return;
	}

	for (integer j = 0; j < this->mesh.n_edges(); j++)
	{
		decimal tmp_length_diff = 0;
		decimal tmp_dihedral_diff = 0;
		if (this->use_sim_or_not == false)
		{
#pragma omp parallel
			{
#pragma omp for reduction(+: tmp_length_diff, tmp_dihedral_diff)
				for (integer k = 0; k < basis.cols(); k++)
				{
					tmp_length_diff += l_bw_unknown(k) * length_basis_vec[k][j];
					tmp_dihedral_diff += d_bw_unknown(k) * dihedral_basis_vec[k][j];
				}
			}
		}
		else
		{
#pragma omp parallel
			{
#pragma omp for reduction(+: tmp_length_diff, tmp_dihedral_diff)
				for (integer k = 0; k < basis.cols(); k++)
				{
					tmp_length_diff += l_bw_unknown(k) * length_basis_vec[k][j];
					tmp_dihedral_diff += d_bw_unknown(k) * dihedral_basis_vec[k][j];
				}
			}
		}
		wnt_length_set[j] = ori_length_set[j] * (1 + tmp_length_diff);
		wnt_dihedral_set[j] = 3.1415926 - ((3.1415926 - ori_dihedral_set[j]) + tmp_dihedral_diff);


	}
	if (alpha_v != 0)
	{
		wnt_volume = ori_volume;
	}
}

void mns::MyDeformationHandler::BuildWntVecByBasis()
{
	ReceiveBasis();
	Eigen::VectorXd basis_weight = Eigen::VectorXd::Zero(this->basis.cols());
	std::ifstream w_ifile("w_unknown.txt", std::ios::in);
	for (integer i = 0; i < this->basis.cols(); i++)
	{
		decimal tmp;
		w_ifile >> tmp;
		basis_weight(i) = tmp;
	}
	w_ifile.close();

	Eigen::VectorXd wanted_vec = this->basis*basis_weight + this->data_avg;

	for (integer i = 0; i < this->mesh.n_edges(); i++)
	{
		wnt_length_set[i] = wanted_vec(2 * i + 0);
		wnt_dihedral_set[i] = wanted_vec(2 * i + 1);
		if (wanted_vec(2 * i + 1) > M_PI)
		{
			wnt_dihedral_set[i] = M_PI;
		}
		if (wanted_vec(2 * i + 1) < 0)
		{
			wnt_dihedral_set[i] = 0;
		}
	}
	if (alpha_v != 0)
	{
		wnt_volume = wanted_vec(2 * this->mesh.n_edges());
	}
}

void mns::MyDeformationHandler::BuildWntVecByBasis2()
{
	if (this->use_basis_or_not == false)
	{
		return;
	}
	std::ifstream w_ifile("setting_weight.txt", std::ios::in);
	int_vector w_idx_vec;
	dec_vector w_vec;
	while (w_ifile.peek() != EOF)
	{
		integer int_tmp;
		w_ifile >> int_tmp;
		w_idx_vec.push_back(int_tmp);
		decimal dec_tmp;
		w_ifile >> dec_tmp;
		w_vec.push_back(dec_tmp);
	}
	/*w_idx_vec.erase(w_idx_vec.begin() + w_idx_vec.size() - 1);
	w_vec.erase(w_vec.begin() + w_vec.size() - 1);*/
	w_ifile.close();

	Eigen::VectorXd wanted_vec = this->data_avg;

	for (integer i = 0; i < w_idx_vec.size(); i++)
	{
		if (w_idx_vec[i] % 2 == 0)
		{
			for (integer j = 0; j < this->mesh.n_edges(); j++)
			{
				wanted_vec[2 * j + 0] += length_basis_vec[w_idx_vec[i] / 2][j] * w_vec[i];
			}
		}
		if (w_idx_vec[i] % 2 == 1)
		{
			for (integer j = 0; j < this->mesh.n_edges(); j++)
			{
				wanted_vec[2 * j + 1] += dihedral_basis_vec[w_idx_vec[i] / 2][j] * w_vec[i];
				if (wanted_vec[2 * j + 1] > M_PI)
				{
					wanted_vec[2 * j + 1] = M_PI;
				}
				if (wanted_vec[2 * j + 1] < 0)
				{
					wanted_vec[2 * j + 1] = 0;
				}
			}
		}
	}

	std::ofstream dl_ofile("DiLenV.txt", std::ios::out);
	for (integer i = 0; i < this->mesh.n_edges(); i++)
	{
		dl_ofile << sign_flag[i] * (M_PI - wanted_vec[2 * i + 1]) << std::endl;
		dl_ofile << wanted_vec[2 * i + 0] << std::endl;
	}
	dl_ofile.close();
}

void mns::MyDeformationHandler::Transform2Example()
{
	std::ifstream nmif("test_nm.txt", std::ios::in);
	std::string ex_nm;
	nmif >> ex_nm;
	nmif.close();
	bool res = OpenMesh::IO::read_mesh(ex_mesh, ex_nm);
	if (res == false)
	{
		std::cout << "read error\n";
	}

	dec_vector eg_length_set;
	dec_vector eg_dihedral_set;

	//length & dihedral init
	for (auto he_it = ex_mesh.halfedges_begin(); he_it != ex_mesh.halfedges_end(); ++he_it)
	{
		auto e = *he_it;
		auto fromVertex = ex_mesh.from_vertex_handle(e);
		auto toVertex = ex_mesh.to_vertex_handle(e);
		auto fromPt = ex_mesh.point(fromVertex);
		auto toPt = ex_mesh.point(toVertex);

		//stretching initialization		
		decimal l = ex_mesh.calc_edge_length(e);
		eg_length_set.push_back(l);

		//bend initialization
		auto f1 = ex_mesh.face_handle(e);
		auto op_he = ex_mesh.opposite_halfedge_handle(e);
		auto f2 = ex_mesh.face_handle(op_he);
		pt_vector pt_set;
		pt_vector t1_set;
		pt_vector t2_set;
		if (f1.idx() == -1 || f2.idx() == -1)
		{
			eg_dihedral_set.push_back(0);

			//edge, not half edge
			++he_it;
			continue;
		}

		pt_set.push_back(fromPt);
		t1_set.push_back(fromPt);
		t2_set.push_back(fromPt);
		pt_set.push_back(toPt);
		t1_set.push_back(toPt);
		t2_set.push_back(toPt);

		for (auto fv_it1 = ex_mesh.fv_begin(f1); fv_it1 != ex_mesh.fv_end(f1); ++fv_it1)
		{
			auto v = *fv_it1;
			auto pt = ex_mesh.point(v);
			if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
			{
				pt_set.push_back(pt);
				t1_set.push_back(pt);
				break;
			}
		}
		integer temp_cnt = 0;
		for (auto fv_it2 = ex_mesh.fv_begin(f2); fv_it2 != ex_mesh.fv_end(f2); ++fv_it2)
		{
			auto v = *fv_it2;
			auto pt = ex_mesh.point(v);
			if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
			{
				t2_set.push_back(pt);
				pt_set.push_back(pt);
				break;
			}
		}

		decimal angle = CalDiheral(pt_set[0], pt_set[1], pt_set[2], pt_set[3]);
		eg_dihedral_set.push_back(angle);

		//edge, not half edge
		++he_it;
	}

	if (alpha_v != 0)
	{
		decimal eg_volume = 0;
		for (auto f_it = ex_mesh.faces_begin(); f_it != ex_mesh.faces_end(); ++f_it)
		{
			auto face = *f_it;
			pt_vector pt_vec;
			for (auto fv_it = ex_mesh.fv_begin(face); fv_it != ex_mesh.fv_end(face); ++fv_it)
			{
				auto v = *fv_it;
				auto pt = ex_mesh.point(v);
				pt_vec.push_back(pt);
			}
			decimal temp_vol = (1.0 / 6.0) * ((pt_vec[0] % pt_vec[1]) | pt_vec[2]);
			eg_volume += temp_vol;
		}

		wnt_volume = eg_volume;
	}

	for (integer i = 0; i < this->mesh.n_edges(); i++)
	{
		wnt_length_set[i] = eg_length_set[i];
		wnt_dihedral_set[i] = eg_dihedral_set[i];
	}
}

void mns::MyDeformationHandler::PrintLengthDihedral()
{
	dec_vector length_set;
	dec_vector my_dihedral_set;
	dec_vector other_dihedral_set;

	for (auto he_it = this->mesh.halfedges_begin(); he_it != this->mesh.halfedges_end(); ++he_it)
	{
		auto e = *he_it;
		auto fromVertex = this->mesh.from_vertex_handle(e);
		auto toVertex = this->mesh.to_vertex_handle(e);
		auto fromPt = this->mesh.point(fromVertex);
		auto toPt = this->mesh.point(toVertex);

		//stretching initialization		
		decimal l = this->mesh.calc_edge_length(e);
		length_set.push_back(l);

		//bend initialization
		auto f1 = this->mesh.face_handle(e);
		auto op_he = this->mesh.opposite_halfedge_handle(e);
		auto f2 = this->mesh.face_handle(op_he);
		pt_vector pt_set;
		if (f1.idx() == -1 || f2.idx() == -1)
		{
			my_dihedral_set.push_back(0);

			sign_flag.push_back(1);
			other_dihedral_set.push_back(0);
			//edge, not half edge
			++he_it;
			continue;
		}

		pt_set.push_back(fromPt);
		pt_set.push_back(toPt);

		for (auto fv_it1 = this->mesh.fv_begin(f1); fv_it1 != this->mesh.fv_end(f1); ++fv_it1)
		{
			auto v = *fv_it1;
			auto pt = this->mesh.point(v);
			if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
			{
				pt_set.push_back(pt);
				break;
			}
		}
		integer temp_cnt = 0;
		for (auto fv_it2 = this->mesh.fv_begin(f2); fv_it2 != this->mesh.fv_end(f2); ++fv_it2)
		{
			auto v = *fv_it2;
			auto pt = this->mesh.point(v);
			if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
			{
				pt_set.push_back(pt);
				break;
			}
		}

		decimal angle = CalDiheral(pt_set[0], pt_set[1], pt_set[2], pt_set[3]);
		my_dihedral_set.push_back(angle);

		decimal phi;
		auto normal1_tmp = this->mesh.normal(f1);
		Eigen::Vector3d normal1;
		normal1 << normal1_tmp.data()[0], normal1_tmp.data()[1], normal1_tmp.data()[2];

		auto normal2_tmp = this->mesh.normal(f2);
		Eigen::Vector3d normal2;
		normal2 << normal2_tmp.data()[0], normal2_tmp.data()[1], normal2_tmp.data()[2];

		auto e_ij_tmp = toPt - fromPt;
		Eigen::Vector3d e_ij;
		e_ij << e_ij_tmp.data()[0], e_ij_tmp.data()[1], e_ij_tmp.data()[2];
		e_ij = e_ij / (e_ij.norm());

		Eigen::Vector3d crossnormal = normal1.cross(normal2);
		decimal dotcross = e_ij.dot(crossnormal);
		if (dotcross >= 0)
		{
			sign_flag.push_back(1);
			double dotNorm1Norm2 = normal1.dot(normal2);
			if (dotNorm1Norm2 < 1 - 1e-6 && dotNorm1Norm2 > -1 + 1e-6)
			{
				phi = std::acos(dotNorm1Norm2);
			}
			else if (dotNorm1Norm2 >= 1 - 1e-6)
			{
				phi = 0;
			}
			else
			{
				phi = 3.1415926;
			}
		}
		else
		{
			sign_flag.push_back(-1);
			double dotNorm1Norm2 = normal1.dot(normal2);
			if (dotNorm1Norm2 < 1 - 1e-6 && dotNorm1Norm2 > -1 + 1e-6)
			{
				phi = -std::acos(dotNorm1Norm2);
			}
			else if (dotNorm1Norm2 >= 1 - 1e-6)
			{
				phi = 0;
			}
			else
			{
				phi = -3.1415926;
			}
		}
		other_dihedral_set.push_back(phi);

		//edge, not half edge
		++he_it;
	}

	std::ofstream my_ofile("my_dile.txt", std::ios::out);
	std::ofstream other_ofile("other_dile.txt", std::ios::out);
	std::ofstream change_ofile("change_dile.txt", std::ios::out);
	for (integer i = 0; i < length_set.size(); i++)
	{
		//my_ofile << length_set[i] << std::endl;
		my_ofile << my_dihedral_set[i] << std::endl << length_set[i] << std::endl;
		other_ofile << other_dihedral_set[i] << std::endl << length_set[i] << std::endl;
		//change_ofile << sign_flag[i] * (3.14159 - my_dihedral_set[i]) << std::endl;
	}
	my_ofile.close();
	other_ofile.close();
	change_ofile.close();
}

void mns::MyDeformationHandler::SetModelMod(integer mod)
{
	this->model_mod = mod;
}

mns::MyPointSet mns::MyDeformationHandler::GetConSetByV(const std::vector<integer>& s_vid_vec, const Eigen::MatrixXd& v_mat)
{
	mns::MyPointSet con_set;

	for (integer i = 0; i < s_vid_vec.size(); i++)
	{
		mns::MyPoint p;
		p.idx = s_vid_vec[i];
		p.x = v_mat(s_vid_vec[i], 0);
		p.y = v_mat(s_vid_vec[i], 1);
		p.z = v_mat(s_vid_vec[i], 2);

		con_set.set.push_back(p);
	}

	return con_set;
}

mns::MyPointSet mns::MyDeformationHandler::GetConSetByInput()
{
	mns::MyPointSet con_set;
	mns::MyPoint p;
	integer int_tmp;
	decimal dec_tmp;
	std::set<integer> idx_set;
	std::vector<integer> idx_vec;

	std::ifstream ifile("conset.txt", std::ios::in);
	while (ifile >> int_tmp)
	{
		if (int_tmp < 0)
		{
			ifile.close();
			return con_set;
		}
		else
		{
			std::set<integer>::iterator idx_set_it = idx_set.find(int_tmp);
			if (idx_set_it != idx_set.end())
			{
				for (integer i = 0; i < con_set.set.size(); i++)
				{
					if (con_set.set[i].idx == int_tmp)
					{
						ifile >> dec_tmp;
						con_set.set[i].x = dec_tmp;
						ifile >> dec_tmp;
						con_set.set[i].y = dec_tmp;
						ifile >> dec_tmp;
						con_set.set[i].z = dec_tmp;
						break;
					}
				}
			}
			else
			{
				idx_set.insert(int_tmp);
				p.idx = int_tmp;
				ifile >> dec_tmp;
				p.x = dec_tmp;
				ifile >> dec_tmp;
				p.y = dec_tmp;
				ifile >> dec_tmp;
				p.z = dec_tmp;
				con_set.set.push_back(p);
			}
		}
	}
}

mns::MyPointSet mns::MyDeformationHandler::GetConSetByInput(string  p_path)
{
	mns::MyPointSet con_set;
	mns::MyPoint p;
	integer int_tmp;
	decimal dec_tmp;
	std::set<integer> idx_set;
	std::vector<integer> idx_vec;

	std::ifstream ifile(p_path, std::ios::in);
	while (ifile >> int_tmp)
	{
		if (int_tmp < 0)
		{
			ifile.close();
			return con_set;
		}
		else
		{
			std::set<integer>::iterator idx_set_it = idx_set.find(int_tmp);
			if (idx_set_it != idx_set.end())
			{
				for (integer i = 0; i < con_set.set.size(); i++)
				{
					if (con_set.set[i].idx == int_tmp)
					{
						ifile >> dec_tmp;
						con_set.set[i].x = dec_tmp;
						ifile >> dec_tmp;
						con_set.set[i].y = dec_tmp;
						ifile >> dec_tmp;
						con_set.set[i].z = dec_tmp;
						break;
					}
				}
			}
			else
			{
				idx_set.insert(int_tmp);
				p.idx = int_tmp;
				ifile >> dec_tmp;
				p.x = dec_tmp;
				ifile >> dec_tmp;
				p.y = dec_tmp;
				ifile >> dec_tmp;
				p.z = dec_tmp;
				con_set.set.push_back(p);
			}
		}
	}
	return con_set;
}

mns::MyPointSet mns::MyDeformationHandler::GetConSetByDif()
{
	std::ifstream ifile("nm.txt", std::ios::in);
	std::string nm;
	ifile >> nm;
	std::string test_nm;
	ifile >> test_nm;
	ifile.close();

	MyMesh target;
	bool res = OpenMesh::IO::read_mesh(target, test_nm);
	if (res == false)
	{
		std::cout << "read error\n";
	}

	mns::MyPointSet con_set;
	mns::MyPoint p;
	auto t_it = target.vertices_begin();
	integer i = 0;
	for (auto it = this->ori_mesh.vertices_begin(); it != this->ori_mesh.vertices_end(); ++it, ++t_it, ++i)
	{
		auto pt = this->ori_mesh.point(*it);
		auto t_pt = target.point(*t_it);
		auto p_dif = pt - t_pt;
		decimal dif = p_dif.length();
		if (dif > 0.1)
		{
			p.idx = i;
			p.x = t_pt.data()[0];
			p.y = t_pt.data()[1];
			p.z = t_pt.data()[2];
			con_set.set.push_back(p);
		}
	}
	return con_set;
}

mns::MyPointSet mns::MyDeformationHandler::GetConSetByMesh()
{
	std::ifstream ifile("conset_nm_test.txt", std::ios::in);
	integer choice;
	ifile >> choice;
	std::string tar_nm;
	ifile >> tar_nm;
	decimal sam_per;
	ifile >> sam_per;
	ifile.close();

	integer NV = this->mesh.n_vertices();
	integer sam_interval = (integer)((decimal)NV*sam_per);
	sam_interval = (sam_interval > 1) ? sam_interval : 2;
	std::cout << sam_interval << std::endl;

	MyMesh target;
	bool res = OpenMesh::IO::read_mesh(target, tar_nm);
	if (res == false)
	{
		std::cout << "read error\n";
	}

	mns::MyPointSet con_set;
	integer i = 0;
	for (auto it = target.vertices_begin(); it != target.vertices_end(); ++it, ++i)
	{
		if (i % sam_interval == 0)
		{
			auto pt = target.point(*it);
			mns::MyPoint p;
			p.idx = i;
			p.x = pt.data()[0];
			p.y = pt.data()[1];
			p.z = pt.data()[2];
			con_set.set.push_back(p);
		}
	}
	target.clear();
	return con_set;
}

mns::MyPointSet mns::MyDeformationHandler::GetConSetByMesh(std::string tar_nm, decimal sam_per)
{
	integer NV = this->mesh.n_vertices();
	integer sam_interval = (integer)((decimal)NV*sam_per);
	sam_interval = (sam_interval > 1) ? sam_interval : 2;

	MyMesh target;
	bool res = OpenMesh::IO::read_mesh(target, tar_nm);
	if (res == false)
	{
		std::cout << "read error\n";
	}

	mns::MyPointSet con_set;
	integer i = 0;
	for (auto it = target.vertices_begin(); it != target.vertices_end(); ++it, ++i)
	{
		if (i % sam_interval == 0)
		{
			auto pt = target.point(*it);
			mns::MyPoint p;
			p.idx = i;
			p.x = pt.data()[0];
			p.y = pt.data()[1];
			p.z = pt.data()[2];
			con_set.set.push_back(p);
		}
	}
	target.clear();
	return con_set;
}

bool setCmp(const mns::MyPoint& a, const mns::MyPoint& b)
{
	if (a.idx > b.idx)
	{
		return false;
	}
	else
	{
		return true;
	}
}

void mns::MyDeformationHandler::CreateConstraintPointSet1(const MyPointSet& con)
{
	integer sz = con.set.size();
	integer lz = cps.set.size();
	cps.set.resize(sz + lz);
	std::copy(con.set.begin(), con.set.end(), &(cps.set[lz]));
	pt_move_or_not.clear();
	pt_move_or_not.resize(sz + lz);
	for (int i = 0; i < sz + lz; i++)
	{
		if (i < sz)
			pt_move_or_not[i] = 0;
		else
			pt_move_or_not[i] = 1;
	}
	this->move_pt_num = sz;
	//std::sort(cps.set.begin(), cps.set.end(), setCmp);
	return;
}

void mns::MyDeformationHandler::CreateConstraintPointSet(const MyPointSet& con)
{
	integer sz = con.set.size();
	cps.set.resize(sz);
	std::copy(con.set.begin(), con.set.end(), cps.set.begin());
	std::sort(cps.set.begin(), cps.set.end(), setCmp);
	pt_move_or_not.clear();
	pt_move_or_not.resize(sz);
	for (int i = 0; i < sz ; i++)
	{
		
			pt_move_or_not[i] = 1;
	}
	this->move_pt_num = sz;
	return;
}

bool mns::MyDeformationHandler::CheckMovedPt(MyMesh::Point m_pt, Eigen::RowVector3d v_pt)
{
	decimal eps = 0.0001;
	Eigen::RowVector3d m2v_pt(m_pt.data()[0], m_pt.data()[1], m_pt.data()[2]);
	decimal err = (m2v_pt - v_pt).norm();
	return (err > eps);
}

decimal mns::MyDeformationHandler::GetConTerm(integer id)
{
	decimal energy = 0;
	if (this->use_sim_or_not == false)
	{
		fc(3 * con_ele_cnt + 0) = swc*(unknown(3 * cps.set[id].idx + 0) - cps.set[id].x);
		fc(3 * con_ele_cnt + 1) = swc*(unknown(3 * cps.set[id].idx + 1) - cps.set[id].y);
		fc(3 * con_ele_cnt + 2) = swc*(unknown(3 * cps.set[id].idx + 2) - cps.set[id].z);

		triple.push_back(Eigen::Triplet<decimal>(3 * con_ele_cnt + 0, 3 * cps.set[id].idx + 0, swc));
		triple.push_back(Eigen::Triplet<decimal>(3 * con_ele_cnt + 1, 3 * cps.set[id].idx + 1, swc));
		triple.push_back(Eigen::Triplet<decimal>(3 * con_ele_cnt + 2, 3 * cps.set[id].idx + 2, swc));
	}

	if (this->use_sim_or_not == true)
	{
		/*for (integer j = 0; j < this->W_idx_set[cps.set[id].idx].size(); j++)
		{
		decimal wij = W_mat(cps.set[id].idx, W_idx_set[cps.set[id].idx][j]);
		triple.push_back(Eigen::Triplet<decimal>(3 * con_ele_cnt + 0, 3 * W_idx_set[cps.set[id].idx][j] + 0, swc * wij));
		triple.push_back(Eigen::Triplet<decimal>(3 * con_ele_cnt + 1, 3 * W_idx_set[cps.set[id].idx][j] + 1, swc * wij));
		triple.push_back(Eigen::Triplet<decimal>(3 * con_ele_cnt + 2, 3 * W_idx_set[cps.set[id].idx][j] + 2, swc * wij));
		}*/
		fc(3 * con_ele_cnt + 0) = swc*(high_v_mat(cps.set[id].idx, 0) - cps.set[id].x);
		fc(3 * con_ele_cnt + 1) = swc*(high_v_mat(cps.set[id].idx, 1) - cps.set[id].y);
		fc(3 * con_ele_cnt + 2) = swc*(high_v_mat(cps.set[id].idx, 2) - cps.set[id].z);

		j_triple.push_back(Eigen::Triplet<decimal>(3 * con_ele_cnt + 0, 3 * cps.set[id].idx + 0, swc));
		j_triple.push_back(Eigen::Triplet<decimal>(3 * con_ele_cnt + 1, 3 * cps.set[id].idx + 1, swc));
		j_triple.push_back(Eigen::Triplet<decimal>(3 * con_ele_cnt + 2, 3 * cps.set[id].idx + 2, swc));
	}

	energy = fc(3 * con_ele_cnt + 0)*fc(3 * con_ele_cnt + 0) + fc(3 * con_ele_cnt + 1)*fc(3 * con_ele_cnt + 1) + fc(3 * con_ele_cnt + 2)*fc(3 * con_ele_cnt + 2);
	++con_ele_cnt;

	return energy;
}

decimal mns::MyDeformationHandler::GetFixConTermNew(integer id, decimal fac)
{
	decimal energy = 0;
	decimal new_swc = swc * fac;
	fc(3 * con_ele_cnt + 0) = 0;
	fc(3 * con_ele_cnt + 1) = 0;
	fc(3 * con_ele_cnt + 2) = 0;
	for (int i = 0; i < fixed_feature_idx[id].size(); i++)
	{
		fc(3 * con_ele_cnt + 0) += fixed_feature_weights[id][i] * (unknown(3 * fixed_feature_idx[id][i] + 0));
		fc(3 * con_ele_cnt + 1) += fixed_feature_weights[id][i] * (unknown(3 * fixed_feature_idx[id][i] + 1));
		fc(3 * con_ele_cnt + 2) += fixed_feature_weights[id][i] * (unknown(3 * fixed_feature_idx[id][i] + 2));
		triple.push_back(Eigen::Triplet<decimal>(3 * con_ele_cnt + 0, 3 * fixed_feature_idx[id][i] + 0, new_swc * fixed_feature_weights[id][i]));
		triple.push_back(Eigen::Triplet<decimal>(3 * con_ele_cnt + 1, 3 * fixed_feature_idx[id][i] + 1, new_swc * fixed_feature_weights[id][i]));
		triple.push_back(Eigen::Triplet<decimal>(3 * con_ele_cnt + 2, 3 * fixed_feature_idx[id][i] + 2, new_swc * fixed_feature_weights[id][i]));
	}
	fc(3 * con_ele_cnt + 0) = new_swc * (fc(3 * con_ele_cnt + 0) - fixed_feature_pos[id](0));
	fc(3 * con_ele_cnt + 1) = new_swc * (fc(3 * con_ele_cnt + 1) - fixed_feature_pos[id](1));
	fc(3 * con_ele_cnt + 2) = new_swc * (fc(3 * con_ele_cnt + 2) - fixed_feature_pos[id](2));

	energy = fc(3 * con_ele_cnt + 0)*fc(3 * con_ele_cnt + 0) + fc(3 * con_ele_cnt + 1)*fc(3 * con_ele_cnt + 1) + fc(3 * con_ele_cnt + 2)*fc(3 * con_ele_cnt + 2);
	++con_ele_cnt;
	return energy;
}

decimal mns::MyDeformationHandler::GetConTermNew(integer id, decimal fac)
{
	decimal energy = 0;
	decimal new_swc = swc * fac;
	fc(3 * con_ele_cnt + 0) = 0;
	fc(3 * con_ele_cnt + 1) = 0;
	fc(3 * con_ele_cnt + 2) = 0;
	for (int i = 0; i < body_feature_idx[id].size(); i++)
	{
		fc(3 * con_ele_cnt + 0) += body_feature_weights[id][i] * (unknown(3 * body_feature_idx[id][i] + 0));
		fc(3 * con_ele_cnt + 1) += body_feature_weights[id][i] * (unknown(3 * body_feature_idx[id][i] + 1));
		fc(3 * con_ele_cnt + 2) += body_feature_weights[id][i] * (unknown(3 * body_feature_idx[id][i] + 2));
		triple.push_back(Eigen::Triplet<decimal>(3 * con_ele_cnt + 0, 3 * body_feature_idx[id][i] + 0, new_swc * body_feature_weights[id][i]));
		triple.push_back(Eigen::Triplet<decimal>(3 * con_ele_cnt + 1, 3 * body_feature_idx[id][i] + 1, new_swc * body_feature_weights[id][i]));
		triple.push_back(Eigen::Triplet<decimal>(3 * con_ele_cnt + 2, 3 * body_feature_idx[id][i] + 2, new_swc * body_feature_weights[id][i]));
	}
	fc(3 * con_ele_cnt + 0) = new_swc * (fc(3 * con_ele_cnt + 0) - body_feature_pos[id](0));
	fc(3 * con_ele_cnt + 1) = new_swc * (fc(3 * con_ele_cnt + 1) - body_feature_pos[id](1));
	fc(3 * con_ele_cnt + 2) = new_swc * (fc(3 * con_ele_cnt + 2) - body_feature_pos[id](2));

	energy = fc(3 * con_ele_cnt + 0)*fc(3 * con_ele_cnt + 0) + fc(3 * con_ele_cnt + 1)*fc(3 * con_ele_cnt + 1) + fc(3 * con_ele_cnt + 2)*fc(3 * con_ele_cnt + 2);
	++con_ele_cnt;
	return energy;
}



decimal mns::MyDeformationHandler::GetConTerm(integer id, decimal fac)
{
	decimal energy = 0;
	decimal new_swc = swc * fac;
	if (this->use_sim_or_not == false)
	{
		fc(3 * con_ele_cnt + 0) = new_swc*(unknown(3 * cps.set[id].idx + 0) - cps.set[id].x);
		fc(3 * con_ele_cnt + 1) = new_swc*(unknown(3 * cps.set[id].idx + 1) - cps.set[id].y);
		fc(3 * con_ele_cnt + 2) = new_swc*(unknown(3 * cps.set[id].idx + 2) - cps.set[id].z);

		triple.push_back(Eigen::Triplet<decimal>(3 * con_ele_cnt + 0, 3 * cps.set[id].idx + 0, new_swc));
		triple.push_back(Eigen::Triplet<decimal>(3 * con_ele_cnt + 1, 3 * cps.set[id].idx + 1, new_swc));
		triple.push_back(Eigen::Triplet<decimal>(3 * con_ele_cnt + 2, 3 * cps.set[id].idx + 2, new_swc));
	}

	if (this->use_sim_or_not == true)
	{
		/*for (integer j = 0; j < this->W_idx_set[cps.set[id].idx].size(); j++)
		{
		decimal wij = W_mat(cps.set[id].idx, W_idx_set[cps.set[id].idx][j]);
		triple.push_back(Eigen::Triplet<decimal>(3 * con_ele_cnt + 0, 3 * W_idx_set[cps.set[id].idx][j] + 0, swc * wij));
		triple.push_back(Eigen::Triplet<decimal>(3 * con_ele_cnt + 1, 3 * W_idx_set[cps.set[id].idx][j] + 1, swc * wij));
		triple.push_back(Eigen::Triplet<decimal>(3 * con_ele_cnt + 2, 3 * W_idx_set[cps.set[id].idx][j] + 2, swc * wij));
		}*/
		fc(3 * con_ele_cnt + 0) = new_swc*(high_v_mat(cps.set[id].idx, 0) - cps.set[id].x);
		fc(3 * con_ele_cnt + 1) = new_swc*(high_v_mat(cps.set[id].idx, 1) - cps.set[id].y);
		fc(3 * con_ele_cnt + 2) = new_swc*(high_v_mat(cps.set[id].idx, 2) - cps.set[id].z);

		j_triple.push_back(Eigen::Triplet<decimal>(3 * con_ele_cnt + 0, 3 * cps.set[id].idx + 0, new_swc));
		j_triple.push_back(Eigen::Triplet<decimal>(3 * con_ele_cnt + 1, 3 * cps.set[id].idx + 1, new_swc));
		j_triple.push_back(Eigen::Triplet<decimal>(3 * con_ele_cnt + 2, 3 * cps.set[id].idx + 2, new_swc));
	}

	energy = fc(3 * con_ele_cnt + 0)*fc(3 * con_ele_cnt + 0) + fc(3 * con_ele_cnt + 1)*fc(3 * con_ele_cnt + 1) + fc(3 * con_ele_cnt + 2)*fc(3 * con_ele_cnt + 2);
	++con_ele_cnt;

	return energy;
}

decimal mns::MyDeformationHandler::EcNew()
{
	integer totalpoints = body_feature_idx.size() + fixed_feature_idx.size();
	integer movepoints = body_feature_idx.size();
	if (fc.rows() != 3 * (body_feature_idx.size() + fixed_feature_idx.size()));
	{
		fc.resize(3 * body_feature_idx.size()+ 3 * fixed_feature_idx.size());
	}
	fc = Eigen::VectorXd::Zero(3 * body_feature_idx.size() + 3 * fixed_feature_idx.size());
	decimal total_energy = 0;
	for (int i = 0; i < body_feature_idx.size(); i++)
	{
		//decimal fac = (pt_move_or_not[i] == 1) ? (std::sqrt(((decimal)cps.set.size()) / ((decimal)move_pt_num))) : (std::sqrt(((decimal)cps.set.size()) / ((decimal)(cps.set.size() - move_pt_num))));
		total_energy += GetConTermNew(i, std::sqrt(((decimal)totalpoints) / ((decimal)movepoints)));
	}
	for (int i = 0; i < fixed_feature_idx.size(); i++)
	{
		//decimal fac = (pt_move_or_not[i] == 1) ? (std::sqrt(((decimal)cps.set.size()) / ((decimal)move_pt_num))) : (std::sqrt(((decimal)cps.set.size()) / ((decimal)(cps.set.size() - move_pt_num))));
		total_energy += GetFixConTermNew(i, 5.0 * std::sqrt(((decimal)totalpoints) / ((decimal)(totalpoints - movepoints))));
	}
	row_cnt += 3 * con_ele_cnt;
	return total_energy;
}


decimal mns::MyDeformationHandler::Ec()
{
	if (fc.rows() != 3 * cps.set.size())
	{
		fc.resize(3 * cps.set.size());
	}
	fc = Eigen::VectorXd::Zero(3 * cps.set.size());
	decimal total_energy = 0;
	cout << "cps_size: " << (decimal)cps.set.size() << "    move_pt_num:" << (decimal)move_pt_num << endl;
	for (int i = 0; i < cps.set.size(); i++)
	{
		//decimal fac = (pt_move_or_not[i] == 1) ? (std::sqrt(1 / ((decimal)move_pt_num))) : (std::sqrt(1 / ((decimal)(cps.set.size() - move_pt_num))));
		decimal fac = (pt_move_or_not[i] == 1) ? ( std::sqrt(((decimal)cps.set.size()) / ((decimal)move_pt_num))) : (std::sqrt(((decimal)cps.set.size()) / ((decimal)(cps.set.size() - move_pt_num))));
		total_energy += GetConTerm(i, fac);
	}

	row_cnt += 3 * con_ele_cnt;
	return total_energy;
}

void mns::MyDeformationHandler::GetFeaturePos(std::vector<Eigen::Vector3d> _body_feature_pos)
{
	this->body_feature_pos = _body_feature_pos;
}

void mns::MyDeformationHandler::setFixed(std::vector<std::vector<int>> _fix_idx, std::vector<std::vector<double>> _fix_weights, std::vector<Eigen::Vector3d> _fix_pos)
{
	this->fixed_feature_idx = _fix_idx;
	this->fixed_feature_weights = _fix_weights;
	this->fixed_feature_pos = _fix_pos;
}

decimal mns::MyDeformationHandler::GetStrTerm1(integer from_id, integer to_id, integer e_id)
{
	//decimal ws = alpha_s / (ori_length_set[e_id] * ori_length_set[e_id]);
	decimal sws = this->sws_set[e_id];
	decimal energy = 0;

	if (this->use_sim_or_not == false)
	{
		//soft con
		decimal x_dif = unknown[3 * to_id + 0] - unknown[3 * from_id + 0];
		decimal y_dif = unknown[3 * to_id + 1] - unknown[3 * from_id + 1];
		decimal z_dif = unknown[3 * to_id + 2] - unknown[3 * from_id + 2];

		//hard con
		/*MyMesh::Point to_point;
		MyMesh::Point from_point;
		if (idx_move[to_id] == -1)
		{
		to_point[0] = cps.set[con_ori_map[to_id]].x;
		to_point[1] = cps.set[con_ori_map[to_id]].y;
		to_point[2] = cps.set[con_ori_map[to_id]].z;
		}
		else
		{
		to_point[0] = unknown[3 * (to_id - idx_move[to_id]) + 0];
		to_point[1] = unknown[3 * (to_id - idx_move[to_id]) + 1];
		to_point[2] = unknown[3 * (to_id - idx_move[to_id]) + 2];
		}
		if (idx_move[from_id] == -1)
		{
		from_point[0] = cps.set[con_ori_map[from_id]].x;
		from_point[1] = cps.set[con_ori_map[from_id]].y;
		from_point[2] = cps.set[con_ori_map[from_id]].z;
		}
		else
		{
		from_point[0] = unknown[3 * (from_id - idx_move[from_id]) + 0];
		from_point[1] = unknown[3 * (from_id - idx_move[from_id]) + 1];
		from_point[2] = unknown[3 * (from_id - idx_move[from_id]) + 2];
		}
		decimal x_dif = to_point.data()[0] - from_point.data()[0];
		decimal y_dif = to_point.data()[1] - from_point.data()[1];
		decimal z_dif = to_point.data()[2] - from_point.data()[2];*/

		decimal cur_sqr_l = x_dif * x_dif + y_dif * y_dif + z_dif * z_dif;
		decimal cur_l = std::sqrt(cur_sqr_l);

		/*fs(str_ele_cnt) = sws*(cur_l - wnt_length_set[e_id]);
		energy = fs(str_ele_cnt)*fs(str_ele_cnt);
		++str_ele_cnt;*/
		fsb(sb_ele_cnt) = sws * (cur_l - wnt_length_set[e_id]);
		energy = fsb(sb_ele_cnt)*fsb(sb_ele_cnt);
		++sb_ele_cnt;

		decimal norm_x = sws * x_dif / cur_l;
		decimal norm_y = sws * y_dif / cur_l;
		decimal norm_z = sws * z_dif / cur_l;

		//soft con
		//from
		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * from_id + 0, -norm_x));
		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * from_id + 1, -norm_y));
		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * from_id + 2, -norm_z));
		//to
		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * to_id + 0, norm_x));
		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * to_id + 1, norm_y));
		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * to_id + 2, norm_z));
		//if eg-driven
		if (this->use_eg_or_not == true && this->use_basis_or_not == false)
		{
			for (integer i = 0; i < this->eg_num; i++)
			{
				decimal tmp = -sws * (eg_length_vec[i][e_id] - ori_length_set[e_id]);
				triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * this->mesh.n_vertices() + i, tmp));
			}
		}
		//if basis exists
		if (this->use_eg_or_not == false && this->use_basis_or_not == true)
		{
			for (integer i = 0; i < this->basis.cols(); i++)
			{
				if (std::abs(/*basis(2 * e_id + 0, i)*/length_basis_vec[i][e_id])>0)
				{
					//decimal tmp = -sws*basis(2 * e_id + 0, i);
					//triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * this->mesh.n_vertices() + i, tmp));
					decimal tmp = -sws * ori_length_set[e_id] * (length_basis_vec[i][e_id] + 1);
					triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * this->mesh.n_vertices() + 2 * i + 0, tmp));
				}
			}
		}
		//

		//hard con
		////from
		//if (idx_move[from_id] != -1)
		//{
		//	triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * (from_id - idx_move[from_id]) + 0, -sws*x_dif / cur_l));
		//	triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * (from_id - idx_move[from_id]) + 1, -sws*y_dif / cur_l));
		//	triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * (from_id - idx_move[from_id]) + 2, -sws*z_dif / cur_l));
		//}
		////to
		//if (idx_move[to_id] != -1)
		//{
		//	triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * (to_id - idx_move[to_id]) + 0, sws*x_dif / cur_l));
		//	triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * (to_id - idx_move[to_id]) + 1, sws*y_dif / cur_l));
		//	triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * (to_id - idx_move[to_id]) + 2, sws*z_dif / cur_l));
		//}	

		////if eg-driven
		//if (this->use_eg_or_not == true)
		//{
		//	for (integer i = 0; i < this->eg_num; i++)
		//	{
		//		decimal tmp = -sws*(eg_length_vec[i][e_id] - ori_length_set[e_id]);
		//		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * (this->mesh.n_vertices() - cps.set.size()) + i, tmp));
		//	}
		//}
		////if basis exists
		//if (this->use_basis_or_not == true)
		//{
		//	for (integer i = 0; i < this->basis.cols(); i++)
		//	{
		//		decimal tmp = -sws*basis(2 * e_id + 0, i);
		//		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * (this->mesh.n_vertices() - cps.set.size()) + i, tmp));
		//	}
		//}
		//
	}
	if (this->use_sim_or_not == true)
	{
		/*decimal from_x = 0;
		decimal from_y = 0;
		decimal from_z = 0;
		decimal to_x = 0;
		decimal to_y = 0;
		decimal to_z = 0;
		for (integer j = 0; j < this->sim_mesh.n_vertices(); j++)
		{
		decimal w_from_j = W_mat(from_id, j);
		from_x += w_from_j * unknown(3 * j + 0);
		from_y += w_from_j * unknown(3 * j + 1);
		from_z += w_from_j * unknown(3 * j + 2);

		decimal w_to_j = W_mat(to_id, j);
		to_x += w_to_j * unknown(3 * j + 0);
		to_y += w_to_j * unknown(3 * j + 1);
		to_z += w_to_j * unknown(3 * j + 2);
		}*/
		//soft con
		decimal x_dif = high_v_mat(to_id, 0) - high_v_mat(from_id, 0);
		decimal y_dif = high_v_mat(to_id, 1) - high_v_mat(from_id, 1);
		decimal z_dif = high_v_mat(to_id, 2) - high_v_mat(from_id, 2);
		//hard con(not ready)

		decimal cur_sqr_l = x_dif * x_dif + y_dif * y_dif + z_dif * z_dif;
		decimal cur_l = std::sqrt(cur_sqr_l);

		/*fs(str_ele_cnt) = sws*(cur_l - wnt_length_set[e_id]);
		energy = fs(str_ele_cnt)*fs(str_ele_cnt);
		++str_ele_cnt;*/
		fsb(sb_ele_cnt) = sws * (cur_l - wnt_length_set[e_id]);
		energy = fsb(sb_ele_cnt)*fsb(sb_ele_cnt);
		++sb_ele_cnt;

		decimal norm_x = sws * x_dif / cur_l;
		decimal norm_y = sws * y_dif / cur_l;
		decimal norm_z = sws * z_dif / cur_l;

		//soft con
		//from
		j_triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * from_id + 0, -norm_x));
		j_triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * from_id + 1, -norm_y));
		j_triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * from_id + 2, -norm_z));
		//to
		j_triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * to_id + 0, norm_x));
		j_triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * to_id + 1, norm_y));
		j_triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * to_id + 2, norm_z));
		/*Eigen::VectorXd triple_tmp = Eigen::VectorXd::Zero(3 * this->sim_mesh.n_vertices());
		for (integer j = 0; j < this->W_idx_set[from_id].size(); j++)
		{
		decimal w_from_j = W_mat(from_id, W_idx_set[from_id][j]);
		triple_tmp(3 * W_idx_set[from_id][j] + 0) += -norm_x * w_from_j;
		triple_tmp(3 * W_idx_set[from_id][j] + 1) += -norm_y * w_from_j;
		triple_tmp(3 * W_idx_set[from_id][j] + 2) += -norm_z * w_from_j;
		}
		for (integer j = 0; j < this->W_idx_set[to_id].size(); j++)
		{
		decimal w_to_j = W_mat(to_id, W_idx_set[to_id][j]);
		triple_tmp(3 * W_idx_set[to_id][j] + 0) += norm_x * w_to_j;
		triple_tmp(3 * W_idx_set[to_id][j] + 1) += norm_y * w_to_j;
		triple_tmp(3 * W_idx_set[to_id][j] + 2) += norm_z * w_to_j;
		}*/
		/*for (integer j = 0; j < this->sim_mesh.n_vertices(); j++)
		{
		decimal w_from_j = W_mat(from_id, j);
		decimal w_to_j = W_mat(to_id, j);

		if (w_from_j > limit)
		{
		triple_tmp(3 * j + 0) += -sws*x_dif / cur_l * w_from_j;
		triple_tmp(3 * j + 1) += -sws*y_dif / cur_l * w_from_j;
		triple_tmp(3 * j + 2) += -sws*z_dif / cur_l * w_from_j;
		}
		if (w_to_j > limit)
		{
		triple_tmp(3 * j + 0) += sws*x_dif / cur_l * w_to_j;
		triple_tmp(3 * j + 1) += sws*y_dif / cur_l * w_to_j;
		triple_tmp(3 * j + 2) += sws*z_dif / cur_l * w_to_j;
		}
		}*/
		/*for (integer j = 0; j < this->W_idx_set[from_id].size(); j++)
		{
		if (triple_tmp(3 * W_idx_set[from_id][j] + 0) > limit)
		{
		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * W_idx_set[from_id][j] + 0, triple_tmp(3 * W_idx_set[from_id][j] + 0)));
		}
		if (triple_tmp(3 * W_idx_set[from_id][j] + 1) > limit)
		{
		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * W_idx_set[from_id][j] + 1, triple_tmp(3 * W_idx_set[from_id][j] + 1)));
		}
		if (triple_tmp(3 * W_idx_set[from_id][j] + 2) > limit)
		{
		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * W_idx_set[from_id][j] + 2, triple_tmp(3 * W_idx_set[from_id][j] + 2)));
		}
		}
		for (integer j = 0; j < this->W_idx_set[to_id].size(); j++)
		{
		if (triple_tmp(3 * W_idx_set[to_id][j] + 0) > limit)
		{
		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * W_idx_set[to_id][j] + 0, triple_tmp(3 * W_idx_set[to_id][j] + 0)));
		}
		if (triple_tmp(3 * W_idx_set[to_id][j] + 1) > limit)
		{
		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * W_idx_set[to_id][j] + 1, triple_tmp(3 * W_idx_set[to_id][j] + 1)));
		}
		if (triple_tmp(3 * W_idx_set[to_id][j] + 2) > limit)
		{
		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * W_idx_set[to_id][j] + 2, triple_tmp(3 * W_idx_set[to_id][j] + 2)));
		}
		}*/
		/*for (integer j = 0; j < this->sim_mesh.n_vertices(); j++)
		{
		if (triple_tmp(3 * j + 0) > limit)
		{
		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * j + 0, triple_tmp(3 * j + 0)));
		}
		if (triple_tmp(3 * j + 1) > limit)
		{
		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * j + 1, triple_tmp(3 * j + 1)));
		}
		if (triple_tmp(3 * j + 2) > limit)
		{
		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * j + 2, triple_tmp(3 * j + 2)));
		}
		}*/

		//if eg-driven
		if (this->use_eg_or_not == true && this->use_basis_or_not == false)
		{
			for (integer i = 0; i < this->eg_num; i++)
			{
				decimal tmp = -sws * (eg_length_vec[i][e_id] - ori_length_set[e_id]);
				triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * this->sim_mesh.n_vertices() + i, tmp));
			}
		}
		//if basis exists
		if (this->use_eg_or_not == false && this->use_basis_or_not == true)
		{
			for (integer i = 0; i < this->basis.cols(); i++)
			{
				if (std::abs(basis(2 * e_id + 0, i)) > 0)
				{
					decimal tmp = -sws * basis(2 * e_id + 0, i) *ori_length_set[e_id];
					triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * this->sim_mesh.n_vertices() + i, tmp));
				}
			}
		}
		//
	}

	++(this->row_cnt);

	return energy;
}

decimal mns::MyDeformationHandler::GetStrTerm(integer from_id, integer to_id, integer e_id)
{
	//decimal ws = alpha_s / (ori_length_set[e_id] * ori_length_set[e_id]);
	decimal sws = this->sws_set[e_id];
	decimal energy = 0;

	if (this->use_sim_or_not == false)
	{
		//soft con
		decimal x_dif = unknown[3 * to_id + 0] - unknown[3 * from_id + 0];
		decimal y_dif = unknown[3 * to_id + 1] - unknown[3 * from_id + 1];
		decimal z_dif = unknown[3 * to_id + 2] - unknown[3 * from_id + 2];

		//hard con
		/*MyMesh::Point to_point;
		MyMesh::Point from_point;
		if (idx_move[to_id] == -1)
		{
		to_point[0] = cps.set[con_ori_map[to_id]].x;
		to_point[1] = cps.set[con_ori_map[to_id]].y;
		to_point[2] = cps.set[con_ori_map[to_id]].z;
		}
		else
		{
		to_point[0] = unknown[3 * (to_id - idx_move[to_id]) + 0];
		to_point[1] = unknown[3 * (to_id - idx_move[to_id]) + 1];
		to_point[2] = unknown[3 * (to_id - idx_move[to_id]) + 2];
		}
		if (idx_move[from_id] == -1)
		{
		from_point[0] = cps.set[con_ori_map[from_id]].x;
		from_point[1] = cps.set[con_ori_map[from_id]].y;
		from_point[2] = cps.set[con_ori_map[from_id]].z;
		}
		else
		{
		from_point[0] = unknown[3 * (from_id - idx_move[from_id]) + 0];
		from_point[1] = unknown[3 * (from_id - idx_move[from_id]) + 1];
		from_point[2] = unknown[3 * (from_id - idx_move[from_id]) + 2];
		}
		decimal x_dif = to_point.data()[0] - from_point.data()[0];
		decimal y_dif = to_point.data()[1] - from_point.data()[1];
		decimal z_dif = to_point.data()[2] - from_point.data()[2];*/

		decimal cur_sqr_l = x_dif*x_dif + y_dif*y_dif + z_dif*z_dif;
		decimal cur_l = std::sqrt(cur_sqr_l);

		/*fs(str_ele_cnt) = sws*(cur_l - wnt_length_set[e_id]);
		energy = fs(str_ele_cnt)*fs(str_ele_cnt);
		++str_ele_cnt;*/
		fsb(sb_ele_cnt) = sws*(cur_l - wnt_length_set[e_id]);
		energy = fsb(sb_ele_cnt)*fsb(sb_ele_cnt);
		++sb_ele_cnt;

		decimal norm_x = sws*x_dif / cur_l;
		decimal norm_y = sws*y_dif / cur_l;
		decimal norm_z = sws*z_dif / cur_l;

		//soft con
		//from
		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * from_id + 0, -norm_x));
		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * from_id + 1, -norm_y));
		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * from_id + 2, -norm_z));
		//to
		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * to_id + 0, norm_x));
		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * to_id + 1, norm_y));
		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * to_id + 2, norm_z));
		//if eg-driven
		if (this->use_eg_or_not == true && this->use_basis_or_not == false)
		{
			for (integer i = 0; i < this->eg_num; i++)
			{
				decimal tmp = -sws*(eg_length_vec[i][e_id] - ori_length_set[e_id]);
				triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * this->mesh.n_vertices() + i, tmp));
			}
		}
		//if basis exists
		if (this->use_eg_or_not == false && this->use_basis_or_not == true)
		{
			for (integer i = 0; i < this->basis.cols(); i++)
			{
				if (std::abs(/*basis(2 * e_id + 0, i)*/length_basis_vec[i][e_id])>0)
				{
					//decimal tmp = -sws*basis(2 * e_id + 0, i);
					//triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * this->mesh.n_vertices() + i, tmp));
					decimal tmp = -sws* ori_length_set[e_id] * length_basis_vec[i][e_id];
					triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * this->mesh.n_vertices() + 2 * i + 0, tmp));
				}
			}
		}
		//

		//hard con
		////from
		//if (idx_move[from_id] != -1)
		//{
		//	triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * (from_id - idx_move[from_id]) + 0, -sws*x_dif / cur_l));
		//	triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * (from_id - idx_move[from_id]) + 1, -sws*y_dif / cur_l));
		//	triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * (from_id - idx_move[from_id]) + 2, -sws*z_dif / cur_l));
		//}
		////to
		//if (idx_move[to_id] != -1)
		//{
		//	triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * (to_id - idx_move[to_id]) + 0, sws*x_dif / cur_l));
		//	triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * (to_id - idx_move[to_id]) + 1, sws*y_dif / cur_l));
		//	triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * (to_id - idx_move[to_id]) + 2, sws*z_dif / cur_l));
		//}	

		////if eg-driven
		//if (this->use_eg_or_not == true)
		//{
		//	for (integer i = 0; i < this->eg_num; i++)
		//	{
		//		decimal tmp = -sws*(eg_length_vec[i][e_id] - ori_length_set[e_id]);
		//		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * (this->mesh.n_vertices() - cps.set.size()) + i, tmp));
		//	}
		//}
		////if basis exists
		//if (this->use_basis_or_not == true)
		//{
		//	for (integer i = 0; i < this->basis.cols(); i++)
		//	{
		//		decimal tmp = -sws*basis(2 * e_id + 0, i);
		//		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * (this->mesh.n_vertices() - cps.set.size()) + i, tmp));
		//	}
		//}
		//
	}
	if (this->use_sim_or_not == true)
	{
		/*decimal from_x = 0;
		decimal from_y = 0;
		decimal from_z = 0;
		decimal to_x = 0;
		decimal to_y = 0;
		decimal to_z = 0;
		for (integer j = 0; j < this->sim_mesh.n_vertices(); j++)
		{
		decimal w_from_j = W_mat(from_id, j);
		from_x += w_from_j * unknown(3 * j + 0);
		from_y += w_from_j * unknown(3 * j + 1);
		from_z += w_from_j * unknown(3 * j + 2);

		decimal w_to_j = W_mat(to_id, j);
		to_x += w_to_j * unknown(3 * j + 0);
		to_y += w_to_j * unknown(3 * j + 1);
		to_z += w_to_j * unknown(3 * j + 2);
		}*/
		//soft con
		decimal x_dif = high_v_mat(to_id, 0) - high_v_mat(from_id, 0);
		decimal y_dif = high_v_mat(to_id, 1) - high_v_mat(from_id, 1);
		decimal z_dif = high_v_mat(to_id, 2) - high_v_mat(from_id, 2);
		//hard con(not ready)

		decimal cur_sqr_l = x_dif*x_dif + y_dif*y_dif + z_dif*z_dif;
		decimal cur_l = std::sqrt(cur_sqr_l);

		/*fs(str_ele_cnt) = sws*(cur_l - wnt_length_set[e_id]);
		energy = fs(str_ele_cnt)*fs(str_ele_cnt);
		++str_ele_cnt;*/
		fsb(sb_ele_cnt) = sws*(cur_l - wnt_length_set[e_id]);
		energy = fsb(sb_ele_cnt)*fsb(sb_ele_cnt);
		++sb_ele_cnt;

		decimal norm_x = sws*x_dif / cur_l;
		decimal norm_y = sws*y_dif / cur_l;
		decimal norm_z = sws*z_dif / cur_l;

		//soft con
		//from
		j_triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * from_id + 0, -norm_x));
		j_triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * from_id + 1, -norm_y));
		j_triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * from_id + 2, -norm_z));
		//to
		j_triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * to_id + 0, norm_x));
		j_triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * to_id + 1, norm_y));
		j_triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * to_id + 2, norm_z));
		/*Eigen::VectorXd triple_tmp = Eigen::VectorXd::Zero(3 * this->sim_mesh.n_vertices());
		for (integer j = 0; j < this->W_idx_set[from_id].size(); j++)
		{
		decimal w_from_j = W_mat(from_id, W_idx_set[from_id][j]);
		triple_tmp(3 * W_idx_set[from_id][j] + 0) += -norm_x * w_from_j;
		triple_tmp(3 * W_idx_set[from_id][j] + 1) += -norm_y * w_from_j;
		triple_tmp(3 * W_idx_set[from_id][j] + 2) += -norm_z * w_from_j;
		}
		for (integer j = 0; j < this->W_idx_set[to_id].size(); j++)
		{
		decimal w_to_j = W_mat(to_id, W_idx_set[to_id][j]);
		triple_tmp(3 * W_idx_set[to_id][j] + 0) += norm_x * w_to_j;
		triple_tmp(3 * W_idx_set[to_id][j] + 1) += norm_y * w_to_j;
		triple_tmp(3 * W_idx_set[to_id][j] + 2) += norm_z * w_to_j;
		}*/
		/*for (integer j = 0; j < this->sim_mesh.n_vertices(); j++)
		{
		decimal w_from_j = W_mat(from_id, j);
		decimal w_to_j = W_mat(to_id, j);

		if (w_from_j > limit)
		{
		triple_tmp(3 * j + 0) += -sws*x_dif / cur_l * w_from_j;
		triple_tmp(3 * j + 1) += -sws*y_dif / cur_l * w_from_j;
		triple_tmp(3 * j + 2) += -sws*z_dif / cur_l * w_from_j;
		}
		if (w_to_j > limit)
		{
		triple_tmp(3 * j + 0) += sws*x_dif / cur_l * w_to_j;
		triple_tmp(3 * j + 1) += sws*y_dif / cur_l * w_to_j;
		triple_tmp(3 * j + 2) += sws*z_dif / cur_l * w_to_j;
		}
		}*/
		/*for (integer j = 0; j < this->W_idx_set[from_id].size(); j++)
		{
		if (triple_tmp(3 * W_idx_set[from_id][j] + 0) > limit)
		{
		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * W_idx_set[from_id][j] + 0, triple_tmp(3 * W_idx_set[from_id][j] + 0)));
		}
		if (triple_tmp(3 * W_idx_set[from_id][j] + 1) > limit)
		{
		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * W_idx_set[from_id][j] + 1, triple_tmp(3 * W_idx_set[from_id][j] + 1)));
		}
		if (triple_tmp(3 * W_idx_set[from_id][j] + 2) > limit)
		{
		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * W_idx_set[from_id][j] + 2, triple_tmp(3 * W_idx_set[from_id][j] + 2)));
		}
		}
		for (integer j = 0; j < this->W_idx_set[to_id].size(); j++)
		{
		if (triple_tmp(3 * W_idx_set[to_id][j] + 0) > limit)
		{
		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * W_idx_set[to_id][j] + 0, triple_tmp(3 * W_idx_set[to_id][j] + 0)));
		}
		if (triple_tmp(3 * W_idx_set[to_id][j] + 1) > limit)
		{
		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * W_idx_set[to_id][j] + 1, triple_tmp(3 * W_idx_set[to_id][j] + 1)));
		}
		if (triple_tmp(3 * W_idx_set[to_id][j] + 2) > limit)
		{
		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * W_idx_set[to_id][j] + 2, triple_tmp(3 * W_idx_set[to_id][j] + 2)));
		}
		}*/
		/*for (integer j = 0; j < this->sim_mesh.n_vertices(); j++)
		{
		if (triple_tmp(3 * j + 0) > limit)
		{
		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * j + 0, triple_tmp(3 * j + 0)));
		}
		if (triple_tmp(3 * j + 1) > limit)
		{
		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * j + 1, triple_tmp(3 * j + 1)));
		}
		if (triple_tmp(3 * j + 2) > limit)
		{
		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * j + 2, triple_tmp(3 * j + 2)));
		}
		}*/

		//if eg-driven
		if (this->use_eg_or_not == true && this->use_basis_or_not == false)
		{
			for (integer i = 0; i < this->eg_num; i++)
			{
				decimal tmp = -sws*(eg_length_vec[i][e_id] - ori_length_set[e_id]);
				triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * this->sim_mesh.n_vertices() + i, tmp));
			}
		}
		//if basis exists
		if (this->use_eg_or_not == false && this->use_basis_or_not == true)
		{
			for (integer i = 0; i < this->basis.cols(); i++)
			{
				if (std::abs(basis(2 * e_id + 0, i)) > 0)
				{
					decimal tmp = -sws*basis(2 * e_id + 0, i) *ori_length_set[e_id];
					triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * this->sim_mesh.n_vertices() + i, tmp));
				}
			}
		}
		//
	}

	++(this->row_cnt);

	return energy;
}

decimal mns::MyDeformationHandler::Es()
{
	if (fs.rows() != this->mesh.n_edges())
	{
		fs.resize(this->mesh.n_edges());
	}
	fs = Eigen::VectorXd::Zero(this->mesh.n_edges());
	decimal total_energy = 0;
	int i = 0;
	for (auto he_it = this->mesh.halfedges_begin(); he_it != this->mesh.halfedges_end(); ++i)
	{
		auto e = *he_it;
		auto fromVertex = this->mesh.from_vertex_handle(e);
		auto toVertex = this->mesh.to_vertex_handle(e);
		total_energy += GetStrTerm(fromVertex.idx(), toVertex.idx(), i);

		++he_it;
		++he_it;
	}
	return total_energy;
}

//decimal mns::MyDeformationHandler::CalDiheral(const MyMesh::Point& p1, const MyMesh::Point& p2, const MyMesh::Point& p3, const MyMesh::Point& p4)
//{
//	MyMesh::Normal e21 = p2 - p1;
//	MyMesh::Normal e31 = p3 - p1;
//	MyMesh::Normal e32 = p3 - p2;
//	MyMesh::Normal e41 = p4 - p1;
//	MyMesh::Normal e42 = p4 - p2;
//
//	MyMesh::Normal n124 = e21 % e41;
//	MyMesh::Normal n123 = e21 % e31;
//	decimal x = n123 | n124;
//	decimal l21 = e21.length();
//	decimal y = l21*(e31 | n124);
//	decimal angle = std::atan2(y, x);
//	return (angle > 0) ? (angle):(angle + M_PI);
//}

decimal mns::MyDeformationHandler::CalDiheral(const MyMesh::Point& rj, const MyMesh::Point& rk, const MyMesh::Point& ri/*left*/, const MyMesh::Point& rl/*right*/)
{
	MyMesh::Normal rij = ri - rj;
	MyMesh::Normal rkj = rk - rj;
	MyMesh::Normal rkl = rk - rl;

	MyMesh::Normal rmj = rij % rkj;
	MyMesh::Normal rnk = rkj % rkl;
	decimal rmj_norm = rmj.length();
	decimal rnk_norm = rnk.length();

	decimal cosine = ((rmj | rnk)) / (rmj_norm * rnk_norm);
	if (cosine > 1)
	{
		cosine = 1;
	}
	if (cosine < -1)
	{
		cosine = -1;
	}
	decimal angle = std::acos(cosine);
	return angle;
}

decimal mns::MyDeformationHandler::CalWLZHDiheral(const Eigen::Vector3d& normal1, const Eigen::Vector3d& normal2, const Eigen::Vector3d& e_ij)
{
	decimal phi = 0;
	Eigen::Vector3d crossnormal = normal1.cross(normal2);
	decimal dotcross = e_ij.dot(crossnormal);
	if (dotcross >= 0)
	{
		double dotNorm1Norm2 = normal1.dot(normal2);
		if (dotNorm1Norm2 < 1 - 1e-6 && dotNorm1Norm2 > -1 + 1e-6)
		{
			phi = std::acos(dotNorm1Norm2);
		}
		else if (dotNorm1Norm2 >= 1 - 1e-6)
		{
			phi = 0;
		}
		else
		{
			phi = 3.1415926;
		}
	}
	else
	{
		double dotNorm1Norm2 = normal1.dot(normal2);
		if (dotNorm1Norm2 < 1 - 1e-6 && dotNorm1Norm2 > -1 + 1e-6)
		{
			phi = -std::acos(dotNorm1Norm2);
		}
		else if (dotNorm1Norm2 >= 1 - 1e-6)
		{
			phi = 0;
		}
		else
		{
			phi = -3.1415926;
		}
	}
	return phi;
}

//decimal mns::MyDeformationHandler::GetBenTerm(integer from_id, integer to_id, integer left_id, integer right_id, integer e_id)
//{
//	decimal wb = (alpha_b * ori_length_set[e_id] * ori_length_set[e_id]) / ori_area_set[e_id];
//	decimal swb = std::sqrt(wb);
//
//	MyMesh::Point p1(unknown(3 * from_id + 0), unknown(3 * from_id + 1), unknown(3 * from_id + 2));
//	MyMesh::Point p2(unknown(3 * to_id + 0), unknown(3 * to_id + 1), unknown(3 * to_id + 2));
//	MyMesh::Point p4(unknown(3 * left_id + 0), unknown(3 * left_id + 1), unknown(3 * left_id + 2));
//	MyMesh::Point p3(unknown(3 * right_id + 0), unknown(3 * right_id + 1), unknown(3 * right_id + 2));
//
//	decimal cur_angle = CalDiheral(p1, p2, p3, p4);
//	fb(ben_ele_cnt) = swb * (cur_angle - wnt_dihedral_set[e_id]);
//	decimal energy = fb(ben_ele_cnt) * fb(ben_ele_cnt);
//	++ben_ele_cnt;
//
//	MyMesh::Normal e = p2 - p1;
//	MyMesh::Normal e31 = p3 - p1;
//	MyMesh::Normal e32 = p3 - p2;
//	MyMesh::Normal e41 = p4 - p1;
//	MyMesh::Normal e42 = p4 - p2;
//	MyMesh::Normal n1 = e32 % e31;
//	MyMesh::Normal n2 = e41 % e42;
//
//	decimal e_norm = e.length();
//	decimal n1_norm = n1.length();
//	decimal n2_norm = n2.length();
//
//	MyMesh::Normal partial_p1 = swb*(((e32 | e) / (e_norm*n1_norm*n1_norm))*n1 + ((e42 | e) / (e_norm*n2_norm*n2_norm))*n2);
//	MyMesh::Normal partial_p2 = swb*(-((e31 | e) / (e_norm*n1_norm*n1_norm))*n1 - ((e41 | e) / (e_norm*n2_norm*n2_norm))*n2);
//	MyMesh::Normal partial_p3 = ((swb*e_norm) / (n1_norm*n1_norm)) * n1;
//	MyMesh::Normal partial_p4 = ((swb*e_norm) / (n2_norm*n2_norm)) * n2;
//
//	MyMesh::Normal partial_sum = partial_p1 + partial_p2 + partial_p3 + partial_p4;
//	//x1
//	triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * from_id + 0, partial_p1.data()[0]));
//	triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * from_id + 1, partial_p1.data()[1]));
//	triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * from_id + 2, partial_p1.data()[2]));
//	//x2
//	triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * to_id + 0, partial_p2.data()[0]));
//	triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * to_id + 1, partial_p2.data()[1]));
//	triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * to_id + 2, partial_p2.data()[2]));
//	//x3
//	triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * right_id + 0, partial_p3.data()[0]));
//	triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * right_id + 1, partial_p3.data()[1]));
//	triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * right_id + 2, partial_p3.data()[2]));
//	//x4
//	triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * left_id + 0, partial_p4.data()[0]));
//	triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * left_id + 1, partial_p4.data()[1]));
//	triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * left_id + 2, partial_p4.data()[2]));
//	++row_cnt;
//
//	return energy;
//}

decimal mns::MyDeformationHandler::GetBenTerm(integer from_id, integer to_id, integer left_id, integer right_id, integer e_id)
{
	//decimal wb = (alpha_b * ori_length_set[e_id] * ori_length_set[e_id]) / ori_area_set[e_id];
	decimal swb = this->swb_set[e_id];
	decimal energy = 0;

	if (this->use_sim_or_not == false)
	{
		//soft con
		MyMesh::Point rj(unknown(3 * from_id + 0), unknown(3 * from_id + 1), unknown(3 * from_id + 2));
		MyMesh::Point rk(unknown(3 * to_id + 0), unknown(3 * to_id + 1), unknown(3 * to_id + 2));
		MyMesh::Point ri(unknown(3 * left_id + 0), unknown(3 * left_id + 1), unknown(3 * left_id + 2));
		MyMesh::Point rl(unknown(3 * right_id + 0), unknown(3 * right_id + 1), unknown(3 * right_id + 2));

		//hard con
		//MyMesh::Point rj;
		//MyMesh::Point rk;
		//MyMesh::Point ri;
		//MyMesh::Point rl;
		////rj
		//if (idx_move[from_id] == -1)
		//{
		//	rj[0] = cps.set[con_ori_map[from_id]].x;
		//	rj[1] = cps.set[con_ori_map[from_id]].y;
		//	rj[2] = cps.set[con_ori_map[from_id]].z;
		//}
		//else
		//{
		//	rj[0] = unknown[3 * (from_id - idx_move[from_id]) + 0];
		//	rj[1] = unknown[3 * (from_id - idx_move[from_id]) + 1];
		//	rj[2] = unknown[3 * (from_id - idx_move[from_id]) + 2];
		//}
		////rk
		//if (idx_move[to_id] == -1)
		//{
		//	rk[0] = cps.set[con_ori_map[to_id]].x;
		//	rk[1] = cps.set[con_ori_map[to_id]].y;
		//	rk[2] = cps.set[con_ori_map[to_id]].z;
		//}
		//else
		//{
		//	rk[0] = unknown[3 * (to_id - idx_move[to_id]) + 0];
		//	rk[1] = unknown[3 * (to_id - idx_move[to_id]) + 1];
		//	rk[2] = unknown[3 * (to_id - idx_move[to_id]) + 2];
		//}
		////ri
		//if (idx_move[left_id] == -1)
		//{
		//	ri[0] = cps.set[con_ori_map[left_id]].x;
		//	ri[1] = cps.set[con_ori_map[left_id]].y;
		//	ri[2] = cps.set[con_ori_map[left_id]].z;
		//}
		//else
		//{
		//	ri[0] = unknown[3 * (left_id - idx_move[left_id]) + 0];
		//	ri[1] = unknown[3 * (left_id - idx_move[left_id]) + 1];
		//	ri[2] = unknown[3 * (left_id - idx_move[left_id]) + 2];
		//}
		////rl
		//if (idx_move[right_id] == -1)
		//{
		//	rl[0] = cps.set[con_ori_map[right_id]].x;
		//	rl[1] = cps.set[con_ori_map[right_id]].y;
		//	rl[2] = cps.set[con_ori_map[right_id]].z;
		//}
		//else
		//{
		//	rl[0] = unknown[3 * (right_id - idx_move[right_id]) + 0];
		//	rl[1] = unknown[3 * (right_id - idx_move[right_id]) + 1];
		//	rl[2] = unknown[3 * (right_id - idx_move[right_id]) + 2];
		//}
		//

		decimal cur_angle = CalDiheral(rj, rk, ri, rl);
		/*fb(ben_ele_cnt) = swb * (cur_angle - wnt_dihedral_set[e_id]);
		energy = fb(ben_ele_cnt) * fb(ben_ele_cnt);
		++ben_ele_cnt;*/
		fsb(sb_ele_cnt) = swb * (cur_angle - wnt_dihedral_set[e_id]);
		energy = fsb(sb_ele_cnt) * fsb(sb_ele_cnt);
		++sb_ele_cnt;

		MyMesh::Normal rij = ri - rj;
		MyMesh::Normal rkj = rk - rj;
		MyMesh::Normal rkl = rk - rl;

		MyMesh::Normal rmj = rij % rkj;
		MyMesh::Normal rnk = rkj % rkl;

		decimal rkj_norm = rkj.length();
		decimal rmj_norm = rmj.length();
		decimal rnk_norm = rnk.length();

		decimal rij_rkj = rij | rkj;
		decimal rkl_rkj = rkl | rkj;

		decimal rkj_sqr_norm = rkj_norm * rkj_norm;

		MyMesh::Normal partial_ri = (rkj_norm / (rmj_norm*rmj_norm)) * rmj;
		MyMesh::Normal partial_rl = -(rkj_norm / (rnk_norm*rnk_norm)) * rnk;
		MyMesh::Normal partial_rj = ((rij_rkj / rkj_sqr_norm - 1) * partial_ri) - ((rkl_rkj / rkj_sqr_norm) * partial_rl);
		MyMesh::Normal partial_rk = ((rkl_rkj / rkj_sqr_norm - 1) * partial_rl) - ((rij_rkj / rkj_sqr_norm) * partial_ri);

		MyMesh::Normal partial_sum = partial_ri + partial_rl + partial_rj + partial_rk;

		//set reverse
		swb = -swb;

		//soft con
		//rj
		triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * from_id + 0, swb*partial_rj.data()[0]));
		triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * from_id + 1, swb*partial_rj.data()[1]));
		triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * from_id + 2, swb*partial_rj.data()[2]));
		//rk
		triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * to_id + 0, swb*partial_rk.data()[0]));
		triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * to_id + 1, swb*partial_rk.data()[1]));
		triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * to_id + 2, swb*partial_rk.data()[2]));
		//rl
		triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * right_id + 0, swb*partial_rl.data()[0]));
		triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * right_id + 1, swb*partial_rl.data()[1]));
		triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * right_id + 2, swb*partial_rl.data()[2]));
		//ri
		triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * left_id + 0, swb*partial_ri.data()[0]));
		triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * left_id + 1, swb*partial_ri.data()[1]));
		triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * left_id + 2, swb*partial_ri.data()[2]));
		//if eg-driven
		if (this->use_eg_or_not == true && this->use_basis_or_not == false)
		{
			for (integer i = 0; i < this->eg_num; i++)
			{
				triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * this->mesh.n_vertices() + i, swb*(eg_dihedral_vec[i][e_id] - ori_dihedral_set[e_id])));
			}
		}
		//if basis exists
		if (this->use_eg_or_not == false && this->use_basis_or_not == true)
		{
			for (integer i = 0; i < this->basis.cols(); i++)
			{
				if (std::abs(/*basis(2 * e_id + 1, i)*/dihedral_basis_vec[i][e_id]) > 0)
				{
					//decimal tmp = swb*basis(2 * e_id + 1, i);
					//triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * this->mesh.n_vertices() + i, tmp));
					decimal tmp = swb*dihedral_basis_vec[i][e_id];
					triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * this->mesh.n_vertices() + 2 * i + 1, tmp));
				}
			}
		}
	}
	if (this->use_sim_or_not == true)
	{
		/*decimal from_x = 0;
		decimal from_y = 0;
		decimal from_z = 0;
		decimal to_x = 0;
		decimal to_y = 0;
		decimal to_z = 0;
		decimal left_x = 0;
		decimal left_y = 0;
		decimal left_z = 0;
		decimal right_x = 0;
		decimal right_y = 0;
		decimal right_z = 0;
		for (integer j = 0; j < this->sim_mesh.n_vertices(); j++)
		{
		decimal w_from_j = W_mat(from_id, j);
		from_x += w_from_j * unknown(3 * j + 0);
		from_y += w_from_j * unknown(3 * j + 1);
		from_z += w_from_j * unknown(3 * j + 2);

		decimal w_to_j = W_mat(to_id, j);
		to_x += w_to_j * unknown(3 * j + 0);
		to_y += w_to_j * unknown(3 * j + 1);
		to_z += w_to_j * unknown(3 * j + 2);

		decimal w_left_j = W_mat(left_id, j);
		left_x += w_left_j * unknown(3 * j + 0);
		left_y += w_left_j * unknown(3 * j + 1);
		left_z += w_left_j * unknown(3 * j + 2);

		decimal w_right_j = W_mat(right_id, j);
		right_x += w_right_j * unknown(3 * j + 0);
		right_y += w_right_j * unknown(3 * j + 1);
		right_z += w_right_j * unknown(3 * j + 2);
		}*/

		MyMesh::Point rj(high_v_mat(from_id, 0), high_v_mat(from_id, 1), high_v_mat(from_id, 2));
		MyMesh::Point rk(high_v_mat(to_id, 0), high_v_mat(to_id, 1), high_v_mat(to_id, 2));
		MyMesh::Point ri(high_v_mat(left_id, 0), high_v_mat(left_id, 1), high_v_mat(left_id, 2));
		MyMesh::Point rl(high_v_mat(right_id, 0), high_v_mat(right_id, 1), high_v_mat(right_id, 2));

		decimal cur_angle = CalDiheral(rj, rk, ri, rl);
		/*fb(ben_ele_cnt) = swb * (cur_angle - wnt_dihedral_set[e_id]);
		energy = fb(ben_ele_cnt) * fb(ben_ele_cnt);
		++ben_ele_cnt;*/
		fsb(sb_ele_cnt) = swb * (cur_angle - wnt_dihedral_set[e_id]);
		energy = fsb(sb_ele_cnt) * fsb(sb_ele_cnt);
		++sb_ele_cnt;

		MyMesh::Normal rij = ri - rj;
		MyMesh::Normal rkj = rk - rj;
		MyMesh::Normal rkl = rk - rl;

		MyMesh::Normal rmj = rij % rkj;
		MyMesh::Normal rnk = rkj % rkl;

		decimal rkj_norm = rkj.length();
		decimal rmj_norm = rmj.length();
		decimal rnk_norm = rnk.length();

		decimal rij_rkj = rij | rkj;
		decimal rkl_rkj = rkl | rkj;

		decimal rkj_sqr_norm = rkj_norm * rkj_norm;

		MyMesh::Normal partial_ri = (rkj_norm / (rmj_norm*rmj_norm)) * rmj;
		MyMesh::Normal partial_rl = -(rkj_norm / (rnk_norm*rnk_norm)) * rnk;
		MyMesh::Normal partial_rj = ((rij_rkj / rkj_sqr_norm - 1) * partial_ri) - ((rkl_rkj / rkj_sqr_norm) * partial_rl);
		MyMesh::Normal partial_rk = ((rkl_rkj / rkj_sqr_norm - 1) * partial_rl) - ((rij_rkj / rkj_sqr_norm) * partial_ri);

		MyMesh::Normal partial_sum = partial_ri + partial_rl + partial_rj + partial_rk;

		//set reverse
		swb = -swb;

		//soft con
		//rj
		j_triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * from_id + 0, swb*partial_rj.data()[0]));
		j_triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * from_id + 1, swb*partial_rj.data()[1]));
		j_triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * from_id + 2, swb*partial_rj.data()[2]));
		//rk
		j_triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * to_id + 0, swb*partial_rk.data()[0]));
		j_triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * to_id + 1, swb*partial_rk.data()[1]));
		j_triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * to_id + 2, swb*partial_rk.data()[2]));
		//rl
		j_triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * right_id + 0, swb*partial_rl.data()[0]));
		j_triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * right_id + 1, swb*partial_rl.data()[1]));
		j_triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * right_id + 2, swb*partial_rl.data()[2]));
		//ri
		j_triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * left_id + 0, swb*partial_ri.data()[0]));
		j_triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * left_id + 1, swb*partial_ri.data()[1]));
		j_triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * left_id + 2, swb*partial_ri.data()[2]));

		/*Eigen::VectorXd triple_tmp = Eigen::VectorXd::Zero(3 * this->sim_mesh.n_vertices());
		for (integer j = 0; j < this->W_idx_set[from_id].size(); j++)
		{
		decimal w_from_j = W_mat(from_id, W_idx_set[from_id][j]);
		triple_tmp(3 * W_idx_set[from_id][j] + 0) += swb*partial_rj.data()[0] * w_from_j;
		triple_tmp(3 * W_idx_set[from_id][j] + 1) += swb*partial_rj.data()[1] * w_from_j;
		triple_tmp(3 * W_idx_set[from_id][j] + 2) += swb*partial_rj.data()[2] * w_from_j;
		}
		for (integer j = 0; j < this->W_idx_set[to_id].size(); j++)
		{
		decimal w_to_j = W_mat(to_id, W_idx_set[to_id][j]);
		triple_tmp(3 * W_idx_set[to_id][j] + 0) += swb*partial_rk.data()[0] * w_to_j;
		triple_tmp(3 * W_idx_set[to_id][j] + 1) += swb*partial_rk.data()[1] * w_to_j;
		triple_tmp(3 * W_idx_set[to_id][j] + 2) += swb*partial_rk.data()[2] * w_to_j;
		}
		for (integer j = 0; j < this->W_idx_set[left_id].size(); j++)
		{
		decimal w_left_j = W_mat(left_id, W_idx_set[left_id][j]);
		triple_tmp(3 * W_idx_set[left_id][j] + 0) += swb*partial_ri.data()[0] * w_left_j;
		triple_tmp(3 * W_idx_set[left_id][j] + 1) += swb*partial_ri.data()[1] * w_left_j;
		triple_tmp(3 * W_idx_set[left_id][j] + 2) += swb*partial_ri.data()[2] * w_left_j;
		}
		for (integer j = 0; j < this->W_idx_set[right_id].size(); j++)
		{
		decimal w_right_j = W_mat(right_id, W_idx_set[right_id][j]);
		triple_tmp(3 * W_idx_set[right_id][j] + 0) += swb*partial_rl.data()[0] * w_right_j;
		triple_tmp(3 * W_idx_set[right_id][j] + 1) += swb*partial_rl.data()[1] * w_right_j;
		triple_tmp(3 * W_idx_set[right_id][j] + 2) += swb*partial_rl.data()[2] * w_right_j;
		}*/
		/*for (integer j = 0; j < this->sim_mesh.n_vertices(); j++)
		{
		decimal w_from_j = W_mat(from_id, j);
		decimal w_to_j = W_mat(to_id, j);
		decimal w_left_j = W_mat(left_id, j);
		decimal w_right_j = W_mat(right_id, j);

		if (w_from_j > limit)
		{
		triple_tmp(3 * j + 0) += swb*partial_rj.data()[0] * w_from_j;
		triple_tmp(3 * j + 1) += swb*partial_rj.data()[1] * w_from_j;
		triple_tmp(3 * j + 2) += swb*partial_rj.data()[2] * w_from_j;
		}
		if (w_to_j > limit)
		{
		triple_tmp(3 * j + 0) += swb*partial_rk.data()[0] * w_to_j;
		triple_tmp(3 * j + 1) += swb*partial_rk.data()[1] * w_to_j;
		triple_tmp(3 * j + 2) += swb*partial_rk.data()[2] * w_to_j;
		}
		if (w_left_j > limit)
		{
		triple_tmp(3 * j + 0) += swb*partial_ri.data()[0] * w_left_j;
		triple_tmp(3 * j + 1) += swb*partial_ri.data()[1] * w_left_j;
		triple_tmp(3 * j + 2) += swb*partial_ri.data()[2] * w_left_j;
		}
		if (w_right_j > limit)
		{
		triple_tmp(3 * j + 0) += swb*partial_rl.data()[0] * w_right_j;
		triple_tmp(3 * j + 1) += swb*partial_rl.data()[1] * w_right_j;
		triple_tmp(3 * j + 2) += swb*partial_rl.data()[2] * w_right_j;
		}
		}*/
		/*for (integer j = 0; j < this->sim_mesh.n_vertices(); j++)
		{
		if (triple_tmp(3 * j + 0) > limit)
		{
		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * j + 0, triple_tmp(3 * j + 0)));
		}
		if (triple_tmp(3 * j + 1) > limit)
		{
		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * j + 1, triple_tmp(3 * j + 1)));
		}
		if (triple_tmp(3 * j + 2) > limit)
		{
		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * j + 2, triple_tmp(3 * j + 2)));
		}
		}*/

		//if eg-driven
		if (this->use_eg_or_not == true && this->use_basis_or_not == false)
		{
			for (integer i = 0; i < this->eg_num; i++)
			{
				triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * this->sim_mesh.n_vertices() + i, swb*(eg_dihedral_vec[i][e_id] - ori_dihedral_set[e_id])));
			}
		}

		//if basis exists
		if (this->use_eg_or_not == false && this->use_basis_or_not == true)
		{
			for (integer i = 0; i < this->basis.cols(); i++)
			{
				if (std::abs(basis(2 * e_id + 1, i)) > 0)
				{
					decimal tmp = swb*basis(2 * e_id + 1, i);
					triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * this->sim_mesh.n_vertices() + i, tmp));
				}
			}
		}
	}

	//hard con
	////rj
	//if (idx_move[from_id] != -1)
	//{
	//	triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * (from_id - idx_move[from_id]) + 0, swb*partial_rj.data()[0]));
	//	triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * (from_id - idx_move[from_id]) + 1, swb*partial_rj.data()[1]));
	//	triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * (from_id - idx_move[from_id]) + 2, swb*partial_rj.data()[2]));
	//}	
	////rk
	//if (idx_move[to_id] != -1)
	//{
	//	triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * (to_id - idx_move[to_id]) + 0, swb*partial_rk.data()[0]));
	//	triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * (to_id - idx_move[to_id]) + 1, swb*partial_rk.data()[1]));
	//	triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * (to_id - idx_move[to_id]) + 2, swb*partial_rk.data()[2]));
	//}	
	////rl
	//if (idx_move[right_id] != -1)
	//{
	//	triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * (right_id - idx_move[right_id]) + 0, swb*partial_rl.data()[0]));
	//	triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * (right_id - idx_move[right_id]) + 1, swb*partial_rl.data()[1]));
	//	triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * (right_id - idx_move[right_id]) + 2, swb*partial_rl.data()[2]));
	//}	
	////ri
	//if (idx_move[left_id] != -1)
	//{
	//	triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * (left_id - idx_move[left_id]) + 0, swb*partial_ri.data()[0]));
	//	triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * (left_id - idx_move[left_id]) + 1, swb*partial_ri.data()[1]));
	//	triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * (left_id - idx_move[left_id]) + 2, swb*partial_ri.data()[2]));
	//}	
	////if eg-driven
	//if (this->use_eg_or_not == true)
	//{
	//	for (integer i = 0; i < this->eg_num; i++)
	//	{
	//		triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * (this->mesh.n_vertices() - cps.set.size()) + i, swb*(eg_dihedral_vec[i][e_id] - ori_dihedral_set[e_id])));
	//	}
	//}
	////if basis exists
	//if (this->use_basis_or_not == true)
	//{
	//	for (integer i = 0; i < this->basis.cols(); i++)
	//	{
	//		decimal tmp = swb*basis(2 * e_id + 1, i);
	//		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * (this->mesh.n_vertices() - cps.set.size()) + i, tmp));
	//	}
	//}
	//
	++row_cnt;

	return energy;
}

decimal mns::MyDeformationHandler::Eb()
{
	if (fb.rows() != this->mesh.n_edges())
	{
		fb.resize(this->mesh.n_edges());
	}
	fb = Eigen::VectorXd::Zero(this->mesh.n_edges());
	decimal total_energy = 0;
	int i = 0;
	for (auto he_it = this->mesh.halfedges_begin(); he_it != this->mesh.halfedges_end(); ++i)
	{
		auto e = *he_it;
		auto f1 = this->mesh.face_handle(e);
		auto op_he = this->mesh.opposite_halfedge_handle(e);
		auto f2 = this->mesh.face_handle(op_he);
		auto fromVertex = this->mesh.from_vertex_handle(e);
		auto toVertex = this->mesh.to_vertex_handle(e);
		if (f1.idx() == -1 || f2.idx() == -1)
		{
			++sb_ele_cnt;
			++row_cnt;

			//edge, not half edge
			++he_it;
			++he_it;
			continue;
		}

		integer leftVertex, rightVertex;
		//integer temp_cnt = 0;
		for (auto fv_it1 = this->mesh.fv_begin(f1); fv_it1 != this->mesh.fv_end(f1); ++fv_it1)
		{
			auto v = *fv_it1;
			//temp_cnt++;
			if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
			{
				leftVertex = v.idx();
				break;
			}
		}
		//temp_cnt = 0;
		for (auto fv_it2 = this->mesh.fv_begin(f2); fv_it2 != this->mesh.fv_end(f2); ++fv_it2)
		{
			auto v = *fv_it2;
			//temp_cnt++;
			if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
			{
				rightVertex = v.idx();
				break;
			}
		}
		total_energy += GetBenTerm(fromVertex.idx(), toVertex.idx(), leftVertex, rightVertex, i);
		if (_isnan(total_energy))
		{
			//std::cout << i;
		}
		++he_it;
		++he_it;
	}

	return total_energy;
}

std::vector<decimal> mns::MyDeformationHandler::Esb()
{
	std::vector<decimal> energy;
	energy.push_back(0);
	energy.push_back(0);

	/*if (fs.rows() != this->mesh.n_edges())
	{
	fs.resize(this->mesh.n_edges());
	}
	fs = Eigen::VectorXd::Zero(this->mesh.n_edges());
	if (fb.rows() != this->mesh.n_edges())
	{
	fb.resize(this->mesh.n_edges());
	}
	fb = Eigen::VectorXd::Zero(this->mesh.n_edges());*/
	if (fsb.rows() != 2 * this->mesh.n_edges())
	{
		fsb.resize(2 * this->mesh.n_edges());
	}
	fsb = Eigen::VectorXd::Zero(2 * this->mesh.n_edges());
	if (alpha_b == 0 || alpha_s == 0)
	{
		fsb.resize(this->mesh.n_edges());
		fsb = Eigen::VectorXd::Zero(this->mesh.n_edges());
	}

	int i = 0;
	for (auto he_it = this->mesh.halfedges_begin(); he_it != this->mesh.halfedges_end(); ++i)
	{
		auto e = *he_it;
		auto f1 = this->mesh.face_handle(e);
		auto op_he = this->mesh.opposite_halfedge_handle(e);
		auto f2 = this->mesh.face_handle(op_he);
		auto fromVertex = this->mesh.from_vertex_handle(e);
		auto toVertex = this->mesh.to_vertex_handle(e);

		if (alpha_s != 0)
		{
			energy[0] += GetStrTerm(fromVertex.idx(), toVertex.idx(), i);
		}

		if (alpha_b != 0)
		{
			if (f1.idx() == -1 || f2.idx() == -1)
			{
				++sb_ele_cnt;
				++row_cnt;

				//edge, not half edge
				++he_it;
				++he_it;
				continue;
			}

			integer leftVertex, rightVertex;
			//integer temp_cnt = 0;
			for (auto fv_it1 = this->mesh.fv_begin(f1); fv_it1 != this->mesh.fv_end(f1); ++fv_it1)
			{
				auto v = *fv_it1;
				//temp_cnt++;
				if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
				{
					leftVertex = v.idx();
					break;
				}
			}
			//temp_cnt = 0;
			for (auto fv_it2 = this->mesh.fv_begin(f2); fv_it2 != this->mesh.fv_end(f2); ++fv_it2)
			{
				auto v = *fv_it2;
				//temp_cnt++;
				if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
				{
					rightVertex = v.idx();
					break;
				}
			}
			energy[1] += GetBenTerm(fromVertex.idx(), toVertex.idx(), leftVertex, rightVertex, i);
		}

		if (_isnan(energy[1]))
		{
			std::cout << i;
		}

		++he_it;
		++he_it;
	}

	return energy;
}
decimal mns::MyDeformationHandler::GetBenTermNew1(integer from_id, integer to_id, integer left_id, integer right_id, integer e_id, const Eigen::Vector3d& normal1, const Eigen::Vector3d& normal2, const Eigen::Vector3d& e_ij)
{
	decimal swb = this->swb_set[e_id];
	decimal energy = 0;

	if (this->use_sim_or_not == false)
	{
		//soft con
		Eigen::Vector3d p3(unknown(3 * from_id + 0), unknown(3 * from_id + 1), unknown(3 * from_id + 2));
		Eigen::Vector3d p4(unknown(3 * to_id + 0), unknown(3 * to_id + 1), unknown(3 * to_id + 2));
		Eigen::Vector3d p1(unknown(3 * left_id + 0), unknown(3 * left_id + 1), unknown(3 * left_id + 2));
		Eigen::Vector3d p2(unknown(3 * right_id + 0), unknown(3 * right_id + 1), unknown(3 * right_id + 2));

		/*MyMesh::Point rj(unknown(3 * from_id + 0), unknown(3 * from_id + 1), unknown(3 * from_id + 2));
		MyMesh::Point rk(unknown(3 * to_id + 0), unknown(3 * to_id + 1), unknown(3 * to_id + 2));
		MyMesh::Point ri(unknown(3 * left_id + 0), unknown(3 * left_id + 1), unknown(3 * left_id + 2));
		MyMesh::Point rl(unknown(3 * right_id + 0), unknown(3 * right_id + 1), unknown(3 * right_id + 2));
		*/
		decimal n1_norm = normal1.norm();
		decimal n2_norm = normal2.norm();
		decimal e_norm = e_ij.norm();
		Eigen::Vector3d n1 = normal1 / n1_norm;
		Eigen::Vector3d n2 = normal2 / n2_norm;
		Eigen::Vector3d e = e_ij / e_norm;

		decimal cur_angle = CalWLZHDiheral(n1, n2, e);
		fsb(sb_ele_cnt) = swb * (cur_angle - (wnt_wlz_dihedral_set[e_id] ) );
		energy = fsb(sb_ele_cnt) * fsb(sb_ele_cnt);
		if (_isnan(energy))
		{
			cout << swb << " " << cur_angle << " " << wnt_wlz_dihedral_set[e_id] << " " << endl;
			std::cout << energy;
		}
		++sb_ele_cnt;

		decimal area1 = CalTriangleArea(p3, p4, p1);
		Eigen::Vector3d p34 = (p3 - p4);
		Eigen::Vector3d p41 = (p4 - p1);
		Eigen::Vector3d p13 = (p1 - p3);
		Eigen::Matrix3d dn1dp1 = -CalTensor(p34) / (2 * area1);
		Eigen::Matrix3d dn1dp2 = Eigen::Matrix3d::Zero();
		Eigen::Matrix3d dn1dp3 = -CalTensor(p41) / (2 * area1);
		Eigen::Matrix3d dn1dp4 = -CalTensor(p13) / (2 * area1);

		decimal area2 = CalTriangleArea(p3, p4, p2); // ?
		Eigen::Vector3d p43 = (p4 - p3);
		Eigen::Vector3d p24 = (p2 - p4);
		Eigen::Vector3d p32 = (p3 - p2);
		Eigen::Matrix3d dn2dp1 = Eigen::Matrix3d::Zero();
		Eigen::Matrix3d dn2dp2 = -CalTensor(p43) / (2 * area2);
		Eigen::Matrix3d dn2dp3 = -CalTensor(p24) / (2 * area2);
		Eigen::Matrix3d dn2dp4 = -CalTensor(p32) / (2 * area2);

		Eigen::Matrix3d dedp1 = Eigen::Matrix3d::Zero();
		Eigen::Matrix3d dedp2 = Eigen::Matrix3d::Zero();
		Eigen::Matrix3d dedp3 = -Eigen::Matrix3d::Identity() / p34.norm();
		Eigen::Matrix3d dedp4 = Eigen::Matrix3d::Identity() / p34.norm();

		Eigen::Matrix3d n1_tensor = CalTensor(n1);
		Eigen::Matrix3d n2_tensor = CalTensor(n2);
		Eigen::Matrix3d dn1crossn2dp1 = n1_tensor * dn2dp1 - n2_tensor * dn1dp1;
		Eigen::Matrix3d dn1crossn2dp2 = n1_tensor * dn2dp2 - n2_tensor * dn1dp2;
		Eigen::Matrix3d dn1crossn2dp3 = n1_tensor * dn2dp3 - n2_tensor * dn1dp3;
		Eigen::Matrix3d dn1crossn2dp4 = n1_tensor * dn2dp4 - n2_tensor * dn1dp4;

		Eigen::Vector3d dsindp1 = dn1crossn2dp1.transpose()*e + dedp1.transpose()*(n1.cross(n2));
		Eigen::Vector3d dsindp2 = dn1crossn2dp2.transpose()*e + dedp2.transpose()*(n1.cross(n2));
		Eigen::Vector3d dsindp3 = dn1crossn2dp3.transpose()*e + dedp3.transpose()*(n1.cross(n2));
		Eigen::Vector3d dsindp4 = dn1crossn2dp4.transpose()*e + dedp4.transpose()*(n1.cross(n2));

		Eigen::Vector3d dcosdp1 = dn1dp1.transpose()*n2 + dn2dp1.transpose()*n1;
		Eigen::Vector3d dcosdp2 = dn1dp2.transpose()*n2 + dn2dp2.transpose()*n1;
		Eigen::Vector3d dcosdp3 = dn1dp3.transpose()*n2 + dn2dp3.transpose()*n1;
		Eigen::Vector3d dcosdp4 = dn1dp4.transpose()*n2 + dn2dp4.transpose()*n1;

		Eigen::Vector3d dthetadp1 = std::cos(cur_angle) * dsindp1 - std::sin(cur_angle) * dcosdp1;
		Eigen::Vector3d dthetadp2 = std::cos(cur_angle) * dsindp2 - std::sin(cur_angle) * dcosdp2;
		Eigen::Vector3d dthetadp3 = std::cos(cur_angle) * dsindp3 - std::sin(cur_angle) * dcosdp3;
		Eigen::Vector3d dthetadp4 = std::cos(cur_angle) * dsindp4 - std::sin(cur_angle) * dcosdp4;

		//soft con
		////p3
		triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * from_id + 0, swb*dthetadp3(0)));
		triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * from_id + 1, swb*dthetadp3(1)));
		triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * from_id + 2, swb*dthetadp3(2)));
		//p4
		triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * to_id + 0, swb*dthetadp4(0)));
		triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * to_id + 1, swb*dthetadp4(1)));
		triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * to_id + 2, swb*dthetadp4(2)));
		//p1
		triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * left_id + 0, swb*dthetadp1(0)));
		triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * left_id + 1, swb*dthetadp1(1)));
		triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * left_id + 2, swb*dthetadp1(2)));
		//p2
		triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * right_id + 0, swb*dthetadp2(0)));
		triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * right_id + 1, swb*dthetadp2(1)));
		triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * right_id + 2, swb*dthetadp2(2)));

		/*MyMesh::Normal rij = ri - rj;
		MyMesh::Normal rkj = rk - rj;
		MyMesh::Normal rkl = rk - rl;

		MyMesh::Normal rmj = rij % rkj;
		MyMesh::Normal rnk = rkj % rkl;

		decimal rkj_norm = rkj.length();
		decimal rmj_norm = rmj.length();
		decimal rnk_norm = rnk.length();

		decimal rij_rkj = rij | rkj;
		decimal rkl_rkj = rkl | rkj;

		decimal rkj_sqr_norm = rkj_norm * rkj_norm;

		MyMesh::Normal partial_ri = (rkj_norm / (rmj_norm*rmj_norm)) * rmj;
		MyMesh::Normal partial_rl = -(rkj_norm / (rnk_norm*rnk_norm)) * rnk;
		MyMesh::Normal partial_rj = ((rij_rkj / rkj_sqr_norm - 1) * partial_ri) - ((rkl_rkj / rkj_sqr_norm) * partial_rl);
		MyMesh::Normal partial_rk = ((rkl_rkj / rkj_sqr_norm - 1) * partial_rl) - ((rij_rkj / rkj_sqr_norm) * partial_ri);

		MyMesh::Normal partial_sum = partial_ri + partial_rl + partial_rj + partial_rk;*/

		//set reverse
		//swb = -swb;

		//soft con
		////rj
		//triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * from_id + 0, swb*partial_rj.data()[0]));
		//triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * from_id + 1, swb*partial_rj.data()[1]));
		//triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * from_id + 2, swb*partial_rj.data()[2]));
		////rk
		//triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * to_id + 0, swb*partial_rk.data()[0]));
		//triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * to_id + 1, swb*partial_rk.data()[1]));
		//triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * to_id + 2, swb*partial_rk.data()[2]));
		////rl
		//triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * right_id + 0, swb*partial_rl.data()[0]));
		//triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * right_id + 1, swb*partial_rl.data()[1]));
		//triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * right_id + 2, swb*partial_rl.data()[2]));
		////ri
		//triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * left_id + 0, swb*partial_ri.data()[0]));
		//triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * left_id + 1, swb*partial_ri.data()[1]));
		//triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * left_id + 2, swb*partial_ri.data()[2]));
		//if eg-driven
		if (this->use_eg_or_not == true && this->use_basis_or_not == false)
		{
			for (integer i = 0; i < this->eg_num; i++)
			{
				triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * this->mesh.n_vertices() + i, swb*(eg_dihedral_vec[i][e_id] - ori_dihedral_set[e_id])));
			}
		}
		//if basis exists
		if (this->use_eg_or_not == false && this->use_basis_or_not == true)
		{
			for (integer i = 0; i < this->basis.cols(); i++)
			{
				if (std::abs(/*basis(2 * e_id + 1, i)*/dihedral_basis_vec[i][e_id]) > 0)
				{
					//decimal tmp = swb*basis(2 * e_id + 1, i);
					//triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * this->mesh.n_vertices() + i, tmp));

					decimal tmp = swb * dihedral_basis_vec[i][e_id];
					triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * this->mesh.n_vertices() + 2 * i + 1, tmp));
				}
			}
		}
	}
	if (this->use_sim_or_not == true)
	{
		MyMesh::Point rj(high_v_mat(from_id, 0), high_v_mat(from_id, 1), high_v_mat(from_id, 2));
		MyMesh::Point rk(high_v_mat(to_id, 0), high_v_mat(to_id, 1), high_v_mat(to_id, 2));
		MyMesh::Point ri(high_v_mat(left_id, 0), high_v_mat(left_id, 1), high_v_mat(left_id, 2));
		MyMesh::Point rl(high_v_mat(right_id, 0), high_v_mat(right_id, 1), high_v_mat(right_id, 2));

		decimal cur_angle = CalDiheral(rj, rk, ri, rl);
		fsb(sb_ele_cnt) = swb * (cur_angle - wnt_dihedral_set[e_id]);
		energy = fsb(sb_ele_cnt) * fsb(sb_ele_cnt);
		++sb_ele_cnt;

		MyMesh::Normal rij = ri - rj;
		MyMesh::Normal rkj = rk - rj;
		MyMesh::Normal rkl = rk - rl;

		MyMesh::Normal rmj = rij % rkj;
		MyMesh::Normal rnk = rkj % rkl;

		decimal rkj_norm = rkj.length();
		decimal rmj_norm = rmj.length();
		decimal rnk_norm = rnk.length();

		decimal rij_rkj = rij | rkj;
		decimal rkl_rkj = rkl | rkj;

		decimal rkj_sqr_norm = rkj_norm * rkj_norm;

		MyMesh::Normal partial_ri = (rkj_norm / (rmj_norm*rmj_norm)) * rmj;
		MyMesh::Normal partial_rl = -(rkj_norm / (rnk_norm*rnk_norm)) * rnk;
		MyMesh::Normal partial_rj = ((rij_rkj / rkj_sqr_norm - 1) * partial_ri) - ((rkl_rkj / rkj_sqr_norm) * partial_rl);
		MyMesh::Normal partial_rk = ((rkl_rkj / rkj_sqr_norm - 1) * partial_rl) - ((rij_rkj / rkj_sqr_norm) * partial_ri);

		MyMesh::Normal partial_sum = partial_ri + partial_rl + partial_rj + partial_rk;

		//set reverse
		swb = -swb;

		//soft con
		//rj
		j_triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * from_id + 0, swb*partial_rj.data()[0]));
		j_triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * from_id + 1, swb*partial_rj.data()[1]));
		j_triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * from_id + 2, swb*partial_rj.data()[2]));
		//rk
		j_triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * to_id + 0, swb*partial_rk.data()[0]));
		j_triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * to_id + 1, swb*partial_rk.data()[1]));
		j_triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * to_id + 2, swb*partial_rk.data()[2]));
		//rl
		j_triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * right_id + 0, swb*partial_rl.data()[0]));
		j_triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * right_id + 1, swb*partial_rl.data()[1]));
		j_triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * right_id + 2, swb*partial_rl.data()[2]));
		//ri
		j_triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * left_id + 0, swb*partial_ri.data()[0]));
		j_triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * left_id + 1, swb*partial_ri.data()[1]));
		j_triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * left_id + 2, swb*partial_ri.data()[2]));

		//if eg-driven
		if (this->use_eg_or_not == true && this->use_basis_or_not == false)
		{
			for (integer i = 0; i < this->eg_num; i++)
			{
				triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * this->sim_mesh.n_vertices() + i, swb*(eg_dihedral_vec[i][e_id] - ori_dihedral_set[e_id])));
			}
		}

		//if basis exists
		if (this->use_eg_or_not == false && this->use_basis_or_not == true)
		{
			for (integer i = 0; i < this->basis.cols(); i++)
			{
				if (std::abs(basis(2 * e_id + 1, i)) > 0)
				{
					decimal tmp = -swb * basis(2 * e_id + 1, i);
					triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * this->sim_mesh.n_vertices() + i, tmp));
				}
			}
		}
	}

	++row_cnt;

	return energy;
}

decimal mns::MyDeformationHandler::GetBenTermNew(integer from_id, integer to_id, integer left_id, integer right_id, integer e_id, const Eigen::Vector3d& normal1, const Eigen::Vector3d& normal2, const Eigen::Vector3d& e_ij)
{
	decimal swb = this->swb_set[e_id];
	decimal energy = 0;

	if (this->use_sim_or_not == false)
	{
		//soft con
		Eigen::Vector3d p3(unknown(3 * from_id + 0), unknown(3 * from_id + 1), unknown(3 * from_id + 2));
		Eigen::Vector3d p4(unknown(3 * to_id + 0), unknown(3 * to_id + 1), unknown(3 * to_id + 2));
		Eigen::Vector3d p1(unknown(3 * left_id + 0), unknown(3 * left_id + 1), unknown(3 * left_id + 2));
		Eigen::Vector3d p2(unknown(3 * right_id + 0), unknown(3 * right_id + 1), unknown(3 * right_id + 2));

		/*MyMesh::Point rj(unknown(3 * from_id + 0), unknown(3 * from_id + 1), unknown(3 * from_id + 2));
		MyMesh::Point rk(unknown(3 * to_id + 0), unknown(3 * to_id + 1), unknown(3 * to_id + 2));
		MyMesh::Point ri(unknown(3 * left_id + 0), unknown(3 * left_id + 1), unknown(3 * left_id + 2));
		MyMesh::Point rl(unknown(3 * right_id + 0), unknown(3 * right_id + 1), unknown(3 * right_id + 2));
*/
		decimal n1_norm = normal1.norm();
		decimal n2_norm = normal2.norm();
		decimal e_norm = e_ij.norm();
		Eigen::Vector3d n1 = normal1 / n1_norm;
		Eigen::Vector3d n2 = normal2 / n2_norm;
		Eigen::Vector3d e = e_ij / e_norm;

		decimal cur_angle = CalWLZHDiheral(n1, n2, e);
		fsb(sb_ele_cnt) = swb * (cur_angle - wnt_wlz_dihedral_set[e_id]);
		energy = fsb(sb_ele_cnt) * fsb(sb_ele_cnt);
		if (_isnan(energy))
		{
			cout << swb<<" "<< cur_angle << " " << wnt_wlz_dihedral_set[e_id] << " " << endl;
			std::cout << energy;
		}
		++sb_ele_cnt;

		decimal area1 = CalTriangleArea(p3, p4, p1);
		Eigen::Vector3d p34 = (p3 - p4);
		Eigen::Vector3d p41 = (p4 - p1);
		Eigen::Vector3d p13 = (p1 - p3);
		Eigen::Matrix3d dn1dp1 = -CalTensor(p34) / (2 * area1);
		Eigen::Matrix3d dn1dp2 = Eigen::Matrix3d::Zero();
		Eigen::Matrix3d dn1dp3 = -CalTensor(p41) / (2 * area1);
		Eigen::Matrix3d dn1dp4 = -CalTensor(p13) / (2 * area1);

		decimal area2 = CalTriangleArea(p3, p4, p2); // ?
		Eigen::Vector3d p43 = (p4 - p3);
		Eigen::Vector3d p24 = (p2 - p4);
		Eigen::Vector3d p32 = (p3 - p2);
		Eigen::Matrix3d dn2dp1 = Eigen::Matrix3d::Zero();
		Eigen::Matrix3d dn2dp2 = -CalTensor(p43) / (2 * area2); 
		Eigen::Matrix3d dn2dp3 = -CalTensor(p24) / (2 * area2);
		Eigen::Matrix3d dn2dp4 = -CalTensor(p32) / (2 * area2);

		Eigen::Matrix3d dedp1 = Eigen::Matrix3d::Zero();
		Eigen::Matrix3d dedp2 = Eigen::Matrix3d::Zero();
		Eigen::Matrix3d dedp3 = -Eigen::Matrix3d::Identity() / p34.norm();
		Eigen::Matrix3d dedp4 = Eigen::Matrix3d::Identity() / p34.norm();

		Eigen::Matrix3d n1_tensor = CalTensor(n1);
		Eigen::Matrix3d n2_tensor = CalTensor(n2);
		Eigen::Matrix3d dn1crossn2dp1 = n1_tensor*dn2dp1 - n2_tensor*dn1dp1;
		Eigen::Matrix3d dn1crossn2dp2 = n1_tensor*dn2dp2 - n2_tensor*dn1dp2;
		Eigen::Matrix3d dn1crossn2dp3 = n1_tensor*dn2dp3 - n2_tensor*dn1dp3;
		Eigen::Matrix3d dn1crossn2dp4 = n1_tensor*dn2dp4 - n2_tensor*dn1dp4;

		Eigen::Vector3d dsindp1 = dn1crossn2dp1.transpose()*e + dedp1.transpose()*(n1.cross(n2));
		Eigen::Vector3d dsindp2 = dn1crossn2dp2.transpose()*e + dedp2.transpose()*(n1.cross(n2));
		Eigen::Vector3d dsindp3 = dn1crossn2dp3.transpose()*e + dedp3.transpose()*(n1.cross(n2));
		Eigen::Vector3d dsindp4 = dn1crossn2dp4.transpose()*e + dedp4.transpose()*(n1.cross(n2));

		Eigen::Vector3d dcosdp1 = dn1dp1.transpose()*n2 + dn2dp1.transpose()*n1;
		Eigen::Vector3d dcosdp2 = dn1dp2.transpose()*n2 + dn2dp2.transpose()*n1;
		Eigen::Vector3d dcosdp3 = dn1dp3.transpose()*n2 + dn2dp3.transpose()*n1;
		Eigen::Vector3d dcosdp4 = dn1dp4.transpose()*n2 + dn2dp4.transpose()*n1;

		Eigen::Vector3d dthetadp1 = std::cos(cur_angle) * dsindp1 - std::sin(cur_angle) * dcosdp1;
		Eigen::Vector3d dthetadp2 = std::cos(cur_angle) * dsindp2 - std::sin(cur_angle) * dcosdp2;
		Eigen::Vector3d dthetadp3 = std::cos(cur_angle) * dsindp3 - std::sin(cur_angle) * dcosdp3;
		Eigen::Vector3d dthetadp4 = std::cos(cur_angle) * dsindp4 - std::sin(cur_angle) * dcosdp4;

		//soft con
		////p3
		triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * from_id + 0, swb*dthetadp3(0)));
		triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * from_id + 1, swb*dthetadp3(1)));
		triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * from_id + 2, swb*dthetadp3(2)));
		//p4
		triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * to_id + 0, swb*dthetadp4(0)));
		triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * to_id + 1, swb*dthetadp4(1)));
		triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * to_id + 2, swb*dthetadp4(2)));
		//p1
		triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * left_id + 0, swb*dthetadp1(0)));
		triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * left_id + 1, swb*dthetadp1(1)));
		triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * left_id + 2, swb*dthetadp1(2)));
		//p2
		triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * right_id + 0, swb*dthetadp2(0)));
		triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * right_id + 1, swb*dthetadp2(1)));
		triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * right_id + 2, swb*dthetadp2(2)));

		/*MyMesh::Normal rij = ri - rj;
		MyMesh::Normal rkj = rk - rj;
		MyMesh::Normal rkl = rk - rl;

		MyMesh::Normal rmj = rij % rkj;
		MyMesh::Normal rnk = rkj % rkl;

		decimal rkj_norm = rkj.length();
		decimal rmj_norm = rmj.length();
		decimal rnk_norm = rnk.length();

		decimal rij_rkj = rij | rkj;
		decimal rkl_rkj = rkl | rkj;

		decimal rkj_sqr_norm = rkj_norm * rkj_norm;

		MyMesh::Normal partial_ri = (rkj_norm / (rmj_norm*rmj_norm)) * rmj;
		MyMesh::Normal partial_rl = -(rkj_norm / (rnk_norm*rnk_norm)) * rnk;
		MyMesh::Normal partial_rj = ((rij_rkj / rkj_sqr_norm - 1) * partial_ri) - ((rkl_rkj / rkj_sqr_norm) * partial_rl);
		MyMesh::Normal partial_rk = ((rkl_rkj / rkj_sqr_norm - 1) * partial_rl) - ((rij_rkj / rkj_sqr_norm) * partial_ri);

		MyMesh::Normal partial_sum = partial_ri + partial_rl + partial_rj + partial_rk;*/

		//set reverse
		//swb = -swb;

		//soft con
		////rj
		//triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * from_id + 0, swb*partial_rj.data()[0]));
		//triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * from_id + 1, swb*partial_rj.data()[1]));
		//triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * from_id + 2, swb*partial_rj.data()[2]));
		////rk
		//triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * to_id + 0, swb*partial_rk.data()[0]));
		//triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * to_id + 1, swb*partial_rk.data()[1]));
		//triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * to_id + 2, swb*partial_rk.data()[2]));
		////rl
		//triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * right_id + 0, swb*partial_rl.data()[0]));
		//triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * right_id + 1, swb*partial_rl.data()[1]));
		//triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * right_id + 2, swb*partial_rl.data()[2]));
		////ri
		//triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * left_id + 0, swb*partial_ri.data()[0]));
		//triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * left_id + 1, swb*partial_ri.data()[1]));
		//triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * left_id + 2, swb*partial_ri.data()[2]));
		//if eg-driven
		if (this->use_eg_or_not == true && this->use_basis_or_not == false)
		{
			for (integer i = 0; i < this->eg_num; i++)
			{
				triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * this->mesh.n_vertices() + i, swb*(eg_dihedral_vec[i][e_id] - ori_dihedral_set[e_id])));
			}
		}
		//if basis exists
		if (this->use_eg_or_not == false && this->use_basis_or_not == true)
		{
			for (integer i = 0; i < this->basis.cols(); i++)
			{
				if (std::abs(/*basis(2 * e_id + 1, i)*/dihedral_basis_vec[i][e_id]) > 0)
				{
					//decimal tmp = swb*basis(2 * e_id + 1, i);
					//triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * this->mesh.n_vertices() + i, tmp));
					
					decimal tmp =   swb*dihedral_basis_vec[i][e_id];
					triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * this->mesh.n_vertices() + 2 * i + 1, tmp));
				}
			}
		}
	}
	if (this->use_sim_or_not == true)
	{
		MyMesh::Point rj(high_v_mat(from_id, 0), high_v_mat(from_id, 1), high_v_mat(from_id, 2));
		MyMesh::Point rk(high_v_mat(to_id, 0), high_v_mat(to_id, 1), high_v_mat(to_id, 2));
		MyMesh::Point ri(high_v_mat(left_id, 0), high_v_mat(left_id, 1), high_v_mat(left_id, 2));
		MyMesh::Point rl(high_v_mat(right_id, 0), high_v_mat(right_id, 1), high_v_mat(right_id, 2));

		decimal cur_angle = CalDiheral(rj, rk, ri, rl);
		fsb(sb_ele_cnt) = swb * (cur_angle - wnt_dihedral_set[e_id]);
		energy = fsb(sb_ele_cnt) * fsb(sb_ele_cnt);
		++sb_ele_cnt;

		MyMesh::Normal rij = ri - rj;
		MyMesh::Normal rkj = rk - rj;
		MyMesh::Normal rkl = rk - rl;

		MyMesh::Normal rmj = rij % rkj;
		MyMesh::Normal rnk = rkj % rkl;

		decimal rkj_norm = rkj.length();
		decimal rmj_norm = rmj.length();
		decimal rnk_norm = rnk.length();

		decimal rij_rkj = rij | rkj;
		decimal rkl_rkj = rkl | rkj;

		decimal rkj_sqr_norm = rkj_norm * rkj_norm;

		MyMesh::Normal partial_ri = (rkj_norm / (rmj_norm*rmj_norm)) * rmj;
		MyMesh::Normal partial_rl = -(rkj_norm / (rnk_norm*rnk_norm)) * rnk;
		MyMesh::Normal partial_rj = ((rij_rkj / rkj_sqr_norm - 1) * partial_ri) - ((rkl_rkj / rkj_sqr_norm) * partial_rl);
		MyMesh::Normal partial_rk = ((rkl_rkj / rkj_sqr_norm - 1) * partial_rl) - ((rij_rkj / rkj_sqr_norm) * partial_ri);

		MyMesh::Normal partial_sum = partial_ri + partial_rl + partial_rj + partial_rk;

		//set reverse
		swb = -swb;

		//soft con
		//rj
		j_triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * from_id + 0, swb*partial_rj.data()[0]));
		j_triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * from_id + 1, swb*partial_rj.data()[1]));
		j_triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * from_id + 2, swb*partial_rj.data()[2]));
		//rk
		j_triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * to_id + 0, swb*partial_rk.data()[0]));
		j_triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * to_id + 1, swb*partial_rk.data()[1]));
		j_triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * to_id + 2, swb*partial_rk.data()[2]));
		//rl
		j_triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * right_id + 0, swb*partial_rl.data()[0]));
		j_triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * right_id + 1, swb*partial_rl.data()[1]));
		j_triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * right_id + 2, swb*partial_rl.data()[2]));
		//ri
		j_triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * left_id + 0, swb*partial_ri.data()[0]));
		j_triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * left_id + 1, swb*partial_ri.data()[1]));
		j_triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * left_id + 2, swb*partial_ri.data()[2]));

		//if eg-driven
		if (this->use_eg_or_not == true && this->use_basis_or_not == false)
		{
			for (integer i = 0; i < this->eg_num; i++)
			{
				triple.push_back(Eigen::Triplet<decimal>(row_cnt, 3 * this->sim_mesh.n_vertices() + i, swb*(eg_dihedral_vec[i][e_id] - ori_dihedral_set[e_id])));
			}
		}

		//if basis exists
		if (this->use_eg_or_not == false && this->use_basis_or_not == true)
		{
			for (integer i = 0; i < this->basis.cols(); i++)
			{
				if (std::abs(basis(2 * e_id + 1, i)) > 0)
				{
					decimal tmp = - swb*basis(2 * e_id + 1, i);
					triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * this->sim_mesh.n_vertices() + i, tmp));
				}
			}
		}
	}

	++row_cnt;

	return energy;
}

decimal mns::MyDeformationHandler::EbNew()
{
	if (fb.rows() != this->mesh.n_edges())
	{
		fb.resize(this->mesh.n_edges());
	}
	fb = Eigen::VectorXd::Zero(this->mesh.n_edges());
	decimal total_energy = 0;
	int i = 0;
	for (auto he_it = this->mesh.halfedges_begin(); he_it != this->mesh.halfedges_end(); ++i)
	{
		auto e = *he_it;
		auto f1 = this->mesh.face_handle(e);
		auto op_he = this->mesh.opposite_halfedge_handle(e);
		auto f2 = this->mesh.face_handle(op_he);
		auto fromVertex = this->mesh.from_vertex_handle(e);
		auto toVertex = this->mesh.to_vertex_handle(e);
		if (f1.idx() == -1 || f2.idx() == -1)
		{
			++sb_ele_cnt;
			++row_cnt;

			//edge, not half edge
			++he_it;
			++he_it;
			continue;
		}

		integer leftVertex, rightVertex;
		//integer temp_cnt = 0;
		for (auto fv_it1 = this->mesh.fv_begin(f1); fv_it1 != this->mesh.fv_end(f1); ++fv_it1)
		{
			auto v = *fv_it1;
			//temp_cnt++;
			if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
			{
				leftVertex = v.idx();
				break;
			}
		}
		//temp_cnt = 0;
		for (auto fv_it2 = this->mesh.fv_begin(f2); fv_it2 != this->mesh.fv_end(f2); ++fv_it2)
		{
			auto v = *fv_it2;
			//temp_cnt++;
			if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
			{
				rightVertex = v.idx();
				break;
			}
		}
		total_energy += GetBenTerm(fromVertex.idx(), toVertex.idx(), leftVertex, rightVertex, i);
		if (_isnan(total_energy))
		{
			std::cout << i;
		}
		++he_it;
		++he_it;
	}

	return total_energy;
}

std::vector<decimal> mns::MyDeformationHandler::EsbNew1()
{
	std::vector<decimal> energy;
	energy.push_back(0);
	energy.push_back(0);

	if (fsb.rows() != 2 * this->mesh.n_edges())
	{
		fsb.resize(2 * this->mesh.n_edges());
	}
	fsb = Eigen::VectorXd::Zero(2 * this->mesh.n_edges());
	if (alpha_b == 0 || alpha_s == 0)
	{
		fsb.resize(this->mesh.n_edges());
		fsb = Eigen::VectorXd::Zero(this->mesh.n_edges());
	}

	int i = 0;
	for (auto he_it = this->mesh.halfedges_begin(); he_it != this->mesh.halfedges_end(); ++i)
	{
		auto e = *he_it;
		auto f1 = this->mesh.face_handle(e);
		auto op_he = this->mesh.opposite_halfedge_handle(e);
		auto f2 = this->mesh.face_handle(op_he);
		auto fromVertex = this->mesh.from_vertex_handle(e);
		auto toVertex = this->mesh.to_vertex_handle(e);
		integer from_idx = fromVertex.idx();
		integer to_idx = toVertex.idx();

		if (alpha_s != 0)
		{
			energy[0] += GetStrTerm1(fromVertex.idx(), toVertex.idx(), i);
		}

		if (alpha_b != 0)
		{
			if (f1.idx() == -1 || f2.idx() == -1)
			{
				++sb_ele_cnt;
				++row_cnt;

				//edge, not half edge
				++he_it;
				++he_it;
				continue;
			}

			integer leftVertex, rightVertex;
			//integer temp_cnt = 0;
			for (auto fv_it1 = this->mesh.fv_begin(f1); fv_it1 != this->mesh.fv_end(f1); ++fv_it1)
			{
				auto v = *fv_it1;
				//temp_cnt++;
				if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
				{
					leftVertex = v.idx();
					break;
				}
			}
			//temp_cnt = 0;
			for (auto fv_it2 = this->mesh.fv_begin(f2); fv_it2 != this->mesh.fv_end(f2); ++fv_it2)
			{
				auto v = *fv_it2;
				//temp_cnt++;
				if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
				{
					rightVertex = v.idx();
					break;
				}
			}

			//normal computation
			Eigen::Vector3d normal1;
			Eigen::Vector3d normal2;
			{
				auto fhe_it1 = this->mesh.fh_begin(f1);
				auto he1 = *fhe_it1;
				auto he1idx = he1.idx();
				auto fromVertex1 = this->mesh.from_vertex_handle(he1);
				auto toVertex1 = this->mesh.to_vertex_handle(he1);
				integer from1_idx = fromVertex1.idx();
				integer to1_idx = toVertex1.idx();
				Eigen::Vector3d fromPt1(unknown(3 * from1_idx + 0), unknown(3 * from1_idx + 1), unknown(3 * from1_idx + 2));
				Eigen::Vector3d toPt1(unknown(3 * to1_idx + 0), unknown(3 * to1_idx + 1), unknown(3 * to1_idx + 2));
				auto e1 = toPt1 - fromPt1;
				Eigen::Vector3d ev1;
				ev1 << e1[0], e1[1], e1[2];
				++fhe_it1;
				auto he2 = *fhe_it1;
				auto he2idx = he2.idx();
				auto fromVertex2 = this->mesh.from_vertex_handle(he2);
				auto toVertex2 = this->mesh.to_vertex_handle(he2);
				integer from2_idx = fromVertex2.idx();
				integer to2_idx = toVertex2.idx();
				Eigen::Vector3d fromPt2(unknown(3 * from2_idx + 0), unknown(3 * from2_idx + 1), unknown(3 * from2_idx + 2));
				Eigen::Vector3d toPt2(unknown(3 * to2_idx + 0), unknown(3 * to2_idx + 1), unknown(3 * to2_idx + 2));
				auto e2 = toPt2 - fromPt2;
				Eigen::Vector3d ev2;
				ev2 << e2[0], e2[1], e2[2];
				normal1 = ev1.cross(ev2);
				//normal1 = normal1 / normal1.norm();
			}
			{
				auto fhe_it2 = this->mesh.fh_begin(f2);
				auto he1 = *fhe_it2;
				auto he1idx = he1.idx();
				auto fromVertex1 = this->mesh.from_vertex_handle(he1);
				auto toVertex1 = this->mesh.to_vertex_handle(he1);
				integer from1_idx = fromVertex1.idx();
				integer to1_idx = toVertex1.idx();
				Eigen::Vector3d fromPt1(unknown(3 * from1_idx + 0), unknown(3 * from1_idx + 1), unknown(3 * from1_idx + 2));
				Eigen::Vector3d toPt1(unknown(3 * to1_idx + 0), unknown(3 * to1_idx + 1), unknown(3 * to1_idx + 2));
				auto e1 = toPt1 - fromPt1;
				Eigen::Vector3d ev1;
				ev1 << e1[0], e1[1], e1[2];
				++fhe_it2;
				auto he2 = *fhe_it2;
				auto he2idx = he2.idx();
				auto fromVertex2 = this->mesh.from_vertex_handle(he2);
				auto toVertex2 = this->mesh.to_vertex_handle(he2);
				integer from2_idx = fromVertex2.idx();
				integer to2_idx = toVertex2.idx();
				Eigen::Vector3d fromPt2(unknown(3 * from2_idx + 0), unknown(3 * from2_idx + 1), unknown(3 * from2_idx + 2));
				Eigen::Vector3d toPt2(unknown(3 * to2_idx + 0), unknown(3 * to2_idx + 1), unknown(3 * to2_idx + 2));
				auto e2 = toPt2 - fromPt2;
				Eigen::Vector3d ev2;
				ev2 << e2[0], e2[1], e2[2];
				normal2 = ev1.cross(ev2);
				//normal2 = normal2 / normal2.norm();
			}

			Eigen::Vector3d fromPt(unknown(3 * from_idx + 0), unknown(3 * from_idx + 1), unknown(3 * from_idx + 2));
			Eigen::Vector3d toPt(unknown(3 * to_idx + 0), unknown(3 * to_idx + 1), unknown(3 * to_idx + 2));
			auto e_ij_tmp = toPt - fromPt;
			Eigen::Vector3d e_ij;
			e_ij << e_ij_tmp[0], e_ij_tmp[1], e_ij_tmp[2];
			//e_ij = e_ij / (e_ij.norm());

			energy[1] += GetBenTermNew(fromVertex.idx(), toVertex.idx(), leftVertex, rightVertex, i, normal1, normal2, e_ij);

		}

		if (_isnan(energy[1]))
		{
			std::cout << i;
		}

		++he_it;
		++he_it;
	}

	return energy;
}

std::vector<decimal> mns::MyDeformationHandler::EsbNew()
{
	std::vector<decimal> energy;
	energy.push_back(0);
	energy.push_back(0);

	if (fsb.rows() != 2 * this->mesh.n_edges())
	{
		fsb.resize(2 * this->mesh.n_edges());
	}
	fsb = Eigen::VectorXd::Zero(2 * this->mesh.n_edges());
	if (alpha_b == 0 || alpha_s == 0)
	{
		fsb.resize(this->mesh.n_edges());
		fsb = Eigen::VectorXd::Zero(this->mesh.n_edges());
	}

	int i = 0;
	for (auto he_it = this->mesh.halfedges_begin(); he_it != this->mesh.halfedges_end(); ++i)
	{
		auto e = *he_it;
		auto f1 = this->mesh.face_handle(e);
		auto op_he = this->mesh.opposite_halfedge_handle(e);
		auto f2 = this->mesh.face_handle(op_he);
		auto fromVertex = this->mesh.from_vertex_handle(e);
		auto toVertex = this->mesh.to_vertex_handle(e);
		integer from_idx = fromVertex.idx();
		integer to_idx = toVertex.idx();

		if (alpha_s != 0)
		{
			energy[0] += GetStrTerm(fromVertex.idx(), toVertex.idx(), i);
		}

		if (alpha_b != 0)
		{
			if (f1.idx() == -1 || f2.idx() == -1)
			{
				++sb_ele_cnt;
				++row_cnt;

				//edge, not half edge
				++he_it;
				++he_it;
				continue;
			}

			integer leftVertex, rightVertex;
			//integer temp_cnt = 0;
			for (auto fv_it1 = this->mesh.fv_begin(f1); fv_it1 != this->mesh.fv_end(f1); ++fv_it1)
			{
				auto v = *fv_it1;
				//temp_cnt++;
				if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
				{
					leftVertex = v.idx();
					break;
				}
			}
			//temp_cnt = 0;
			for (auto fv_it2 = this->mesh.fv_begin(f2); fv_it2 != this->mesh.fv_end(f2); ++fv_it2)
			{
				auto v = *fv_it2;
				//temp_cnt++;
				if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
				{
					rightVertex = v.idx();
					break;
				}
			}

			//normal computation
			Eigen::Vector3d normal1;
			Eigen::Vector3d normal2;
			{
				auto fhe_it1 = this->mesh.fh_begin(f1);
				auto he1 = *fhe_it1;
				auto he1idx = he1.idx();
				auto fromVertex1 = this->mesh.from_vertex_handle(he1);
				auto toVertex1 = this->mesh.to_vertex_handle(he1);
				integer from1_idx = fromVertex1.idx();
				integer to1_idx = toVertex1.idx();
				Eigen::Vector3d fromPt1(unknown(3 * from1_idx + 0), unknown(3 * from1_idx + 1), unknown(3 * from1_idx + 2));
				Eigen::Vector3d toPt1(unknown(3 * to1_idx + 0), unknown(3 * to1_idx + 1), unknown(3 * to1_idx + 2));
				auto e1 = toPt1 - fromPt1;
				Eigen::Vector3d ev1;
				ev1 << e1[0], e1[1], e1[2];
				++fhe_it1;
				auto he2 = *fhe_it1;
				auto he2idx = he2.idx();
				auto fromVertex2 = this->mesh.from_vertex_handle(he2);
				auto toVertex2 = this->mesh.to_vertex_handle(he2);
				integer from2_idx = fromVertex2.idx();
				integer to2_idx = toVertex2.idx();
				Eigen::Vector3d fromPt2(unknown(3 * from2_idx + 0), unknown(3 * from2_idx + 1), unknown(3 * from2_idx + 2));
				Eigen::Vector3d toPt2(unknown(3 * to2_idx + 0), unknown(3 * to2_idx + 1), unknown(3 * to2_idx + 2));
				auto e2 = toPt2 - fromPt2;
				Eigen::Vector3d ev2;
				ev2 << e2[0], e2[1], e2[2];
				normal1 = ev1.cross(ev2);
				//normal1 = normal1 / normal1.norm();
			}
			{
				auto fhe_it2 = this->mesh.fh_begin(f2);
				auto he1 = *fhe_it2;
				auto he1idx = he1.idx();
				auto fromVertex1 = this->mesh.from_vertex_handle(he1);
				auto toVertex1 = this->mesh.to_vertex_handle(he1);
				integer from1_idx = fromVertex1.idx();
				integer to1_idx = toVertex1.idx();
				Eigen::Vector3d fromPt1(unknown(3 * from1_idx + 0), unknown(3 * from1_idx + 1), unknown(3 * from1_idx + 2));
				Eigen::Vector3d toPt1(unknown(3 * to1_idx + 0), unknown(3 * to1_idx + 1), unknown(3 * to1_idx + 2));
				auto e1 = toPt1 - fromPt1;
				Eigen::Vector3d ev1;
				ev1 << e1[0], e1[1], e1[2];
				++fhe_it2;
				auto he2 = *fhe_it2;
				auto he2idx = he2.idx();
				auto fromVertex2 = this->mesh.from_vertex_handle(he2);
				auto toVertex2 = this->mesh.to_vertex_handle(he2);
				integer from2_idx = fromVertex2.idx();
				integer to2_idx = toVertex2.idx();
				Eigen::Vector3d fromPt2(unknown(3 * from2_idx + 0), unknown(3 * from2_idx + 1), unknown(3 * from2_idx + 2));
				Eigen::Vector3d toPt2(unknown(3 * to2_idx + 0), unknown(3 * to2_idx + 1), unknown(3 * to2_idx + 2));
				auto e2 = toPt2 - fromPt2;
				Eigen::Vector3d ev2;
				ev2 << e2[0], e2[1], e2[2];
				normal2 = ev1.cross(ev2);
				//normal2 = normal2 / normal2.norm();
			}

			Eigen::Vector3d fromPt(unknown(3 * from_idx + 0), unknown(3 * from_idx + 1), unknown(3 * from_idx + 2));
			Eigen::Vector3d toPt(unknown(3 * to_idx + 0), unknown(3 * to_idx + 1), unknown(3 * to_idx + 2));
			auto e_ij_tmp = toPt - fromPt;
			Eigen::Vector3d e_ij;
			e_ij << e_ij_tmp[0], e_ij_tmp[1], e_ij_tmp[2];
			//e_ij = e_ij / (e_ij.norm());

			energy[1] += GetBenTermNew(fromVertex.idx(), toVertex.idx(), leftVertex, rightVertex, i, normal1, normal2, e_ij);
			
		}

		if (_isnan(energy[1]))
		{
			std::cout << i;
		}

		++he_it;
		++he_it;
	}

	return energy;
}

decimal mns::MyDeformationHandler::Ev()
{
	decimal total_energy = 0;
	decimal total_tmp_vol = 0;
	fv = 0;
	//u_vec = Eigen::VectorXd::Zero(unknown.size());
	if (this->use_sim_or_not == false)
	{
		u_vec = Eigen::VectorXd::Zero(unknown.size());
		for (auto f_it = this->mesh.faces_begin(); f_it != this->mesh.faces_end(); ++f_it)
		{
			auto face = *f_it;
			pt_vector pt_vec;
			int_vector idx_vec;
			for (auto fv_it = this->mesh.fv_begin(face); fv_it != this->mesh.fv_end(face); ++fv_it)
			{
				auto v = *fv_it;
				//soft con
				MyMesh::Point pt(unknown(3 * v.idx() + 0), unknown(3 * v.idx() + 1), unknown(3 * v.idx() + 2));

				//hard con
				/*MyMesh::Point pt;
				if (idx_move[v.idx()] == -1)
				{
				pt[0] = cps.set[con_ori_map[v.idx()]].x;
				pt[1] = cps.set[con_ori_map[v.idx()]].y;
				pt[2] = cps.set[con_ori_map[v.idx()]].z;
				}
				else
				{
				pt[0] = unknown[3 * (v.idx() - idx_move[v.idx()]) + 0];
				pt[1] = unknown[3 * (v.idx() - idx_move[v.idx()]) + 1];
				pt[2] = unknown[3 * (v.idx() - idx_move[v.idx()]) + 2];
				}*/
				idx_vec.push_back(v.idx());
				pt_vec.push_back(pt);
			}
			decimal temp_vol = (1.0 / 6.0) * ((pt_vec[0] % pt_vec[1]) | pt_vec[2]);
			total_tmp_vol += temp_vol;

			//soft con
			//xi
			u_vec(3 * idx_vec[0] + 0) += -(1.0 / 6.0) * swv*(pt_vec[2].data()[1] * pt_vec[1].data()[2] - pt_vec[2].data()[2] * pt_vec[1].data()[1]);
			u_vec(3 * idx_vec[0] + 1) += (1.0 / 6.0) * swv*(pt_vec[2].data()[0] * pt_vec[1].data()[2] - pt_vec[2].data()[2] * pt_vec[1].data()[0]);
			u_vec(3 * idx_vec[0] + 2) += -(1.0 / 6.0) * swv*(pt_vec[2].data()[0] * pt_vec[1].data()[1] - pt_vec[2].data()[1] * pt_vec[1].data()[0]);
			//xj
			u_vec(3 * idx_vec[1] + 0) += (1.0 / 6.0) * swv*(pt_vec[2].data()[1] * pt_vec[0].data()[2] - pt_vec[2].data()[2] * pt_vec[0].data()[1]);
			u_vec(3 * idx_vec[1] + 1) += -(1.0 / 6.0) * swv*(pt_vec[2].data()[0] * pt_vec[0].data()[2] - pt_vec[2].data()[2] * pt_vec[0].data()[0]);
			u_vec(3 * idx_vec[1] + 2) += (1.0 / 6.0) * swv*(pt_vec[2].data()[0] * pt_vec[0].data()[1] - pt_vec[2].data()[1] * pt_vec[0].data()[0]);
			//xk
			u_vec(3 * idx_vec[2] + 0) += (1.0 / 6.0) * swv*(pt_vec[0].data()[1] * pt_vec[1].data()[2] - pt_vec[1].data()[1] * pt_vec[0].data()[2]);
			u_vec(3 * idx_vec[2] + 1) += -(1.0 / 6.0) * swv*(pt_vec[0].data()[0] * pt_vec[1].data()[2] - pt_vec[1].data()[0] * pt_vec[0].data()[2]);
			u_vec(3 * idx_vec[2] + 2) += (1.0 / 6.0) * swv*(pt_vec[0].data()[0] * pt_vec[1].data()[1] - pt_vec[1].data()[0] * pt_vec[0].data()[1]);

			//hard con
			////xi
			//if (idx_move[idx_vec[0]] != -1)
			//{
			//	u_vec(3 * (idx_vec[0] - idx_move[idx_vec[0]]) + 0) += -(1.0 / 6.0) * swv*(pt_vec[2].data()[1] * pt_vec[1].data()[2] - pt_vec[2].data()[2] * pt_vec[1].data()[1]);
			//	u_vec(3 * (idx_vec[0] - idx_move[idx_vec[0]]) + 1) += (1.0 / 6.0) * swv*(pt_vec[2].data()[0] * pt_vec[1].data()[2] - pt_vec[2].data()[2] * pt_vec[1].data()[0]);
			//	u_vec(3 * (idx_vec[0] - idx_move[idx_vec[0]]) + 2) += -(1.0 / 6.0) * swv*(pt_vec[2].data()[0] * pt_vec[1].data()[1] - pt_vec[2].data()[1] * pt_vec[1].data()[0]);
			//}		

			////xj
			//if (idx_move[idx_vec[1]] != -1)
			//{
			//	u_vec(3 * (idx_vec[1] - idx_move[idx_vec[1]]) + 0) += (1.0 / 6.0) * swv*(pt_vec[2].data()[1] * pt_vec[0].data()[2] - pt_vec[2].data()[2] * pt_vec[0].data()[1]);
			//	u_vec(3 * (idx_vec[1] - idx_move[idx_vec[1]]) + 1) += -(1.0 / 6.0) * swv*(pt_vec[2].data()[0] * pt_vec[0].data()[2] - pt_vec[2].data()[2] * pt_vec[0].data()[0]);
			//	u_vec(3 * (idx_vec[1] - idx_move[idx_vec[1]]) + 2) += (1.0 / 6.0) * swv*(pt_vec[2].data()[0] * pt_vec[0].data()[1] - pt_vec[2].data()[1] * pt_vec[0].data()[0]);
			//}

			////xk
			//if (idx_move[idx_vec[2]] != -1)
			//{
			//	u_vec(3 * (idx_vec[2] - idx_move[idx_vec[2]]) + 0) += (1.0 / 6.0) * swv*(pt_vec[0].data()[1] * pt_vec[1].data()[2] - pt_vec[1].data()[1] * pt_vec[0].data()[2]);
			//	u_vec(3 * (idx_vec[2] - idx_move[idx_vec[2]]) + 1) += -(1.0 / 6.0) * swv*(pt_vec[0].data()[0] * pt_vec[1].data()[2] - pt_vec[1].data()[0] * pt_vec[0].data()[2]);
			//	u_vec(3 * (idx_vec[2] - idx_move[idx_vec[2]]) + 2) += (1.0 / 6.0) * swv*(pt_vec[0].data()[0] * pt_vec[1].data()[1] - pt_vec[1].data()[0] * pt_vec[0].data()[1]);
			//}
			//
		}
	}
	if (this->use_sim_or_not == true)
	{
		if (this->use_eg_or_not == true)
		{
			u_vec = Eigen::VectorXd::Zero(3 * this->mesh.n_vertices() + this->eg_num);
		}
		else if (this->use_basis_or_not == true)
		{
			u_vec = Eigen::VectorXd::Zero(3 * this->mesh.n_vertices() + 2 * this->basis.cols());
		}
		else
		{
			u_vec = Eigen::VectorXd::Zero(3 * this->mesh.n_vertices());
		}

		for (auto f_it = this->mesh.faces_begin(); f_it != this->mesh.faces_end(); ++f_it)
		{
			auto face = *f_it;
			pt_vector pt_vec;
			int_vector idx_vec;
			for (auto fv_it = this->mesh.fv_begin(face); fv_it != this->mesh.fv_end(face); ++fv_it)
			{
				auto v = *fv_it;
				//soft con
				/*decimal tmp_x = 0;
				decimal tmp_y = 0;
				decimal tmp_z = 0;
				for (integer j = 0; j < this->sim_mesh.n_vertices(); j++)
				{
				decimal wij = W_mat(v.idx(), j);
				tmp_x += wij * unknown(3 * j + 0);
				tmp_y += wij * unknown(3 * j + 1);
				tmp_z += wij * unknown(3 * j + 2);
				}*/
				MyMesh::Point pt(high_v_mat(v.idx(), 0), high_v_mat(v.idx(), 1), high_v_mat(v.idx(), 2));

				//hard con
				/*MyMesh::Point pt;
				if (idx_move[v.idx()] == -1)
				{
				pt[0] = cps.set[con_ori_map[v.idx()]].x;
				pt[1] = cps.set[con_ori_map[v.idx()]].y;
				pt[2] = cps.set[con_ori_map[v.idx()]].z;
				}
				else
				{
				pt[0] = unknown[3 * (v.idx() - idx_move[v.idx()]) + 0];
				pt[1] = unknown[3 * (v.idx() - idx_move[v.idx()]) + 1];
				pt[2] = unknown[3 * (v.idx() - idx_move[v.idx()]) + 2];
				}*/
				idx_vec.push_back(v.idx());
				pt_vec.push_back(pt);
			}
			decimal temp_vol = (1.0 / 6.0) * ((pt_vec[0] % pt_vec[1]) | pt_vec[2]);
			total_tmp_vol += temp_vol;

			//soft con
			//xi
			u_vec(3 * idx_vec[0] + 0) += -(1.0 / 6.0) * swv*(pt_vec[2].data()[1] * pt_vec[1].data()[2] - pt_vec[2].data()[2] * pt_vec[1].data()[1]);
			u_vec(3 * idx_vec[0] + 1) += (1.0 / 6.0) * swv*(pt_vec[2].data()[0] * pt_vec[1].data()[2] - pt_vec[2].data()[2] * pt_vec[1].data()[0]);
			u_vec(3 * idx_vec[0] + 2) += -(1.0 / 6.0) * swv*(pt_vec[2].data()[0] * pt_vec[1].data()[1] - pt_vec[2].data()[1] * pt_vec[1].data()[0]);
			//xj
			u_vec(3 * idx_vec[1] + 0) += (1.0 / 6.0) * swv*(pt_vec[2].data()[1] * pt_vec[0].data()[2] - pt_vec[2].data()[2] * pt_vec[0].data()[1]);
			u_vec(3 * idx_vec[1] + 1) += -(1.0 / 6.0) * swv*(pt_vec[2].data()[0] * pt_vec[0].data()[2] - pt_vec[2].data()[2] * pt_vec[0].data()[0]);
			u_vec(3 * idx_vec[1] + 2) += (1.0 / 6.0) * swv*(pt_vec[2].data()[0] * pt_vec[0].data()[1] - pt_vec[2].data()[1] * pt_vec[0].data()[0]);
			//xk
			u_vec(3 * idx_vec[2] + 0) += (1.0 / 6.0) * swv*(pt_vec[0].data()[1] * pt_vec[1].data()[2] - pt_vec[1].data()[1] * pt_vec[0].data()[2]);
			u_vec(3 * idx_vec[2] + 1) += -(1.0 / 6.0) * swv*(pt_vec[0].data()[0] * pt_vec[1].data()[2] - pt_vec[1].data()[0] * pt_vec[0].data()[2]);
			u_vec(3 * idx_vec[2] + 2) += (1.0 / 6.0) * swv*(pt_vec[0].data()[0] * pt_vec[1].data()[1] - pt_vec[1].data()[0] * pt_vec[0].data()[1]);

			//xi
			/*decimal i0 = -(1.0 / 6.0) * swv*(pt_vec[2].data()[1] * pt_vec[1].data()[2] - pt_vec[2].data()[2] * pt_vec[1].data()[1]);
			decimal i1 = (1.0 / 6.0) * swv*(pt_vec[2].data()[0] * pt_vec[1].data()[2] - pt_vec[2].data()[2] * pt_vec[1].data()[0]);
			decimal i2 = -(1.0 / 6.0) * swv*(pt_vec[2].data()[0] * pt_vec[1].data()[1] - pt_vec[2].data()[1] * pt_vec[1].data()[0]);*/
			//xj
			/*decimal j0 = (1.0 / 6.0) * swv*(pt_vec[2].data()[1] * pt_vec[0].data()[2] - pt_vec[2].data()[2] * pt_vec[0].data()[1]);
			decimal j1 = -(1.0 / 6.0) * swv*(pt_vec[2].data()[0] * pt_vec[0].data()[2] - pt_vec[2].data()[2] * pt_vec[0].data()[0]);
			decimal j2 = (1.0 / 6.0) * swv*(pt_vec[2].data()[0] * pt_vec[0].data()[1] - pt_vec[2].data()[1] * pt_vec[0].data()[0]);*/
			//xk
			/*decimal k0 = (1.0 / 6.0) * swv*(pt_vec[0].data()[1] * pt_vec[1].data()[2] - pt_vec[1].data()[1] * pt_vec[0].data()[2]);
			decimal k1 = -(1.0 / 6.0) * swv*(pt_vec[0].data()[0] * pt_vec[1].data()[2] - pt_vec[1].data()[0] * pt_vec[0].data()[2]);
			decimal k2 = (1.0 / 6.0) * swv*(pt_vec[0].data()[0] * pt_vec[1].data()[1] - pt_vec[1].data()[0] * pt_vec[0].data()[1]);*/

			/*for (integer j = 0; j < this->W_idx_set[idx_vec[0]].size(); j++)
			{
			decimal w_0_j = W_mat(idx_vec[0], W_idx_set[idx_vec[0]][j]);
			u_vec(3 * W_idx_set[idx_vec[0]][j] + 0) += i0 * w_0_j;
			u_vec(3 * W_idx_set[idx_vec[0]][j] + 1) += i1 * w_0_j;
			u_vec(3 * W_idx_set[idx_vec[0]][j] + 2) += i2 * w_0_j;
			}
			for (integer j = 0; j < this->W_idx_set[idx_vec[1]].size(); j++)
			{
			decimal w_1_j = W_mat(idx_vec[1], W_idx_set[idx_vec[1]][j]);
			u_vec(3 * W_idx_set[idx_vec[1]][j] + 0) += j0 * w_1_j;
			u_vec(3 * W_idx_set[idx_vec[1]][j] + 1) += j1 * w_1_j;
			u_vec(3 * W_idx_set[idx_vec[1]][j] + 2) += j2 * w_1_j;
			}
			for (integer j = 0; j < this->W_idx_set[idx_vec[2]].size(); j++)
			{
			decimal w_2_j = W_mat(idx_vec[2], W_idx_set[idx_vec[2]][j]);
			u_vec(3 * W_idx_set[idx_vec[2]][j] + 0) += k0 * w_2_j;
			u_vec(3 * W_idx_set[idx_vec[2]][j] + 1) += k1 * w_2_j;
			u_vec(3 * W_idx_set[idx_vec[2]][j] + 2) += k2 * w_2_j;
			}*/
			/*for (integer j = 0; j < this->sim_mesh.n_vertices(); j++)
			{
			decimal w_0_j = W_mat(idx_vec[0], j);
			decimal w_1_j = W_mat(idx_vec[1], j);
			decimal w_2_j = W_mat(idx_vec[2], j);

			if (w_0_j > limit)
			{
			u_vec(3 * j + 0) += i0 * w_0_j;
			u_vec(3 * j + 1) += i1 * w_0_j;
			u_vec(3 * j + 2) += i2 * w_0_j;
			}
			if (w_1_j > limit)
			{
			u_vec(3 * j + 0) += j0 * w_1_j;
			u_vec(3 * j + 1) += j1 * w_1_j;
			u_vec(3 * j + 2) += j2 * w_1_j;
			}
			if (w_2_j > limit)
			{
			u_vec(3 * j + 0) += k0 * w_2_j;
			u_vec(3 * j + 1) += k1 * w_2_j;
			u_vec(3 * j + 2) += k2 * w_2_j;
			}
			}*/
		}
	}

	++row_cnt;
	fv = swv*(total_tmp_vol - wnt_volume);
	total_energy = fv*fv;

	return total_energy;
}

decimal mns::MyDeformationHandler::GetWWTerm(integer idx)
{
	decimal energy = 0;
	if (this->use_sim_or_not == false)
	{
		const integer vn = this->mesh.n_vertices();
		const integer vb = this->basis.cols();
		decimal uli_minus_wli = unknown(3 * vn + 2 * vb + 2 * idx + 0) - unknown(3 * vn + 2 * idx + 0);
		decimal uli_plus_wli = unknown(3 * vn + 2 * vb + 2 * idx + 0) + unknown(3 * vn + 2 * idx + 0);
		decimal udi_minus_wdi = unknown(3 * vn + 2 * vb + 2 * idx + 1) - unknown(3 * vn + 2 * idx + 1);
		decimal udi_plus_wdi = unknown(3 * vn + 2 * vb + 2 * idx + 1) + unknown(3 * vn + 2 * idx + 1);

		decimal sqrt_lsum = alpha_w*unknown(3 * vn + 2 * vb + 2 * idx + 0) - log(uli_minus_wli) - log(uli_plus_wli);
		sqrt_lsum = sqrt(sqrt_lsum);

		decimal sqrt_dsum = alpha_w*unknown(3 * vn + 2 * vb + 2 * idx + 1) - log(udi_minus_wdi) - log(udi_plus_wdi);
		sqrt_dsum = sqrt(sqrt_dsum);

		fw(ww_ele_cnt) = sqrt_lsum;
		ww_ele_cnt++;
		fw(ww_ele_cnt) = sqrt_dsum;
		ww_ele_cnt++;
		energy = sqrt_lsum*sqrt_lsum + sqrt_dsum*sqrt_dsum;

		decimal diff_f_uli = -(1.0 / uli_minus_wli + 1.0 / uli_plus_wli - alpha_w) / (2 * sqrt_lsum);
		decimal diff_f_wli = (1.0 / uli_minus_wli - 1.0 / uli_plus_wli) / (2 * sqrt_lsum);

		decimal diff_f_udi = -(1.0 / udi_minus_wdi + 1.0 / udi_plus_wdi - alpha_w) / (2 * sqrt_dsum);
		decimal diff_f_wdi = (1.0 / udi_minus_wdi - 1.0 / udi_plus_wdi) / (2 * sqrt_dsum);

		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * vn + 2 * vb + 2 * idx + 0, diff_f_uli));
		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * vn + 2 * idx + 0, diff_f_wli));
		++(this->row_cnt);

		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * vn + 2 * vb + 2 * idx + 1, diff_f_udi));
		triple.push_back(Eigen::Triplet<decimal>(this->row_cnt, 3 * vn + 2 * idx + 1, diff_f_wdi));
		++(this->row_cnt);
	}
	else
	{

	}
	return energy;
}

decimal mns::MyDeformationHandler::Ew()
{
	if (fw.rows() != 2 * this->basis.cols())
	{
		fw.resize(2 * this->basis.cols());
	}
	fw = Eigen::VectorXd::Zero(2 * this->basis.cols());

	decimal total_energy = 0;

	for (integer i = 0; i < this->basis.cols(); i++)
	{
		total_energy += GetWWTerm(i);
	}

	return total_energy;
}

void mns::MyDeformationHandler::ColorWeightedEdge(Eigen::MatrixXd& vc_mat, decimal fac)
{
	std::ifstream w_ifile("w_unknown.txt", std::ios::in);
	dec_vector w_vec;
	while (w_ifile.peek() != EOF)
	{
		decimal tmp;
		w_ifile >> tmp;
		w_vec.push_back(tmp);
	}
	w_vec.erase(w_vec.begin() + w_vec.size() - 1);
	w_ifile.close();
	integer l_max_col = -1;
	integer d_max_col = -1;
	decimal l_max = 0;
	decimal d_max = 0;
	integer cols = (integer)(w_vec.size() / 2);

	//#pragma omp parallel for
	for (integer i = 0; i < cols; i++)
	{
		//#pragma omp critical
		{
			if (std::abs(w_vec[2 * i + 0]) > l_max)
			{
				l_max = std::abs(w_vec[2 * i + 0]);
				l_max_col = i;
			}
			if (std::abs(w_vec[2 * i + 1]) > d_max)
			{
				d_max = std::abs(w_vec[2 * i + 1]);
				d_max_col = i;
			}
		}
	}

	int_vector l_idx_vec;
	int_vector d_idx_vec;
	dec_vector l_com_vec;
	dec_vector d_com_vec;
	decimal l_max_col_sum = 0;
	decimal d_max_col_sum = 0;
#pragma omp parallel for
	for (integer i = 0; i < this->mesh.n_edges(); i++)
	{
#pragma omp critical
		{
			if (length_basis_vec[l_max_col][i] != 0)
			{
				l_idx_vec.push_back(i);
				l_com_vec.push_back(std::abs(length_basis_vec[l_max_col][i]));
				l_max_col_sum += std::abs(length_basis_vec[l_max_col][i]);
			}
			if (dihedral_basis_vec[d_max_col][i] != 0)
			{
				d_idx_vec.push_back(i);
				d_com_vec.push_back(std::abs(dihedral_basis_vec[d_max_col][i]));
				d_max_col_sum += std::abs(dihedral_basis_vec[d_max_col][i]);
			}
		}
	}
	decimal scale_fac = fac;

	integer j = 0;
	integer l_cnt = 0;
	integer d_cnt = 0;
	integer l_flag = 1;
	integer d_flag = 1;
	for (auto he_it = this->mesh.halfedges_begin(); he_it != this->mesh.halfedges_end(); ++he_it, ++he_it, ++j)
	{
		if (/*l_cnt >= l_idx_vec.size() && */d_cnt >= d_idx_vec.size())
		{
			break;
		}
		auto e = *he_it;
		auto fromVertex = this->mesh.from_vertex_handle(e);
		auto toVertex = this->mesh.to_vertex_handle(e);
		/*if (l_flag && j == l_idx_vec[l_cnt])
		{
		if (l_cnt < l_idx_vec.size())
		{
		vc_mat.row(fromVertex.idx()) << vc_mat.row(fromVertex.idx()) + Eigen::RowVector3d(0, scale_fac * l_com_vec[l_cnt] / l_max_col_sum, 0);
		vc_mat.row(toVertex.idx()) <<  vc_mat.row(toVertex.idx()) + Eigen::RowVector3d(0, scale_fac * l_com_vec[l_cnt] / l_max_col_sum, 0);
		l_cnt++;
		if (l_cnt == l_idx_vec.size())
		{
		l_flag = 0;
		l_cnt--;
		}
		}
		}*/
		if (d_flag && j == d_idx_vec[d_cnt])
		{
			if (d_cnt < d_idx_vec.size())
			{
				vc_mat.row(fromVertex.idx()) << vc_mat.row(fromVertex.idx()) + Eigen::RowVector3d(0, scale_fac * d_com_vec[d_cnt] / d_max_col_sum, 0);
				vc_mat.row(toVertex.idx()) << vc_mat.row(toVertex.idx()) + Eigen::RowVector3d(0, scale_fac * d_com_vec[d_cnt] / d_max_col_sum, 0);
				d_cnt++;
				if (d_cnt == d_idx_vec.size())
				{
					d_flag = 0;
					d_cnt--;
					break;
				}
			}
		}
	}

	//j = 0;
	//for (auto he_it = this->mesh.halfedges_begin(); he_it != this->mesh.halfedges_end(); ++he_it, ++he_it, ++j)
	//{
	//	if (l_cnt >= l_idx_vec.size()/* && d_cnt >= d_idx_vec.size()*/)
	//	{
	//		break;
	//	}
	//	auto e = *he_it;
	//	auto fromVertex = this->mesh.from_vertex_handle(e);
	//	auto toVertex = this->mesh.to_vertex_handle(e);
	//	if (l_flag && j == l_idx_vec[l_cnt])
	//	{
	//		if (l_cnt < l_idx_vec.size())
	//		{
	//			vc_mat.row(fromVertex.idx()) << vc_mat.row(fromVertex.idx()) + Eigen::RowVector3d(0, 0, scale_fac * l_com_vec[l_cnt] / l_max_col_sum);
	//			vc_mat.row(toVertex.idx()) << vc_mat.row(toVertex.idx()) + Eigen::RowVector3d(0, 0, scale_fac * l_com_vec[l_cnt] / l_max_col_sum);
	//			l_cnt++;
	//			if (l_cnt == l_idx_vec.size())
	//			{
	//				l_flag = 0;
	//				l_cnt--;
	//				break;
	//			}
	//		}
	//	}
	//	/*if (d_flag && j == d_idx_vec[d_cnt])
	//	{
	//		if (d_cnt < d_idx_vec.size())
	//		{
	//			vc_mat.row(fromVertex.idx()) << vc_mat.row(fromVertex.idx()) + Eigen::RowVector3d(0, 0, scale_fac * d_com_vec[l_cnt] / d_max_col_sum);
	//			vc_mat.row(toVertex.idx()) << vc_mat.row(toVertex.idx()) + Eigen::RowVector3d(0, 0, scale_fac * d_com_vec[l_cnt] / d_max_col_sum);
	//			d_cnt++;
	//			if (d_cnt == d_idx_vec.size())
	//			{
	//				d_flag = 0;
	//				d_cnt--;
	//			}
	//		}
	//	}*/
	//}
}

void mns::MyDeformationHandler::ColorWeightedEdge1(Eigen::MatrixXd& vc_mat, decimal fac)
{
	std::cout << "set move pt num? (now " << this->move_pt_num << ")\n";
	int decision;
	std::cin >> decision;
	this->move_pt_num = (decision > 0) ? decision : this->move_pt_num;
	std::ifstream w_ifile("w_unknown.txt", std::ios::in);
	dec_vector w_vec;
	while (w_ifile.peek() != EOF)
	{
		decimal tmp;
		w_ifile >> tmp;
		w_vec.push_back(tmp);
	}
	w_vec.erase(w_vec.begin() + w_vec.size() - 1);
	w_ifile.close();
	std::vector<integer> top_w_max_col_vec;//integer d_max_col = -1;
	std::vector<decimal> top_w_max_vec;//decimal d_max = 0;
	std::vector<integer> top_w_max_ctr_vec;
	integer cols = w_vec.size() / 2;
	integer close_err = (integer)(this->mesh.n_edges() / cols);

#pragma omp parallel for
	for (integer i = 0; i < this->move_pt_num; i++)
	{
#pragma omp critical
		{
			top_w_max_col_vec.push_back(-1);
			top_w_max_vec.push_back(0);
			top_w_max_ctr_vec.push_back(-close_err);
		}
	}

	//#pragma omp parallel for
	for (integer i = 0; i < cols; i++)
	{
		//#pragma omp critical
		{
			for (int j = 0; j < move_pt_num; j++)
			{
				if (std::abs(w_vec[2 * i + 1]) > top_w_max_vec[j])
				{
					//check if the found component is close to the existing components
					decimal max_edge_com = 0;
					integer max_edge_idx = -1;
					for (int c = 0; c < this->mesh.n_edges(); c++)
					{
						if (dihedral_basis_vec[i][c] > max_edge_com)
						{
							max_edge_com = dihedral_basis_vec[i][c];
							max_edge_idx = c;
						}
					}
					bool check_close = false;
					integer collision_d = -1;
					for (integer d = 0; d < top_w_max_ctr_vec.size(); d++)
					{
						if (std::abs(max_edge_idx - top_w_max_ctr_vec[d]) < close_err)
						{
							check_close = true;
							collision_d = d;
							break;
						}
					}
					if (!check_close)
					{
						for (int k = j + 1; k < move_pt_num; k++)
						{
							top_w_max_vec[k] = top_w_max_vec[k - 1];
							top_w_max_col_vec[k] = top_w_max_col_vec[k - 1];
							top_w_max_ctr_vec[k] = top_w_max_ctr_vec[k - 1];
						}
						top_w_max_vec[j] = std::abs(w_vec[2 * i + 1]);
						top_w_max_col_vec[j] = i;
						top_w_max_ctr_vec[j] = max_edge_idx;
						break;
					}
					else
					{
						if (std::abs(w_vec[2 * i + 1]) < top_w_max_vec[collision_d])
						{
							break;
						}
						else
						{
							top_w_max_vec[collision_d] = std::abs(w_vec[2 * i + 1]);
							top_w_max_col_vec[collision_d] = i;
							top_w_max_ctr_vec[collision_d] = max_edge_idx;
							if ((collision_d > 0) && top_w_max_vec[collision_d - 1] < top_w_max_vec[collision_d])
							{
								for (int e = collision_d; e > 0; e--)
								{
									if (top_w_max_vec[e - 1] > top_w_max_vec[e])
									{
										break;
									}
									decimal w_tmp = top_w_max_vec[e];
									integer w_col_tmp = top_w_max_col_vec[e];
									decimal w_ctr_tmp = top_w_max_ctr_vec[e];
									top_w_max_vec[e] = top_w_max_vec[e - 1];
									top_w_max_col_vec[e] = top_w_max_col_vec[e - 1];
									top_w_max_ctr_vec[e] = top_w_max_ctr_vec[e - 1];
									top_w_max_vec[e - 1] = w_tmp;
									top_w_max_col_vec[e - 1] = w_col_tmp;
									top_w_max_ctr_vec[e - 1] = w_ctr_tmp;
								}
							}
							else if ((collision_d < (this->move_pt_num - 1)) && top_w_max_vec[collision_d + 1] > top_w_max_vec[collision_d])
							{
								for (int e = collision_d; e < (this->move_pt_num - 1); e++)
								{
									if (top_w_max_vec[e + 1] < top_w_max_vec[e])
									{
										break;
									}
									decimal w_tmp = top_w_max_vec[e];
									integer w_col_tmp = top_w_max_col_vec[e];
									decimal w_ctr_tmp = top_w_max_ctr_vec[e];
									top_w_max_vec[e] = top_w_max_vec[e + 1];
									top_w_max_col_vec[e] = top_w_max_col_vec[e + 1];
									top_w_max_ctr_vec[e] = top_w_max_ctr_vec[e + 1];
									top_w_max_vec[e + 1] = w_tmp;
									top_w_max_col_vec[e + 1] = w_col_tmp;
									top_w_max_ctr_vec[e + 1] = w_ctr_tmp;
								}
							}
						}
					}
				}
			}
		}
	}

	for (int c = 0; c < move_pt_num; c++)
	{
		int_vector d_idx_vec;
		dec_vector d_com_vec;
		decimal d_max_col_sum = 0;
#pragma omp parallel for
		for (integer i = 0; i < this->mesh.n_edges(); i++)
		{
#pragma omp critical
			{
				if (dihedral_basis_vec[top_w_max_col_vec[c]][i] != 0)
				{
					d_idx_vec.push_back(i);
					d_com_vec.push_back(std::abs(dihedral_basis_vec[top_w_max_col_vec[c]][i]));
					d_max_col_sum += std::abs(dihedral_basis_vec[top_w_max_col_vec[c]][i]);
				}
			}
		}
		decimal scale_fac = fac;

		integer j = 0;
		integer d_cnt = 0;
		integer d_flag = 1;
		for (auto he_it = this->mesh.halfedges_begin(); he_it != this->mesh.halfedges_end(); ++he_it, ++he_it, ++j)
		{
			if (/*l_cnt >= l_idx_vec.size() && */d_cnt >= d_idx_vec.size())
			{
				break;
			}
			auto e = *he_it;
			auto fromVertex = this->mesh.from_vertex_handle(e);
			auto toVertex = this->mesh.to_vertex_handle(e);
			/*if (l_flag && j == l_idx_vec[l_cnt])
			{
			if (l_cnt < l_idx_vec.size())
			{
			vc_mat.row(fromVertex.idx()) << vc_mat.row(fromVertex.idx()) + Eigen::RowVector3d(0, scale_fac * l_com_vec[l_cnt] / l_max_col_sum, 0);
			vc_mat.row(toVertex.idx()) <<  vc_mat.row(toVertex.idx()) + Eigen::RowVector3d(0, scale_fac * l_com_vec[l_cnt] / l_max_col_sum, 0);
			l_cnt++;
			if (l_cnt == l_idx_vec.size())
			{
			l_flag = 0;
			l_cnt--;
			}
			}
			}*/
			if (d_flag && j == d_idx_vec[d_cnt])
			{
				if (d_cnt < d_idx_vec.size())
				{
					vc_mat.row(fromVertex.idx()) << vc_mat.row(fromVertex.idx()) + Eigen::RowVector3d(0, scale_fac * d_com_vec[d_cnt] / d_max_col_sum, 0);
					vc_mat.row(toVertex.idx()) << vc_mat.row(toVertex.idx()) + Eigen::RowVector3d(0, scale_fac * d_com_vec[d_cnt] / d_max_col_sum, 0);
					d_cnt++;
					if (d_cnt == d_idx_vec.size())
					{
						d_flag = 0;
						d_cnt--;
						break;
					}
				}
			}
		}

		//j = 0;
		//for (auto he_it = this->mesh.halfedges_begin(); he_it != this->mesh.halfedges_end(); ++he_it, ++he_it, ++j)
		//{
		//	if (l_cnt >= l_idx_vec.size()/* && d_cnt >= d_idx_vec.size()*/)
		//	{
		//		break;
		//	}
		//	auto e = *he_it;
		//	auto fromVertex = this->mesh.from_vertex_handle(e);
		//	auto toVertex = this->mesh.to_vertex_handle(e);
		//	if (l_flag && j == l_idx_vec[l_cnt])
		//	{
		//		if (l_cnt < l_idx_vec.size())
		//		{
		//			vc_mat.row(fromVertex.idx()) << vc_mat.row(fromVertex.idx()) + Eigen::RowVector3d(0, 0, scale_fac * l_com_vec[l_cnt] / l_max_col_sum);
		//			vc_mat.row(toVertex.idx()) << vc_mat.row(toVertex.idx()) + Eigen::RowVector3d(0, 0, scale_fac * l_com_vec[l_cnt] / l_max_col_sum);
		//			l_cnt++;
		//			if (l_cnt == l_idx_vec.size())
		//			{
		//				l_flag = 0;
		//				l_cnt--;
		//				break;
		//			}
		//		}
		//	}
		//	/*if (d_flag && j == d_idx_vec[d_cnt])
		//	{
		//		if (d_cnt < d_idx_vec.size())
		//		{
		//			vc_mat.row(fromVertex.idx()) << vc_mat.row(fromVertex.idx()) + Eigen::RowVector3d(0, 0, scale_fac * d_com_vec[l_cnt] / d_max_col_sum);
		//			vc_mat.row(toVertex.idx()) << vc_mat.row(toVertex.idx()) + Eigen::RowVector3d(0, 0, scale_fac * d_com_vec[l_cnt] / d_max_col_sum);
		//			d_cnt++;
		//			if (d_cnt == d_idx_vec.size())
		//			{
		//				d_flag = 0;
		//				d_cnt--;
		//			}
		//		}
		//	}*/
		//}
	}
}


void mns::MyDeformationHandler::ColorWeightedEdge2(Eigen::MatrixXd& vc_mat, decimal fac)
{
	std::ifstream w_ifile("d_bw_unknown.txt", std::ios::in);
	dec_vector w_vec;
	while (w_ifile.peek() != EOF)
	{
		decimal tmp;
		w_ifile >> tmp;
		w_vec.push_back(tmp);
	}
	w_vec.erase(w_vec.begin() + w_vec.size() - 1);
	w_ifile.close();
	integer d_max_col = -1;
	decimal d_max = 0;
	integer cols = w_vec.size();

	//#pragma omp parallel for
	for (integer i = 0; i < cols; i++)
	{
		//#pragma omp critical
		{
			if (std::abs(w_vec[i]) > d_max)
			{
				d_max = std::abs(w_vec[i]);
				d_max_col = i;
			}
		}
	}

	int_vector d_idx_vec;
	dec_vector d_com_vec;
	decimal d_max_col_sum = 0;
#pragma omp parallel for
	for (integer i = 0; i < this->mesh.n_edges(); i++)
	{
#pragma omp critical
		{
			if (dihedral_basis_vec[d_max_col][i] != 0)
			{
				d_idx_vec.push_back(i);
				d_com_vec.push_back(std::abs(dihedral_basis_vec[d_max_col][i]));
				d_max_col_sum += std::abs(dihedral_basis_vec[d_max_col][i]);
			}
		}
	}
	decimal scale_fac = fac;

	integer j = 0;
	integer d_cnt = 0;
	integer d_flag = 1;
	for (auto he_it = this->mesh.halfedges_begin(); he_it != this->mesh.halfedges_end(); ++he_it, ++he_it, ++j)
	{
		if (/*l_cnt >= l_idx_vec.size() && */d_cnt >= d_idx_vec.size())
		{
			break;
		}
		auto e = *he_it;
		auto fromVertex = this->mesh.from_vertex_handle(e);
		auto toVertex = this->mesh.to_vertex_handle(e);
		/*if (l_flag && j == l_idx_vec[l_cnt])
		{
		if (l_cnt < l_idx_vec.size())
		{
		vc_mat.row(fromVertex.idx()) << vc_mat.row(fromVertex.idx()) + Eigen::RowVector3d(0, scale_fac * l_com_vec[l_cnt] / l_max_col_sum, 0);
		vc_mat.row(toVertex.idx()) <<  vc_mat.row(toVertex.idx()) + Eigen::RowVector3d(0, scale_fac * l_com_vec[l_cnt] / l_max_col_sum, 0);
		l_cnt++;
		if (l_cnt == l_idx_vec.size())
		{
		l_flag = 0;
		l_cnt--;
		}
		}
		}*/
		if (d_flag && j == d_idx_vec[d_cnt])
		{
			if (d_cnt < d_idx_vec.size())
			{
				vc_mat.row(fromVertex.idx()) << vc_mat.row(fromVertex.idx()) + Eigen::RowVector3d(0, scale_fac * d_com_vec[d_cnt] / d_max_col_sum, 0);
				vc_mat.row(toVertex.idx()) << vc_mat.row(toVertex.idx()) + Eigen::RowVector3d(0, scale_fac * d_com_vec[d_cnt] / d_max_col_sum, 0);
				d_cnt++;
				if (d_cnt == d_idx_vec.size())
				{
					d_flag = 0;
					d_cnt--;
					break;
				}
			}
		}
	}

	//j = 0;
	//for (auto he_it = this->mesh.halfedges_begin(); he_it != this->mesh.halfedges_end(); ++he_it, ++he_it, ++j)
	//{
	//	if (l_cnt >= l_idx_vec.size()/* && d_cnt >= d_idx_vec.size()*/)
	//	{
	//		break;
	//	}
	//	auto e = *he_it;
	//	auto fromVertex = this->mesh.from_vertex_handle(e);
	//	auto toVertex = this->mesh.to_vertex_handle(e);
	//	if (l_flag && j == l_idx_vec[l_cnt])
	//	{
	//		if (l_cnt < l_idx_vec.size())
	//		{
	//			vc_mat.row(fromVertex.idx()) << vc_mat.row(fromVertex.idx()) + Eigen::RowVector3d(0, 0, scale_fac * l_com_vec[l_cnt] / l_max_col_sum);
	//			vc_mat.row(toVertex.idx()) << vc_mat.row(toVertex.idx()) + Eigen::RowVector3d(0, 0, scale_fac * l_com_vec[l_cnt] / l_max_col_sum);
	//			l_cnt++;
	//			if (l_cnt == l_idx_vec.size())
	//			{
	//				l_flag = 0;
	//				l_cnt--;
	//				break;
	//			}
	//		}
	//	}
	//	/*if (d_flag && j == d_idx_vec[d_cnt])
	//	{
	//		if (d_cnt < d_idx_vec.size())
	//		{
	//			vc_mat.row(fromVertex.idx()) << vc_mat.row(fromVertex.idx()) + Eigen::RowVector3d(0, 0, scale_fac * d_com_vec[l_cnt] / d_max_col_sum);
	//			vc_mat.row(toVertex.idx()) << vc_mat.row(toVertex.idx()) + Eigen::RowVector3d(0, 0, scale_fac * d_com_vec[l_cnt] / d_max_col_sum);
	//			d_cnt++;
	//			if (d_cnt == d_idx_vec.size())
	//			{
	//				d_flag = 0;
	//				d_cnt--;
	//			}
	//		}
	//	}*/
	//}
}

void mns::MyDeformationHandler::ColorWeightedEdge3(Eigen::MatrixXd& vc_mat, decimal fac)
{
	std::cout << "set move pt num? (now " << this->move_pt_num << ")\n";
	int decision;
	std::cin >> decision;
	this->move_pt_num = (decision > 0) ? decision : this->move_pt_num;
	std::ifstream w_ifile("d_bw_unknown.txt", std::ios::in);
	dec_vector w_vec;
	while (w_ifile.peek() != EOF)
	{
		decimal tmp;
		w_ifile >> tmp;
		w_vec.push_back(tmp);
	}
	w_vec.erase(w_vec.begin() + w_vec.size() - 1);
	w_ifile.close();
	std::vector<integer> top_w_max_col_vec;//integer d_max_col = -1;
	std::vector<decimal> top_w_max_vec;//decimal d_max = 0;
	std::vector<integer> top_w_max_ctr_vec;
	integer cols = w_vec.size();
	integer close_err = (integer)(this->mesh.n_edges() / cols);

#pragma omp parallel for
	for (integer i = 0; i < this->move_pt_num; i++)
	{
#pragma omp critical
		{
			top_w_max_col_vec.push_back(-1);
			top_w_max_vec.push_back(0);
			top_w_max_ctr_vec.push_back(-close_err);
		}
	}

	//#pragma omp parallel for
	for (integer i = 0; i < cols; i++)
	{
		//#pragma omp critical
		{
			for (int j = 0; j < move_pt_num; j++)
			{
				if (std::abs(w_vec[i]) > top_w_max_vec[j])
				{
					//check if the found component is close to the existing components
					decimal max_edge_com = 0;
					integer max_edge_idx = -1;
					for (int c = 0; c < this->mesh.n_edges(); c++)
					{
						if (dihedral_basis_vec[i][c] > max_edge_com)
						{
							max_edge_com = dihedral_basis_vec[i][c];
							max_edge_idx = c;
						}
					}
					bool check_close = false;
					integer collision_d = -1;
					for (integer d = 0; d < top_w_max_ctr_vec.size(); d++)
					{
						if (std::abs(max_edge_idx - top_w_max_ctr_vec[d]) < close_err)
						{
							check_close = true;
							collision_d = d;
							break;
						}
					}
					if (!check_close)
					{
						for (int k = j + 1; k < move_pt_num; k++)
						{
							top_w_max_vec[k] = top_w_max_vec[k - 1];
							top_w_max_col_vec[k] = top_w_max_col_vec[k - 1];
							top_w_max_ctr_vec[k] = top_w_max_ctr_vec[k - 1];
						}
						top_w_max_vec[j] = std::abs(w_vec[i]);
						top_w_max_col_vec[j] = i;
						top_w_max_ctr_vec[j] = max_edge_idx;
						break;
					}
					else
					{
						if (std::abs(w_vec[i]) < top_w_max_vec[collision_d])
						{
							break;
						}
						else
						{
							top_w_max_vec[collision_d] = std::abs(w_vec[i]);
							top_w_max_col_vec[collision_d] = i;
							top_w_max_ctr_vec[collision_d] = max_edge_idx;
							if ((collision_d > 0) && top_w_max_vec[collision_d - 1] < top_w_max_vec[collision_d])
							{
								for (int e = collision_d; e > 0; e--)
								{
									if (top_w_max_vec[e - 1] > top_w_max_vec[e])
									{
										break;
									}
									decimal w_tmp = top_w_max_vec[e];
									integer w_col_tmp = top_w_max_col_vec[e];
									decimal w_ctr_tmp = top_w_max_ctr_vec[e];
									top_w_max_vec[e] = top_w_max_vec[e - 1];
									top_w_max_col_vec[e] = top_w_max_col_vec[e - 1];
									top_w_max_ctr_vec[e] = top_w_max_ctr_vec[e - 1];
									top_w_max_vec[e - 1] = w_tmp;
									top_w_max_col_vec[e - 1] = w_col_tmp;
									top_w_max_ctr_vec[e - 1] = w_ctr_tmp;
								}
							}
							else if ((collision_d < (this->move_pt_num - 1)) && top_w_max_vec[collision_d + 1] > top_w_max_vec[collision_d])
							{
								for (int e = collision_d; e < (this->move_pt_num - 1); e++)
								{
									if (top_w_max_vec[e + 1] < top_w_max_vec[e])
									{
										break;
									}
									decimal w_tmp = top_w_max_vec[e];
									integer w_col_tmp = top_w_max_col_vec[e];
									decimal w_ctr_tmp = top_w_max_ctr_vec[e];
									top_w_max_vec[e] = top_w_max_vec[e + 1];
									top_w_max_col_vec[e] = top_w_max_col_vec[e + 1];
									top_w_max_ctr_vec[e] = top_w_max_ctr_vec[e + 1];
									top_w_max_vec[e + 1] = w_tmp;
									top_w_max_col_vec[e + 1] = w_col_tmp;
									top_w_max_ctr_vec[e + 1] = w_ctr_tmp;
								}
							}
						}
					}
				}
			}
		}
	}

	for (int c = 0; c < move_pt_num; c++)
	{
		int_vector d_idx_vec;
		dec_vector d_com_vec;
		decimal d_max_col_sum = 0;
#pragma omp parallel for
		for (integer i = 0; i < this->mesh.n_edges(); i++)
		{
#pragma omp critical
			{
				if (dihedral_basis_vec[top_w_max_col_vec[c]][i] != 0)
				{
					d_idx_vec.push_back(i);
					d_com_vec.push_back(std::abs(dihedral_basis_vec[top_w_max_col_vec[c]][i]));
					d_max_col_sum += std::abs(dihedral_basis_vec[top_w_max_col_vec[c]][i]);
				}
			}
		}
		decimal scale_fac = fac;

		integer j = 0;
		integer d_cnt = 0;
		integer d_flag = 1;
		for (auto he_it = this->mesh.halfedges_begin(); he_it != this->mesh.halfedges_end(); ++he_it, ++he_it, ++j)
		{
			if (/*l_cnt >= l_idx_vec.size() && */d_cnt >= d_idx_vec.size())
			{
				break;
			}
			auto e = *he_it;
			auto fromVertex = this->mesh.from_vertex_handle(e);
			auto toVertex = this->mesh.to_vertex_handle(e);
			/*if (l_flag && j == l_idx_vec[l_cnt])
			{
			if (l_cnt < l_idx_vec.size())
			{
			vc_mat.row(fromVertex.idx()) << vc_mat.row(fromVertex.idx()) + Eigen::RowVector3d(0, scale_fac * l_com_vec[l_cnt] / l_max_col_sum, 0);
			vc_mat.row(toVertex.idx()) <<  vc_mat.row(toVertex.idx()) + Eigen::RowVector3d(0, scale_fac * l_com_vec[l_cnt] / l_max_col_sum, 0);
			l_cnt++;
			if (l_cnt == l_idx_vec.size())
			{
			l_flag = 0;
			l_cnt--;
			}
			}
			}*/
			if (d_flag && j == d_idx_vec[d_cnt])
			{
				if (d_cnt < d_idx_vec.size())
				{
					vc_mat.row(fromVertex.idx()) << vc_mat.row(fromVertex.idx()) + Eigen::RowVector3d(0, scale_fac * d_com_vec[d_cnt] / d_max_col_sum, 0);
					vc_mat.row(toVertex.idx()) << vc_mat.row(toVertex.idx()) + Eigen::RowVector3d(0, scale_fac * d_com_vec[d_cnt] / d_max_col_sum, 0);
					d_cnt++;
					if (d_cnt == d_idx_vec.size())
					{
						d_flag = 0;
						d_cnt--;
						break;
					}
				}
			}
		}

		//j = 0;
		//for (auto he_it = this->mesh.halfedges_begin(); he_it != this->mesh.halfedges_end(); ++he_it, ++he_it, ++j)
		//{
		//	if (l_cnt >= l_idx_vec.size()/* && d_cnt >= d_idx_vec.size()*/)
		//	{
		//		break;
		//	}
		//	auto e = *he_it;
		//	auto fromVertex = this->mesh.from_vertex_handle(e);
		//	auto toVertex = this->mesh.to_vertex_handle(e);
		//	if (l_flag && j == l_idx_vec[l_cnt])
		//	{
		//		if (l_cnt < l_idx_vec.size())
		//		{
		//			vc_mat.row(fromVertex.idx()) << vc_mat.row(fromVertex.idx()) + Eigen::RowVector3d(0, 0, scale_fac * l_com_vec[l_cnt] / l_max_col_sum);
		//			vc_mat.row(toVertex.idx()) << vc_mat.row(toVertex.idx()) + Eigen::RowVector3d(0, 0, scale_fac * l_com_vec[l_cnt] / l_max_col_sum);
		//			l_cnt++;
		//			if (l_cnt == l_idx_vec.size())
		//			{
		//				l_flag = 0;
		//				l_cnt--;
		//				break;
		//			}
		//		}
		//	}
		//	/*if (d_flag && j == d_idx_vec[d_cnt])
		//	{
		//		if (d_cnt < d_idx_vec.size())
		//		{
		//			vc_mat.row(fromVertex.idx()) << vc_mat.row(fromVertex.idx()) + Eigen::RowVector3d(0, 0, scale_fac * d_com_vec[l_cnt] / d_max_col_sum);
		//			vc_mat.row(toVertex.idx()) << vc_mat.row(toVertex.idx()) + Eigen::RowVector3d(0, 0, scale_fac * d_com_vec[l_cnt] / d_max_col_sum);
		//			d_cnt++;
		//			if (d_cnt == d_idx_vec.size())
		//			{
		//				d_flag = 0;
		//				d_cnt--;
		//			}
		//		}
		//	}*/
		//}
	}
}

void mns::MyDeformationHandler::ColorWeightedEdge4(Eigen::MatrixXd& vc_mat, decimal fac)
{
	std::cout << "set move pt num? (now " << this->move_pt_num << ")\n";
	int decision;
	std::cin >> decision;
	this->move_pt_num = (decision > 0) ? decision : this->move_pt_num;
	std::ifstream w_ifile("l_bw_unknown.txt", std::ios::in);
	dec_vector w_vec;
	while (w_ifile.peek() != EOF)
	{
		decimal tmp;
		w_ifile >> tmp;
		w_vec.push_back(tmp);
	}
	w_vec.erase(w_vec.begin() + w_vec.size() - 1);
	w_ifile.close();
	std::vector<integer> top_w_max_col_vec;//integer d_max_col = -1;
	std::vector<decimal> top_w_max_vec;//decimal d_max = 0;
	std::vector<integer> top_w_max_ctr_vec;
	integer cols = w_vec.size();
	integer close_err = (integer)(this->mesh.n_edges() / cols);

#pragma omp parallel for
	for (integer i = 0; i < this->move_pt_num; i++)
	{
#pragma omp critical
		{
			top_w_max_col_vec.push_back(-1);
			top_w_max_vec.push_back(0);
			top_w_max_ctr_vec.push_back(-close_err);
		}
	}

	//#pragma omp parallel for
	for (integer i = 0; i < cols; i++)
	{
		//#pragma omp critical
		{
			for (int j = 0; j < move_pt_num; j++)
			{
				if (std::abs(w_vec[i]) > top_w_max_vec[j])
				{
					//check if the found component is close to the existing components
					decimal max_edge_com = 0;
					integer max_edge_idx = -1;
					for (int c = 0; c < this->mesh.n_edges(); c++)
					{
						if (dihedral_basis_vec[i][c] > max_edge_com)
						{
							max_edge_com = dihedral_basis_vec[i][c];
							max_edge_idx = c;
						}
					}
					bool check_close = false;
					integer collision_d = -1;
					for (integer d = 0; d < top_w_max_ctr_vec.size(); d++)
					{
						if (std::abs(max_edge_idx - top_w_max_ctr_vec[d]) < close_err)
						{
							check_close = true;
							collision_d = d;
							break;
						}
					}
					if (!check_close)
					{
						for (int k = j + 1; k < move_pt_num; k++)
						{
							top_w_max_vec[k] = top_w_max_vec[k - 1];
							top_w_max_col_vec[k] = top_w_max_col_vec[k - 1];
							top_w_max_ctr_vec[k] = top_w_max_ctr_vec[k - 1];
						}
						top_w_max_vec[j] = std::abs(w_vec[i]);
						top_w_max_col_vec[j] = i;
						top_w_max_ctr_vec[j] = max_edge_idx;
						break;
					}
					else
					{
						if (std::abs(w_vec[i]) < top_w_max_vec[collision_d])
						{
							break;
						}
						else
						{
							top_w_max_vec[collision_d] = std::abs(w_vec[i]);
							top_w_max_col_vec[collision_d] = i;
							top_w_max_ctr_vec[collision_d] = max_edge_idx;
							if ((collision_d > 0) && top_w_max_vec[collision_d - 1] < top_w_max_vec[collision_d])
							{
								for (int e = collision_d; e > 0; e--)
								{
									if (top_w_max_vec[e - 1] > top_w_max_vec[e])
									{
										break;
									}
									decimal w_tmp = top_w_max_vec[e];
									integer w_col_tmp = top_w_max_col_vec[e];
									decimal w_ctr_tmp = top_w_max_ctr_vec[e];
									top_w_max_vec[e] = top_w_max_vec[e - 1];
									top_w_max_col_vec[e] = top_w_max_col_vec[e - 1];
									top_w_max_ctr_vec[e] = top_w_max_ctr_vec[e - 1];
									top_w_max_vec[e - 1] = w_tmp;
									top_w_max_col_vec[e - 1] = w_col_tmp;
									top_w_max_ctr_vec[e - 1] = w_ctr_tmp;
								}
							}
							else if ((collision_d < (this->move_pt_num - 1)) && top_w_max_vec[collision_d + 1] > top_w_max_vec[collision_d])
							{
								for (int e = collision_d; e < (this->move_pt_num - 1); e++)
								{
									if (top_w_max_vec[e + 1] < top_w_max_vec[e])
									{
										break;
									}
									decimal w_tmp = top_w_max_vec[e];
									integer w_col_tmp = top_w_max_col_vec[e];
									decimal w_ctr_tmp = top_w_max_ctr_vec[e];
									top_w_max_vec[e] = top_w_max_vec[e + 1];
									top_w_max_col_vec[e] = top_w_max_col_vec[e + 1];
									top_w_max_ctr_vec[e] = top_w_max_ctr_vec[e + 1];
									top_w_max_vec[e + 1] = w_tmp;
									top_w_max_col_vec[e + 1] = w_col_tmp;
									top_w_max_ctr_vec[e + 1] = w_ctr_tmp;
								}
							}
						}
					}
				}
			}
		}
	}

	for (int c = 0; c < move_pt_num; c++)
	{
		int_vector l_idx_vec;
		dec_vector l_com_vec;
		decimal l_max_col_sum = 0;
#pragma omp parallel for
		for (integer i = 0; i < this->mesh.n_edges(); i++)
		{
#pragma omp critical
			{
				if (length_basis_vec[top_w_max_col_vec[c]][i] != 0)
				{
					l_idx_vec.push_back(i);
					l_com_vec.push_back(std::abs(length_basis_vec[top_w_max_col_vec[c]][i]));
					l_max_col_sum += std::abs(length_basis_vec[top_w_max_col_vec[c]][i]);
				}
			}
		}
		decimal scale_fac = fac;

		integer j = 0;
		integer l_cnt = 0;
		integer l_flag = 1;
		for (auto he_it = this->mesh.halfedges_begin(); he_it != this->mesh.halfedges_end(); ++he_it, ++he_it, ++j)
		{
			if (/*l_cnt >= l_idx_vec.size() && */l_cnt >= l_idx_vec.size())
			{
				break;
			}
			auto e = *he_it;
			auto fromVertex = this->mesh.from_vertex_handle(e);
			auto toVertex = this->mesh.to_vertex_handle(e);
			if (l_flag && j == l_idx_vec[l_cnt])
			{
				if (l_cnt < l_idx_vec.size())
				{
					vc_mat.row(fromVertex.idx()) << vc_mat.row(fromVertex.idx()) + Eigen::RowVector3d(0, scale_fac * l_com_vec[l_cnt] / l_max_col_sum, 0);
					vc_mat.row(toVertex.idx()) << vc_mat.row(toVertex.idx()) + Eigen::RowVector3d(0, scale_fac * l_com_vec[l_cnt] / l_max_col_sum, 0);
					l_cnt++;
					if (l_cnt == l_idx_vec.size())
					{
						l_flag = 0;
						l_cnt--;
						break;
					}
				}
			}
		}
	}
}

bool mns::MyDeformationHandler::ColorBwDecision()
{
	return (this->use_eg_or_not == true && this->use_basis_or_not == true);
}

void mns::MyDeformationHandler::ColorSpecificComponent(Eigen::MatrixXd& vc_mat, integer choice, decimal fac)
{
	int_vector l_idx_vec;
	int_vector d_idx_vec;
	dec_vector l_com_vec;
	dec_vector d_com_vec;
	decimal l_max_col_sum = 0;
	decimal d_max_col_sum = 0;
#pragma omp parallel for
	for (integer i = 0; i < this->mesh.n_edges(); i++)
	{
#pragma omp critical
		{
			if (length_basis_vec[choice][i] != 0)
			{
				l_idx_vec.push_back(i);
				l_com_vec.push_back(std::abs(length_basis_vec[choice][i]));
				l_max_col_sum += std::abs(length_basis_vec[choice][i]);
			}
			if (dihedral_basis_vec[choice][i] != 0)
			{
				d_idx_vec.push_back(i);
				d_com_vec.push_back(std::abs(dihedral_basis_vec[choice][i]));
				d_max_col_sum += std::abs(dihedral_basis_vec[choice][i]);
			}
		}
	}
	decimal scale_fac = fac;

	integer j = 0;
	integer l_cnt = 0;
	integer d_cnt = 0;
	integer l_flag = 1;
	integer d_flag = 1;
	for (auto he_it = this->mesh.halfedges_begin(); he_it != this->mesh.halfedges_end(); ++he_it, ++he_it, ++j)
	{
		if (/*l_cnt >= l_idx_vec.size() && */d_cnt >= d_idx_vec.size())
		{
			break;
		}
		auto e = *he_it;
		auto fromVertex = this->mesh.from_vertex_handle(e);
		auto toVertex = this->mesh.to_vertex_handle(e);
		/*if (l_flag && j == l_idx_vec[l_cnt])
		{
		if (l_cnt < l_idx_vec.size())
		{
		vc_mat.row(fromVertex.idx()) << vc_mat.row(fromVertex.idx()) + Eigen::RowVector3d(0, scale_fac * l_com_vec[l_cnt] / l_max_col_sum, 0);
		vc_mat.row(toVertex.idx()) <<  vc_mat.row(toVertex.idx()) + Eigen::RowVector3d(0, scale_fac * l_com_vec[l_cnt] / l_max_col_sum, 0);
		l_cnt++;
		if (l_cnt == l_idx_vec.size())
		{
		l_flag = 0;
		l_cnt--;
		}
		}
		}*/
		if (d_flag && j == d_idx_vec[d_cnt])
		{
			if (d_cnt < d_idx_vec.size())
			{
				vc_mat.row(fromVertex.idx()) << vc_mat.row(fromVertex.idx()) + Eigen::RowVector3d(0, scale_fac * d_com_vec[d_cnt] / d_max_col_sum, 0);
				vc_mat.row(toVertex.idx()) << vc_mat.row(toVertex.idx()) + Eigen::RowVector3d(0, scale_fac * d_com_vec[d_cnt] / d_max_col_sum, 0);
				d_cnt++;
				if (d_cnt == d_idx_vec.size())
				{
					d_flag = 0;
					d_cnt--;
					break;
				}
			}
		}
	}
}

void mns::MyDeformationHandler::ColorSpecificComponent2(Eigen::MatrixXd& vc_mat, integer choice, decimal fac)
{
	int_vector l_idx_vec;
	dec_vector l_com_vec;
	decimal l_max_col_sum = 0;
#pragma omp parallel for
	for (integer i = 0; i < this->mesh.n_edges(); i++)
	{
#pragma omp critical
		{
			if (length_basis_vec[choice][i] != 0)
			{
				l_idx_vec.push_back(i);
				l_com_vec.push_back(std::abs(length_basis_vec[choice][i]));
				l_max_col_sum += std::abs(length_basis_vec[choice][i]);
			}
		}
	}
	decimal scale_fac = fac;

	integer j = 0;
	integer l_cnt = 0;
	integer l_flag = 1;
	for (auto he_it = this->mesh.halfedges_begin(); he_it != this->mesh.halfedges_end(); ++he_it, ++he_it, ++j)
	{
		if (l_cnt >= l_idx_vec.size())
		{
			break;
		}
		auto e = *he_it;
		auto fromVertex = this->mesh.from_vertex_handle(e);
		auto toVertex = this->mesh.to_vertex_handle(e);

		if (l_flag && j == l_idx_vec[l_cnt])
		{
			if (l_cnt < l_idx_vec.size())
			{
				vc_mat.row(fromVertex.idx()) << vc_mat.row(fromVertex.idx()) + Eigen::RowVector3d(0, scale_fac * l_com_vec[l_cnt] / l_max_col_sum, 0);
				vc_mat.row(toVertex.idx()) << vc_mat.row(toVertex.idx()) + Eigen::RowVector3d(0, scale_fac * l_com_vec[l_cnt] / l_max_col_sum, 0);
				l_cnt++;
				if (l_cnt == l_idx_vec.size())
				{
					l_flag = 0;
					l_cnt--;
					break;
				}
			}
		}
	}
}

void mns::MyDeformationHandler::ColorComponent(Eigen::MatrixXd& vc_mat, decimal fac)
{
	std::ifstream r_ifile("total_R.txt", std::ios::in);
	dec_vector r_vec;
	while (r_ifile.peek() != EOF)
	{
		decimal tmp;
		r_ifile >> tmp;
		r_vec.push_back(tmp);
	}
	r_vec.erase(r_vec.begin() + r_vec.size() - 1);
	r_ifile.close();

	Eigen::VectorXd r_mat = Eigen::VectorXd::Zero(r_vec.size());
#pragma omp parallel for
	for (integer i = 0; i < r_vec.size(); i++)
	{
#pragma omp critical
		{
			r_mat[i] = sqrt(abs(r_vec[i]));
		}
	}
	r_mat = r_mat / r_mat.maxCoeff();

	decimal scale_fac = fac;
	integer j = 0;
	for (auto he_it = this->mesh.halfedges_begin(); he_it != this->mesh.halfedges_end(); ++he_it, ++he_it, ++j)
	{
		auto e = *he_it;
		auto fromVertex = this->mesh.from_vertex_handle(e);
		auto toVertex = this->mesh.to_vertex_handle(e);

		vc_mat.row(fromVertex.idx()) << vc_mat.row(fromVertex.idx()) + Eigen::RowVector3d(0, r_mat[j] * scale_fac, 0);
		vc_mat.row(toVertex.idx()) << vc_mat.row(toVertex.idx()) + Eigen::RowVector3d(0, r_mat[j] * scale_fac, 0);
	}
}


void mns::MyDeformationHandler::BuildJmat(SMatXd& Jmat, SMatXd& Jbar, bool vol)
{
	if (f.rows() != row_cnt)
	{
		f.resize(row_cnt);
	}
	f = Eigen::VectorXd::Zero(row_cnt);

	int cnt = 0;
	if (alpha_c != 0)
	{
#pragma omp parallel for
		for (int i = cnt; i < cnt + fc.size(); i++)
		{
#pragma omp critical
			{
				f(i) = fc(i - cnt);
			}
		}
		cnt += fc.size();
	}
	if (alpha_s != 0 || alpha_b != 0)
	{
#pragma omp parallel for
		for (int i = cnt; i < cnt + fsb.size(); i++)
		{
#pragma omp critical
			{
				f(i) = fsb(i - cnt);
			}
		}
		cnt += fsb.size();
	}
	if (alpha_w != 0 && this->use_ww_or_not == true)
	{
#pragma omp parallel for
		for (int i = cnt; i < cnt + fw.size(); i++)
		{
#pragma omp critical
			{
				f(i) = fw(i - cnt);
			}
		}
		cnt += fw.size();
	}
	/*if (alpha_s != 0)
	{
	for (int i = cnt; i < cnt + fs.size(); i++)
	{
	f(i) = fs(i - cnt);
	}
	cnt += fs.size();
	}
	if (alpha_b != 0)
	{
	for (int i = cnt; i < cnt + fb.size(); i++)
	{
	f(i) = fb(i - cnt);
	}
	cnt += fb.size();
	}*/
	if (vol)
	{
		f(cnt) = fv;
		Jbar = SMatXd(f.size() - 1, unknown.size());
		if (this->use_sim_or_not == false)
		{
			//Jbar = SMatXd(f.size() - 1, unknown.size());
			Jbar.setFromTriplets(triple.begin(), triple.end());
#pragma omp parallel for
			for (integer i = 0; i < this->mesh.n_vertices(); i++)
			{
#pragma omp critical
				{
					//soft con
					triple.push_back(Eigen::Triplet<decimal>(row_cnt - 1, 3 * i + 0, u_vec(3 * i + 0)));
					triple.push_back(Eigen::Triplet<decimal>(row_cnt - 1, 3 * i + 1, u_vec(3 * i + 1)));
					triple.push_back(Eigen::Triplet<decimal>(row_cnt - 1, 3 * i + 2, u_vec(3 * i + 2)));

					//hard con
					/*if (idx_move[i] != -1)
					{
					triple.push_back(Eigen::Triplet<decimal>(row_cnt - 1, 3 * (i - idx_move[i]) + 0, u_vec(3 * (i - idx_move[i]) + 0)));
					triple.push_back(Eigen::Triplet<decimal>(row_cnt - 1, 3 * (i - idx_move[i]) + 1, u_vec(3 * (i - idx_move[i]) + 1)));
					triple.push_back(Eigen::Triplet<decimal>(row_cnt - 1, 3 * (i - idx_move[i]) + 2, u_vec(3 * (i - idx_move[i]) + 2)));
					}*/
				}
			}
			if (this->use_eg_or_not == true && this->use_basis_or_not == false)
			{
#pragma omp parallel for
				for (integer i = 0; i < this->eg_num; i++)
				{
#pragma omp critical
					{
						triple.push_back(Eigen::Triplet<decimal>(row_cnt - 1, 3 * this->mesh.n_vertices() + i, -swv*(eg_volume_vec[i] - ori_volume)));
						u_vec(3 * this->mesh.n_vertices() + i) = -swv*(eg_volume_vec[i] - ori_volume);
					}
				}
			}
			if (this->use_eg_or_not == false && this->use_basis_or_not == true)
			{
#pragma omp parallel for
				for (integer i = 0; i < this->basis.cols(); i++)
				{
#pragma omp critical
					{
						triple.push_back(Eigen::Triplet<decimal>(row_cnt - 1, 3 * this->mesh.n_vertices() + i, -swv*basis(2 * this->mesh.n_edges(), i)));
						u_vec(3 * this->mesh.n_vertices() + i) = -swv*basis(2 * this->mesh.n_edges(), i);
					}
				}
			}
			//hard con
			/*if (this->use_eg_or_not == true)
			{
			for (integer i = 0; i < this->eg_num; i++)
			{
			triple.push_back(Eigen::Triplet<decimal>(row_cnt - 1, 3 * (this->mesh.n_vertices() - cps.set.size()) + i, -swv*(eg_volume_vec[i] - ori_volume)));
			}
			}*/
		}

		if (this->use_sim_or_not == true)
		{
			/*SMatXd Jbar_tmp;
			if (this->use_eg_or_not == true)
			{
			Jbar_tmp = SMatXd(f.size() - 1, 3 * this->mesh.n_vertices() + this->eg_num);
			}
			if (this->use_basis_or_not == true)
			{
			Jbar_tmp = SMatXd(f.size() - 1, 3 * this->mesh.n_vertices() + this->basis.cols());
			}
			if (this->use_eg_or_not == false && this->use_basis_or_not == false)
			{
			Jbar_tmp = SMatXd(f.size() - 1, 3 * this->mesh.n_vertices());
			}*/

			SMatXd Jbar_tmp1 = SMatXd(f.size() - 1, 3 * this->mesh.n_vertices());
			Jbar_tmp1.setFromTriplets(j_triple.begin(), j_triple.end());
			SMatXd Jbar_tmp2 = Jbar_tmp1 * huge_W_mat;

			Jbar.setFromTriplets(triple.begin(), triple.end());
			Jbar.leftCols(3 * this->sim_mesh.n_vertices()) = Jbar_tmp2;

#pragma omp parallel for
			for (integer i = 0; i < this->mesh.n_vertices(); i++)
			{
#pragma omp critical
				{
					//soft con
					j_triple.push_back(Eigen::Triplet<decimal>(row_cnt - 1, 3 * i + 0, u_vec(3 * i + 0)));
					j_triple.push_back(Eigen::Triplet<decimal>(row_cnt - 1, 3 * i + 1, u_vec(3 * i + 1)));
					j_triple.push_back(Eigen::Triplet<decimal>(row_cnt - 1, 3 * i + 2, u_vec(3 * i + 2)));

					//hard con(not ready)
				}
			}
			if (this->use_eg_or_not == true && this->use_basis_or_not == false)
			{
#pragma omp parallel for
				for (integer i = 0; i < this->eg_num; i++)
				{
#pragma omp critical
					{
						triple.push_back(Eigen::Triplet<decimal>(row_cnt - 1, 3 * this->sim_mesh.n_vertices() + i, -swv*(eg_volume_vec[i] - ori_volume)));
					}
				}
			}
			if (this->use_eg_or_not == false && this->use_basis_or_not == true)
			{
#pragma omp parallel for
				for (integer i = 0; i < this->basis.cols(); i++)
				{
#pragma omp critical
					{
						triple.push_back(Eigen::Triplet<decimal>(row_cnt - 1, 3 * this->sim_mesh.n_vertices() + i, -swv*basis(2 * this->mesh.n_edges(), i)));
					}
				}
			}
		}

		++cnt;
	}

	Jmat = SMatXd(f.size(), unknown.size());
	Jmat.setFromTriplets(triple.begin(), triple.end());
	if (this->use_sim_or_not == true)
	{
		SMatXd tmpJmat = SMatXd(f.size(), 3 * this->mesh.n_vertices());
		tmpJmat.setFromTriplets(j_triple.begin(), j_triple.end());
		SMatXd left_Jmat = tmpJmat * huge_W_mat;
		Jmat.leftCols(3 * this->sim_mesh.n_vertices()) = left_Jmat;
		Jmat.makeCompressed();
	}

	triple.clear();
	j_triple.clear();
	/*
	std::ofstream j_ofile("J.txt", std::ios::out);
	std::ofstream f_ofile("fx.txt", std::ios::out);
	if (j_ofile.is_open())
	{
	j_ofile << Jmat;
	}
	j_ofile.close();
	if (f_ofile.is_open())
	{
	f_ofile << f;
	}
	f_ofile.close();
	*/
	std::cout << "pause" << std::endl;
}

void mns::MyDeformationHandler::Deform_new(std::vector<integer>& s_vid_vec, Eigen::MatrixXd& v_mat, Eigen::MatrixXd& con_v_mat)
{
	std::ofstream energy_ofile("energy.txt", std::ios::out);
	energy_ofile.close();

	std::ofstream unk_ofile("unknowndif.txt", std::ios::out);
	unk_ofile.close();

	Eigen::SimplicialLLT<SMatXd> llt;
	clock_t start = clock();
	while (1)
	{
		alpha_c = alpha_c * 1;
		alpha_s = alpha_s * 1;
		alpha_b = alpha_b * 1;
		//mod unknown by basis
		/*if (this->cal_bw_or_not == true && this->use_eg_or_not == false && this->use_basis_or_not == true)
		{
			CalWunknownByBasis();
		}
		if (this->cal_bw_or_not == true && this->use_eg_or_not == true && this->use_basis_or_not == true)
		{
			CalWunknownByBasis2();
		}*/
		//

		if (this->use_eg_or_not == false && this->use_basis_or_not == true)
		{
			UpdateWntVecWithBasis2();
		}
		

		decimal Econ = 0;
		if (alpha_c != 0)
		{
			Econ = this->Ec();
		}
		decimal Estr = 0;
		/*if (alpha_s != 0)
		{
		Estr = this->Es();
		}*/
		decimal Eben = 0;
		/*if (alpha_b != 0)
		{
		Eben = this->Eb();
		}*/
		if (alpha_s != 0 || alpha_b != 0)
		{
			//std::vector<decimal> esb = Esb();
			std::vector<decimal> esb = EsbNew();
			Estr = esb[0];
			Eben = esb[1];
		}
		decimal Ewwc = 0;
		if (alpha_w != 0 && this->use_ww_or_not == true)
		{
			Ewwc = this->Ew();
		}
		decimal Evol = 0;
		if (alpha_v != 0)
		{
			Evol = Ev();
		}
		if (isnan(Ewwc))
		{
			unknown = pre_unknown;
			break;
		}

		decimal Energy = Econ + Estr + Eben + Evol + Ewwc;
		//decimal Energy = Econ + Estr;
		std::ofstream energy_ofile("energy.txt", std::ios::app);
		energy_ofile << iter << "-th iteration: " << Econ << " " << Estr << " " << Eben << " " << Energy << std::endl;
		energy_ofile.close();
		cout << iter << "-th iteration: " << Econ << " " << Estr << " " << Eben << " " << Energy << std::endl;
		if (Econ < 100)
		{
			break;
		}
		//YP
		if (this->pre_energy < Energy)
		{
			std::ofstream energy_ofile("energy.txt", std::ios::app);
			energy_ofile << "iteration ends by 'pre_energy < energy'" << std::endl;
			energy_ofile.close();
			cout << "iteration ends by 'pre_energy < energy'" << std::endl;
			unknown = pre_unknown;

			/*if (this->cal_bw_or_not == true && this->use_eg_or_not == false && this->use_basis_or_not == true)
			{
				CalWunknownByBasis();
			}
			if (this->cal_bw_or_not == true && this->use_eg_or_not == true && this->use_basis_or_not == true)
			{
				CalWunknownByBasis2();
			}
			break;*/
		}

		if (fitmode == 0)
		{
			epsilon2 = 5;
		}
		if (fitmode == 1)
		{
			epsilon2 = 0.01;
		}
		//YP
		if ((abs(Energy - pre_energy) <= (epsilon1*pre_energy) || (abs(Energy)<epsilon2)) && (this->ch_con == 0))
		{
			std::ofstream energy_ofile("energy.txt", std::ios::app);
			energy_ofile << "iteration ends by meet relative error threshold." << std::endl;
			energy_ofile.close();
			break;
		}

		if (iter > max_iter)
		{
			break;
		}
		decimal pre_e = this->pre_energy;
		this->pre_energy = Energy;

		SMatXd Jmat;
		SMatXd Jbar;
		Eigen::VectorXd hg;
		BuildJmat(Jmat, Jbar, (alpha_v != 0));

		//Jmat.makeCompressed();
		if (alpha_v == 0)
		{
			SMatXd Jt = Jmat.transpose();
			SMatXd JtJprime = (Jt*Jmat).pruned();

			triple.clear();
			SMatXd sI = SMatXd(JtJprime.rows(), JtJprime.cols());
			for (int i = 0; i<sI.rows(); i++)
			{
				triple.push_back(Eigen::Triplet<decimal>(i, i, 0.00001));
			}
			sI.setFromTriplets(triple.begin(), triple.end());
			triple.clear();

			SMatXd JtJ = (JtJprime + sI).pruned();

			/*std::ofstream j_ofile("Jtj.txt", std::ios::out);
			if (j_ofile.is_open())
			{
			j_ofile << JtJ;
			}
			j_ofile.close();*/

			llt.compute(JtJ);
			//jtj = JtJ;
			if (llt.info() != Eigen::Success)
			{
				std::cout << llt.info();
				return;
			}
			hg = llt.solve((-1)*Jt*f);
			/*std::ofstream jf_ofile("jtf.txt", std::ios::out);
			if (jf_ofile.is_open())
			{
			jf_ofile << Jt*f;
			}
			jf_ofile.close();
			std::ofstream h_ofile("hg.txt", std::ios::out);
			if (h_ofile.is_open())
			{
			h_ofile << hg;
			}
			h_ofile.close();*/
			if (llt.info() != Eigen::Success)
			{
				return;
			}
		}
		else
		{
			SMatXd Jt = Jmat.transpose();
			SMatXd Jbt = Jbar.transpose();
			SMatXd tmp_A = (Jbt*Jbar).pruned();
			Eigen::MatrixXd b = -Jt * f;

			triple.clear();
			SMatXd sI = SMatXd(tmp_A.rows(), tmp_A.cols());
#pragma omp parallel for 
			for (int i = 0; i<sI.rows(); i++)
			{
#pragma omp critical
				{
					triple.push_back(Eigen::Triplet<decimal>(i, i, 0.00001));
				}
			}

			sI.setFromTriplets(triple.begin(), triple.end());
			triple.clear();

			SMatXd A = tmp_A + sI;

			llt.compute(A);
			//jtj = JtJ;
			if (llt.info() != Eigen::Success)
			{
				//std::cout << llt.info();
				return;
			}
			Eigen::VectorXd y = llt.solve(b);
			if (llt.info() != Eigen::Success)
			{
				return;
			}
			if (this->use_sim_or_not == false)
			{
				Eigen::VectorXd z = llt.solve(u_vec);
				if (llt.info() != Eigen::Success)
				{
					return;
				}
				hg = y - ((z*(u_vec.transpose()*y)) / (1 + u_vec.transpose()*z));
			}
			if (this->use_sim_or_not == true)
			{
				Eigen::RowVectorXd u_rowvec_tmp1 = u_vec.transpose();
				Eigen::RowVectorXd u_rowvec_left_tmp = u_rowvec_tmp1.block(0, 0, 1, 3 * this->mesh.n_vertices());
				Eigen::RowVectorXd u_rowvec_right_tmp = u_rowvec_tmp1.block(0, 3 * this->mesh.n_vertices(), 1, u_vec.size() - 3 * this->mesh.n_vertices());
				Eigen::RowVectorXd u_rowvec = Eigen::RowVectorXd::Zero(unknown.size());
				u_rowvec.block(0, 0, 1, 3 * this->sim_mesh.n_vertices()) = u_rowvec_left_tmp * huge_W_mat;
				u_rowvec.block(0, 3 * this->sim_mesh.n_vertices(), 1, unknown.size() - 3 * this->sim_mesh.n_vertices()) = u_rowvec_right_tmp;
				Eigen::VectorXd new_u_vec = u_rowvec.transpose();
				Eigen::VectorXd z = llt.solve(new_u_vec);
				if (llt.info() != Eigen::Success)
				{
					return;
				}
				hg = y - ((z*(new_u_vec.transpose()*y)) / (1 + new_u_vec.transpose()*z));
			}
		}

		//calculate optimal step length
		//decimal step_len = (iter == 0) ? ori_step : CalStepLength(pre_e, Jmat, this->f, hg, this->unknown);
		decimal step_len = ori_step;
		//

		pre_unknown = unknown;
#pragma omp parallel
		{
#pragma omp for
			for (int i = 0; i < 3 *  this->mesh.n_vertices(); i++)
			{
#pragma omp critical
				{
					unknown(i) += step_len * hg(i);
				}
			}
		}
		std::ofstream unk_ofile("unknowndif.txt", std::ios::app);
		unk_ofile << "unknown difference iteration between " << iter << " and " << iter + 1 << ": " << unknown.norm() - pre_unknown.norm() << std::endl;
		unk_ofile.close();

		//decrease step length
		//step = step*0.75;

		/*if (step < epsilon2)
		{
		break;
		}*/
		if (this->use_sim_or_not == true)
		{
			Eigen::MatrixXd low_v_mat = Eigen::MatrixXd::Zero(sim_mesh.n_vertices(), 3);
#pragma omp parallel for 
			for (integer i = 0; i < sim_mesh.n_vertices(); i++)
			{
#pragma omp critical
				{
					low_v_mat.row(i) << unknown(3 * i + 0), unknown(3 * i + 1), unknown(3 * i + 2);
				}
			}
			high_v_mat = W_mat * low_v_mat;
		}

		con_ele_cnt = 0;
		str_ele_cnt = 0;
		ben_ele_cnt = 0;
		sb_ele_cnt = 0;
		ww_ele_cnt = 0;
		row_cnt = 0;

		++iter;
	}
	clock_t ends = clock();
	std::ofstream time_ofile("running_time.txt");
	time_ofile << (double)start << " " << (double)ends << " " << (double)(ends - start) / CLOCKS_PER_SEC;
	time_ofile.close();

	/*if (u_ofile.is_open())
	{
	u_ofile << unknown;
	}*/
	if (this->use_eg_or_not == false && this->use_basis_or_not == true)
	{
		std::ofstream u_ofile("w_unknown.txt", std::ios::out);
		for (integer i = 0; i < 2 * basis.cols(); i++)
		{
			decimal tmp = unknown(3 * mesh.n_vertices() + i);
			//decimal tmp = (unknown(3 * mesh.n_vertices() + i) >(0.1 / basis.cols())) ? unknown(3 * mesh.n_vertices() + i) : 0;
			u_ofile << tmp << std::endl;
		}
		u_ofile.close();
		//CreateDiLenMat();
	}
	
}
void mns::MyDeformationHandler::DeformLD(std::vector<integer>& s_vid_vec, Eigen::MatrixXd& v_mat, Eigen::MatrixXd& con_v_mat)
{
	DeformInitLD(s_vid_vec, v_mat, con_v_mat);
	//BuildWntVecByBasis(		);
	//Transform2Example();

	std::ofstream energy_ofile("energy.txt", std::ios::out);
	energy_ofile.close();

	std::ofstream unk_ofile("unknowndif.txt", std::ios::out);
	unk_ofile.close();

	Eigen::SimplicialLLT<SMatXd> llt;
	clock_t start = clock();
	while (1)
	{
		decimal Econ = 0;
		if (alpha_c != 0)
		{
			Econ = this->EcNew();
		}
		decimal Estr = 0;
		decimal Eben = 0;
		if (alpha_s != 0 || alpha_b != 0)
		{
			std::vector<decimal> esb = EsbNew();
			Estr = esb[0];
			Eben = esb[1];
		}
		decimal Ewwc = 0;
		if (alpha_w != 0 && this->use_ww_or_not == true)
		{
			Ewwc = this->Ew();
		}
		decimal Evol = 0;
		if (alpha_v != 0)
		{
			Evol = Ev();
		}
		if (isnan(Ewwc))
		{
			unknown = pre_unknown;
			break;
		}

		decimal Energy = Econ;
		std::ofstream energy_ofile("energy.txt", std::ios::app);
		energy_ofile << iter << "-th iteration: " << Econ << " " << Estr << " " << Eben << " " << Energy << std::endl;
		cout << iter << "-th iteration: " << Econ << " " << Estr << " " << Eben << " " << Energy << std::endl;
		energy_ofile.close();

		if (this->pre_energy < Energy)
		{
			std::ofstream energy_ofile("energy.txt", std::ios::app);
			energy_ofile << "iteration ends by 'pre_energy < energy'" << std::endl;
			energy_ofile.close();
			unknown = pre_unknown;
			break;
		}
		if (fitmode == 0)
			epsilon2 = 0.1;
		if (fitmode == 1)
		{
			epsilon2 = 0.001;
			epsilon1 = 0.001;
		}
		//if ((abs(Energy - pre_energy) <= (epsilon1*pre_energy) || (abs(Energy)<epsilon2)) && (this->ch_con == 0))
		if(Econ < 0.2)
		{
			std::ofstream energy_ofile("energy.txt", std::ios::app);
			energy_ofile << "iteration ends by meet relative error threshold." << std::endl;
			energy_ofile.close();
			break;
		}
		if (iter > max_iter)
		{
			break;
		}
		decimal pre_e = this->pre_energy;
		this->pre_energy = Energy;

		SMatXd Jmat;
		SMatXd Jbar;
		Eigen::VectorXd hg;
		BuildJmat(Jmat, Jbar, (alpha_v != 0));

		//Jmat.makeCompressed();
		if (alpha_v == 0)
		{
			SMatXd Jt = Jmat.transpose();
			SMatXd JtJprime = (Jt*Jmat).pruned();

			triple.clear();
			SMatXd sI = SMatXd(JtJprime.rows(), JtJprime.cols());
			for (int i = 0; i<sI.rows(); i++)
			{
				triple.push_back(Eigen::Triplet<decimal>(i, i, 0.00001));
			}
			sI.setFromTriplets(triple.begin(), triple.end());
			triple.clear();

			SMatXd JtJ = (JtJprime + sI).pruned();

			llt.compute(JtJ);
			//jtj = JtJ;
			if (llt.info() != Eigen::Success)
			{
				std::cout << llt.info();
				return;
			}
			hg = llt.solve((-1)*Jt*f);
			if (llt.info() != Eigen::Success)
			{
				return;
			}
		}
		else
		{
			SMatXd Jt = Jmat.transpose();
			SMatXd Jbt = Jbar.transpose();
			SMatXd tmp_A = (Jbt*Jbar).pruned();
			Eigen::MatrixXd b = -Jt * f;

			triple.clear();
			SMatXd sI = SMatXd(tmp_A.rows(), tmp_A.cols());
#pragma omp parallel for 
			for (int i = 0; i<sI.rows(); i++)
			{
#pragma omp critical
				{
					triple.push_back(Eigen::Triplet<decimal>(i, i, 0.00001));
				}
			}

			sI.setFromTriplets(triple.begin(), triple.end());
			triple.clear();

			SMatXd A = tmp_A + sI;

			llt.compute(A);
			if (llt.info() != Eigen::Success)
			{
				return;
			}
			Eigen::VectorXd y = llt.solve(b);
			if (llt.info() != Eigen::Success)
			{
				return;
			}
			if (this->use_sim_or_not == false)
			{
				Eigen::VectorXd z = llt.solve(u_vec);
				if (llt.info() != Eigen::Success)
				{
					return;
				}
				hg = y - ((z*(u_vec.transpose()*y)) / (1 + u_vec.transpose()*z));
			}
			if (this->use_sim_or_not == true)
			{
				Eigen::RowVectorXd u_rowvec_tmp1 = u_vec.transpose();
				Eigen::RowVectorXd u_rowvec_left_tmp = u_rowvec_tmp1.block(0, 0, 1, 3 * this->mesh.n_vertices());
				Eigen::RowVectorXd u_rowvec_right_tmp = u_rowvec_tmp1.block(0, 3 * this->mesh.n_vertices(), 1, u_vec.size() - 3 * this->mesh.n_vertices());
				Eigen::RowVectorXd u_rowvec = Eigen::RowVectorXd::Zero(unknown.size());
				u_rowvec.block(0, 0, 1, 3 * this->sim_mesh.n_vertices()) = u_rowvec_left_tmp * huge_W_mat;
				u_rowvec.block(0, 3 * this->sim_mesh.n_vertices(), 1, unknown.size() - 3 * this->sim_mesh.n_vertices()) = u_rowvec_right_tmp;
				Eigen::VectorXd new_u_vec = u_rowvec.transpose();
				Eigen::VectorXd z = llt.solve(new_u_vec);
				if (llt.info() != Eigen::Success)
				{
					return;
				}
				hg = y - ((z*(new_u_vec.transpose()*y)) / (1 + new_u_vec.transpose()*z));
			}
		}

		decimal step_len = ori_step;
		//

		pre_unknown = unknown;
#pragma omp parallel
		{
#pragma omp for
			for (int i = 0; i < unknown.size(); i++)
			{
#pragma omp critical
				{
					unknown(i) += step_len * hg(i);
				}
			}
		}
		std::ofstream unk_ofile("unknowndif.txt", std::ios::app);
		unk_ofile << "unknown difference iteration between " << iter << " and " << iter + 1 << ": " << unknown.norm() - pre_unknown.norm() << std::endl;
		unk_ofile.close();

		if (this->use_sim_or_not == true)
		{
			Eigen::MatrixXd low_v_mat = Eigen::MatrixXd::Zero(sim_mesh.n_vertices(), 3);
#pragma omp parallel for 
			for (integer i = 0; i < sim_mesh.n_vertices(); i++)
			{
#pragma omp critical
				{
					low_v_mat.row(i) << unknown(3 * i + 0), unknown(3 * i + 1), unknown(3 * i + 2);
				}
			}
			high_v_mat = W_mat * low_v_mat;
		}

		con_ele_cnt = 0;
		str_ele_cnt = 0;
		ben_ele_cnt = 0;
		sb_ele_cnt = 0;
		ww_ele_cnt = 0;
		row_cnt = 0;

		++iter;
	}
	clock_t ends = clock();
	std::ofstream time_ofile("running_time.txt");
	time_ofile << (double)start << " " << (double)ends << " " << (double)(ends - start) / CLOCKS_PER_SEC;
	time_ofile.close();
}

void mns::MyDeformationHandler::DeformInitLD(std::vector<integer>& s_vid_vec, Eigen::MatrixXd& v_mat, Eigen::MatrixXd& con_v_mat)
{
	this->con_ele_cnt = 0;
	this->str_ele_cnt = 0;
	this->ben_ele_cnt = 0;
	this->sb_ele_cnt = 0;
	this->row_cnt = 0;
	this->pre_energy = 100000000000000;
	this->iter = 0;

	//std::ofstream mydi_ofile("mydihedral.txt");
	//std::ofstream wdi_ofile("wdihedral.txt");

	for (auto he_it = this->mesh.halfedges_begin(); he_it != this->mesh.halfedges_end(); ++he_it)
	{
		auto e = *he_it;
		auto fromVertex = this->mesh.from_vertex_handle(e);
		auto toVertex = this->mesh.to_vertex_handle(e);
		auto fromPt = this->mesh.point(fromVertex);
		auto toPt = this->mesh.point(toVertex);

		//stretching initialization		
		decimal l = this->mesh.calc_edge_length(e);
		this->ori_length_set.push_back(l);
		this->wnt_length_set.push_back(l);

		//bend initialization
		auto f1 = this->mesh.face_handle(e);
		auto op_he = this->mesh.opposite_halfedge_handle(e);
		auto f2 = this->mesh.face_handle(op_he);
		pt_vector pt_set;
		pt_vector t1_set;
		pt_vector t2_set;

		if (f1.idx() == -1 || f2.idx() == -1)
		{
			this->ori_dihedral_set.push_back(0);
			this->wnt_dihedral_set.push_back(0);
			ori_area_set.push_back(0);

			this->ori_wlz_dihedral_set.push_back(0);
			this->wnt_wlz_dihedral_set.push_back(0);
			sign_flag.push_back(1);

			//edge, not half edge
			++he_it;
			continue;
		}

		pt_set.push_back(fromPt);
		t1_set.push_back(fromPt);
		t2_set.push_back(fromPt);
		pt_set.push_back(toPt);
		t1_set.push_back(toPt);
		t2_set.push_back(toPt);

		for (auto fv_it1 = this->mesh.fv_begin(f1); fv_it1 != this->mesh.fv_end(f1); ++fv_it1)
		{
			auto v = *fv_it1;
			auto pt = this->mesh.point(v);
			if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
			{
				pt_set.push_back(pt);
				t1_set.push_back(pt);
				break;
			}
		}
		integer temp_cnt = 0;
		for (auto fv_it2 = this->mesh.fv_begin(f2); fv_it2 != this->mesh.fv_end(f2); ++fv_it2)
		{
			auto v = *fv_it2;
			auto pt = this->mesh.point(v);
			if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
			{
				t2_set.push_back(pt);
				pt_set.push_back(pt);
				break;
			}
		}
		decimal t1_side[3];
		t1_side[0] = (t1_set[0] - t1_set[1]).length();
		t1_side[1] = (t1_set[2] - t1_set[1]).length();
		t1_side[2] = (t1_set[0] - t1_set[2]).length();

		decimal p1 = (t1_side[0] + t1_side[1] + t1_side[2]) / 2; //半周长;  
		decimal area1 = std::sqrt(p1*(p1 - t1_side[0])*(p1 - t1_side[1])*(p1 - t1_side[2]));

		decimal t2_side[3];
		t2_side[0] = (t2_set[0] - t2_set[1]).length();
		t2_side[1] = (t2_set[2] - t2_set[1]).length();
		t2_side[2] = (t2_set[0] - t2_set[2]).length();

		decimal p2 = (t2_side[0] + t2_side[1] + t2_side[2]) / 2; //半周长;  
		decimal area2 = std::sqrt(p2*(p2 - t2_side[0])*(p2 - t2_side[1])*(p2 - t2_side[2]));
		ori_area_set.push_back(area1 + area2);

		decimal angle = CalDiheral(pt_set[0], pt_set[1], pt_set[2], pt_set[3]);
		this->ori_dihedral_set.push_back(angle);
		this->wnt_dihedral_set.push_back(angle);

		//mydi_ofile << angle << std::endl;

		//sign_flag computation
		Eigen::Vector3d normal1;
		Eigen::Vector3d normal2;
		{
			auto fhe_it1 = this->mesh.fh_begin(f1);
			auto he1 = *fhe_it1;
			auto he1idx = he1.idx();
			auto fromVertex1 = this->mesh.from_vertex_handle(he1);
			auto toVertex1 = this->mesh.to_vertex_handle(he1);
			auto fromPt1 = this->mesh.point(fromVertex1);
			auto toPt1 = this->mesh.point(toVertex1);
			auto e1 = toPt1 - fromPt1;
			Eigen::Vector3d ev1;
			ev1 << e1.data()[0], e1.data()[1], e1.data()[2];
			++fhe_it1;
			auto he2 = *fhe_it1;
			auto he2idx = he2.idx();
			auto fromVertex2 = this->mesh.from_vertex_handle(he2);
			auto toVertex2 = this->mesh.to_vertex_handle(he2);
			auto fromPt2 = this->mesh.point(fromVertex2);
			auto toPt2 = this->mesh.point(toVertex2);
			auto e2 = toPt2 - fromPt2;
			Eigen::Vector3d ev2;
			ev2 << e2.data()[0], e2.data()[1], e2.data()[2];
			normal1 = ev1.cross(ev2);
			normal1 = normal1 / normal1.norm();
		}
		{
			auto fhe_it2 = this->mesh.fh_begin(f2);
			auto he1 = *fhe_it2;
			auto he1idx = he1.idx();
			auto fromVertex1 = this->mesh.from_vertex_handle(he1);
			auto toVertex1 = this->mesh.to_vertex_handle(he1);
			auto fromPt1 = this->mesh.point(fromVertex1);
			auto toPt1 = this->mesh.point(toVertex1);
			auto e1 = toPt1 - fromPt1;
			Eigen::Vector3d ev1;
			ev1 << e1.data()[0], e1.data()[1], e1.data()[2];
			++fhe_it2;
			auto he2 = *fhe_it2;
			auto he2idx = he2.idx();
			auto fromVertex2 = this->mesh.from_vertex_handle(he2);
			auto toVertex2 = this->mesh.to_vertex_handle(he2);
			auto fromPt2 = this->mesh.point(fromVertex2);
			auto toPt2 = this->mesh.point(toVertex2);
			auto e2 = toPt2 - fromPt2;
			Eigen::Vector3d ev2;
			ev2 << e2.data()[0], e2.data()[1], e2.data()[2];
			normal2 = ev1.cross(ev2);
			normal2 = normal2 / normal2.norm();
		}

		auto e_ij_tmp = toPt - fromPt;
		Eigen::Vector3d e_ij;
		e_ij << e_ij_tmp.data()[0], e_ij_tmp.data()[1], e_ij_tmp.data()[2];
		e_ij = e_ij / (e_ij.norm());

		decimal phi = CalWLZHDiheral(normal1, normal2, e_ij);
		this->ori_wlz_dihedral_set.push_back(phi);
		this->wnt_wlz_dihedral_set.push_back(phi);
		/*this->ori_dihedral_set.push_back(phi);
		this->wnt_dihedral_set.push_back(phi);*/

		if (phi > 0)
		{
			sign_flag.push_back(1);
		}
		else
		{
			sign_flag.push_back(-1);
		}

		//wdi_ofile << phi << std::endl;

		//edge, not half edge
		++he_it;
	}
	//mydi_ofile.close();
	//wdi_ofile.close();

	for (integer i = 0; i < this->ori_dihedral_set.size(); i++)
	{
		decimal ws = alpha_s / (ori_length_set[i] * ori_length_set[i]);
		decimal sqr_ws = std::sqrt(ws);
		this->sws_set.push_back(sqr_ws);
		decimal wb = (alpha_b * ori_length_set[i] * ori_length_set[i]) / ori_area_set[i];
		decimal sqr_wb = std::sqrt(wb);
		this->swb_set.push_back(sqr_wb);
		//cout << sqr_wb << endl;
	}

	//volume initialization
	ori_volume = 0.0;
	if (alpha_v != 0)
	{
		for (auto f_it = this->mesh.faces_begin(); f_it != this->mesh.faces_end(); ++f_it)
		{
			auto face = *f_it;
			pt_vector pt_vec;
			for (auto fv_it = this->mesh.fv_begin(face); fv_it != this->mesh.fv_end(face); ++fv_it)
			{
				auto v = *fv_it;
				auto pt = this->mesh.point(v);
				pt_vec.push_back(pt);
			}
			decimal temp_vol = (1.0 / 6.0) * ((pt_vec[0] % pt_vec[1]) | pt_vec[2]);
			ori_volume += temp_vol;
		}
		wnt_volume = ori_volume;
		decimal wv = alpha_v / (ori_volume * ori_volume);
		swv = std::sqrt(wv);
	}

	//average
	Eigen::VectorXd avg = Eigen::VectorXd::Zero(2 * this->mesh.n_edges());
#pragma omp parallel
	{
#pragma omp for
		for (integer i = 0; i < this->mesh.n_edges(); i++)
		{
			avg(2 * i + 1) = ori_dihedral_set[i];
			avg(2 * i + 0) = ori_length_set[i];
		}
	}

	this->data_avg = avg;
	Eigen::VectorXd tmp_avg = Eigen::VectorXd::Zero(2 * this->mesh.n_edges() + 1);
	if (alpha_v != 0)
	{
		tmp_avg.block(0, 0, avg.rows(), 1) = avg;
		tmp_avg(avg.rows()) = ori_volume;
		this->data_avg.resize(2 * this->mesh.n_edges() + 1);
		this->data_avg = tmp_avg;
	}


	//constraint initialization
	
	//PrintPtControl();
	decimal con_set_num = ((decimal)body_feature_idx.size());
	decimal wc = alpha_c / con_set_num;
	swc = std::sqrt(wc);
	integer cnt = 0;
	bool cnt_flag = true;
	this->move_pt_num = 0;

	//build unknown	
	unknown = Eigen::VectorXd::Zero(3 * this->mesh.n_vertices());

	int i = 0;
	for (auto it = this->mesh.vertices_begin(); it != this->mesh.vertices_end(); ++it, ++i)
	{
		auto point = this->mesh.point(*it);


		//soft con
		unknown(3 * i + 0) = point.data()[0];
		unknown(3 * i + 1) = point.data()[1];
		unknown(3 * i + 2) = point.data()[2];

	}

	pre_unknown = unknown;
}
void mns::MyDeformationHandler::Deform1(std::vector<integer>& s_vid_vec, Eigen::MatrixXd& v_mat, Eigen::MatrixXd& con_v_mat)
{
	//DeformInit(s_vid_vec, v_mat, con_v_mat);
	//BuildWntVecByBasis();
	//Transform2Example();

	std::ofstream energy_ofile("`11energy.txt", std::ios::out);
	energy_ofile.close();

	std::ofstream unk_ofile("unknowndif.txt", std::ios::out);
	unk_ofile.close();

	Eigen::SimplicialLLT<SMatXd> llt;
	clock_t start = clock();
	while (1)
	{
		alpha_c = alpha_c * 1;
		alpha_s = alpha_s * 1;
		alpha_b = alpha_b * 1;
		//mod unknown by basis
		if (this->cal_bw_or_not == true && this->use_eg_or_not == false && this->use_basis_or_not == true)
		{
			CalWunknownByBasis();
		}
		if (this->cal_bw_or_not == true && this->use_eg_or_not == true && this->use_basis_or_not == true)
		{
			CalWunknownByBasis2();
		}
		//

		if (this->use_eg_or_not == true && this->use_basis_or_not == false)
		{
			UpdateWntVec();
		}
		if (this->use_eg_or_not == false && this->use_basis_or_not == true)
		{
			UpdateWntVecWithBasis7();
		}
		if (this->use_eg_or_not == true && this->use_basis_or_not == true)
		{
			UpdateWntVecWithBasis3();
		}

		decimal Econ = 0;
		if (alpha_c != 0)
		{
			Econ = this->EcNew();
		}
		decimal Estr = 0;
		/*if (alpha_s != 0)
		{
		Estr = this->Es();
		}*/
		decimal Eben = 0;
		/*if (alpha_b != 0)
		{
		Eben = this->Eb();
		}*/
		if (alpha_s != 0 || alpha_b != 0)
		{
			//std::vector<decimal> esb = Esb();
			std::vector<decimal> esb = EsbNew1();
			Estr = esb[0];
			Eben = esb[1];
		}
		decimal Ewwc = 0;
		if (alpha_w != 0 && this->use_ww_or_not == true)
		{
			Ewwc = this->Ew();
		}
		decimal Evol = 0;
		if (alpha_v != 0)
		{
			Evol = Ev();
		}
		if (isnan(Ewwc))
		{
			unknown = pre_unknown;
			break;
		}

		decimal Energy = Econ + Evol + Ewwc;
		//decimal Energy = Econ + Estr;
		std::ofstream energy_ofile("energy.txt", std::ios::app);
		energy_ofile << iter << "-th iteration: " << Econ << " " << Estr << " " << Eben << " " << Energy << std::endl;
		energy_ofile.close();
		cout << iter << "-th iteration: " << Econ << " " << Estr << " " << Eben << " " << Energy << std::endl;

		//YP
		if (this->pre_energy < Energy)
		{
			std::ofstream energy_ofile("energy.txt", std::ios::app);
			energy_ofile << "iteration ends by 'pre_energy < energy'" << std::endl;
			energy_ofile.close();
			cout << "iteration ends by 'pre_energy < energy'" << std::endl;
			unknown = pre_unknown;

			if (this->cal_bw_or_not == true && this->use_eg_or_not == false && this->use_basis_or_not == true)
			{
				CalWunknownByBasis();
			}
			if (this->cal_bw_or_not == true && this->use_eg_or_not == true && this->use_basis_or_not == true)
			{
				CalWunknownByBasis2();
			}
			break;
		}

		/*if (this->pre_energy < Energy)
		{
		unknown = pre_unknown;

		if (this->ch_con > 0)
		{
		this->ch_con = 0;
		}
		else
		{
		if ((alpha_c_ch == false) && (need_ch[0] == true))
		{
		Alpha_C_Change(0.5 / need_fac);
		Clear4ReCal();
		this->ch_con += 1;
		}
		if ((alpha_s_ch == false) && (need_ch[1] == true))
		{
		Alpha_S_Change(need_fac);
		Clear4ReCal();
		this->ch_con += 1;
		}
		if ((alpha_b_ch == false) && (need_ch[2] == true))
		{
		Alpha_B_Change(need_fac);
		Clear4ReCal();
		this->ch_con += 1;
		}
		if ((alpha_v_ch == false) && (need_ch[3] == true))
		{
		Alpha_V_Change(need_fac);
		Clear4ReCal();
		this->ch_con += 1;
		}
		if ((alpha_w_ch == false) && (need_ch[4] == true))
		{
		Alpha_W_Change(need_fac);
		Clear4ReCal();
		this->ch_con += 1;
		}
		if (Is_Alpha_Change() && (this->ch_con > 0))
		{
		continue;
		}
		if (Is_Alpha_Change() && (this->ch_con == 0))
		{
		break;
		}
		}
		}*/

		//YP
		if ((abs(Energy - pre_energy) <= (epsilon1*pre_energy) || (abs(Energy)<epsilon2)) && (this->ch_con == 0) || (Econ < econ_threshold))
		{
			std::ofstream energy_ofile("energy.txt", std::ios::app);
			energy_ofile << "iteration ends by meet relative error threshold." << std::endl;
			energy_ofile.close();
			break;
		}

		if (iter > max_iter)
		{
			break;
		}
		decimal pre_e = this->pre_energy;
		this->pre_energy = Energy;

		SMatXd Jmat;
		SMatXd Jbar;
		Eigen::VectorXd hg;
		BuildJmat(Jmat, Jbar, (alpha_v != 0));

		//Jmat.makeCompressed();
		if (alpha_v == 0)
		{
			SMatXd Jt = Jmat.transpose();
			SMatXd JtJprime = (Jt*Jmat).pruned();

			triple.clear();
			SMatXd sI = SMatXd(JtJprime.rows(), JtJprime.cols());
			for (int i = 0; i<sI.rows(); i++)
			{
				triple.push_back(Eigen::Triplet<decimal>(i, i, 0.00001));
			}
			sI.setFromTriplets(triple.begin(), triple.end());
			triple.clear();

			SMatXd JtJ = (JtJprime + sI).pruned();

			/*std::ofstream j_ofile("Jtj.txt", std::ios::out);
			if (j_ofile.is_open())
			{
			j_ofile << JtJ;
			}
			j_ofile.close();*/

			llt.compute(JtJ);
			//jtj = JtJ;
			if (llt.info() != Eigen::Success)
			{
				std::cout << llt.info();
				return;
			}
			hg = llt.solve((-1)*Jt*f);
			/*std::ofstream jf_ofile("jtf.txt", std::ios::out);
			if (jf_ofile.is_open())
			{
			jf_ofile << Jt*f;
			}
			jf_ofile.close();
			std::ofstream h_ofile("hg.txt", std::ios::out);
			if (h_ofile.is_open())
			{
			h_ofile << hg;
			}
			h_ofile.close();*/
			if (llt.info() != Eigen::Success)
			{
				return;
			}
		}
		else
		{
			SMatXd Jt = Jmat.transpose();
			SMatXd Jbt = Jbar.transpose();
			SMatXd tmp_A = (Jbt*Jbar).pruned();
			Eigen::MatrixXd b = -Jt * f;

			triple.clear();
			SMatXd sI = SMatXd(tmp_A.rows(), tmp_A.cols());
#pragma omp parallel for 
			for (int i = 0; i<sI.rows(); i++)
			{
#pragma omp critical
				{
					triple.push_back(Eigen::Triplet<decimal>(i, i, 0.00001));
				}
			}

			sI.setFromTriplets(triple.begin(), triple.end());
			triple.clear();

			SMatXd A = tmp_A + sI;

			llt.compute(A);
			//jtj = JtJ;
			if (llt.info() != Eigen::Success)
			{
				//std::cout << llt.info();
				return;
			}
			Eigen::VectorXd y = llt.solve(b);
			if (llt.info() != Eigen::Success)
			{
				return;
			}
			if (this->use_sim_or_not == false)
			{
				Eigen::VectorXd z = llt.solve(u_vec);
				if (llt.info() != Eigen::Success)
				{
					return;
				}
				hg = y - ((z*(u_vec.transpose()*y)) / (1 + u_vec.transpose()*z));
			}
			if (this->use_sim_or_not == true)
			{
				Eigen::RowVectorXd u_rowvec_tmp1 = u_vec.transpose();
				Eigen::RowVectorXd u_rowvec_left_tmp = u_rowvec_tmp1.block(0, 0, 1, 3 * this->mesh.n_vertices());
				Eigen::RowVectorXd u_rowvec_right_tmp = u_rowvec_tmp1.block(0, 3 * this->mesh.n_vertices(), 1, u_vec.size() - 3 * this->mesh.n_vertices());
				Eigen::RowVectorXd u_rowvec = Eigen::RowVectorXd::Zero(unknown.size());
				u_rowvec.block(0, 0, 1, 3 * this->sim_mesh.n_vertices()) = u_rowvec_left_tmp * huge_W_mat;
				u_rowvec.block(0, 3 * this->sim_mesh.n_vertices(), 1, unknown.size() - 3 * this->sim_mesh.n_vertices()) = u_rowvec_right_tmp;
				Eigen::VectorXd new_u_vec = u_rowvec.transpose();
				Eigen::VectorXd z = llt.solve(new_u_vec);
				if (llt.info() != Eigen::Success)
				{
					return;
				}
				hg = y - ((z*(new_u_vec.transpose()*y)) / (1 + new_u_vec.transpose()*z));
			}
		}

		//calculate optimal step length
		//decimal step_len = (iter == 0) ? ori_step : CalStepLength(pre_e, Jmat, this->f, hg, this->unknown);
		decimal step_len = ori_step;
		//

		pre_unknown = unknown;
#pragma omp parallel
		{
#pragma omp for
			for (int i = 0; i < unknown.size(); i++)
			{
#pragma omp critical
				{
					unknown(i) += step_len * hg(i);
				}
			}
		}
		std::ofstream unk_ofile("unknowndif.txt", std::ios::app);
		unk_ofile << "unknown difference iteration between " << iter << " and " << iter + 1 << ": " << unknown.norm() - pre_unknown.norm() << std::endl;
		unk_ofile.close();

		//decrease step length
		//step = step*0.75;

		/*if (step < epsilon2)
		{
		break;
		}*/
		if (this->use_sim_or_not == true)
		{
			Eigen::MatrixXd low_v_mat = Eigen::MatrixXd::Zero(sim_mesh.n_vertices(), 3);
#pragma omp parallel for 
			for (integer i = 0; i < sim_mesh.n_vertices(); i++)
			{
#pragma omp critical
				{
					low_v_mat.row(i) << unknown(3 * i + 0), unknown(3 * i + 1), unknown(3 * i + 2);
				}
			}
			high_v_mat = W_mat * low_v_mat;
		}

		con_ele_cnt = 0;
		str_ele_cnt = 0;
		ben_ele_cnt = 0;
		sb_ele_cnt = 0;
		ww_ele_cnt = 0;
		row_cnt = 0;

		++iter;
	}
	clock_t ends = clock();
	std::ofstream time_ofile("running_time.txt");
	time_ofile << (double)start << " " << (double)ends << " " << (double)(ends - start) / CLOCKS_PER_SEC;
	time_ofile.close();

	/*if (u_ofile.is_open())
	{
	u_ofile << unknown;
	}*/
	if (this->use_eg_or_not == false && this->use_basis_or_not == true)
	{
		std::ofstream u_ofile("w_unknown.txt", std::ios::out);
		for (integer i = 0; i < 2 * basis.cols(); i++)
		{
			decimal tmp = unknown(3 * mesh.n_vertices() + i);
			//decimal tmp = (unknown(3 * mesh.n_vertices() + i) >(0.1 / basis.cols())) ? unknown(3 * mesh.n_vertices() + i) : 0;
			u_ofile << tmp << std::endl;
		}
		u_ofile.close();
		//CreateDiLenMat();
	}
	if (this->use_eg_or_not == true && this->use_basis_or_not == true)
	{
		std::ofstream l_ofile("l_bw_unknown.txt", std::ios::out);
		std::ofstream d_ofile("d_bw_unknown.txt", std::ios::out);
		for (integer i = 0; i < basis.cols(); i++)
		{
			//decimal tmp = (unknown(3 * mesh.n_vertices() + i) >(0.1 / basis.cols())) ? unknown(3 * mesh.n_vertices() + i) : 0;
			l_ofile << l_bw_unknown(i) << std::endl;
			d_ofile << d_bw_unknown(i) << std::endl;
		}
		l_ofile.close();
		d_ofile.close();
		//CreateDiLenMat2();
	}
	if (this->use_eg_or_not == true && this->use_basis_or_not == false)
	{
		std::ofstream u_ofile("eg_w_unknown.txt", std::ios::out);
		for (integer i = 0; i < this->eg_num; i++)
		{
			u_ofile << unknown(3 * mesh.n_vertices() + i) << " ";
		}
		u_ofile.close();
	}

}


void mns::MyDeformationHandler::Deform(std::vector<integer>& s_vid_vec, Eigen::MatrixXd& v_mat, Eigen::MatrixXd& con_v_mat)
{
	//DeformInit(s_vid_vec, v_mat, con_v_mat);
	//BuildWntVecByBasis();
	//Transform2Example();

	std::ofstream energy_ofile("`11energy.txt", std::ios::out);
	energy_ofile.close();

	std::ofstream unk_ofile("unknowndif.txt", std::ios::out);
	unk_ofile.close();

	Eigen::SimplicialLLT<SMatXd> llt;
	clock_t start = clock();
	while (1)
	{
		alpha_c = alpha_c * 1;
		alpha_s = alpha_s * 1;
		alpha_b = alpha_b * 1;
		//mod unknown by basis
		if (this->cal_bw_or_not == true && this->use_eg_or_not == false && this->use_basis_or_not == true)
		{
			CalWunknownByBasis();
		}
		if (this->cal_bw_or_not == true && this->use_eg_or_not == true && this->use_basis_or_not == true)
		{
			CalWunknownByBasis2();
		}
		//

		if (this->use_eg_or_not == true && this->use_basis_or_not == false)
		{
			UpdateWntVec();
		}
		if (this->use_eg_or_not == false && this->use_basis_or_not == true)
		{
			UpdateWntVecWithBasis2();
		}
		if (this->use_eg_or_not == true && this->use_basis_or_not == true)
		{
			UpdateWntVecWithBasis3();
		}

		decimal Econ = 0;
		if (alpha_c != 0)
		{
			Econ = this->EcNew();
		}
		decimal Estr = 0;
		/*if (alpha_s != 0)
		{
		Estr = this->Es();
		}*/
		decimal Eben = 0;
		/*if (alpha_b != 0)
		{
		Eben = this->Eb();
		}*/
		if (alpha_s != 0 || alpha_b != 0)
		{
			//std::vector<decimal> esb = Esb();
			std::vector<decimal> esb = EsbNew();
			Estr = esb[0];
			Eben = esb[1];
		}
		decimal Ewwc = 0;
		if (alpha_w != 0 && this->use_ww_or_not == true)
		{
			Ewwc = this->Ew();
		}
		decimal Evol = 0;
		if (alpha_v != 0)
		{
			Evol = Ev();
		}
		if (isnan(Ewwc))
		{
			unknown = pre_unknown;
			break;
		}

		//decimal Energy =  Econ + Estr + Evol + Ewwc;
		decimal Energy = Econ + 0.01 * Estr + 0.0001 * Eben;
		std::ofstream energy_ofile("energy.txt", std::ios::app);
		energy_ofile << iter << "-th iteration: " << Econ << " " << Estr << " " << Eben << " " << Energy << std::endl;
		energy_ofile.close();
		cout << iter << "-th iteration: " << Econ << " " << Estr << " " << Eben << " " << Energy << std::endl;
		
		//YP
		if (this->pre_energy < Energy)
		{
			std::ofstream energy_ofile("energy.txt", std::ios::app);
			energy_ofile << "iteration ends by 'pre_energy < energy'" << std::endl;
			energy_ofile.close();
			cout << "iteration ends by 'pre_energy < energy'" << std::endl;
			unknown = pre_unknown;
			
			if (this->cal_bw_or_not == true && this->use_eg_or_not == false && this->use_basis_or_not == true)
			{
				CalWunknownByBasis();
			}
			if (this->cal_bw_or_not == true && this->use_eg_or_not == true && this->use_basis_or_not == true)
			{
				CalWunknownByBasis2();
			}
			break;
		}

		/*if (this->pre_energy < Energy)
		{
		unknown = pre_unknown;

		if (this->ch_con > 0)
		{
		this->ch_con = 0;
		}
		else
		{
		if ((alpha_c_ch == false) && (need_ch[0] == true))
		{
		Alpha_C_Change(0.5 / need_fac);
		Clear4ReCal();
		this->ch_con += 1;
		}
		if ((alpha_s_ch == false) && (need_ch[1] == true))
		{
		Alpha_S_Change(need_fac);
		Clear4ReCal();
		this->ch_con += 1;
		}
		if ((alpha_b_ch == false) && (need_ch[2] == true))
		{
		Alpha_B_Change(need_fac);
		Clear4ReCal();
		this->ch_con += 1;
		}
		if ((alpha_v_ch == false) && (need_ch[3] == true))
		{
		Alpha_V_Change(need_fac);
		Clear4ReCal();
		this->ch_con += 1;
		}
		if ((alpha_w_ch == false) && (need_ch[4] == true))
		{
		Alpha_W_Change(need_fac);
		Clear4ReCal();
		this->ch_con += 1;
		}
		if (Is_Alpha_Change() && (this->ch_con > 0))
		{
		continue;
		}
		if (Is_Alpha_Change() && (this->ch_con == 0))
		{
		break;
		}
		}
		}*/
		
		//YP
		epsilon1 = 0.01;
		if ((abs(Energy - pre_energy) <= (epsilon1*pre_energy) || (abs(Energy)<epsilon2)) && (this->ch_con == 0) || (Econ < econ_threshold))
		{
			std::ofstream energy_ofile("energy.txt", std::ios::app);
			energy_ofile << "iteration ends by meet relative error threshold." << std::endl;
			energy_ofile.close();
			break;
		}

		if (iter > max_iter)
		{
			break;
		}
		decimal pre_e = this->pre_energy;
		this->pre_energy = Energy;

		SMatXd Jmat;
		SMatXd Jbar;
		Eigen::VectorXd hg;
		BuildJmat(Jmat, Jbar, (alpha_v != 0));

		//Jmat.makeCompressed();
		if (alpha_v == 0)
		{
			SMatXd Jt = Jmat.transpose();
			SMatXd JtJprime = (Jt*Jmat).pruned();

			triple.clear();
			SMatXd sI = SMatXd(JtJprime.rows(), JtJprime.cols());
			for (int i = 0; i<sI.rows(); i++)
			{
				triple.push_back(Eigen::Triplet<decimal>(i, i, 0.00001));
			}
			sI.setFromTriplets(triple.begin(), triple.end());
			triple.clear();

			SMatXd JtJ = (JtJprime + sI).pruned();

			/*std::ofstream j_ofile("Jtj.txt", std::ios::out);
			if (j_ofile.is_open())
			{
			j_ofile << JtJ;
			}
			j_ofile.close();*/

			llt.compute(JtJ);
			//jtj = JtJ;
			if (llt.info() != Eigen::Success)
			{
				std::cout << llt.info();
				return;
			}
			hg = llt.solve((-1) * Jt * f);
			/*std::ofstream jf_ofile("jtf.txt", std::ios::out);
			if (jf_ofile.is_open())
			{
			jf_ofile << Jt*f;
			}
			jf_ofile.close();
			std::ofstream h_ofile("hg.txt", std::ios::out);
			if (h_ofile.is_open())
			{
			h_ofile << hg;
			}
			h_ofile.close();*/
			if (llt.info() != Eigen::Success)
			{
				return;
			}
		}
		else
		{
			SMatXd Jt = Jmat.transpose();
			SMatXd Jbt = Jbar.transpose();
			SMatXd tmp_A = (Jbt*Jbar).pruned();
			Eigen::MatrixXd b = -Jt*f;

			triple.clear();
			SMatXd sI = SMatXd(tmp_A.rows(), tmp_A.cols());
#pragma omp parallel for 
			for (int i = 0; i<sI.rows(); i++)
			{
#pragma omp critical
				{
					triple.push_back(Eigen::Triplet<decimal>(i, i, 0.00001));
				}
			}

			sI.setFromTriplets(triple.begin(), triple.end());
			triple.clear();

			SMatXd A = tmp_A + sI;

			llt.compute(A);
			//jtj = JtJ;
			if (llt.info() != Eigen::Success)
			{
				//std::cout << llt.info();
				return;
			}
			Eigen::VectorXd y = llt.solve(b);
			if (llt.info() != Eigen::Success)
			{
				return;
			}
			if (this->use_sim_or_not == false)
			{
				Eigen::VectorXd z = llt.solve(u_vec);
				if (llt.info() != Eigen::Success)
				{
					return;
				}
				hg = y - ((z*(u_vec.transpose()*y)) / (1 + u_vec.transpose()*z));
			}
			if (this->use_sim_or_not == true)
			{
				Eigen::RowVectorXd u_rowvec_tmp1 = u_vec.transpose();
				Eigen::RowVectorXd u_rowvec_left_tmp = u_rowvec_tmp1.block(0, 0, 1, 3 * this->mesh.n_vertices());
				Eigen::RowVectorXd u_rowvec_right_tmp = u_rowvec_tmp1.block(0, 3 * this->mesh.n_vertices(), 1, u_vec.size() - 3 * this->mesh.n_vertices());
				Eigen::RowVectorXd u_rowvec = Eigen::RowVectorXd::Zero(unknown.size());
				u_rowvec.block(0, 0, 1, 3 * this->sim_mesh.n_vertices()) = u_rowvec_left_tmp * huge_W_mat;
				u_rowvec.block(0, 3 * this->sim_mesh.n_vertices(), 1, unknown.size() - 3 * this->sim_mesh.n_vertices()) = u_rowvec_right_tmp;
				Eigen::VectorXd new_u_vec = u_rowvec.transpose();
				Eigen::VectorXd z = llt.solve(new_u_vec);
				if (llt.info() != Eigen::Success)
				{
					return;
				}
				hg = y - ((z*(new_u_vec.transpose()*y)) / (1 + new_u_vec.transpose()*z));
			}
		}

		//calculate optimal step length
		//decimal step_len = (iter == 0) ? ori_step : CalStepLength(pre_e, Jmat, this->f, hg, this->unknown);
		decimal step_len = ori_step;
		//

		pre_unknown = unknown;
#pragma omp parallel
		{
#pragma omp for
			for (int i = 0; i < unknown.size(); i++)
			{
#pragma omp critical
				{
					unknown(i) += step_len*hg(i);
				}
			}
		}
		std::ofstream unk_ofile("unknowndif.txt", std::ios::app);
		unk_ofile << "unknown difference iteration between " << iter << " and " << iter + 1 << ": " << unknown.norm() - pre_unknown.norm() << std::endl;
		unk_ofile.close();

		//decrease step length
		//step = step*0.75;

		/*if (step < epsilon2)
		{
		break;
		}*/
		if (this->use_sim_or_not == true)
		{
			Eigen::MatrixXd low_v_mat = Eigen::MatrixXd::Zero(sim_mesh.n_vertices(), 3);
#pragma omp parallel for 
			for (integer i = 0; i < sim_mesh.n_vertices(); i++)
			{
#pragma omp critical
				{
					low_v_mat.row(i) << unknown(3 * i + 0), unknown(3 * i + 1), unknown(3 * i + 2);
				}
			}
			high_v_mat = W_mat * low_v_mat;
		}

		con_ele_cnt = 0;
		str_ele_cnt = 0;
		ben_ele_cnt = 0;
		sb_ele_cnt = 0;
		ww_ele_cnt = 0;
		row_cnt = 0;

		++iter;
	}
	clock_t ends = clock();
	std::ofstream time_ofile("running_time.txt");
	time_ofile << (double)start << " " << (double)ends << " " << (double)(ends - start) / CLOCKS_PER_SEC;
	time_ofile.close();

	/*if (u_ofile.is_open())
	{
	u_ofile << unknown;
	}*/
	if (this->use_eg_or_not == false && this->use_basis_or_not == true)
	{
		std::ofstream u_ofile("w_unknown.txt", std::ios::out);
		for (integer i = 0; i < 2 * basis.cols(); i++)
		{
			decimal tmp = unknown(3 * mesh.n_vertices() + i);
			//decimal tmp = (unknown(3 * mesh.n_vertices() + i) >(0.1 / basis.cols())) ? unknown(3 * mesh.n_vertices() + i) : 0;
			u_ofile << tmp << std::endl;
		}
		u_ofile.close();
		//CreateDiLenMat();
	}
	if (this->use_eg_or_not == true && this->use_basis_or_not == true)
	{
		std::ofstream l_ofile("l_bw_unknown.txt", std::ios::out);
		std::ofstream d_ofile("d_bw_unknown.txt", std::ios::out);
		for (integer i = 0; i < basis.cols(); i++)
		{
			//decimal tmp = (unknown(3 * mesh.n_vertices() + i) >(0.1 / basis.cols())) ? unknown(3 * mesh.n_vertices() + i) : 0;
			l_ofile << l_bw_unknown(i) << std::endl;
			d_ofile << d_bw_unknown(i) << std::endl;
		}
		l_ofile.close();
		d_ofile.close();
		//CreateDiLenMat2();
	}
	if (this->use_eg_or_not == true && this->use_basis_or_not == false)
	{
		std::ofstream u_ofile("eg_w_unknown.txt", std::ios::out);
		for (integer i = 0; i < this->eg_num; i++)
		{
			u_ofile << unknown(3 * mesh.n_vertices() + i) << " ";
		}
		u_ofile.close();
	}
	
}


void mns::MyDeformationHandler::ModVMat(Eigen::MatrixXd& v_mat)
{
	for (integer i = 0; i < v_mat.rows(); i++)
	{
		//soft con
		v_mat.row(i) << unknown(3 * i + 0), unknown(3 * i + 1), unknown(3 * i + 2);

		//hard con
		/*if (idx_move[i] != -1)
		{
		v_mat.row(i) << unknown(3 * (i - idx_move[i]) + 0), unknown(3 * (i - idx_move[i]) + 1), unknown(3 * (i - idx_move[i]) + 2);
		}*/
	}
}

void mns::MyDeformationHandler::VmatModAnotherMesh(const Eigen::MatrixXd& v_mat)
{
	int i = 0;
	if (this->use_sim_or_not == true)
	{
		if (!this->mesh.has_vertex_status())
			this->mesh.request_vertex_status();
		if (!this->mesh.has_face_status())
			this->mesh.request_face_status();
		if (!this->mesh.has_edge_status())
			this->mesh.request_edge_status();
		if (!this->mesh.has_halfedge_status())
			this->mesh.has_halfedge_status();

		for (auto v_it = this->mesh.vertices_begin(); v_it != this->mesh.vertices_end(); ++v_it, ++i)
		{
			//soft con
			MyMesh::Point new_v(v_mat(i, 0), v_mat(i, 1), v_mat(i, 2));
			auto vertex = *v_it;
			this->mesh.set_point(vertex, new_v);

			//hard con
			/*if (idx_move[i] != -1)
			{
			MyMesh::Point new_v(unknown(3 * (i - idx_move[i]) + 0), unknown(3 * (i - idx_move[i]) + 1), unknown(3 * (i - idx_move[i]) + 2));
			auto vertex = *v_it;
			this->mesh.set_point(vertex, new_v);
			}*/
		}

		if (this->mesh.has_vertex_status())
			this->mesh.release_vertex_status();/**/
		if (this->mesh.has_face_status())
			this->mesh.release_face_status();
		if (this->mesh.has_edge_status())
			this->mesh.release_edge_status();
		if (this->mesh.has_halfedge_status())
			this->mesh.release_halfedge_status();
	}
}

void mns::MyDeformationHandler::ModMesh(/*MyMesh& mesh*/)
{
	int i = 0;
	if (this->use_sim_or_not == false)
	{
		if (!this->mesh.has_vertex_status())
			this->mesh.request_vertex_status();
		if (!this->mesh.has_face_status())
			this->mesh.request_face_status();
		if (!this->mesh.has_edge_status())
			this->mesh.request_edge_status();
		if (!this->mesh.has_halfedge_status())
			this->mesh.has_halfedge_status();

		for (auto v_it = this->mesh.vertices_begin(); v_it != this->mesh.vertices_end(); ++v_it, ++i)
		{
			//soft con
			MyMesh::Point new_v(unknown(3 * i + 0), unknown(3 * i + 1), unknown(3 * i + 2));
			auto vertex = *v_it;
			this->mesh.set_point(vertex, new_v);

			//hard con
			/*if (idx_move[i] != -1)
			{
			MyMesh::Point new_v(unknown(3 * (i - idx_move[i]) + 0), unknown(3 * (i - idx_move[i]) + 1), unknown(3 * (i - idx_move[i]) + 2));
			auto vertex = *v_it;
			this->mesh.set_point(vertex, new_v);
			}*/
		}

		if (this->mesh.has_vertex_status())
			this->mesh.release_vertex_status();/**/
		if (this->mesh.has_face_status())
			this->mesh.release_face_status();
		if (this->mesh.has_edge_status())
			this->mesh.release_edge_status();
		if (this->mesh.has_halfedge_status())
			this->mesh.release_halfedge_status();
	}
	if (this->use_sim_or_not == true)
	{
		if (!this->sim_mesh.has_vertex_status())
			this->sim_mesh.request_vertex_status();
		if (!this->sim_mesh.has_face_status())
			this->sim_mesh.request_face_status();
		if (!this->sim_mesh.has_edge_status())
			this->sim_mesh.request_edge_status();
		if (!this->sim_mesh.has_halfedge_status())
			this->sim_mesh.has_halfedge_status();

		for (auto v_it = this->sim_mesh.vertices_begin(); v_it != this->sim_mesh.vertices_end(); ++v_it, ++i)
		{
			//soft con
			MyMesh::Point new_v(unknown(3 * i + 0), unknown(3 * i + 1), unknown(3 * i + 2));
			auto vertex = *v_it;
			this->sim_mesh.set_point(vertex, new_v);

			//hard con
			/*if (idx_move[i] != -1)
			{
			MyMesh::Point new_v(unknown(3 * (i - idx_move[i]) + 0), unknown(3 * (i - idx_move[i]) + 1), unknown(3 * (i - idx_move[i]) + 2));
			auto vertex = *v_it;
			this->mesh.set_point(vertex, new_v);
			}*/
		}

		if (this->sim_mesh.has_vertex_status())
			this->sim_mesh.release_vertex_status();/**/
		if (this->sim_mesh.has_face_status())
			this->sim_mesh.release_face_status();
		if (this->sim_mesh.has_edge_status())
			this->sim_mesh.release_edge_status();
		if (this->sim_mesh.has_halfedge_status())
			this->sim_mesh.release_halfedge_status();
	}

	//SaveDeformedMesh();
}

void mns::MyDeformationHandler::ModOriginalMesh(MyMesh& mesh)
{
	if (!mesh.has_vertex_status())
		mesh.request_vertex_status();
	if (!mesh.has_face_status())
		mesh.request_face_status();
	if (!mesh.has_edge_status())
		mesh.request_edge_status();
	if (!mesh.has_halfedge_status())
		mesh.has_halfedge_status();

	int i = 0;
	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it, ++i)
	{
		//soft con
		MyMesh::Point new_v(unknown(3 * i + 0), unknown(3 * i + 1), unknown(3 * i + 2));
		auto vertex = *v_it;
		mesh.set_point(vertex, new_v);

		//hard con
		/*if (idx_move[i] != -1)
		{
		MyMesh::Point new_v(unknown(3 * (i - idx_move[i]) + 0), unknown(3 * (i - idx_move[i]) + 1), unknown(3 * (i - idx_move[i]) + 2));
		auto vertex = *v_it;
		mesh.set_point(vertex, new_v);
		}*/
	}

	if (mesh.has_vertex_status())
		mesh.release_vertex_status();/**/
	if (mesh.has_face_status())
		mesh.release_face_status();
	if (mesh.has_edge_status())
		mesh.release_edge_status();
	if (mesh.has_halfedge_status())
		mesh.release_halfedge_status();

	bool res = OpenMesh::IO::write_mesh(mesh, "ori_Deformed.ply");
}

bool mns::MyDeformationHandler::SaveDeformedMesh()
{
	bool res = OpenMesh::IO::write_mesh(this->mesh, "Deformed.ply");
	return res;
}

void mns::MyDeformationHandler::Subdivision()
{
	mesh.add_property(vhandle, "The new handle on each vertex of old mesh");

	for (MyMesh::VertexIter vit = mesh.vertices_begin(); vit != mesh.vertices_end(); ++vit)
	{
		MyMesh::Point p = mesh.point(*vit);
		mesh.property(vhandle, *vit) = subdivMesh.add_vertex(p);
	}

	mesh.add_property(ehandle, "the new vertex handle on each edge");
	for (MyMesh::EdgeIter eit = mesh.edges_begin(); eit != mesh.edges_end(); ++eit)
	{
		MyMesh::HalfedgeHandle heit = mesh.halfedge_handle(*eit, 0);
		MyMesh::VertexHandle vit0 = mesh.from_vertex_handle(heit);
		MyMesh::VertexHandle vit1 = mesh.to_vertex_handle(heit);
		MyMesh::Point p0 = mesh.point(vit0);
		MyMesh::Point p1 = mesh.point(vit1);
		MyMesh::Point p((p0 + p1) / 2.0);
		mesh.property(ehandle, *eit) = subdivMesh.add_vertex(p);
	}
	for (MyMesh::FaceIter fit = mesh.faces_begin(); fit != mesh.faces_end(); ++fit)
	{
		MyMesh::FaceHalfedgeIter fhit = mesh.fh_begin(*fit);
		MyMesh::HalfedgeHandle h0 = *fhit;
		++fhit;
		MyMesh::HalfedgeHandle h1 = *fhit;
		++fhit;
		MyMesh::HalfedgeHandle h2 = *fhit;
		MyMesh::VertexHandle v0 = mesh.from_vertex_handle(h0);
		MyMesh::VertexHandle v1 = mesh.to_vertex_handle(h0);
		MyMesh::VertexHandle v2 = mesh.to_vertex_handle(h1);
		MyMesh::VertexHandle p0 = mesh.property(vhandle, v0);
		MyMesh::VertexHandle p1 = mesh.property(vhandle, v1);
		MyMesh::VertexHandle p2 = mesh.property(vhandle, v2);
		MyMesh::VertexHandle p3 = mesh.property(ehandle, mesh.edge_handle(h0));
		MyMesh::VertexHandle p4 = mesh.property(ehandle, mesh.edge_handle(h1));
		MyMesh::VertexHandle p5 = mesh.property(ehandle, mesh.edge_handle(h2));
		std::vector< MyMesh::VertexHandle > face_vHandle(3);
		face_vHandle.clear();
		face_vHandle.push_back(p0);
		face_vHandle.push_back(p3);
		face_vHandle.push_back(p5);
		subdivMesh.add_face(face_vHandle);

		face_vHandle.clear();
		face_vHandle.push_back(p3);
		face_vHandle.push_back(p1);
		face_vHandle.push_back(p4);
		subdivMesh.add_face(face_vHandle);

		face_vHandle.clear();
		face_vHandle.push_back(p4);
		face_vHandle.push_back(p2);
		face_vHandle.push_back(p5);
		subdivMesh.add_face(face_vHandle);

		face_vHandle.clear();
		face_vHandle.push_back(p3);
		face_vHandle.push_back(p4);
		face_vHandle.push_back(p5);
		subdivMesh.add_face(face_vHandle);
	}
	//OpenMesh::IO::write_mesh(subdivMesh, "output_subdiv.ply");
}

void mns::MyDeformationHandler::LoadData()
{
	//传细分网格面的坐标参数
	std::ofstream facestxt("faces.txt", std::ios::out);
	for (MyMesh::FaceIter fit = subdivMesh.faces_begin(); fit != subdivMesh.faces_end(); ++fit)
	{

		MyMesh::FVIter fvit = subdivMesh.fv_begin(*fit);
		facestxt << (*fvit).idx() << " ";
		++fvit;
		facestxt << (*fvit).idx() << " ";
		++fvit;
		facestxt << (*fvit).idx() << std::endl;
	}
	facestxt.close();

	std::ofstream verttxt("vert.txt", std::ios::out);
	for (int i = 0; i < subdivMesh.n_vertices(); i++)
	{
		MyMesh::Point p = subdivMesh.point(subdivMesh.vertex_handle(i));
		for (int j = 0; j < 3; j++)
		{
			verttxt << p[j] << " ";
		}
		verttxt << std::endl;
	}
	verttxt.close();

	//传原来网格的边在细分网格中点的index
	std::ofstream indextxt("index.txt", std::ios::out);
	for (MyMesh::EdgeIter eit = mesh.edges_begin(); eit != mesh.edges_end(); ++eit)
	{
		MyMesh::VertexHandle v = mesh.property(ehandle, *eit);
		indextxt << v.idx() << std::endl;
	}
	indextxt.close();
}

void mns::MyDeformationHandler::SetRBFlag(bool flag)
{
	this->rb_flag = flag;
}

// first dihedral then length
void mns::MyDeformationHandler::ReceiveBasis() // first dihedral then length
{
	std::ifstream comptxt("C.txt", std::ios::in);
	this->basis = Eigen::MatrixXd::Zero(2 * this->mesh.n_edges(), basis_col_num);
	for (integer j = 0; j < basis_col_num; j++)
	{
		for (integer i = 0; i < this->mesh.n_edges(); i++)
		{
			decimal tmp;
			comptxt >> tmp;
			basis(2 * i + 0, j) = tmp;
			comptxt >> tmp;
			basis(2 * i + 1, j) = tmp;
		}
	}
	comptxt.close();
}

void mns::MyDeformationHandler::ReceiveBasis2(string p_path)
{
	std::ifstream comptxt(p_path, std::ios::in);
	this->basis = Eigen::MatrixXd::Zero(2 * this->mesh.n_edges(), basis_col_num);
	for (integer j = 0; j < basis_col_num; j++)
	{
		for (integer i = 0; i < this->mesh.n_edges(); i++)
		{
			decimal tmp;
			comptxt >> tmp;
			basis(2 * i + 1, j) = tmp;
			comptxt >> tmp;
			basis(2 * i + 0, j) = tmp;
		}
	}
	comptxt.close();

	//basis modification
	Eigen::MatrixXd verts = Eigen::MatrixXd::Zero(2 * this->mesh.n_edges(), eg_num);
	std::ifstream vertstxt("verts.txt", std::ios::in);
	for (integer j = 0; j < eg_num; j++)
	{
		for (integer i = 0; i < this->mesh.n_edges(); i++)
		{
			decimal tmp;
			vertstxt >> tmp;
			verts(2 * i + 1, j) = tmp;
			vertstxt >> tmp;
			verts(2 * i + 0, j) = tmp;
		}
	}
	vertstxt.close();



	//sparate basis into two parts
	decimal small_limit = 0.0000001;
#pragma omp parallel for
	for (integer k = 0; k < basis.cols(); k++)
	{
#pragma omp critical
		{
			dec_vector length_vec;
			dec_vector dihedral_vec;
			Eigen::VectorXd l_tmp = Eigen::VectorXd::Zero(this->mesh.n_edges());
			Eigen::VectorXd d_tmp = Eigen::VectorXd::Zero(this->mesh.n_edges());
			for (integer j = 0; j < basis.rows(); j++)
			{
				//basis(j, k) = basis(j, k) / row_sum(k);
				if (j % 2 == 0)
				{
					l_tmp(j / 2) = basis(j, k);
				}
				else
				{
					d_tmp(j / 2) = basis(j, k);
				}
			}

			l_tmp = l_tmp.normalized();
			d_tmp = d_tmp.normalized();
			for (integer j = 0; j < basis.rows(); j++)
			{
				//basis(j, k) = basis(j, k) / row_sum(k);
				if (j % 2 == 0)
				{
					//length_vec.push_back(basis(j, k));
					length_vec.push_back(l_tmp(j / 2));
				}
				else
				{
					//dihedral_vec.push_back(basis(j, k));
					dihedral_vec.push_back(d_tmp(j / 2));
				}
			}
			this->length_basis_vec.push_back(length_vec);
			this->dihedral_basis_vec.push_back(dihedral_vec);
		}
	}

	Eigen::MatrixXd tmp_basis = Eigen::MatrixXd::Zero(2 * this->mesh.n_edges() + 1, basis_col_num);
	if (alpha_v != 0)
	{
		tmp_basis.block(0, 0, basis.rows(), basis.cols()) = basis;
		decimal tmp = -0.1;
#pragma omp parallel for
		for (integer i = 0; i < basis.cols(); i++)
		{
#pragma omp critical 
			{
				tmp_basis(2 * this->mesh.n_edges(), i) = -tmp;
				tmp = -1 * tmp;
			}
		}
		this->basis.resizeLike(tmp_basis);
		this->basis = tmp_basis;
	}
}


void mns::MyDeformationHandler::ReceiveBasis2()
{
	std::ifstream comptxt("F:\\scan\\male_expressions\\DE\\C.txt", std::ios::in);
	this->basis = Eigen::MatrixXd::Zero(2 * this->mesh.n_edges(), basis_col_num);
	for (integer j = 0; j < basis_col_num; j++)
	{
		for (integer i = 0; i < this->mesh.n_edges(); i++)
		{
			decimal tmp;
			comptxt >> tmp;
			basis(2 * i + 1, j) = tmp;
			comptxt >> tmp;
			basis(2 * i + 0, j) = tmp;
		}
	}
	comptxt.close();

	//basis modification
	Eigen::MatrixXd verts = Eigen::MatrixXd::Zero(2 * this->mesh.n_edges(), eg_num);
	std::ifstream vertstxt("verts.txt", std::ios::in);
	for (integer j = 0; j < eg_num; j++)
	{
		for (integer i = 0; i < this->mesh.n_edges(); i++)
		{
			decimal tmp;
			vertstxt >> tmp;
			verts(2 * i + 1, j) = tmp;
			vertstxt >> tmp;
			verts(2 * i + 0, j) = tmp;
		}
	}
	vertstxt.close();

	

	//sparate basis into two parts
	decimal small_limit = 0.0000001;
#pragma omp parallel for
	for (integer k = 0; k < basis.cols(); k++)
	{
#pragma omp critical
		{
			dec_vector length_vec;
			dec_vector dihedral_vec;
			Eigen::VectorXd l_tmp = Eigen::VectorXd::Zero(this->mesh.n_edges());
			Eigen::VectorXd d_tmp = Eigen::VectorXd::Zero(this->mesh.n_edges());
			for (integer j = 0; j < basis.rows(); j++)
			{
				//basis(j, k) = basis(j, k) / row_sum(k);
				if (j % 2 == 0)
				{
					l_tmp(j / 2) = basis(j, k);
				}
				else
				{
					d_tmp(j / 2) = basis(j, k);
				}
			}

			l_tmp = l_tmp.normalized();
			d_tmp = d_tmp.normalized();
			for (integer j = 0; j < basis.rows(); j++)
			{
				//basis(j, k) = basis(j, k) / row_sum(k);
				if (j % 2 == 0)
				{
					//length_vec.push_back(basis(j, k));
					length_vec.push_back(l_tmp(j / 2));
				}
				else
				{
					//dihedral_vec.push_back(basis(j, k));
					dihedral_vec.push_back(d_tmp(j / 2));
				}
			}
			this->length_basis_vec.push_back(length_vec);
			this->dihedral_basis_vec.push_back(dihedral_vec);
		}
	}

	Eigen::MatrixXd tmp_basis = Eigen::MatrixXd::Zero(2 * this->mesh.n_edges() + 1, basis_col_num);
	if (alpha_v != 0)
	{
		tmp_basis.block(0, 0, basis.rows(), basis.cols()) = basis;
		decimal tmp = -0.1;
#pragma omp parallel for
		for (integer i = 0; i < basis.cols(); i++)
		{
#pragma omp critical 
			{
				tmp_basis(2 * this->mesh.n_edges(), i) = -tmp;
				tmp = -1 * tmp;
			}
		}
		this->basis.resizeLike(tmp_basis);
		this->basis = tmp_basis;
	}
}

void mns::MyDeformationHandler::ReceiveBasis3()
{
	std::ifstream lcomptxt("lC.txt", std::ios::in);
	std::ifstream dcomptxt("dC.txt", std::ios::in);
	this->basis = Eigen::MatrixXd::Zero(2 * this->mesh.n_edges(), basis_col_num);
	for (integer j = 0; j < basis_col_num; j++)
	{
		for (integer i = 0; i < this->mesh.n_edges(); i++)
		{
			decimal tmp;
			lcomptxt >> tmp;
			lcomptxt >> tmp;
			basis(2 * i + 0, j) = tmp;
			dcomptxt >> tmp;
			basis(2 * i + 1, j) = tmp;
			dcomptxt >> tmp;
		}
	}
	lcomptxt.close();
	dcomptxt.close();
	/*
	std::ofstream ba_txt("output_basis.txt");
	ba_txt << basis;
	ba_txt.close();
	*/
	//sparate basis into two parts
	decimal small_limit = 0.0000001;
#pragma omp parallel for
	for (integer k = 0; k < basis.cols(); k++)
	{
#pragma omp critical
		{
			dec_vector length_vec;
			dec_vector dihedral_vec;
			Eigen::VectorXd l_tmp = Eigen::VectorXd::Zero(this->mesh.n_edges());
			Eigen::VectorXd d_tmp = Eigen::VectorXd::Zero(this->mesh.n_edges());
			for (integer j = 0; j < basis.rows(); j++)
			{
				//basis(j, k) = basis(j, k) / row_sum(k);
				if (j % 2 == 0)
				{
					l_tmp(j / 2) = basis(j, k);
				}
				else
				{
					d_tmp(j / 2) = basis(j, k);
				}
			}

			l_tmp = l_tmp.normalized();
			d_tmp = d_tmp.normalized();
			for (integer j = 0; j < basis.rows(); j++)
			{
				//basis(j, k) = basis(j, k) / row_sum(k);
				if (j % 2 == 0)
				{
					//length_vec.push_back(basis(j, k));
					length_vec.push_back(l_tmp(j / 2));
				}
				else
				{
					//dihedral_vec.push_back(basis(j, k));
					dihedral_vec.push_back(d_tmp(j / 2));
				}
			}
			this->length_basis_vec.push_back(length_vec);
			this->dihedral_basis_vec.push_back(dihedral_vec);
		}
	}

	Eigen::MatrixXd tmp_basis = Eigen::MatrixXd::Zero(2 * this->mesh.n_edges() + 1, basis_col_num);
	if (alpha_v != 0)
	{
		tmp_basis.block(0, 0, basis.rows(), basis.cols()) = basis;
		decimal tmp = -0.1;
#pragma omp parallel for
		for (integer i = 0; i < basis.cols(); i++)
		{
#pragma omp critical 
			{
				tmp_basis(2 * this->mesh.n_edges(), i) = -tmp;
				tmp = -1 * tmp;
			}
		}
		this->basis.resizeLike(tmp_basis);
		this->basis = tmp_basis;
	}
}

void mns::MyDeformationHandler::ReceiveBasis4()
{
	std::ifstream lcomptxt("lC.txt", std::ios::in);
	std::ifstream dcomptxt("dC.txt", std::ios::in);
	this->basis = Eigen::MatrixXd::Zero(2 * this->mesh.n_edges(), basis_col_num);
	for (integer j = 0; j < basis_col_num; j++)
	{
		for (integer i = 0; i < this->mesh.n_edges(); i++)
		{
			decimal tmp;
			lcomptxt >> tmp;
			basis(2 * i + 0, j) = tmp;
			dcomptxt >> tmp;
			basis(2 * i + 1, j) = tmp;
		}
	}
	lcomptxt.close();
	dcomptxt.close();
	/*
	std::ofstream ba_txt("output_basis.txt");
	ba_txt << basis;
	ba_txt.close();
	*/
	//sparate basis into two parts
	decimal small_limit = 0.0000001;
#pragma omp parallel for
	for (integer k = 0; k < basis.cols(); k++)
	{
#pragma omp critical
		{
			dec_vector length_vec;
			dec_vector dihedral_vec;
			Eigen::VectorXd l_tmp = Eigen::VectorXd::Zero(this->mesh.n_edges());
			Eigen::VectorXd d_tmp = Eigen::VectorXd::Zero(this->mesh.n_edges());
			for (integer j = 0; j < basis.rows(); j++)
			{
				//basis(j, k) = basis(j, k) / row_sum(k);
				if (j % 2 == 0)
				{
					l_tmp(j / 2) = basis(j, k);
				}
				else
				{
					d_tmp(j / 2) = basis(j, k);
				}
			}

			l_tmp = l_tmp.normalized();
			d_tmp = d_tmp.normalized();
			for (integer j = 0; j < basis.rows(); j++)
			{
				//basis(j, k) = basis(j, k) / row_sum(k);
				if (j % 2 == 0)
				{
					//length_vec.push_back(basis(j, k));
					length_vec.push_back(l_tmp(j / 2));
				}
				else
				{
					//dihedral_vec.push_back(basis(j, k));
					dihedral_vec.push_back(d_tmp(j / 2));
				}
			}
			this->length_basis_vec.push_back(length_vec);
			this->dihedral_basis_vec.push_back(dihedral_vec);
		}
	}

	Eigen::MatrixXd tmp_basis = Eigen::MatrixXd::Zero(2 * this->mesh.n_edges() + 1, basis_col_num);
	if (alpha_v != 0)
	{
		tmp_basis.block(0, 0, basis.rows(), basis.cols()) = basis;
		decimal tmp = -0.1;
#pragma omp parallel for
		for (integer i = 0; i < basis.cols(); i++)
		{
#pragma omp critical 
			{
				tmp_basis(2 * this->mesh.n_edges(), i) = -tmp;
				tmp = -1 * tmp;
			}
		}
		this->basis.resizeLike(tmp_basis);
		this->basis = tmp_basis;
	}
}

// first dihedral then length
void mns::MyDeformationHandler::ReceiveBasisWithName(std::string nm, integer NB)
{
	SetRBFlag(1);
	std::ifstream comptxt(nm, std::ios::in);
	this->basis = Eigen::MatrixXd::Zero(2 * this->mesh.n_edges(), NB);
	for (integer j = 0; j < basis_col_num; j++)
	{
		for (integer i = 0; i < this->mesh.n_edges(); i++)
		{
			decimal tmp;
			comptxt >> tmp;
			basis(2 * i + 0, j) = tmp;
			comptxt >> tmp;
			basis(2 * i + 1, j) = tmp;
		}
	}
	comptxt.close();
}

void mns::MyDeformationHandler::SeparateVerts(const Eigen::MatrixXd& verts)
{
	Eigen::MatrixXd l_verts = Eigen::MatrixXd::Zero(2 * this->mesh.n_edges(), eg_num);
	Eigen::MatrixXd d_verts = Eigen::MatrixXd::Zero(2 * this->mesh.n_edges(), eg_num);

#pragma omp parallel for
	for (integer j = 0; j < eg_num; j++)
	{
#pragma omp parallel for
		for (integer i = 0; i < this->mesh.n_edges(); i++)
		{
#pragma omp critical 
			{
				l_verts(2 * i + 0, j) = verts(2 * i + 0, j);
				d_verts(2 * i + 1, j) = verts(2 * i + 1, j);
			}
		}
	}

	std::ofstream lvertstxt("lverts.txt", std::ios::out);
	std::ofstream dvertstxt("dverts.txt", std::ios::out);
	auto lvt = l_verts.transpose();
	auto dvt = d_verts.transpose();
	for (integer j = 0; j < eg_num; j++)
	{
		for (integer i = 0; i < 2 * this->mesh.n_edges(); i++)
		{
			lvertstxt << lvt(j, i) << " ";
			dvertstxt << dvt(j, i) << " ";
		}
		lvertstxt << std::endl;
		dvertstxt << std::endl;
	}
	lvertstxt.close();
	dvertstxt.close();
}

integer mns::MyDeformationHandler::GetBasisCols()
{
	return this->basis.cols();
}

void mns::MyDeformationHandler::CreateDiLenMat2()
{
	if (use_basis_or_not == false)
	{
		return;
	}

	dec_vector length_set;
	dec_vector dihedral_set;

	for (integer j = 0; j < this->mesh.n_edges(); j++)
	{
		decimal tmp_length_diff = 0;
		decimal tmp_dihedral_diff = 0;
		for (integer k = 0; k < basis.cols(); k++)
		{
			tmp_length_diff += l_bw_unknown(k) * length_basis_vec[k][j];
			tmp_dihedral_diff += d_bw_unknown(k) * dihedral_basis_vec[k][j];
		}
		length_set.push_back(ori_length_set[j] + tmp_length_diff);
		dihedral_set.push_back(ori_dihedral_set[j] + tmp_dihedral_diff);
	}

	std::ofstream ofile("DiLenVec.txt", std::ios::out);
	for (integer j = 0; j < this->mesh.n_edges(); j++)
	{
		ofile << sign_flag[j] * (3.14159 - dihedral_set[j]) << std::endl;
		ofile << length_set[j] << std::endl;
	}
	ofile.close();
}

void mns::MyDeformationHandler::CreateDiLenMat()
{
	if (use_basis_or_not == false)
	{
		return;
	}

	dec_vector length_set;
	dec_vector dihedral_set;

	for (integer j = 0; j < this->mesh.n_edges(); j++)
	{
		decimal tmp_length_diff = 0;
		decimal tmp_dihedral_diff = 0;
		for (integer k = 0; k < basis.cols(); k++)
		{
			tmp_length_diff += unknown(3 * this->mesh.n_vertices() + 2 * k + 0) * length_basis_vec[k][j];
			tmp_dihedral_diff += unknown(3 * this->mesh.n_vertices() + 2 * k + 1) * dihedral_basis_vec[k][j];
		}
		length_set.push_back(ori_length_set[j] + tmp_length_diff);
		dihedral_set.push_back(ori_dihedral_set[j] + tmp_dihedral_diff);
	}

	std::ofstream ofile("DiLenVec.txt", std::ios::out);
	for (integer j = 0; j < this->mesh.n_edges(); j++)
	{
		ofile << sign_flag[j] * (3.14159 - dihedral_set[j]) << std::endl;
		ofile << length_set[j] << std::endl;
	}
	ofile.close();
}

bool mns::MyDeformationHandler::GetRBF()
{
	return this->rb_flag;
}

void mns::MyDeformationHandler::CalWeightByBasis()
{
	MyMesh mesh;
	bool res = OpenMesh::IO::read_mesh(mesh, "high.obj");
	if (res == false)
	{
		std::cout << "read error\n";
	}
	dec_vector length_set;
	dec_vector dihedral_set;
	for (auto he_it = mesh.halfedges_begin(); he_it != mesh.halfedges_end(); ++he_it)
	{
		auto e = *he_it;
		auto fromVertex = mesh.from_vertex_handle(e);
		auto toVertex = mesh.to_vertex_handle(e);
		auto fromPt = mesh.point(fromVertex);
		auto toPt = mesh.point(toVertex);

		//stretching initialization		
		decimal l = mesh.calc_edge_length(e);
		length_set.push_back(l);

		//bend initialization
		auto f1 = mesh.face_handle(e);
		auto op_he = mesh.opposite_halfedge_handle(e);
		auto f2 = mesh.face_handle(op_he);
		pt_vector pt_set;
		pt_vector t1_set;
		pt_vector t2_set;

		if (f1.idx() == -1 || f2.idx() == -1)
		{
			dihedral_set.push_back(0);

			//edge, not half edge
			++he_it;
			continue;
		}

		pt_set.push_back(fromPt);
		t1_set.push_back(fromPt);
		t2_set.push_back(fromPt);
		pt_set.push_back(toPt);
		t1_set.push_back(toPt);
		t2_set.push_back(toPt);

		for (auto fv_it1 = mesh.fv_begin(f1); fv_it1 != mesh.fv_end(f1); ++fv_it1)
		{
			auto v = *fv_it1;
			auto pt = mesh.point(v);
			if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
			{
				pt_set.push_back(pt);
				t1_set.push_back(pt);
				break;
			}
		}
		integer temp_cnt = 0;
		for (auto fv_it2 = mesh.fv_begin(f2); fv_it2 != mesh.fv_end(f2); ++fv_it2)
		{
			auto v = *fv_it2;
			auto pt = mesh.point(v);
			if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
			{
				t2_set.push_back(pt);
				pt_set.push_back(pt);
				break;
			}
		}

		decimal angle = CalDiheral(pt_set[0], pt_set[1], pt_set[2], pt_set[3]);
		dihedral_set.push_back(angle);

		//edge, not half edge
		++he_it;
	}

	Eigen::VectorXd b_dih = Eigen::VectorXd::Zero(dihedral_set.size());
	Eigen::MatrixXd A_dih = Eigen::MatrixXd::Zero(dihedral_set.size(), dihedral_basis_vec.size());

	for (integer k = 0; k <= dihedral_basis_vec.size(); k++)
	{
		for (integer i = 0; i < dihedral_set.size(); i++)
		{
			if (k < dihedral_basis_vec.size())
			{
				A_dih(i, k) = dihedral_basis_vec[k][i];
			}
			else
			{
				b_dih(i) = dihedral_set[i];
			}
		}
	}

	std::ofstream x_dih_ofile("cal_w_dih.txt", std::ios::out);
	Eigen::VectorXd x_dih = A_dih.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b_dih);
	x_dih_ofile << x_dih;
	x_dih_ofile.close();

	Eigen::VectorXd b_len = Eigen::VectorXd::Zero(length_set.size());
	Eigen::MatrixXd A_len = Eigen::MatrixXd::Zero(length_set.size(), length_basis_vec.size());

	for (integer k = 0; k <= length_basis_vec.size(); k++)
	{
		for (integer i = 0; i < length_set.size(); i++)
		{
			if (k < length_basis_vec.size())
			{
				A_len(i, k) = length_basis_vec[k][i];
			}
			else
			{
				b_len(i) = length_set[i];
			}
		}
	}

	std::ofstream x_len_ofile("cal_w_len.txt", std::ios::out);
	Eigen::VectorXd x_len = A_len.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b_len);
	x_len_ofile << x_len;
	x_len_ofile.close();
}

void mns::MyDeformationHandler::CalWunknownByBasis()
{
	integer nE = this->mesh.n_edges();
	integer nV = this->mesh.n_vertices();
	integer nB = this->basis_col_num;

	dec_vector length_set;
	dec_vector dihedral_set;

	for (auto he_it = mesh.halfedges_begin(); he_it != mesh.halfedges_end(); ++he_it)
	{
		auto e = *he_it;
		auto fromVertex = mesh.from_vertex_handle(e);
		auto toVertex = mesh.to_vertex_handle(e);
		integer fromIdx = fromVertex.idx();
		auto fromPt = MyMesh::Point(unknown(3 * fromIdx + 0), unknown(3 * fromIdx + 1), unknown(3 * fromIdx + 2));
		integer toIdx = toVertex.idx();
		auto toPt = MyMesh::Point(unknown(3 * toIdx + 0), unknown(3 * toIdx + 1), unknown(3 * toIdx + 2));
		/*auto fromPt = mesh.point(fromVertex);
		auto toPt = mesh.point(toVertex);*/

		//stretching initialization		
		decimal l = (toPt - fromPt).norm();
		length_set.push_back(l);

		//bend initialization
		auto f1 = mesh.face_handle(e);
		auto op_he = mesh.opposite_halfedge_handle(e);
		auto f2 = mesh.face_handle(op_he);
		pt_vector pt_set;
		pt_vector t1_set;
		pt_vector t2_set;

		if (f1.idx() == -1 || f2.idx() == -1)
		{
			dihedral_set.push_back(0);

			//edge, not half edge
			++he_it;
			continue;
		}

		pt_set.push_back(fromPt);
		t1_set.push_back(fromPt);
		t2_set.push_back(fromPt);
		pt_set.push_back(toPt);
		t1_set.push_back(toPt);
		t2_set.push_back(toPt);

		for (auto fv_it1 = mesh.fv_begin(f1); fv_it1 != mesh.fv_end(f1); ++fv_it1)
		{
			auto v = *fv_it1;
			integer vIdx = v.idx();
			auto pt = MyMesh::Point(unknown(3 * vIdx + 0), unknown(3 * vIdx + 1), unknown(3 * vIdx + 2));
			if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
			{
				pt_set.push_back(pt);
				t1_set.push_back(pt);
				break;
			}
		}
		integer temp_cnt = 0;
		for (auto fv_it2 = mesh.fv_begin(f2); fv_it2 != mesh.fv_end(f2); ++fv_it2)
		{
			auto v = *fv_it2;
			integer vIdx = v.idx();
			auto pt = MyMesh::Point(unknown(3 * vIdx + 0), unknown(3 * vIdx + 1), unknown(3 * vIdx + 2));
			if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
			{
				t2_set.push_back(pt);
				pt_set.push_back(pt);
				break;
			}
		}

		decimal angle = CalDiheral(pt_set[0], pt_set[1], pt_set[2], pt_set[3]);
		dihedral_set.push_back(angle);

		//edge, not half edge
		++he_it;
	}

	Eigen::VectorXd b_dih = Eigen::VectorXd::Zero(dihedral_set.size());
	Eigen::MatrixXd A_dih = Eigen::MatrixXd::Zero(dihedral_set.size(), dihedral_basis_vec.size());

#pragma omp parallel for 
	for (integer k = 0; k <= dihedral_basis_vec.size(); k++)
	{
#pragma omp critical
		{
			for (integer i = 0; i < dihedral_set.size(); i++)
			{
				if (k < dihedral_basis_vec.size())
				{
					A_dih(i, k) = dihedral_basis_vec[k][i];
				}
				else
				{
					b_dih(i) = dihedral_set[i] - this->data_avg[2 * i + 1];
				}
			}
		}
	}

	Eigen::VectorXd b_len = Eigen::VectorXd::Zero(length_set.size());
	Eigen::MatrixXd A_len = Eigen::MatrixXd::Zero(length_set.size(), length_basis_vec.size());

#pragma omp parallel for
	for (integer k = 0; k <= length_basis_vec.size(); k++)
	{
#pragma omp critical
		{
			for (integer i = 0; i < length_set.size(); i++)
			{
				if (k < length_basis_vec.size())
				{
					A_len(i, k) = length_basis_vec[k][i];
				}
				else
				{
					b_len(i) = length_set[i] - this->data_avg[2 * i + 0];
				}
			}
		}
	}

	Eigen::VectorXd x_dih = A_dih.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b_dih);
	Eigen::VectorXd x_len = A_len.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b_len);

#pragma omp parallel for
	for (integer i = 0; i < nB; i++)
	{
#pragma omp critical
		{

			unknown(3 * nV + 2 * i + 0) = x_len(i);
			unknown(3 * nV + 2 * i + 1) = x_dih(i);
		}
	}

	/*std::ofstream x_dih_ofile("cal_wu_dih.txt", std::ios::out);
	x_dih_ofile << x_dih;
	x_dih_ofile.close();

	std::ofstream x_len_ofile("cal_wu_len.txt", std::ios::out);
	x_len_ofile << x_len;
	x_len_ofile.close();*/
}

void mns::MyDeformationHandler::CalWunknownByBasis2()
{
	integer nE = this->mesh.n_edges();
	integer nV = this->mesh.n_vertices();
	integer nB = this->basis_col_num;

	dec_vector length_set;
	dec_vector dihedral_set;

	for (auto he_it = mesh.halfedges_begin(); he_it != mesh.halfedges_end(); ++he_it)
	{
		auto e = *he_it;
		auto fromVertex = mesh.from_vertex_handle(e);
		auto toVertex = mesh.to_vertex_handle(e);
		integer fromIdx = fromVertex.idx();
		auto fromPt = MyMesh::Point(unknown(3 * fromIdx + 0), unknown(3 * fromIdx + 1), unknown(3 * fromIdx + 2));
		integer toIdx = toVertex.idx();
		auto toPt = MyMesh::Point(unknown(3 * toIdx + 0), unknown(3 * toIdx + 1), unknown(3 * toIdx + 2));
		/*auto fromPt = mesh.point(fromVertex);
		auto toPt = mesh.point(toVertex);*/

		//stretching initialization		
		decimal l = (toPt - fromPt).norm();
		length_set.push_back(l);

		//bend initialization
		auto f1 = mesh.face_handle(e);
		auto op_he = mesh.opposite_halfedge_handle(e);
		auto f2 = mesh.face_handle(op_he);
		pt_vector pt_set;
		pt_vector t1_set;
		pt_vector t2_set;

		if (f1.idx() == -1 || f2.idx() == -1)
		{
			dihedral_set.push_back(0);

			//edge, not half edge
			++he_it;
			continue;
		}

		pt_set.push_back(fromPt);
		t1_set.push_back(fromPt);
		t2_set.push_back(fromPt);
		pt_set.push_back(toPt);
		t1_set.push_back(toPt);
		t2_set.push_back(toPt);

		for (auto fv_it1 = mesh.fv_begin(f1); fv_it1 != mesh.fv_end(f1); ++fv_it1)
		{
			auto v = *fv_it1;
			integer vIdx = v.idx();
			auto pt = MyMesh::Point(unknown(3 * vIdx + 0), unknown(3 * vIdx + 1), unknown(3 * vIdx + 2));
			if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
			{
				pt_set.push_back(pt);
				t1_set.push_back(pt);
				break;
			}
		}
		integer temp_cnt = 0;
		for (auto fv_it2 = mesh.fv_begin(f2); fv_it2 != mesh.fv_end(f2); ++fv_it2)
		{
			auto v = *fv_it2;
			integer vIdx = v.idx();
			auto pt = MyMesh::Point(unknown(3 * vIdx + 0), unknown(3 * vIdx + 1), unknown(3 * vIdx + 2));
			if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
			{
				t2_set.push_back(pt);
				pt_set.push_back(pt);
				break;
			}
		}

		decimal angle = CalDiheral(pt_set[0], pt_set[1], pt_set[2], pt_set[3]);
		dihedral_set.push_back(angle);

		//edge, not half edge
		++he_it;
	}

	Eigen::VectorXd b_dih = Eigen::VectorXd::Zero(dihedral_set.size());
	Eigen::MatrixXd A_dih = Eigen::MatrixXd::Zero(dihedral_set.size(), dihedral_basis_vec.size());

#pragma omp parallel for 
	for (integer k = 0; k <= dihedral_basis_vec.size(); k++)
	{
#pragma omp critical
		{
			for (integer i = 0; i < dihedral_set.size(); i++)
			{
				if (k < dihedral_basis_vec.size())
				{
					A_dih(i, k) = dihedral_basis_vec[k][i];
				}
				else
				{
					b_dih(i) = dihedral_set[i] - this->data_avg[2 * i + 1];
				}
			}
		}
	}

	Eigen::VectorXd b_len = Eigen::VectorXd::Zero(length_set.size());
	Eigen::MatrixXd A_len = Eigen::MatrixXd::Zero(length_set.size(), length_basis_vec.size());

#pragma omp parallel for
	for (integer k = 0; k <= length_basis_vec.size(); k++)
	{
#pragma omp critical
		{
			for (integer i = 0; i < length_set.size(); i++)
			{
				if (k < length_basis_vec.size())
				{
					A_len(i, k) = length_basis_vec[k][i];
				}
				else
				{
					b_len(i) = length_set[i] - this->data_avg[2 * i + 0];
				}
			}
		}
	}

	d_bw_unknown = A_dih.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b_dih);
	l_bw_unknown = A_len.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b_len);

}

void mns::MyDeformationHandler::Clear4ReCal()
{
	con_ele_cnt = 0;
	str_ele_cnt = 0;
	ben_ele_cnt = 0;
	sb_ele_cnt = 0;
	ww_ele_cnt = 0;
	row_cnt = 0;

	triple.clear();
	j_triple.clear();
}

bool mns::MyDeformationHandler::Is_Alpha_Change()
{
	return (alpha_c_ch || alpha_v_ch || alpha_s_ch || alpha_b_ch || alpha_w_ch);
}

void mns::MyDeformationHandler::Alpha_C_Change(decimal fac)
{
	this->swc = this->swc * fac;
	this->alpha_c_ch = true;
}

void mns::MyDeformationHandler::Alpha_V_Change(decimal fac)
{
	this->swv = this->swv * fac;
	this->alpha_v_ch = true;
}

void mns::MyDeformationHandler::Alpha_S_Change(decimal fac)
{
	for (integer i = 0; i < sws_set.size(); i++)
	{
		sws_set[i] = sws_set[i] * fac;
	}
	this->alpha_s_ch = true;
}

void mns::MyDeformationHandler::Alpha_B_Change(decimal fac)
{
	for (integer i = 0; i < swb_set.size(); i++)
	{
		swb_set[i] = swb_set[i] * fac;
	}
	this->alpha_b_ch = true;
}

void mns::MyDeformationHandler::Alpha_W_Change(decimal fac)
{
	this->alpha_w = this->alpha_w * fac;
	this->alpha_w_ch = true;
}

decimal mns::MyDeformationHandler::CalStepLength(decimal pre_en, SMatXd& J, Eigen::VectorXd& f, Eigen::VectorXd& hg, Eigen::VectorXd& un)
{
	SMatXd Jt = J.transpose();
	Eigen::VectorXd xk_plus = this->unknown + ori_step * hg;
	Eigen::VectorXd gk = Jt*f;
	decimal phi_prime_zero = gk.transpose()*hg;
	decimal phi_zero = pre_en;

	decimal tmp_step = 0;
	decimal Econ = 0;
	if (alpha_c != 0)
	{
		Econ = this->CSL_Ec(xk_plus);
	}
	decimal Estr = 0;
	decimal Eben = 0;
	if (alpha_s != 0 || alpha_b != 0)
	{
		std::vector<decimal> esb = CSL_Esb(xk_plus);
		Estr = esb[0];
		Eben = esb[1];
	}
	decimal Ewwc = 0;
	if (alpha_w != 0 && this->use_ww_or_not == true)
	{
		Ewwc = this->CSL_Ew(xk_plus);
	}
	decimal Evol = 0;
	if (alpha_v != 0)
	{
		Evol = CSL_Ev(xk_plus);
	}
	if (isnan(Ewwc))
	{
		unknown = pre_unknown;
		return 0;
	}

	decimal phi_k = Econ + Estr + Eben + Evol + Ewwc;
	decimal denom = 2 * (phi_k - phi_zero - phi_prime_zero*ori_step);
	decimal numer = phi_prime_zero*ori_step*ori_step - 2 * ori_step*(phi_k - phi_zero);
	tmp_step = -numer / denom;

	return tmp_step;
}

decimal mns::MyDeformationHandler::CSL_GetConTerm(integer id, Eigen::VectorXd& un)
{
	decimal energy = 0;
	Eigen::Vector3d f = Eigen::Vector3d::Zero(3);
	if (this->use_sim_or_not == false)
	{
		f(0) = swc*(un(3 * cps.set[id].idx + 0) - cps.set[id].x);
		f(1) = swc*(un(3 * cps.set[id].idx + 1) - cps.set[id].y);
		f(2) = swc*(un(3 * cps.set[id].idx + 2) - cps.set[id].z);
	}

	energy = f(0)*f(0) + f(1)*f(1) + f(2)*f(2);

	return energy;
}

decimal mns::MyDeformationHandler::CSL_Ec(Eigen::VectorXd& un)
{
	decimal total_energy = 0;
	for (int i = 0; i < cps.set.size(); i++)
	{
		total_energy += CSL_GetConTerm(i, un);
	}

	return total_energy;
}

decimal mns::MyDeformationHandler::CSL_GetStrTerm(integer from_id, integer to_id, integer e_id, Eigen::VectorXd& un)
{
	decimal sws = this->sws_set[e_id];
	decimal energy = 0;
	decimal f = 0;

	if (this->use_sim_or_not == false)
	{
		//soft con
		decimal x_dif = un[3 * to_id + 0] - un[3 * from_id + 0];
		decimal y_dif = un[3 * to_id + 1] - un[3 * from_id + 1];
		decimal z_dif = un[3 * to_id + 2] - un[3 * from_id + 2];

		//hard con
		/*MyMesh::Point to_point;
		MyMesh::Point from_point;
		if (idx_move[to_id] == -1)
		{
		to_point[0] = cps.set[con_ori_map[to_id]].x;
		to_point[1] = cps.set[con_ori_map[to_id]].y;
		to_point[2] = cps.set[con_ori_map[to_id]].z;
		}
		else
		{
		to_point[0] = unknown[3 * (to_id - idx_move[to_id]) + 0];
		to_point[1] = unknown[3 * (to_id - idx_move[to_id]) + 1];
		to_point[2] = unknown[3 * (to_id - idx_move[to_id]) + 2];
		}
		if (idx_move[from_id] == -1)
		{
		from_point[0] = cps.set[con_ori_map[from_id]].x;
		from_point[1] = cps.set[con_ori_map[from_id]].y;
		from_point[2] = cps.set[con_ori_map[from_id]].z;
		}
		else
		{
		from_point[0] = unknown[3 * (from_id - idx_move[from_id]) + 0];
		from_point[1] = unknown[3 * (from_id - idx_move[from_id]) + 1];
		from_point[2] = unknown[3 * (from_id - idx_move[from_id]) + 2];
		}
		decimal x_dif = to_point.data()[0] - from_point.data()[0];
		decimal y_dif = to_point.data()[1] - from_point.data()[1];
		decimal z_dif = to_point.data()[2] - from_point.data()[2];*/

		decimal cur_sqr_l = x_dif*x_dif + y_dif*y_dif + z_dif*z_dif;
		decimal cur_l = std::sqrt(cur_sqr_l);

		/*fs(str_ele_cnt) = sws*(cur_l - wnt_length_set[e_id]);
		energy = fs(str_ele_cnt)*fs(str_ele_cnt);
		++str_ele_cnt;*/
		f = sws*(cur_l - wnt_length_set[e_id]);
		energy = f*f;
	}
	//if (this->use_sim_or_not == true)
	//{
	//	/*decimal from_x = 0;
	//	decimal from_y = 0;
	//	decimal from_z = 0;
	//	decimal to_x = 0;
	//	decimal to_y = 0;
	//	decimal to_z = 0;
	//	for (integer j = 0; j < this->sim_mesh.n_vertices(); j++)
	//	{
	//	decimal w_from_j = W_mat(from_id, j);
	//	from_x += w_from_j * unknown(3 * j + 0);
	//	from_y += w_from_j * unknown(3 * j + 1);
	//	from_z += w_from_j * unknown(3 * j + 2);

	//	decimal w_to_j = W_mat(to_id, j);
	//	to_x += w_to_j * unknown(3 * j + 0);
	//	to_y += w_to_j * unknown(3 * j + 1);
	//	to_z += w_to_j * unknown(3 * j + 2);
	//	}*/
	//	//soft con
	//	decimal x_dif = high_v_mat(to_id, 0) - high_v_mat(from_id, 0);
	//	decimal y_dif = high_v_mat(to_id, 1) - high_v_mat(from_id, 1);
	//	decimal z_dif = high_v_mat(to_id, 2) - high_v_mat(from_id, 2);
	//	//hard con(not ready)

	//	decimal cur_sqr_l = x_dif*x_dif + y_dif*y_dif + z_dif*z_dif;
	//	decimal cur_l = std::sqrt(cur_sqr_l);

	//	/*fs(str_ele_cnt) = sws*(cur_l - wnt_length_set[e_id]);
	//	energy = fs(str_ele_cnt)*fs(str_ele_cnt);
	//	++str_ele_cnt;*/
	//	f = sws*(cur_l - wnt_length_set[e_id]);
	//	energy = f*f;
	//}

	return energy;
}

decimal mns::MyDeformationHandler::CSL_GetBenTerm(integer from_id, integer to_id, integer left_id, integer right_id, integer e_id, Eigen::VectorXd& un)
{
	decimal swb = this->swb_set[e_id];
	decimal energy = 0;
	decimal f = 0;

	if (this->use_sim_or_not == false)
	{
		//soft con
		MyMesh::Point rj(un(3 * from_id + 0), un(3 * from_id + 1), un(3 * from_id + 2));
		MyMesh::Point rk(un(3 * to_id + 0), un(3 * to_id + 1), un(3 * to_id + 2));
		MyMesh::Point ri(un(3 * left_id + 0), un(3 * left_id + 1), un(3 * left_id + 2));
		MyMesh::Point rl(un(3 * right_id + 0), un(3 * right_id + 1), un(3 * right_id + 2));

		//hard con
		//MyMesh::Point rj;
		//MyMesh::Point rk;
		//MyMesh::Point ri;
		//MyMesh::Point rl;
		////rj
		//if (idx_move[from_id] == -1)
		//{
		//	rj[0] = cps.set[con_ori_map[from_id]].x;
		//	rj[1] = cps.set[con_ori_map[from_id]].y;
		//	rj[2] = cps.set[con_ori_map[from_id]].z;
		//}
		//else
		//{
		//	rj[0] = unknown[3 * (from_id - idx_move[from_id]) + 0];
		//	rj[1] = unknown[3 * (from_id - idx_move[from_id]) + 1];
		//	rj[2] = unknown[3 * (from_id - idx_move[from_id]) + 2];
		//}
		////rk
		//if (idx_move[to_id] == -1)
		//{
		//	rk[0] = cps.set[con_ori_map[to_id]].x;
		//	rk[1] = cps.set[con_ori_map[to_id]].y;
		//	rk[2] = cps.set[con_ori_map[to_id]].z;
		//}
		//else
		//{
		//	rk[0] = unknown[3 * (to_id - idx_move[to_id]) + 0];
		//	rk[1] = unknown[3 * (to_id - idx_move[to_id]) + 1];
		//	rk[2] = unknown[3 * (to_id - idx_move[to_id]) + 2];
		//}
		////ri
		//if (idx_move[left_id] == -1)
		//{
		//	ri[0] = cps.set[con_ori_map[left_id]].x;
		//	ri[1] = cps.set[con_ori_map[left_id]].y;
		//	ri[2] = cps.set[con_ori_map[left_id]].z;
		//}
		//else
		//{
		//	ri[0] = unknown[3 * (left_id - idx_move[left_id]) + 0];
		//	ri[1] = unknown[3 * (left_id - idx_move[left_id]) + 1];
		//	ri[2] = unknown[3 * (left_id - idx_move[left_id]) + 2];
		//}
		////rl
		//if (idx_move[right_id] == -1)
		//{
		//	rl[0] = cps.set[con_ori_map[right_id]].x;
		//	rl[1] = cps.set[con_ori_map[right_id]].y;
		//	rl[2] = cps.set[con_ori_map[right_id]].z;
		//}
		//else
		//{
		//	rl[0] = unknown[3 * (right_id - idx_move[right_id]) + 0];
		//	rl[1] = unknown[3 * (right_id - idx_move[right_id]) + 1];
		//	rl[2] = unknown[3 * (right_id - idx_move[right_id]) + 2];
		//}
		//

		decimal cur_angle = CalDiheral(rj, rk, ri, rl);
		/*fb(ben_ele_cnt) = swb * (cur_angle - wnt_dihedral_set[e_id]);
		energy = fb(ben_ele_cnt) * fb(ben_ele_cnt);
		++ben_ele_cnt;*/
		f = swb * (cur_angle - wnt_dihedral_set[e_id]);
		//if(_isnan(energy)
		energy = f * f;
	}
	if (this->use_sim_or_not == true)
	{
		/*decimal from_x = 0;
		decimal from_y = 0;
		decimal from_z = 0;
		decimal to_x = 0;
		decimal to_y = 0;
		decimal to_z = 0;
		decimal left_x = 0;
		decimal left_y = 0;
		decimal left_z = 0;
		decimal right_x = 0;
		decimal right_y = 0;
		decimal right_z = 0;
		for (integer j = 0; j < this->sim_mesh.n_vertices(); j++)
		{
		decimal w_from_j = W_mat(from_id, j);
		from_x += w_from_j * unknown(3 * j + 0);
		from_y += w_from_j * unknown(3 * j + 1);
		from_z += w_from_j * unknown(3 * j + 2);

		decimal w_to_j = W_mat(to_id, j);
		to_x += w_to_j * unknown(3 * j + 0);
		to_y += w_to_j * unknown(3 * j + 1);
		to_z += w_to_j * unknown(3 * j + 2);

		decimal w_left_j = W_mat(left_id, j);
		left_x += w_left_j * unknown(3 * j + 0);
		left_y += w_left_j * unknown(3 * j + 1);
		left_z += w_left_j * unknown(3 * j + 2);

		decimal w_right_j = W_mat(right_id, j);
		right_x += w_right_j * unknown(3 * j + 0);
		right_y += w_right_j * unknown(3 * j + 1);
		right_z += w_right_j * unknown(3 * j + 2);
		}*/

		MyMesh::Point rj(high_v_mat(from_id, 0), high_v_mat(from_id, 1), high_v_mat(from_id, 2));
		MyMesh::Point rk(high_v_mat(to_id, 0), high_v_mat(to_id, 1), high_v_mat(to_id, 2));
		MyMesh::Point ri(high_v_mat(left_id, 0), high_v_mat(left_id, 1), high_v_mat(left_id, 2));
		MyMesh::Point rl(high_v_mat(right_id, 0), high_v_mat(right_id, 1), high_v_mat(right_id, 2));

		decimal cur_angle = CalDiheral(rj, rk, ri, rl);
		/*fb(ben_ele_cnt) = swb * (cur_angle - wnt_dihedral_set[e_id]);
		energy = fb(ben_ele_cnt) * fb(ben_ele_cnt);
		++ben_ele_cnt;*/
		f = swb * (cur_angle - wnt_dihedral_set[e_id]);
		energy = f * f;
	}

	return energy;
}

std::vector<decimal> mns::MyDeformationHandler::CSL_Esb(Eigen::VectorXd& un)
{
	std::vector<decimal> energy;
	energy.push_back(0);
	energy.push_back(0);

	int i = 0;
	for (auto he_it = this->mesh.halfedges_begin(); he_it != this->mesh.halfedges_end(); ++i)
	{
		auto e = *he_it;
		auto f1 = this->mesh.face_handle(e);
		auto op_he = this->mesh.opposite_halfedge_handle(e);
		auto f2 = this->mesh.face_handle(op_he);
		auto fromVertex = this->mesh.from_vertex_handle(e);
		auto toVertex = this->mesh.to_vertex_handle(e);

		if (alpha_s != 0)
		{
			energy[0] += CSL_GetStrTerm(fromVertex.idx(), toVertex.idx(), i, un);
		}

		if (alpha_b != 0)
		{
			if (f1.idx() == -1 || f2.idx() == -1)
			{
				//edge, not half edge
				++he_it;
				++he_it;
				continue;
			}
			integer leftVertex, rightVertex;
			//integer temp_cnt = 0;
			for (auto fv_it1 = this->mesh.fv_begin(f1); fv_it1 != this->mesh.fv_end(f1); ++fv_it1)
			{
				auto v = *fv_it1;
				//temp_cnt++;
				if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
				{
					leftVertex = v.idx();
					break;
				}
			}
			//temp_cnt = 0;
			for (auto fv_it2 = this->mesh.fv_begin(f2); fv_it2 != this->mesh.fv_end(f2); ++fv_it2)
			{
				auto v = *fv_it2;
				//temp_cnt++;
				if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
				{
					rightVertex = v.idx();
					break;
				}
			}
			energy[1] += CSL_GetBenTerm(fromVertex.idx(), toVertex.idx(), leftVertex, rightVertex, i, un);
		}

		if (_isnan(energy[1]))
		{
			std::cout << i;
		}

		++he_it;
		++he_it;
	}

	return energy;
}

decimal mns::MyDeformationHandler::CSL_Ev(Eigen::VectorXd& un)
{
	decimal total_energy = 0;
	decimal total_tmp_vol = 0;
	decimal f = 0;
	if (this->use_sim_or_not == false)
	{
		for (auto f_it = this->mesh.faces_begin(); f_it != this->mesh.faces_end(); ++f_it)
		{
			auto face = *f_it;
			pt_vector pt_vec;
			int_vector idx_vec;
			for (auto fv_it = this->mesh.fv_begin(face); fv_it != this->mesh.fv_end(face); ++fv_it)
			{
				auto v = *fv_it;
				MyMesh::Point pt(un(3 * v.idx() + 0), un(3 * v.idx() + 1), un(3 * v.idx() + 2));

				idx_vec.push_back(v.idx());
				pt_vec.push_back(pt);
			}
			decimal temp_vol = (1.0 / 6.0) * ((pt_vec[0] % pt_vec[1]) | pt_vec[2]);
			total_tmp_vol += temp_vol;
		}
	}
	if (this->use_sim_or_not == true)
	{
		for (auto f_it = this->mesh.faces_begin(); f_it != this->mesh.faces_end(); ++f_it)
		{
			auto face = *f_it;
			pt_vector pt_vec;
			int_vector idx_vec;
			for (auto fv_it = this->mesh.fv_begin(face); fv_it != this->mesh.fv_end(face); ++fv_it)
			{
				auto v = *fv_it;
				MyMesh::Point pt(high_v_mat(v.idx(), 0), high_v_mat(v.idx(), 1), high_v_mat(v.idx(), 2));

				idx_vec.push_back(v.idx());
				pt_vec.push_back(pt);
			}
			decimal temp_vol = (1.0 / 6.0) * ((pt_vec[0] % pt_vec[1]) | pt_vec[2]);
			total_tmp_vol += temp_vol;
		}
	}

	f = swv*(total_tmp_vol - wnt_volume);
	total_energy = f*f;

	return total_energy;
}

decimal mns::MyDeformationHandler::CSL_GetWWTerm(integer idx, Eigen::VectorXd& un)
{
	decimal energy = 0;
	decimal f = 0;
	if (this->use_sim_or_not == false)
	{
		const integer vn = this->mesh.n_vertices();
		const integer vb = this->basis.cols();
		decimal uli_minus_wli = un(3 * vn + 2 * vb + 2 * idx + 0) - un(3 * vn + 2 * idx + 0);
		decimal uli_plus_wli = un(3 * vn + 2 * vb + 2 * idx + 0) + un(3 * vn + 2 * idx + 0);
		decimal udi_minus_wdi = un(3 * vn + 2 * vb + 2 * idx + 1) - un(3 * vn + 2 * idx + 1);
		decimal udi_plus_wdi = un(3 * vn + 2 * vb + 2 * idx + 1) + un(3 * vn + 2 * idx + 1);

		decimal sqrt_lsum = alpha_w*un(3 * vn + 2 * vb + 2 * idx + 0) - log(uli_minus_wli) - log(uli_plus_wli);
		sqrt_lsum = sqrt(sqrt_lsum);

		decimal sqrt_dsum = alpha_w*un(3 * vn + 2 * vb + 2 * idx + 1) - log(udi_minus_wdi) - log(udi_plus_wdi);
		sqrt_dsum = sqrt(sqrt_dsum);

		energy = sqrt_lsum*sqrt_lsum + sqrt_dsum*sqrt_dsum;
	}
	else
	{

	}
	return energy;
}

decimal mns::MyDeformationHandler::CSL_Ew(Eigen::VectorXd& un)
{
	decimal total_energy = 0;

	for (integer i = 0; i < this->basis.cols(); i++)
	{
		total_energy += CSL_GetWWTerm(i, un);
	}

	return total_energy;
}

void mns::MyDeformationHandler::SetSubspaceW(const Eigen::MatrixXd& w_mat, decimal lim)
{
	this->W_mat = w_mat;
	this->limit = lim;
	this->huge_W_mat = SMatXd(3 * W_mat.rows(), 3 * W_mat.cols());
	std::vector<T> w_triple;
	//std::vector<decimal> row_sum_vec;
	for (integer i = 0; i < W_mat.rows(); i++)
	{
		std::vector<integer> idx_vec;
		//decimal row_sum = 0;
		for (integer j = 0; j < W_mat.cols(); j++)
		{
			decimal wij = W_mat(i, j);
			if (std::abs(wij) > limit)
			{
				w_triple.push_back(Eigen::Triplet<decimal>(3 * i + 0, 3 * j + 0, wij));
				w_triple.push_back(Eigen::Triplet<decimal>(3 * i + 1, 3 * j + 1, wij));
				w_triple.push_back(Eigen::Triplet<decimal>(3 * i + 2, 3 * j + 2, wij));
				idx_vec.push_back(j);
				//row_sum += wij;
			}
			else
			{
				W_mat(i, j) = 0;
			}
		}
		//row_sum_vec.push_back(row_sum);
		this->W_idx_set.push_back(idx_vec);
	}
	/*for (integer i = 0; i < W_mat.rows(); i++)
	{
	decimal row_sum = row_sum_vec[i];
	for (integer j = 0; j < W_mat.cols(); j++)
	{
	if (std::abs(W_mat(i, j)) > 0)
	{
	W_mat(i, j) = W_mat(i, j) / row_sum;
	decimal wij = W_mat(i, j);
	w_triple.push_back(Eigen::Triplet<decimal>(3 * i + 0, 3 * j + 0, wij));
	w_triple.push_back(Eigen::Triplet<decimal>(3 * i + 1, 3 * j + 1, wij));
	w_triple.push_back(Eigen::Triplet<decimal>(3 * i + 2, 3 * j + 2, wij));
	}
	}
	}*/
	huge_W_mat.setFromTriplets(w_triple.begin(), w_triple.end());

	w_triple.clear();
}

bool mns::MyDeformationHandler::GetUseSim()
{
	return this->use_sim_or_not;
}

//Mesh weight estimation // first dihedral then length
void mns::MyDeformationHandler::MeshWeightEst()
{
	std::ifstream po_ifile("pose_nm.txt", std::ios::in);
	std::string pose_nm, basis_loc;
	std::string site;
	po_ifile >> site;
	std::string m_nm;
	po_ifile >> m_nm;
	std::string append;
	po_ifile >> append;
	po_ifile >> basis_loc;
	integer NB = 0;
	po_ifile >> NB;
	po_ifile.close();
	pose_nm = site + m_nm + append;

	if (this->GetRBF() == false)
	{
		ReceiveBasisWithName(basis_loc, NB);
	}
	Eigen::MatrixXd A = this->basis;

	dec_vector m_length_set;
	dec_vector m_dihedral_set;
	for (auto he_it = this->mesh.halfedges_begin(); he_it != this->mesh.halfedges_end(); ++he_it)
	{
		auto e = *he_it;
		auto fromVertex = this->mesh.from_vertex_handle(e);
		auto toVertex = this->mesh.to_vertex_handle(e);
		auto fromPt = this->mesh.point(fromVertex);
		auto toPt = this->mesh.point(toVertex);

		//stretching initialization		
		decimal l = this->mesh.calc_edge_length(e);
		m_length_set.push_back(l);

		//bend initialization
		auto f1 = this->mesh.face_handle(e);
		auto op_he = this->mesh.opposite_halfedge_handle(e);
		auto f2 = this->mesh.face_handle(op_he);
		pt_vector pt_set;
		pt_vector t1_set;
		pt_vector t2_set;

		if (f1.idx() == -1 || f2.idx() == -1)
		{
			m_dihedral_set.push_back(0);

			//edge, not half edge
			++he_it;
			continue;
		}

		pt_set.push_back(fromPt);
		t1_set.push_back(fromPt);
		t2_set.push_back(fromPt);
		pt_set.push_back(toPt);
		t1_set.push_back(toPt);
		t2_set.push_back(toPt);

		for (auto fv_it1 = this->mesh.fv_begin(f1); fv_it1 != this->mesh.fv_end(f1); ++fv_it1)
		{
			auto v = *fv_it1;
			auto pt = this->mesh.point(v);
			if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
			{
				pt_set.push_back(pt);
				t1_set.push_back(pt);
				break;
			}
		}
		integer temp_cnt = 0;
		for (auto fv_it2 = this->mesh.fv_begin(f2); fv_it2 != this->mesh.fv_end(f2); ++fv_it2)
		{
			auto v = *fv_it2;
			auto pt = this->mesh.point(v);
			if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
			{
				t2_set.push_back(pt);
				pt_set.push_back(pt);
				break;
			}
		}

		decimal angle = CalDiheral(pt_set[0], pt_set[1], pt_set[2], pt_set[3]);
		m_dihedral_set.push_back(angle);

		//edge, not half edge
		++he_it;
	}

	Eigen::VectorXd avg = Eigen::VectorXd::Zero(2 * this->mesh.n_edges());
	for (int i = 0; i < this->mesh.n_edges(); i++)
	{
		avg(2 * i + 0) = m_dihedral_set[i];
		avg(2 * i + 1) = m_length_set[i];
	}

	MyMesh pose;
	bool res = OpenMesh::IO::read_mesh(pose, pose_nm);
	if (res == false)
	{
		std::cout << "read error\n";
	}

	dec_vector po_length_set;
	dec_vector po_dihedral_set;

	//length & dihedral init
	for (auto he_it = pose.halfedges_begin(); he_it != pose.halfedges_end(); ++he_it)
	{
		auto e = *he_it;
		auto fromVertex = pose.from_vertex_handle(e);
		auto toVertex = pose.to_vertex_handle(e);
		auto fromPt = pose.point(fromVertex);
		auto toPt = pose.point(toVertex);

		//stretching initialization		
		decimal l = pose.calc_edge_length(e);
		po_length_set.push_back(l);

		//bend initialization
		auto f1 = pose.face_handle(e);
		auto op_he = pose.opposite_halfedge_handle(e);
		auto f2 = pose.face_handle(op_he);
		pt_vector pt_set;
		pt_vector t1_set;
		pt_vector t2_set;

		if (f1.idx() == -1 || f2.idx() == -1)
		{
			po_dihedral_set.push_back(0);

			//edge, not half edge
			++he_it;
			continue;
		}

		pt_set.push_back(fromPt);
		t1_set.push_back(fromPt);
		t2_set.push_back(fromPt);
		pt_set.push_back(toPt);
		t1_set.push_back(toPt);
		t2_set.push_back(toPt);

		for (auto fv_it1 = pose.fv_begin(f1); fv_it1 != pose.fv_end(f1); ++fv_it1)
		{
			auto v = *fv_it1;
			auto pt = pose.point(v);
			if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
			{
				pt_set.push_back(pt);
				t1_set.push_back(pt);
				break;
			}
		}
		integer temp_cnt = 0;
		for (auto fv_it2 = pose.fv_begin(f2); fv_it2 != pose.fv_end(f2); ++fv_it2)
		{
			auto v = *fv_it2;
			auto pt = pose.point(v);
			if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
			{
				t2_set.push_back(pt);
				pt_set.push_back(pt);
				break;
			}
		}

		decimal angle = CalDiheral(pt_set[0], pt_set[1], pt_set[2], pt_set[3]);
		po_dihedral_set.push_back(angle);

		//edge, not half edge
		++he_it;
	}

	Eigen::VectorXd eg = Eigen::VectorXd::Zero(2 * this->mesh.n_edges());
	for (int i = 0; i < this->mesh.n_edges(); i++)
	{
		eg(2 * i + 0) = po_dihedral_set[i];
		eg(2 * i + 1) = po_length_set[i];
	}

	Eigen::VectorXd b = eg - avg;
	//Eigen::VectorXd w = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
	Eigen::VectorXd w = A.colPivHouseholderQr().solve(b);

	std::string txt_head = "edit1_";
	std::string save_nm = site + txt_head + m_nm + ".txt";
	std::ofstream w_ofile(save_nm, std::ios::out);
	w_ofile << w.rows() << std::endl;
	for (int r = 0; r < w.rows(); r++)
	{
		w_ofile << r << " " << w[r] << std::endl;
	}
	w_ofile.close();

	Eigen::VectorXd aw = A*w;
	std::ofstream re_ofile("aw_cal.txt", std::ios::out);
	for (integer i = 0; i < aw.rows(); i++)
	{
		re_ofile << aw[i] << std::endl;
	}
	re_ofile.close();
	std::cout << "error: " << (aw - b).norm() << ", relative error: " << (A*w - b).norm() / b.norm();
}

//Pose weight estimation // first dihedral then length
void mns::MyDeformationHandler::PoseWeightEst()
{
	if (this->GetRBF() == false)
	{
		ReceiveBasis();
	}
	Eigen::MatrixXd A = this->basis;

	dec_vector m_length_set;
	dec_vector m_dihedral_set;
	for (auto he_it = this->mesh.halfedges_begin(); he_it != this->mesh.halfedges_end(); ++he_it)
	{
		auto e = *he_it;
		auto fromVertex = this->mesh.from_vertex_handle(e);
		auto toVertex = this->mesh.to_vertex_handle(e);
		auto fromPt = this->mesh.point(fromVertex);
		auto toPt = this->mesh.point(toVertex);

		//stretching initialization		
		decimal l = this->mesh.calc_edge_length(e);
		m_length_set.push_back(l);

		//bend initialization
		auto f1 = this->mesh.face_handle(e);
		auto op_he = this->mesh.opposite_halfedge_handle(e);
		auto f2 = this->mesh.face_handle(op_he);
		pt_vector pt_set;
		pt_vector t1_set;
		pt_vector t2_set;

		if (f1.idx() == -1 || f2.idx() == -1)
		{
			m_dihedral_set.push_back(0);

			//edge, not half edge
			++he_it;
			continue;
		}

		pt_set.push_back(fromPt);
		t1_set.push_back(fromPt);
		t2_set.push_back(fromPt);
		pt_set.push_back(toPt);
		t1_set.push_back(toPt);
		t2_set.push_back(toPt);

		for (auto fv_it1 = this->mesh.fv_begin(f1); fv_it1 != this->mesh.fv_end(f1); ++fv_it1)
		{
			auto v = *fv_it1;
			auto pt = this->mesh.point(v);
			if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
			{
				pt_set.push_back(pt);
				t1_set.push_back(pt);
				break;
			}
		}
		integer temp_cnt = 0;
		for (auto fv_it2 = this->mesh.fv_begin(f2); fv_it2 != this->mesh.fv_end(f2); ++fv_it2)
		{
			auto v = *fv_it2;
			auto pt = this->mesh.point(v);
			if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
			{
				t2_set.push_back(pt);
				pt_set.push_back(pt);
				break;
			}
		}

		decimal angle = CalDiheral(pt_set[0], pt_set[1], pt_set[2], pt_set[3]);
		m_dihedral_set.push_back(angle);

		//edge, not half edge
		++he_it;
	}

	Eigen::VectorXd avg = Eigen::VectorXd::Zero(2 * this->mesh.n_edges());
	for (int i = 0; i < this->mesh.n_edges(); i++)
	{
		avg(2 * i + 0) = m_dihedral_set[i];
		avg(2 * i + 1) = m_length_set[i];
	}

	std::ifstream ex_ifile("ex_nm.txt", std::ios::in);
	std::string example_str;
	ex_ifile >> example_str;
	integer eg_num = 0;
	ex_ifile >> eg_num;
	decimal sam_per = 0.01;
	ex_ifile >> sam_per;
	ex_ifile.close();

	std::string eg_begin = example_str;
	std::string eg_end;
	std::string eg_idx0 = "1";
	std::string eg_idx1 = "0";
	std::string eg_idx2 = "0";
	std::string eg_idx3 = "0";
	if (model_mod == 1)
	{
		eg_end = ".ply";
	}
	if (model_mod == 2)
	{
		eg_end = ".off";
	}
	if (model_mod == 3)
	{
		eg_end = ".obj";
	}

	for (integer k = 0; k < eg_num; k++)
	{
		MyMesh eg_mesh;
		std::string eg_idx;
		if ((k + 1) % 10 == 0 && (k + 1) % 100 != 0)
		{
			eg_idx0 = "0";
			eg_idx1[0]++;
			//eg_idx = eg_idx1 + eg_idx0;
		}
		if ((k + 1) % 100 == 0 && (k + 1) % 1000 != 0)
		{
			eg_idx0 = "0";
			eg_idx1 = "0";
			eg_idx2[0]++;
		}
		eg_idx = eg_idx2 + eg_idx1 + eg_idx0;
		std::string eg_nm = eg_begin + eg_idx + eg_end;

		bool res = OpenMesh::IO::read_mesh(eg_mesh, eg_nm);
		if (res == false)
		{
			std::cout << "read error\n";
		}

		dec_vector eg_length_set;
		dec_vector eg_dihedral_set;

		//length & dihedral init
		for (auto he_it = eg_mesh.halfedges_begin(); he_it != eg_mesh.halfedges_end(); ++he_it)
		{
			auto e = *he_it;
			auto fromVertex = eg_mesh.from_vertex_handle(e);
			auto toVertex = eg_mesh.to_vertex_handle(e);
			auto fromPt = eg_mesh.point(fromVertex);
			auto toPt = eg_mesh.point(toVertex);

			//stretching initialization		
			decimal l = eg_mesh.calc_edge_length(e);
			eg_length_set.push_back(l);

			//bend initialization
			auto f1 = eg_mesh.face_handle(e);
			auto op_he = eg_mesh.opposite_halfedge_handle(e);
			auto f2 = eg_mesh.face_handle(op_he);
			pt_vector pt_set;
			pt_vector t1_set;
			pt_vector t2_set;

			if (f1.idx() == -1 || f2.idx() == -1)
			{
				eg_dihedral_set.push_back(0);

				//edge, not half edge
				++he_it;
				continue;
			}

			pt_set.push_back(fromPt);
			t1_set.push_back(fromPt);
			t2_set.push_back(fromPt);
			pt_set.push_back(toPt);
			t1_set.push_back(toPt);
			t2_set.push_back(toPt);

			for (auto fv_it1 = eg_mesh.fv_begin(f1); fv_it1 != eg_mesh.fv_end(f1); ++fv_it1)
			{
				auto v = *fv_it1;
				auto pt = eg_mesh.point(v);
				if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
				{
					pt_set.push_back(pt);
					t1_set.push_back(pt);
					break;
				}
			}
			integer temp_cnt = 0;
			for (auto fv_it2 = eg_mesh.fv_begin(f2); fv_it2 != eg_mesh.fv_end(f2); ++fv_it2)
			{
				auto v = *fv_it2;
				auto pt = eg_mesh.point(v);
				if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
				{
					t2_set.push_back(pt);
					pt_set.push_back(pt);
					break;
				}
			}

			decimal angle = CalDiheral(pt_set[0], pt_set[1], pt_set[2], pt_set[3]);
			eg_dihedral_set.push_back(angle);

			//edge, not half edge
			++he_it;
		}

		Eigen::VectorXd eg = Eigen::VectorXd::Zero(2 * this->mesh.n_edges());
		for (int i = 0; i < this->mesh.n_edges(); i++)
		{
			eg(2 * i + 0) = eg_dihedral_set[i];
			eg(2 * i + 1) = eg_length_set[i];
		}

		Eigen::VectorXd b = eg - avg;
		Eigen::VectorXd w = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);

		std::string txt_head = "edit1_";
		std::string save_nm = txt_head + eg_idx + ".txt";
		std::ofstream w_ofile(save_nm, std::ios::out);
		w_ofile << w.rows() << std::endl;
		for (int r = 0; r < w.rows(); r++)
		{
			w_ofile << r << " " << w[r] << std::endl;
		}
		w_ofile.close();
	}
}

// first dihedral then length
void mns::MyDeformationHandler::PoseAndWeightCompare(std::string pose_nm, std::string wei_nm, integer NB, integer selected_frame, std::string type)
{
	if (this->GetRBF() == false)
	{
		return;
	}
	Eigen::MatrixXd A = this->basis;

	dec_vector m_length_set;
	dec_vector m_dihedral_set;
	for (auto he_it = this->mesh.halfedges_begin(); he_it != this->mesh.halfedges_end(); ++he_it)
	{
		auto e = *he_it;
		auto fromVertex = this->mesh.from_vertex_handle(e);
		auto toVertex = this->mesh.to_vertex_handle(e);
		auto fromPt = this->mesh.point(fromVertex);
		auto toPt = this->mesh.point(toVertex);

		//stretching initialization		
		decimal l = this->mesh.calc_edge_length(e);
		m_length_set.push_back(l);

		//bend initialization
		auto f1 = this->mesh.face_handle(e);
		auto op_he = this->mesh.opposite_halfedge_handle(e);
		auto f2 = this->mesh.face_handle(op_he);
		pt_vector pt_set;
		pt_vector t1_set;
		pt_vector t2_set;

		if (f1.idx() == -1 || f2.idx() == -1)
		{
			m_dihedral_set.push_back(0);

			//edge, not half edge
			++he_it;
			continue;
		}

		pt_set.push_back(fromPt);
		t1_set.push_back(fromPt);
		t2_set.push_back(fromPt);
		pt_set.push_back(toPt);
		t1_set.push_back(toPt);
		t2_set.push_back(toPt);

		for (auto fv_it1 = this->mesh.fv_begin(f1); fv_it1 != this->mesh.fv_end(f1); ++fv_it1)
		{
			auto v = *fv_it1;
			auto pt = this->mesh.point(v);
			if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
			{
				pt_set.push_back(pt);
				t1_set.push_back(pt);
				break;
			}
		}
		integer temp_cnt = 0;
		for (auto fv_it2 = this->mesh.fv_begin(f2); fv_it2 != this->mesh.fv_end(f2); ++fv_it2)
		{
			auto v = *fv_it2;
			auto pt = this->mesh.point(v);
			if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
			{
				t2_set.push_back(pt);
				pt_set.push_back(pt);
				break;
			}
		}

		decimal angle = CalDiheral(pt_set[0], pt_set[1], pt_set[2], pt_set[3]);
		m_dihedral_set.push_back(angle);

		//edge, not half edge
		++he_it;
	}

	Eigen::VectorXd avg = Eigen::VectorXd::Zero(2 * this->mesh.n_edges());
	for (int i = 0; i < this->mesh.n_edges(); i++)
	{
		avg(2 * i + 0) = m_dihedral_set[i];
		avg(2 * i + 1) = m_length_set[i];
	}

	MyMesh pose;
	bool res = OpenMesh::IO::read_mesh(pose, pose_nm);
	if (res == false)
	{
		std::cout << "read error\n";
	}

	dec_vector po_length_set;
	dec_vector po_dihedral_set;

	//length & dihedral init
	for (auto he_it = pose.halfedges_begin(); he_it != pose.halfedges_end(); ++he_it)
	{
		auto e = *he_it;
		auto fromVertex = pose.from_vertex_handle(e);
		auto toVertex = pose.to_vertex_handle(e);
		auto fromPt = pose.point(fromVertex);
		auto toPt = pose.point(toVertex);

		//stretching initialization		
		decimal l = pose.calc_edge_length(e);
		po_length_set.push_back(l);

		//bend initialization
		auto f1 = pose.face_handle(e);
		auto op_he = pose.opposite_halfedge_handle(e);
		auto f2 = pose.face_handle(op_he);
		pt_vector pt_set;
		pt_vector t1_set;
		pt_vector t2_set;

		if (f1.idx() == -1 || f2.idx() == -1)
		{
			po_dihedral_set.push_back(0);

			//edge, not half edge
			++he_it;
			continue;
		}

		pt_set.push_back(fromPt);
		t1_set.push_back(fromPt);
		t2_set.push_back(fromPt);
		pt_set.push_back(toPt);
		t1_set.push_back(toPt);
		t2_set.push_back(toPt);

		for (auto fv_it1 = pose.fv_begin(f1); fv_it1 != pose.fv_end(f1); ++fv_it1)
		{
			auto v = *fv_it1;
			auto pt = pose.point(v);
			if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
			{
				pt_set.push_back(pt);
				t1_set.push_back(pt);
				break;
			}
		}
		integer temp_cnt = 0;
		for (auto fv_it2 = pose.fv_begin(f2); fv_it2 != pose.fv_end(f2); ++fv_it2)
		{
			auto v = *fv_it2;
			auto pt = pose.point(v);
			if (v.idx() != fromVertex.idx() && v.idx() != toVertex.idx())
			{
				t2_set.push_back(pt);
				pt_set.push_back(pt);
				break;
			}
		}

		decimal angle = CalDiheral(pt_set[0], pt_set[1], pt_set[2], pt_set[3]);
		po_dihedral_set.push_back(angle);

		//edge, not half edge
		++he_it;
	}

	Eigen::VectorXd eg = Eigen::VectorXd::Zero(2 * this->mesh.n_edges());
	for (int i = 0; i < this->mesh.n_edges(); i++)
	{
		eg(2 * i + 0) = po_dihedral_set[i];
		eg(2 * i + 1) = po_length_set[i];
	}

	Eigen::VectorXd b = eg - avg;

	Eigen::VectorXd W = Eigen::VectorXd::Zero(A.cols());
	/*integer cnt = 0;
	std::ifstream in_file(wei_nm);
	in_file >> cnt;
	for (int i = 0; i < cnt; i++)
	{
	integer itmp;
	decimal dtmp;
	in_file >> itmp >> dtmp;
	W(itmp) = dtmp;
	}
	in_file.close();*/
	std::ifstream w_ifile(wei_nm, std::ios::in);
	dec_vector w_vec;
	while (w_ifile.peek() != EOF)
	{
		decimal dec_tmp;
		w_ifile >> dec_tmp;
		w_vec.push_back(dec_tmp);
	}
	w_vec.erase(w_vec.begin() + w_vec.size() - 1);
	w_ifile.close();

	std::stringstream newstr;
	newstr << selected_frame;

	std::string output_nm = "edit1_" + type + "_sf_" + newstr.str() + "_fromW.txt";
	std::ofstream ofile(output_nm);
	ofile << NB << std::endl;
	for (int i = 0; i < NB; i++)
	{
		W[i] = w_vec[selected_frame * NB + i];
		ofile << i << " " << W[i] << std::endl;
	}
	ofile.close();

	Eigen::VectorXd aw = A*W;
	std::ofstream re_ofile("aw_w.txt", std::ios::out);
	for (integer i = 0; i < aw.rows(); i++)
	{
		re_ofile << aw[i] << std::endl;
	}
	re_ofile.close();
	decimal error = (aw - b).norm();
	decimal bnorm = b.norm();
	std::cout << "error: " << error << ", relative error: " << error / bnorm;
}


/** 根据点坐标导出网格*/
