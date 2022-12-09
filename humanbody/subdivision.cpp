#include "subdivision.h"
using namespace std;
using namespace Eigen;

AdaptiveSubdivision::AdaptiveSubdivision(PGMesh * _mesh)
{
	subdividedMesh = new PGMesh();
	mesh_ = _mesh;
	nV = mesh_->n_vertices();
	nE = mesh_->n_edges();
	nF = mesh_->n_faces();
	mesh_->add_property(is_head_edge);
	mesh_->add_property(is_head_face);
	mesh_->add_property(is_head_vertice);
	mesh_->add_property(edge_vertice_idx);
	for (auto vit = mesh_->vertices_begin(); vit != mesh_->vertices_end(); ++vit)
	{
		mesh_->property(is_head_vertice, *vit) = false;
	}
	for (auto eit = mesh_->edges_begin(); eit != mesh_->edges_end(); ++eit)
	{
		mesh_->property(is_head_edge, *eit) = false;
		mesh_->property(edge_vertice_idx, *eit) = -1;
	}
	for (auto fit = mesh_->faces_begin(); fit != mesh_->faces_end(); ++fit)
	{
		mesh_->property(is_head_face, *fit) = false;
	}
}

void AdaptiveSubdivision::getHeadVertices(char* h_path)
{
	ifstream in(h_path);
	affected_vertices.resize(0);
	int iter = 0;
	while (!in.eof() ) {
		int tmp;
		in >> tmp;
		affected_vertices.push_back(tmp);
		iter++;
		auto vit = mesh_->vertex_handle(tmp);
		mesh_->property(is_head_vertice, vit) = true;
		//cout << i << " " <<j << endl;
	}
	n_affected_vertices = iter;
	in.close();
}

void AdaptiveSubdivision::setSubdividedArea()
{
	int e_idx = 0;
	affected_edges.resize(0);
	for (auto eit = mesh_->edges_begin(); eit != mesh_->edges_end(); ++eit)
	{
		auto heit = mesh_->halfedge_handle(*eit, 0);
		auto v0 = mesh_->from_vertex_handle(heit);
		auto v1 = mesh_->to_vertex_handle(heit);
		if (mesh_->property(is_head_vertice, v0) == true && mesh_->property(is_head_vertice, v1) == true)
		{
			affected_edges.push_back(eit->idx());
			mesh_->property(edge_vertice_idx, *eit) = nV + e_idx;
			e_idx++;
		}
	}
	n_affected_edges = affected_edges.size();
	for (auto fit = mesh_->faces_begin(); fit != mesh_->faces_end(); ++fit)
	{
		auto vit = mesh_->fv_begin(*fit);
		auto v0 = *vit;
		++vit;
		auto v1 = *vit;
		++vit;
		auto v2 = *vit;
		if ((mesh_->property(is_head_vertice, v0))  || (mesh_->property(is_head_vertice, v1) ) || (mesh_->property(is_head_vertice, v2)) )
		{
			mesh_->property(is_head_face, *fit) = true;
		}
	}

}

void AdaptiveSubdivision::subdivision()
{
	/*add vertices*/
	std::vector<PGMesh::VertexHandle> vhandle(nV + n_affected_edges);
	for (auto vit = mesh_->vertices_begin(); vit != mesh_->vertices_end(); ++vit)
	{
		PGMesh::Point p = mesh_->point(*vit);
		vhandle[vit->idx()] = subdividedMesh->add_vertex(p);
	}
	for (int i = 0; i < affected_edges.size(); i ++)
	{
		auto eit = mesh_->edge_handle(affected_edges[i]);
		auto heit = mesh_->halfedge_handle(eit, 0);
		auto v0 = mesh_->from_vertex_handle(heit);
		auto v1 = mesh_->to_vertex_handle(heit);
		PGMesh::Point p0 = mesh_->point(v0);
		PGMesh::Point p1 = mesh_->point(v1);
		PGMesh::Point p = (p0 + p1) / 2;
		vhandle[nV + i] = subdividedMesh->add_vertex(p);
	}
	/*add faces*/
	std::vector<PGMesh::VertexHandle>  face_vhandles;
	for (auto fit = mesh_->faces_begin(); fit != mesh_->faces_end(); ++fit)
	{
		
		auto heit = mesh_->fh_begin(*fit);
		auto he0 = *heit;
		auto e0 = mesh_->edge_handle(he0);
		auto v0 = mesh_->from_vertex_handle(he0);
		auto v1 = mesh_->to_vertex_handle(he0);
		int v0_idx = v0.idx();
		int v1_idx = v1.idx();
		int ve0_idx = mesh_->property(edge_vertice_idx, e0);
		++heit;
		auto he1 = *heit;
		auto e1 = mesh_->edge_handle(he1);
		auto ve1_idx = mesh_->property(edge_vertice_idx, e1);
		auto v2 = mesh_->to_vertex_handle(he1);
		int v2_idx = v2.idx();
		++heit;
		auto he2 = *heit;
		auto e2 = mesh_->edge_handle(he2);
		auto ve2_idx = mesh_->property(edge_vertice_idx, e2);
		if (!(mesh_->property(is_head_face, *fit)))
		{
			face_vhandles.clear();
			face_vhandles.push_back(vhandle[v0_idx]);
			face_vhandles.push_back(vhandle[v1_idx]);
			face_vhandles.push_back(vhandle[v2_idx]);
			subdividedMesh->add_face(face_vhandles);
		}
		
		else if (ve0_idx > -1)
		{
			if (ve1_idx > -1)
			{
				if (ve2_idx > -1)
				{
					face_vhandles.clear();
					face_vhandles.push_back(vhandle[v0_idx]);
					face_vhandles.push_back(vhandle[ve0_idx]);
					face_vhandles.push_back(vhandle[ve2_idx]);
					subdividedMesh->add_face(face_vhandles);

					face_vhandles.clear();
					face_vhandles.push_back(vhandle[ve0_idx]);
					face_vhandles.push_back(vhandle[v1_idx]);
					face_vhandles.push_back(vhandle[ve1_idx]);
					subdividedMesh->add_face(face_vhandles);

					face_vhandles.clear();
					face_vhandles.push_back(vhandle[ve1_idx]);
					face_vhandles.push_back(vhandle[v2_idx]);
					face_vhandles.push_back(vhandle[ve2_idx]);
					subdividedMesh->add_face(face_vhandles);

					face_vhandles.clear();
					face_vhandles.push_back(vhandle[ve0_idx]);
					face_vhandles.push_back(vhandle[ve1_idx]);
					face_vhandles.push_back(vhandle[ve2_idx]);
					subdividedMesh->add_face(face_vhandles);
				}
				else
				{
					face_vhandles.clear();
					face_vhandles.push_back(vhandle[v0_idx]);
					face_vhandles.push_back(vhandle[ve0_idx]);
					face_vhandles.push_back(vhandle[ve1_idx]);
					subdividedMesh->add_face(face_vhandles);

					face_vhandles.clear();
					face_vhandles.push_back(vhandle[ve0_idx]);
					face_vhandles.push_back(vhandle[v1_idx]);
					face_vhandles.push_back(vhandle[ve1_idx]);
					subdividedMesh->add_face(face_vhandles);

					face_vhandles.clear();
					face_vhandles.push_back(vhandle[ve1_idx]);
					face_vhandles.push_back(vhandle[v2_idx]);
					face_vhandles.push_back(vhandle[v0_idx]);
					subdividedMesh->add_face(face_vhandles);
				}
			}
			else
			{
				if (ve2_idx > -1)
				{
					face_vhandles.clear();
					face_vhandles.push_back(vhandle[v0_idx]);
					face_vhandles.push_back(vhandle[ve0_idx]);
					face_vhandles.push_back(vhandle[ve2_idx]);
					subdividedMesh->add_face(face_vhandles);

					face_vhandles.clear();
					face_vhandles.push_back(vhandle[ve0_idx]);
					face_vhandles.push_back(vhandle[v1_idx]);
					face_vhandles.push_back(vhandle[ve2_idx]);
					subdividedMesh->add_face(face_vhandles);

					face_vhandles.clear();
					face_vhandles.push_back(vhandle[v1_idx]);
					face_vhandles.push_back(vhandle[v2_idx]);
					face_vhandles.push_back(vhandle[ve2_idx]);
					subdividedMesh->add_face(face_vhandles);
				}

				else
				{
					face_vhandles.clear();
					face_vhandles.push_back(vhandle[v0_idx]);
					face_vhandles.push_back(vhandle[ve0_idx]);
					face_vhandles.push_back(vhandle[v2_idx]);
					subdividedMesh->add_face(face_vhandles);

					face_vhandles.clear();
					face_vhandles.push_back(vhandle[ve0_idx]);
					face_vhandles.push_back(vhandle[v1_idx]);
					face_vhandles.push_back(vhandle[v2_idx]);
					subdividedMesh->add_face(face_vhandles);
				}
				
			}
		}
		else if (ve1_idx > -1)
		{
			if (ve2_idx > -1)
			{
				face_vhandles.clear();
				face_vhandles.push_back(vhandle[v0_idx]);
				face_vhandles.push_back(vhandle[v1_idx]);
				face_vhandles.push_back(vhandle[ve1_idx]);
				subdividedMesh->add_face(face_vhandles);

				face_vhandles.clear();
				face_vhandles.push_back(vhandle[v0_idx]);
				face_vhandles.push_back(vhandle[ve1_idx]);
				face_vhandles.push_back(vhandle[ve2_idx]);
				subdividedMesh->add_face(face_vhandles);

				face_vhandles.clear();
				face_vhandles.push_back(vhandle[ve2_idx]);
				face_vhandles.push_back(vhandle[ve1_idx]);
				face_vhandles.push_back(vhandle[v2_idx]);
				subdividedMesh->add_face(face_vhandles);
			}
			else
			{
				face_vhandles.clear();
				face_vhandles.push_back(vhandle[v0_idx]);
				face_vhandles.push_back(vhandle[v1_idx]);
				face_vhandles.push_back(vhandle[ve1_idx]);
				subdividedMesh->add_face(face_vhandles);

				face_vhandles.clear();
				face_vhandles.push_back(vhandle[ve1_idx]);
				face_vhandles.push_back(vhandle[v2_idx]);
				face_vhandles.push_back(vhandle[v0_idx]);
				subdividedMesh->add_face(face_vhandles);
			}
		}
		else if (ve2_idx > -1)
		{
			face_vhandles.clear();
			face_vhandles.push_back(vhandle[v0_idx]);
			face_vhandles.push_back(vhandle[v1_idx]);
			face_vhandles.push_back(vhandle[ve2_idx]);
			subdividedMesh->add_face(face_vhandles);

			face_vhandles.clear();
			face_vhandles.push_back(vhandle[ve2_idx]);
			face_vhandles.push_back(vhandle[v1_idx]);
			face_vhandles.push_back(vhandle[v2_idx]);
			subdividedMesh->add_face(face_vhandles);

		}
		else
		{
			face_vhandles.clear();
			face_vhandles.push_back(vhandle[v0_idx]);
			face_vhandles.push_back(vhandle[v1_idx]);
			face_vhandles.push_back(vhandle[v2_idx]);
			subdividedMesh->add_face(face_vhandles);
		}
	}

	OpenMesh::IO::write_mesh(*subdividedMesh, "subdivide.obj");
}

void AdaptiveSubdivision::setSubdivisionCoeffs()
{
	subdivisionCoeff = Eigen::SparseMatrix<double>(nV + n_affected_edges, nV);
	std::vector<Eigen::Triplet<double> > triple;
	for (int i = 0; i < nV; i++)
	{
		triple.push_back(Eigen::Triplet<double>(i, i, 1));
		subdivisionCoeff.insert(i, i) = 1;
	}
	/*std::vector<int> index;
	std::vector<double> coef;*/
	for (int i = 0; i < affected_edges.size(); i++)
	{
		if (nV + i == 11159)
		{
			cout << i << endl;
		}
		auto e = mesh_->edge_handle(affected_edges[i]);
		auto he = mesh_->halfedge_handle(e, 0);
		auto he1 = mesh_->halfedge_handle(e, 1);
		auto v0 = mesh_->from_vertex_handle(he);
		auto v1 = mesh_->to_vertex_handle(he);
		std::vector<int> l0 = det(v0);
		std::vector<int> l1 = det(v1);
		int det_v0 = l0.size();
		int det_v1 = l1.size();
		if (mesh_->is_boundary(e))
		{
			PGMesh::VertexHandle v0b;
			for (auto vvit = mesh_->vv_begin(v0); vvit != mesh_->vv_end(v0); ++vvit)
			{
				if (*vvit != v1 && mesh_->is_boundary(*vvit))
				{
					v0b = *vvit;
				}
			}
			PGMesh::VertexHandle v1b;
			for (auto vvit = mesh_->vv_begin(v1); vvit != mesh_->vv_end(v1); ++vvit)
			{
				if (*vvit != v0 && mesh_->is_boundary(*vvit))
				{
					v1b = *vvit;
				}
			}
			triple.push_back(Eigen::Triplet<double>(nV + i, v0.idx(), 9.0 / 16));
			triple.push_back(Eigen::Triplet<double>(nV + i, v1.idx(), 9.0 / 16));
			triple.push_back(Eigen::Triplet<double>(nV + i, v0b.idx(), -1.0 / 16));
			triple.push_back(Eigen::Triplet<double>(nV + i, v1b.idx(), -1.0 / 16));

		}
		else
		{
			//if (det_v0 == 6 && det_v1 == 6)
			if(1)
			{
				triple.push_back(Eigen::Triplet<double>(nV + i, v0.idx(), 0.5));
				triple.push_back(Eigen::Triplet<double>(nV + i, v1.idx(), 0.5));
				/*subdivisionCoeff.insert(nV + i, v0.idx()) = 0.5;
				subdivisionCoeff.insert(nV + i, v1.idx()) = 0.5;*/
				auto f0 = mesh_->face_handle(he);
				auto f1 = mesh_->face_handle(he1);
				auto fv0 = mesh_->fv_begin(f0);
				PGMesh::VertexHandle b0;
				if (*fv0 != v0 && *fv0 != v1)
					b0 = *fv0;
				++fv0;
				if (*fv0 != v0 && *fv0 != v1)
					b0 = *fv0;
				++fv0;
				if (*fv0 != v0 && *fv0 != v1)
					b0 = *fv0;

				auto fv1 = mesh_->fv_begin(f1);
				PGMesh::VertexHandle b1;
				if (*fv1 != v0 && *fv1 != v1)
					b1 = *fv1;
				++fv1;
				if (*fv1 != v0 && *fv1 != v1)
					b1 = *fv1;
				++fv1;
				if (*fv1 != v0 && *fv1 != v1)
					b1 = *fv1;
				triple.push_back(Eigen::Triplet<double>(nV + i, b0.idx(), 0.125));
				triple.push_back(Eigen::Triplet<double>(nV + i, b1.idx(), 0.125));
				/*subdivisionCoeff.insert(nV + i, b0.idx()) = 0.125;
				subdivisionCoeff.insert(nV + i, b1.idx()) = 0.125;*/
				auto fh0 = mesh_->fh_begin(f0);
				PGMesh::HalfedgeHandle h0[2];
				int idx = 0;
				if (*fh0 != he)
				{
					h0[idx] = *fh0;
					idx++;
				}
				++fh0;
				if (*fh0 != he)
				{
					h0[idx] = *fh0;
					idx++;
				}
				++fh0;
				if (*fh0 != he)
				{
					h0[idx] = *fh0;
					idx++;
				}

				auto fh1 = mesh_->fh_begin(f1);
				PGMesh::HalfedgeHandle h1[2];
				idx = 0;
				if (*fh1 != he1)
				{
					h1[idx] = *fh1;
					idx++;
				}
				++fh1;
				if (*fh1 != he1)
				{
					h1[idx] = *fh1;
					idx++;
				}
				++fh1;
				if (*fh1 != he1)
				{
					h1[idx] = *fh1;
					idx++;
				}

				auto fc0 = mesh_->opposite_face_handle(h0[0]);
				auto fc1 = mesh_->opposite_face_handle(h0[1]);
				auto fc2 = mesh_->opposite_face_handle(h1[0]);
				auto fc3 = mesh_->opposite_face_handle(h1[1]);
				PGMesh::VertexHandle vc[4];
				for (auto fvit = mesh_->fv_begin(fc0); fvit != mesh_->fv_end(fc0); ++fvit)
				{
					if (*fvit != mesh_->from_vertex_handle(h0[0]) && *fvit != mesh_->to_vertex_handle(h0[0]))
					{
						vc[0] = *fvit;
					}
				}

				for (auto fvit = mesh_->fv_begin(fc1); fvit != mesh_->fv_end(fc1); ++fvit)
				{
					if (*fvit != mesh_->from_vertex_handle(h0[1]) && *fvit != mesh_->to_vertex_handle(h0[1]))
					{
						vc[1] = *fvit;
					}
				}

				for (auto fvit = mesh_->fv_begin(fc2); fvit != mesh_->fv_end(fc2); ++fvit)
				{
					if (*fvit != mesh_->from_vertex_handle(h1[0]) && *fvit != mesh_->to_vertex_handle(h1[0]))
					{
						vc[2] = *fvit;
					}
				}

				for (auto fvit = mesh_->fv_begin(fc3); fvit != mesh_->fv_end(fc3); ++fvit)
				{
					if (*fvit != mesh_->from_vertex_handle(h1[1]) && *fvit != mesh_->to_vertex_handle(h1[1]))
					{
						vc[3] = *fvit;
					}
				}

				for (int w = 0; w < 4; w++)
				{
					triple.push_back(Eigen::Triplet<double>(nV + i, vc[w].idx(), -0.125 / 2));
					//subdivisionCoeff.insert(nV + i, vc[w].idx())= 0.125 / 2;
				}
			}

			else if (det_v0 == 6 && det_v1 != 6)
			{
				std::vector<int> adjust_l1 = adjust_seq(l1, v0.idx());
				if (det_v1 == 3)
				{
					triple.push_back(Eigen::Triplet<double>(nV + i, adjust_l1[0], 5.0 / 12));
					triple.push_back(Eigen::Triplet<double>(nV + i, adjust_l1[1], -1.0 / 12));
					triple.push_back(Eigen::Triplet<double>(nV + i, adjust_l1[2], -1.0 / 12));
					triple.push_back(Eigen::Triplet<double>(nV + i, v1.idx(), 3.0/4));
					/*subdivisionCoeff.insert(nV + i, adjust_l1[0]) = 5.0 / 12;
					subdivisionCoeff.insert(nV + i, adjust_l1[0]) = -1.0 / 12;
					subdivisionCoeff.insert(nV + i, adjust_l1[0]) = -1.0 / 12;*/
				}
				if (det_v1 == 4)
				{
					triple.push_back(Eigen::Triplet<double>(nV + i, adjust_l1[0], 3.0 / 8));
					triple.push_back(Eigen::Triplet<double>(nV + i, adjust_l1[2], -1.0 / 8));
					triple.push_back(Eigen::Triplet<double>(nV + i, v1.idx(), 3.0 / 4));
					/*subdivisionCoeff.insert(nV + i, adjust_l1[0]) = 3.0/8 ;
					subdivisionCoeff.insert(nV + i, adjust_l1[2]) = -1.0 / 8;*/

				}
				if (det_v1 >= 5)
				{
					double cc1 = 0;
					for (int hhh = 0; hhh < det_v1; hhh++)
					{
						triple.push_back(Eigen::Triplet<double>(nV + i, adjust_l1[hhh], (0.25 + std::cos(2 * 3.141592653 * hhh / det_v1) + 0.5 * std::cos(4 * 3.141592653 * hhh / det_v1)) / det_v1));
						cc1 += (0.25 + std::cos(2 * 3.141592653 * hhh / det_v1) + 0.5 * std::cos(4 * 3.141592653 * hhh / det_v1)) / det_v1;
						//subdivisionCoeff.insert(nV + i, adjust_l1[hhh]) = (0.25 + std::cos(2 * 3.141592653 * hhh / det_v1) + 0.5 * std::cos(4 * 3.141592653 * hhh / det_v1)) / det_v1;
					}
					triple.push_back(Eigen::Triplet<double>(nV + i, v1.idx(), 1 - cc1));
				}
			}

			else if (det_v0 != 6 && det_v1 == 6)
			{
				std::vector<int> adjust_l0 = adjust_seq(l0, v1.idx());
				if (det_v0 == 3)
				{
					triple.push_back(Eigen::Triplet<double>(nV + i, adjust_l0[0], 5.0 / 12));
					triple.push_back(Eigen::Triplet<double>(nV + i, adjust_l0[1], -1.0 / 12));
					triple.push_back(Eigen::Triplet<double>(nV + i, adjust_l0[2], -1.0 / 12));
					triple.push_back(Eigen::Triplet<double>(nV + i, v0.idx(), 3.0 / 4));
					/*subdivisionCoeff.insert(nV + i, adjust_l0[0]) = 5.0/12;
					subdivisionCoeff.insert(nV + i, adjust_l0[1]) = -1.0/12;
					subdivisionCoeff.insert(nV + i, adjust_l0[2]) = -1.0/12;*/
				}
				if (det_v0 == 4)
				{
					triple.push_back(Eigen::Triplet<double>(nV + i, adjust_l0[0], 3.0 / 8));
					triple.push_back(Eigen::Triplet<double>(nV + i, adjust_l0[2], -1.0 / 8));
					triple.push_back(Eigen::Triplet<double>(nV + i, v0.idx(), 3.0 / 4));
					/*subdivisionCoeff.insert(nV + i, adjust_l0[0]) = 3.0 / 8;
					subdivisionCoeff.insert(nV + i, adjust_l0[2]) = -1.0 / 8;*/
				}
				if (det_v0 >= 5)
				{
					double ccl = 0;
					for (int hhh = 0; hhh < det_v0; hhh++)
					{
						triple.push_back(Eigen::Triplet<double>(nV + i, adjust_l0[hhh], (0.25 + std::cos(2 * 3.141592653 * hhh / det_v1) + 0.5 * std::cos(4 * 3.141592653 * hhh / det_v1)) / det_v1));
						
						ccl += (0.25 + std::cos(2 * 3.141592653 * hhh / det_v1) + 0.5 * std::cos(4 * 3.141592653 * hhh / det_v1)) / det_v1;
					}
					triple.push_back(Eigen::Triplet<double>(nV + i, v1.idx(), 1 - ccl));
				}
			}

			else
			{
				std::vector<int> adjust_l1 = adjust_seq(l1, v0.idx());

				if (det_v1 == 3)
				{
					/*triple.push_back(Eigen::Triplet<double>(nV + i, adjust_l1[0], 5.0 / 24));
					triple.push_back(Eigen::Triplet<double>(nV + i, adjust_l1[1], -1.0 / 24));
					triple.push_back(Eigen::Triplet<double>(nV + i, adjust_l1[2], -1.0 / 24));*/
					index.push_back(adjust_l1[0]);
					coef.push_back(5.0 / 24);
					index.push_back(adjust_l1[1]);
					coef.push_back(-1.0 / 24);
					index.push_back(adjust_l1[2]);
					coef.push_back(-1.0 / 24);
					index.push_back(v1.idx());
					coef.push_back(3.0/8);
				}
				if (det_v1 == 4)
				{
					/*triple.push_back(Eigen::Triplet<double>(nV + i, adjust_l1[0], 3.0 / 16));
					triple.push_back(Eigen::Triplet<double>(nV + i, adjust_l1[2], -1.0 / 16));*/
					index.push_back(adjust_l1[0]);
					coef.push_back(3.0 / 16);
					index.push_back(adjust_l1[2]);
					coef.push_back(-1.0 / 16);
					index.push_back(v1.idx());
					coef.push_back(3.0 / 8);
				}
				if (det_v1 >= 5)
				{
					double ccl = 0;
					for (int hhh = 0; hhh < det_v1; hhh++)
					{
						//triple.push_back(Eigen::Triplet<double>(nV + i, adjust_l1[hhh], (0.25 + std::cos(2 * 3.141592653 * hhh / det_v1) + 0.5 * std::cos(4 * 3.141592653 * hhh / det_v1)) / (2 * det_v1)));
						index.push_back(adjust_l1[hhh]);
						coef.push_back((0.25 + std::cos(2 * 3.141592653 * hhh / det_v1) + 0.5 * std::cos(4 * 3.141592653 * hhh / det_v1)) / (2 * det_v1));
						ccl += (0.25 + std::cos(2 * 3.141592653 * hhh / det_v1) + 0.5 * std::cos(4 * 3.141592653 * hhh / det_v1)) / (2 * det_v1);
					}
					index.push_back(v1.idx());
					coef.push_back(0.5 - ccl);
				}

				std::vector<int> adjust_l0 = adjust_seq(l0, v1.idx());
				if (det_v0 == 3)
				{
					if (isIn(index, adjust_l0[0]) > -1)
					{
						coef[isIn(index, adjust_l0[0])] += 5.0 / 24;
					}
					else
					{
						index.push_back(adjust_l0[0]);
						coef.push_back(5.0 / 24);
					}

					if (isIn(index, adjust_l0[1]) > -1)
					{
						coef[isIn(index, adjust_l0[1])] += -1.0 / 24;
					}
					else
					{
						index.push_back(adjust_l0[1]);
						coef.push_back(-1.0 / 24);
					}

					if (isIn(index, adjust_l0[2]) > -1)
					{
						coef[isIn(index, adjust_l0[2])] += -1.0 / 24;
					}
					else
					{
						index.push_back(adjust_l0[2]);
						coef.push_back(-1.0 / 24);
					}

					if (isIn(index, v0.idx()) > -1)
					{
						coef[isIn(index, v0.idx())] += 3.0 / 8;
					}
					else
					{
						index.push_back(v0.idx());
						coef.push_back(3.0/8);
					}

				}
				if (det_v0 == 4)
				{
					if (isIn(index, adjust_l0[0]) > -1)
					{
						coef[isIn(index, adjust_l0[0])] += 3.0 / 16;
					}
					else
					{
						index.push_back(adjust_l0[0]);
						coef.push_back(3.0 / 16);
					}

					if (isIn(index, adjust_l0[2]) > -1)
					{
						coef[isIn(index, adjust_l0[2])] += (-1.0 / 16);
					}
					else
					{
						index.push_back(adjust_l0[2]);
						coef.push_back(-1.0 / 16);
					}

					if (isIn(index, v0.idx()) > -1)
					{
						coef[isIn(index, v0.idx())] += 3.0 / 8;
					}
					else
					{
						index.push_back(v0.idx());
						coef.push_back(3.0 / 8);
					}

				}
				if (det_v0 >= 5)
				{
					double ccl= 0;
					for (int hhh = 0; hhh < det_v0; hhh++)
					{
						if (isIn(index, adjust_l0[hhh]) > -1)
						{
							ccl += (0.25 + std::cos(2 * 3.141592653 * hhh / det_v0) + 0.5 * std::cos(4 * 3.141592653 * hhh / det_v1)) / (2 * det_v0);
							coef[isIn(index, adjust_l0[hhh])] += (0.25 + std::cos(2 * 3.141592653 * hhh / det_v0) + 0.5 * std::cos(4 * 3.141592653 * hhh / det_v1)) / (2 * det_v0);
						}
						else
						{
							index.push_back(adjust_l0[hhh]);
							ccl += (0.25 + std::cos(2 * 3.141592653 * hhh / det_v0) + 0.5 * std::cos(4 * 3.141592653 * hhh / det_v1)) / (2 * det_v0);
							coef.push_back((0.25 + std::cos(2 * 3.141592653 * hhh / det_v0) + 0.5 * std::cos(4 * 3.141592653 * hhh / det_v1)) / (2 * det_v0));
						}

					}
					if (isIn(index, v0.idx()) > -1)
					{
						coef[isIn(index, v0.idx())] += (0.5 - ccl);
					}
					else
					{
						index.push_back(v0.idx());
						coef.push_back(0.5 - ccl);
					}
				}

				for (int hhh = 0; hhh < index.size(); hhh++)
				{
					triple.push_back(Eigen::Triplet<double>(nV + i, index[hhh], coef[hhh]));
					//subdivisionCoeff.insert(nV + i, index[hhh]) = coef[hhh];
				}

				index.clear();
				coef.clear();
			}
		}
	}

	ofstream file1("coef.txt");
	for (auto x = triple.begin(); x != triple.end(); ++x)
	{
		//auto y = &x;
		file1 << x->row() << "  "<< x->col() << "   "<< x->value() << endl;
	}
	file1.close();
	subdivisionCoeff.setFromTriplets(triple.begin(), triple.end());
	cout << "ok" << endl;

	
}

std::vector<int> AdaptiveSubdivision::det(PGMesh::VertexHandle v_)
{
	std::vector<int> l;
	auto v = v_;
	for (auto vvit = mesh_->vv_begin(v); vvit != mesh_->vv_end(v); ++vvit)
	{
		l.push_back(vvit->idx());
	}
	return l;
}

std::vector<int> AdaptiveSubdivision::adjust_seq(std::vector<int> seq0, int h)
{
	std::vector<int> seq(seq0.size());
	int start_idx = 0;
	int i = 0;
	while (i < seq0.size())
	{
		if (seq0[i] == h)
		{
			start_idx = i;
			break;
		}
		i++;
	}
	int j = start_idx;
	int x = start_idx;
	int l;
	while (j < seq0.size())
	{
		l = j - x;
		seq[l] = seq0[j];
		j = j + 1;
	}
	j = 0;
	while ( j < x)
	{
		if (j == x)
		{
			break;
		}
		seq[seq0.size() - x + j] = seq0[j];
		j++;
	}
	/*for (int xxx = 0; xxx < seq0.size(); xxx++)
	{
		cout << seq0[xxx] << "   ";
	}
	
	cout << endl << h << endl;
	for (int xxx = 0; xxx < seq0.size(); xxx++)
	{
		cout << seq[xxx] << "   ";
	}
	cout << endl;*/
	return seq;
}

int AdaptiveSubdivision::isIn(std::vector<int> seq, int h)
{
	for (int i = 0; i < seq.size(); i++)
	{
		if (seq[i] == h)
		{
			return i;
		}
	}
	return -1;
}

Eigen::MatrixXd AdaptiveSubdivision::getVertexMatrix(PGMesh * _mesh)
{
	Eigen::MatrixXd V(nV, 3);
	for (auto vit = _mesh->vertices_begin(); vit != _mesh->vertices_end(); ++vit)
	{
		PGMesh::Point p = _mesh->point(*vit);
		V(vit->idx(), 0) = p[0];
		V(vit->idx(), 1) = p[1];
		V(vit->idx(), 2) = p[2];
	}
	return V;
}

Eigen::VectorXd AdaptiveSubdivision::getVertexVec(PGMesh *_mesh)
{
	Eigen::VectorXd V(nV*3);
	for (auto vit = _mesh->vertices_begin(); vit != _mesh->vertices_end(); ++vit)
	{
		PGMesh::Point p = _mesh->point(*vit);
		V(3 * vit->idx()+ 0) = p[0];
		V(3 * vit->idx() + 1) = p[1];
		V(3 * vit->idx()+ 2) = p[2];
	}
	return V;
}

Eigen::MatrixXd AdaptiveSubdivision::getSubdivisionPosition(Eigen::MatrixXd V0)
{
	return subdivisionCoeff*V0; 
}

void AdaptiveSubdivision::outMesh(Eigen::MatrixXd cordComponent, string path)
{
	ofstream offfile(path);
	//offfile << "COFF"<<endl;
	//offfile << nV << " " << nF<< " "<< 0 <<endl;

	for (int i = 0; i< nV + n_affected_edges; i++)
	{
		offfile << "v" << " " << cordComponent(i , 0) << " " << cordComponent(i,  1) << " " << cordComponent(i , 2) << endl;
	}
	for (PGMesh::FaceIter fit = subdividedMesh->faces_begin(); fit != subdividedMesh->faces_end(); ++fit)
	{
		offfile << "f" << " ";
		for (PGMesh::FVIter fvit = subdividedMesh->fv_begin(*fit); fvit != subdividedMesh->fv_end(*fit); ++fvit)
		{
			offfile << fvit->idx() + 1 << " ";
		}
		offfile << endl;
		/*PGMesh::FVIter fvit = mesh_->fv_begin(fit.handle());
		PGMesh::VertexHandle vv0 = fvit.handle();
		++fvit;
		PGMesh::VertexHandle vv1 = fvit.handle();
		++fvit;
		PGMesh::VertexHandle vv2 = fvit.handle();
		offfile << vv2.idx() + 1 << " "<< vv1.idx() + 1 << " " << vv0.idx() + 1 << " "  ;
		offfile << endl;*/
	}
	offfile.close();
}