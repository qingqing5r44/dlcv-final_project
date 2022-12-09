#include"SparsePCAForShapeGrad.h"

using namespace std;
using namespace Eigen;


		unsigned char colormap6[300] = 
		{
			153, 	153, 	153, 	
			147, 	157, 	157, 	
			142, 	161, 	161, 	
			136, 	164, 	164, 	
			130, 	168, 	168, 	
			125, 	172, 	172, 	
			119, 	176, 	176, 	
			113, 	179, 	179, 	
			108, 	183, 	183, 	
			102, 	187, 	187, 	
			96, 	191, 	191, 	
			91, 	195, 	195, 	
			85, 	198, 	198, 	
			79, 	202, 	202, 	
			74, 	206, 	206, 	
			68, 	210, 	210, 	
			62, 	213, 	213, 	
			57, 	217, 	217, 	
			51, 	221, 	221, 	
			45, 	225, 	225, 	
			40, 	229, 	229, 	
			34, 	232, 	232, 	
			28, 	236, 	236, 	
			23, 	240, 	240, 	
			17, 	244, 	244, 	
			11, 	247, 	247, 	
			6, 	251, 	251, 	
			0, 	255, 	255, 	
			0, 	251, 	255, 	
			0, 	248, 	255, 	
			0, 	244, 	255, 	
			0, 	241, 	255, 	
			0, 	237, 	255, 	
			0, 	234, 	255, 	
			0, 	230, 	255, 	
			0, 	227, 	255, 	
			0, 	223, 	255, 	
			0, 	220, 	255, 	
			0, 	216, 	255, 	
			0, 	213, 	255, 	
			0, 	209, 	255, 	
			0, 	206, 	255, 	
			0, 	202, 	255, 	
			0, 	199, 	255, 	
			0, 	195, 	255, 	
			0, 	192, 	255, 	
			0, 	188, 	255, 	
			0, 	185, 	255, 	
			0, 	181, 	255, 	
			0, 	178, 	255, 	
			0, 	174, 	255, 	
			0, 	171, 	255, 	
			0, 	167, 	255, 	
			0, 	164, 	255, 	
			0, 	160, 	255, 	
			0, 	157, 	255, 	
			0, 	153, 	255, 	
			0, 	151, 	255, 	
			0, 	148, 	255, 	
			0, 	146, 	255, 	
			0, 	144, 	255, 	
			0, 	141, 	255, 	
			0, 	139, 	255, 	
			0, 	136, 	255, 	
			0, 	134, 	255, 	
			0, 	132, 	255, 	
			0, 	129, 	255, 	
			0, 	127, 	255, 	
			0, 	125, 	255, 	
			0, 	122, 	255, 	
			0, 	120, 	255, 	
			0, 	117, 	255, 	
			0, 	115, 	255, 	
			0, 	113, 	255, 	
			0, 	110, 	255, 	
			0, 	108, 	255, 	
			0, 	106, 	255, 	
			0, 	103, 	255, 	
			0, 	101, 	255, 	
			0, 	98, 	255, 	
			0, 	96, 	255, 	
			0, 	94, 	255, 	
			0, 	91, 	255, 	
			0, 	89, 	255, 	
			0, 	87, 	255, 	
			0, 	84, 	255, 	
			0, 	82, 	255, 	
			0, 	79, 	255, 	
			0, 	77, 	255, 	
			0, 	75, 	255, 	
			0, 	72, 	255, 	
			0, 	70, 	255, 	
			0, 	68, 	255, 	
			0, 	65, 	255, 	
			0, 	63, 	255, 	
			0, 	60, 	255, 	
			0, 	58, 	255, 	
			0, 	56, 	255, 	
			0, 	53, 	255, 	
			0, 	51, 	255, 	

		};
		

		int SparsePCAForShapeGrad::fixFacesNum = 10;

		bool SparsePCAForShapeGrad::setDeformedMesh(std::vector<Eigen::VectorXd> & sd)
		{
			if (sd.size() != nV)
			{
				return false;
			}
#pragma omp parallel for
			for (int i = 0; i < nV; i++)
			{
				SpatialTemData[i][1] = sd[i];
			}
			return true;
		}

		void SparsePCAForShapeGrad::setDeformedMesh(PGMesh *_mesh)
		{
#pragma omp parallel for
			for (int i = 0; i < nV; i++)
			{
				SpatialTemData[i].resize(2);
				PGMesh::VertexHandle vh = _mesh->vertex_handle(i);
				PGMesh::Point &p = _mesh->point(vh);
				Eigen::VectorXd x(3);
				x << p[0], p[1], p[2];
				SpatialTemData[i][1] = x;
			}
		}

		void SparsePCAForShapeGrad::computeDG()
		{
#pragma omp parallel for
			for (int i = 0; i < nF; ++i){
					/*找到每个三角形在第0帧和第j帧的3个顶点坐标v0, v1, v2, v0t, v1t, v2t, 计算V和Vt*/
					/*
					v = [v1 - v0, v2 - v0]
					v' = [v1' - v0', v2' - v0']
					*/
					Eigen::VectorXd v0;
					Eigen::VectorXd v1;
					Eigen::VectorXd v2;
					Eigen::VectorXd v0t;
					Eigen::VectorXd v1t;
					Eigen::VectorXd v2t;
					Eigen::MatrixXd V(3, 3);		//3*2的面矩阵V和V‘
					Eigen::MatrixXd Vt(3, 3);

					PGMesh::FaceHandle fh = mesh_->face_handle(i);
					PGMesh::FVIter fv_it = mesh_->fv_begin(fh);
					int v0No = fv_it->idx();
					v0 = SpatialTemData[v0No][0];
					v0t = SpatialTemData[v0No][1];
					++fv_it;
					int v1No = fv_it->idx();
					v1 = SpatialTemData[v1No][0];
					v1t = SpatialTemData[v1No][1];
					++fv_it;
					int v2No = fv_it->idx();
					v2 = SpatialTemData[v2No][0];
					v2t = SpatialTemData[v2No][1];
					Eigen::Vector3d v1v0(3);
					Eigen::Vector3d v2v0(3);
					Eigen::Vector3d v1v0t(3);
					Eigen::Vector3d v2v0t(3);
					v1v0 = v1 - v0;
					v2v0 = v2 - v0;
					v1v0t = v1t - v0t;
					v2v0t = v2t - v0t;
					Eigen::Vector3d N0 = v1v0.cross(v2v0);
					Eigen::Vector3d N1 = v1v0t.cross(v2v0t);
					V.col(0) = v1v0;
					V.col(1) = v2v0;
					V.col(2) = N0;
					Vt.col(0) = v1v0t;
					Vt.col(1) = v2v0t;
					Vt.col(2) = N1;
					//对V进行QR分解
					/*HouseholderQR< MatrixXd > qr;
					qr.compute(V);
					MatrixXd R_tmp = qr.matrixQR().triangularView<Upper>();
					MatrixXd Q_tmp = qr.householderQ();
					MatrixXd R = R_tmp.block(0, 0, 2, 2);
					MatrixXd Q = Q_tmp.block(0, 0, 3, 2);
					MatrixXd L = R.inverse() * Q.transpose();
					MatrixXd J = Vt * L;*/
					V_inv[i] = V.inverse();
					MatrixXd J = Vt * V_inv[i];
					DGs[i] = J;
					//用JacobiSVD进行极坐标分解
					//U = u * v' P = v * s * v'
					JacobiSVD<MatrixXd> svd(J, ComputeThinU | ComputeThinV);
					MatrixXd s(3, 3);
					s.setZero();
					s.diagonal() = svd.singularValues();
					MatrixXd u = svd.matrixU();
					MatrixXd v = svd.matrixV();
					MatrixXd U = u * v.transpose();
					MatrixXd P = v * s * v.transpose();
					Rs[i] = U;
					Ss[i] = P;
					MatrixXd log_R_i = U.log();
					/*op << r_vec[0]<<" "<< r_vec[1]<<" "<< r_vec[2] << endl;
					op << endl;*/
					VectorXd cm(9);	//单个面的Deformation gradient
					cm << log_R_i(0, 1), log_R_i(0, 2), log_R_i(1, 2),
						P(0, 0) - 1, P(1, 0), P(1, 1) - 1, P(2, 0), P(2, 1), P(2, 2) - 1;
					DGs[i] = cm;
					/*CMComp.segment(9 * i + 0, 9) = cm;
					R_ijs[i] = R_ij;*/

				}
		}

		void SparsePCAForShapeGrad::computeCM()
		{
			
			Eigen::VectorXd CMComp(9 * nE);
			R_ijs.resize(nE);
#pragma omp parallel for
			for (int i = 0; i < nE; i++)
			{
				PGMesh::EdgeHandle eh = mesh_->edge_handle(i);
				PGMesh::HalfedgeHandle heh0 = mesh_->halfedge_handle(eh, 0);
				PGMesh::HalfedgeHandle heh1 = mesh_->halfedge_handle(eh, 1);
				PGMesh::FaceHandle fh0 = mesh_->face_handle(heh0);
				PGMesh::FaceHandle fh1 = mesh_->face_handle(heh1);
				int from;
				int to;
				if (!fh1.is_valid())
				{
					from = fh0.idx();
					to = from;
				}
				else if (!fh0.is_valid())
				{
					to = fh1.idx();
					from = to;
				}
				else
				{
					from = fh0.idx();
					to = fh1.idx();
				}
				Eigen::MatrixXd fromR = Rs[from];
				Eigen::MatrixXd toR = Rs[to];
				Eigen::MatrixXd fromS = Ss[from];
				Eigen::MatrixXd toS = Ss[to];
				Eigen::MatrixXd R_ij = toR * fromR.transpose();
				//HouseholderQR< MatrixXd > qr;
				//qr.compute(fromR);
				//MatrixXd R_tmp = qr.matrixQR().triangularView<Upper>();
				//MatrixXd Q_tmp = qr.householderQ();
				//MatrixXd R = R_tmp.block(0, 0, 3, 3);
				//MatrixXd Q = Q_tmp.block(0, 0, 3, 3);
				//MatrixXd L = R.inverse() * Q.transpose();
				//MatrixXd R_ij = toR * L;
				Eigen::MatrixXd S_ij = (fromS +toS) / 2;

				
				/*double R_matrix[9] = { R_ij(0,0), R_ij(0,1) , R_ij(0,2) , R_ij(1,0) , R_ij(1,1) , R_ij(1,2) , R_ij(2,0) , R_ij(2,1) , R_ij(2,2) };
				double r_vec[3];
				CvMat pr_vec;
				CvMat pR_matrix;
				cvInitMatHeader(&pr_vec, 1, 3, CV_64FC1, r_vec, CV_AUTOSTEP);
				cvInitMatHeader(&pR_matrix, 3, 3, CV_64FC1, R_matrix, CV_AUTOSTEP);

				cvRodrigues2(&pR_matrix, &pr_vec, 0);*/
				
				MatrixXd log_R_ij = R_ij.log();
				/*op << r_vec[0]<<" "<< r_vec[1]<<" "<< r_vec[2] << endl;
				op << endl;*/
				VectorXd cm(9);	//单个面的Deformation gradient
				cm << log_R_ij(0,1), log_R_ij(0,2), log_R_ij(1,2),
					S_ij(0, 0) - 1, S_ij(1, 0), S_ij(1, 1) - 1, S_ij(2, 0), S_ij(2, 1), S_ij(2, 2) - 1;
				CMs[i] = cm;
				CMComp.segment(9 * i + 0, 9) = cm;
				R_ijs[i] = R_ij;
			}
			VecCM = CMComp;
		}
		
		void SparsePCAForShapeGrad::writeCM(char* path)
		{
			fstream op(path, ios::out);
			op << VecCM << endl;
			op.close();
		}

		void SparsePCAForShapeGrad::build_matrix()
		{
			//const int nF = _mesh->n_faces();
			//const int nV = 9;			//每个向量中的元素个数，这里为9；
			const int vector_length = 9;	//每个向量中的元素个数，这里为9；


			//用来存储每个面上各个节点的编号
			struct Node_idx{
				int v0_idx;
				int v1_idx;
				int v2_idx;
			};

			//直接存V1-V0的值，2个为一组，代表一个面的V
			double *x = new double[2 * nF];
			double *y = new double[2 * nF];
			double *z = new double[2 * nF];
			Node_idx* node_idx = new Node_idx[nF];		//每个面上对应三个点对应一个node_idx用于存储各个节点的编号

			int vCounter = 0;
			int tmp_faceCount = 0;


			//遍历每个面，并读取每个面上的所有点的位置信息
			int fNumber = 0;
			for (PGMesh::FaceIter f_it = mesh_->faces_begin(); f_it != mesh_->faces_end(); ++f_it, ++fNumber){
				double x_tmp, y_tmp, z_tmp;
				int i = 0;
				tmp_faceCount++;
				for (PGMesh::FaceVertexIter v_it = mesh_->fv_begin(*f_it); v_it != mesh_->fv_end(*f_it); ++v_it, ++i){
					auto point = mesh_->point(*v_it);
					if (i == 0){
						x_tmp = point.data()[0];
						y_tmp = point.data()[1];
						z_tmp = point.data()[2];
						node_idx[fNumber].v0_idx = v_it->idx();
					}
					else{
						x[vCounter] = point.data()[0] - x_tmp;
						y[vCounter] = point.data()[1] - y_tmp;
						z[vCounter] = point.data()[2] - z_tmp;
						if (i == 1)
							node_idx[fNumber].v1_idx = v_it->idx();
						else node_idx[fNumber].v2_idx = v_it->idx();
						++vCounter;
					}
				}
			}
			
			Eigen::SparseMatrix<double> SCoef(9 * nF + 3 * fixPointIdx.size(), 3 * nV);
			std::vector<Eigen::Triplet<double> > triple;
//#pragma omp parallel for
			for (int i = 0; i < nF; ++i){
				Matrix< double, 3, 2 > V_input;
				HouseholderQR< Matrix< double, 3, 2 > > QR;
				int v0_start = node_idx[i].v0_idx * 3;
				int v1_start = node_idx[i].v1_idx * 3;
				int v2_start = node_idx[i].v2_idx * 3;
				int startLine = i * vector_length;

				V_input(0, 0) = x[i * 2];
				V_input(0, 1) = x[i * 2 + 1];
				V_input(1, 0) = y[i * 2];
				V_input(1, 1) = y[i * 2 + 1];
				V_input(2, 0) = z[i * 2];
				V_input(2, 1) = z[i * 2 + 1];
				QR.compute(V_input);
				Matrix< double, 3, 2 > R_tmp = QR.matrixQR().triangularView< Upper>();
				Matrix3d Q_tmp = QR.householderQ();
				Matrix< double, 3, 2 > Q = Q_tmp.block(0, 0, 3, 2);
				Matrix2d R = R_tmp.block(0, 0, 2, 2);
				Matrix< double, 2, 3 > tmp = R.inverse() * Q.transpose();
				double a11 = tmp(0, 0), a21 = tmp(1, 0), a12 = tmp(0, 1),
					a22 = tmp(1, 1), a13 = tmp(0, 2), a23 = tmp(1, 2);
				double tmp_1 = 0.0 - a11 - a21;
				double tmp_2 = 0.0 - a12 - a22;
				double tmp_3 = 0.0 - a13 - a23;

				//v0
				triple.push_back(Eigen::Triplet<double>(startLine, v0_start, tmp_1));
				triple.push_back(Eigen::Triplet<double>(startLine + 1, v0_start, tmp_2));
				triple.push_back(Eigen::Triplet<double>(startLine + 2, v0_start++, tmp_3));

				triple.push_back(Eigen::Triplet<double>(startLine + 3, v0_start, tmp_1));
				triple.push_back(Eigen::Triplet<double>(startLine + 4, v0_start, tmp_2));
				triple.push_back(Eigen::Triplet<double>(startLine + 5, v0_start++, tmp_3));

				triple.push_back(Eigen::Triplet<double>(startLine + 6, v0_start, tmp_1));
				triple.push_back(Eigen::Triplet<double>(startLine + 7, v0_start, tmp_2));
				triple.push_back(Eigen::Triplet<double>(startLine + 8, v0_start, tmp_3));

				//v1
				triple.push_back(Eigen::Triplet<double>(startLine, v1_start, a11));
				triple.push_back(Eigen::Triplet<double>(startLine + 1, v1_start, a12));
				triple.push_back(Eigen::Triplet<double>(startLine + 2, v1_start++, a13));

				triple.push_back(Eigen::Triplet<double>(startLine + 3, v1_start, a11));
				triple.push_back(Eigen::Triplet<double>(startLine + 4, v1_start, a12));
				triple.push_back(Eigen::Triplet<double>(startLine + 5, v1_start++, a13));

				triple.push_back(Eigen::Triplet<double>(startLine + 6, v1_start, a11));
				triple.push_back(Eigen::Triplet<double>(startLine + 7, v1_start, a12));
				triple.push_back(Eigen::Triplet<double>(startLine + 8, v1_start, a13));

				//v2
				triple.push_back(Eigen::Triplet<double>(startLine, v2_start, a21));
				triple.push_back(Eigen::Triplet<double>(startLine + 1, v2_start, a22));
				triple.push_back(Eigen::Triplet<double>(startLine + 2, v2_start++, a23));

				triple.push_back(Eigen::Triplet<double>(startLine + 3, v2_start, a21));
				triple.push_back(Eigen::Triplet<double>(startLine + 4, v2_start, a22));
				triple.push_back(Eigen::Triplet<double>(startLine + 5, v2_start++, a23));

				triple.push_back(Eigen::Triplet<double>(startLine + 6, v2_start, a21));
				triple.push_back(Eigen::Triplet<double>(startLine + 7, v2_start, a22));
				triple.push_back(Eigen::Triplet<double>(startLine + 8, v2_start, a23));
			}
			for (int i = 0; i < fixPointIdx.size(); i++)
			{
				triple.push_back(Eigen::Triplet<double>(9 * nF + 3 * i + 0, 3 * fixPointIdx[i] + 0, 1));
				triple.push_back(Eigen::Triplet<double>(9 * nF + 3 * i + 1, 3 * fixPointIdx[i] + 1, 1));
				triple.push_back(Eigen::Triplet<double>(9 * nF + 3 * i + 2, 3 * fixPointIdx[i] + 2, 1));
			}
			SCoef.setFromTriplets(triple.begin(), triple.end());
			Vsolver.compute(SCoef.transpose() * SCoef);
			DGCoef = SCoef;
		}

		void SparsePCAForShapeGrad::computeS_matrix()
		{
			SSCoef = Eigen::SparseMatrix<double> (6 * nE, 6 * nF);
			std::vector<Eigen::Triplet<double> > triple;
//#pragma omp parallel for
			for (int i = 0; i < nE; i++)
			{
				PGMesh::EdgeHandle eh = mesh_->edge_handle(i);
				PGMesh::HalfedgeHandle heh0 = mesh_->halfedge_handle(eh, 0);
				PGMesh::HalfedgeHandle heh1 = mesh_->halfedge_handle(eh, 1);
				PGMesh::FaceHandle fh0 = mesh_->face_handle(heh0);
				PGMesh::FaceHandle fh1 = mesh_->face_handle(heh1);
				int from;
				int to;
				if (!fh1.is_valid())
				{
					from = fh0.idx();
					to = from;
				}
				else if (!fh0.is_valid())
				{
					to = fh1.idx();
					from = to;
				}
				else
				{
					from = fh0.idx();
					to = fh1.idx();
				}
				for (int j = 0; j < 6; j++)
				{
					if (from != to)
					{
						triple.push_back(Eigen::Triplet<double>(6 * i + j, 6 * from + j, 0.5));
						triple.push_back(Eigen::Triplet<double>(6 * i + j, 6 * to + j, 0.5));
					}
					else
					{
						triple.push_back(Eigen::Triplet<double>(6 * i + j, 6 * to + j, 1));
					}
				}
			}
			SSCoef.setFromTriplets(triple.begin(), triple.end());
			Ssolver.compute(SSCoef.transpose() * SSCoef);
		}

		void SparsePCAForShapeGrad::compute_V_inv()
		{
			for (int i = 0; i < nF; ++i){
				/*找到每个三角形在第0帧和第j帧的3个顶点坐标v0, v1, v2, v0t, v1t, v2t, 计算V和Vt*/
				/*
				v = [v1 - v0, v2 - v0]
				v' = [v1' - v0', v2' - v0']
				*/
				Eigen::VectorXd v0;
				Eigen::VectorXd v1;
				Eigen::VectorXd v2;
				Eigen::VectorXd v0t;
				Eigen::VectorXd v1t;
				Eigen::VectorXd v2t;
				Eigen::MatrixXd V(3, 2);		//3*2的面矩阵V和V‘

				PGMesh::FaceHandle fh = mesh_->face_handle(i);
				PGMesh::FVIter fv_it = mesh_->fv_begin(fh);
				int v0No = fv_it->idx();
				v0 = SpatialTemData[v0No][0];
				++fv_it;
				int v1No = fv_it->idx();
				v1 = SpatialTemData[v1No][0];
				++fv_it;
				int v2No = fv_it->idx();
				v2 = SpatialTemData[v2No][0];
				Eigen::VectorXd v1v0(3);
				Eigen::VectorXd v2v0(3);
				v1v0 = v1 - v0;
				v2v0 = v2 - v0;

				V.col(0) = v1v0;
				V.col(1) = v2v0;
				//对V进行QR分解
				HouseholderQR< MatrixXd > qr;
				qr.compute(V);
				MatrixXd R_tmp = qr.matrixQR().triangularView<Upper>();
				MatrixXd Q_tmp = qr.householderQ();
				MatrixXd R = R_tmp.block(0, 0, 2, 2);
				MatrixXd Q = Q_tmp.block(0, 0, 3, 2);
				MatrixXd L = R.inverse() * Q.transpose();

			}
		}

		void SparsePCAForShapeGrad::vecToMat(Eigen::VectorXd CMComponent, std::vector<Eigen::MatrixXd> &R_ij, std::vector<Eigen::MatrixXd> &S_i)
		{
			R_ij.resize(nE);
			S_i.resize(nF);
			Eigen::VectorXd y(6 * nE);
			Eigen::VectorXd x;
#pragma omp parallel for
			for (int i = 0; i < nE; i++)
			{
				MatrixXd U(3, 3);			// counterpart of U;
				/*double R_matrix[9];
				double r_vec[3] = { CMComponent(9 * i + 0), CMComponent(9 * i + 1), CMComponent(9 * i + 2) };
				CvMat pr_vec;
				CvMat pR_matrix;
				cvInitMatHeader(&pr_vec, 1, 3, CV_64FC1, r_vec, CV_AUTOSTEP);
				cvInitMatHeader(&pR_matrix, 3, 3, CV_64FC1, R_matrix, CV_AUTOSTEP);

				cvRodrigues2(&pr_vec, &pR_matrix, 0);
				U << R_matrix[0], R_matrix[1], R_matrix[2],
					R_matrix[3], R_matrix[4], R_matrix[5],
					R_matrix[6], R_matrix[7], R_matrix[8];
				op << "the " << i << "-tyh R_ij:" << endl;*/
				MatrixXd U_(3, 3);
				U_.setZero();
				U_(0, 1) = CMComponent(9 * i + 0);
				U_(0, 2) = CMComponent(9 * i + 1);
				U_(1, 2) = CMComponent(9 * i + 2);
				U_(1, 0) = -U_(0, 1);
				U_(2, 0) = -U_(0, 2);
				U_(2, 1) = -U_(1, 2);
				U = U_.exp();
				//addend
				
				/*op << endl;
				op << r_vec[0]<<" "<< r_vec[1]<<" "<< r_vec[2] << endl;*/
				R_ij[i] = U;
				y.segment(6 * i, 6) = CMComponent.segment(9 * i + 3, 6);
			}
			x = Ssolver.solve(SSCoef.transpose() * y);
#pragma omp parallel for
			for (int i = 0; i < nF; i++)
			{
				MatrixXd U_(3, 3);
				U_ << x(6 * i + 0) + 1, 0, 0,
					x(6 * i + 1), x(6 * i + 2) + 1, 0,
					x(6 * i + 3), x(6 * i + 4), x(6 * i + 5 ) + 1;
				U_(0, 1) = U_(1, 0);
				U_(0, 2) = U_(2, 0);
				U_(1, 2) = U_(2, 1);
				S_i[i] = U_;
			}
		}

		Eigen::VectorXd SparsePCAForShapeGrad::reconstructFromRS(PGMesh * _mesh, std::vector<Eigen::MatrixXd> R_ij, std::vector<Eigen::MatrixXd> S_i)
		{
			Eigen::SparseMatrix<double> RCoef(9 * nE + 9 * fixFaceIdx.size(), 9 * nF);
			RCoef.setZero();
			std::vector<Eigen::Triplet<double> > triple;
			Eigen::VectorXd y(9 * nE + 9 * fixFaceIdx.size());
			y.setZero();
//#pragma omp parallel for
			for (int i = 0; i < nE; i++)
			{
				int edgeNo = i;
				MatrixXd rot = R_ij[i];
				PGMesh::EdgeHandle eh = mesh_->edge_handle(i);
				PGMesh::HalfedgeHandle heh0 = mesh_->halfedge_handle(eh, 0);
				PGMesh::HalfedgeHandle heh1 = mesh_->halfedge_handle(eh, 1);
				PGMesh::FaceHandle fh0 = mesh_->face_handle(heh0);
				PGMesh::FaceHandle fh1 = mesh_->face_handle(heh1);
				int fromFaceNo;
				int toFaceNo;
				if (!fh1.is_valid())
				{
					fromFaceNo = fh0.idx();
					toFaceNo = fromFaceNo;
				}
				else if (!fh0.is_valid())
				{
					toFaceNo = fh1.idx();
					fromFaceNo = toFaceNo;
				}
				else
				{
					fromFaceNo = fh0.idx();
					toFaceNo = fh1.idx();
				}
				if (fromFaceNo != toFaceNo)
				{
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 0, 9 * fromFaceNo + 0, rot(0, 0)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 0, 9 * fromFaceNo + 3, rot(0, 1)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 0, 9 * fromFaceNo + 6, rot(0, 2)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 0, 9 * toFaceNo + 0, -1));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 1, 9 * fromFaceNo + 1, rot(0, 0)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 1, 9 * fromFaceNo + 4, rot(0, 1)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 1, 9 * fromFaceNo + 7, rot(0, 2)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 1, 9 * toFaceNo + 1, -1));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 2, 9 * fromFaceNo + 2, rot(0, 0)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 2, 9 * fromFaceNo + 5, rot(0, 1)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 2, 9 * fromFaceNo + 8, rot(0, 2)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 2, 9 * toFaceNo + 2, -1));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 3, 9 * fromFaceNo + 0, rot(1, 0)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 3, 9 * fromFaceNo + 3, rot(1, 1)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 3, 9 * fromFaceNo + 6, rot(1, 2)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 3, 9 * toFaceNo + 3, -1));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 4, 9 * fromFaceNo + 1, rot(1, 0)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 4, 9 * fromFaceNo + 4, rot(1, 1)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 4, 9 * fromFaceNo + 7, rot(1, 2)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 4, 9 * toFaceNo + 4, -1));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 5, 9 * fromFaceNo + 2, rot(1, 0)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 5, 9 * fromFaceNo + 5, rot(1, 1)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 5, 9 * fromFaceNo + 8, rot(1, 2)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 5, 9 * toFaceNo + 5, -1));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 6, 9 * fromFaceNo + 0, rot(2, 0)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 6, 9 * fromFaceNo + 3, rot(2, 1)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 6, 9 * fromFaceNo + 6, rot(2, 2)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 6, 9 * toFaceNo + 6, -1));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 7, 9 * fromFaceNo + 1, rot(2, 0)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 7, 9 * fromFaceNo + 4, rot(2, 1)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 7, 9 * fromFaceNo + 7, rot(2, 2)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 7, 9 * toFaceNo + 7, -1));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 8, 9 * fromFaceNo + 2, rot(2, 0)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 8, 9 * fromFaceNo + 5, rot(2, 1)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 8, 9 * fromFaceNo + 8, rot(2, 2)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 8, 9 * toFaceNo + 8, -1));
					/*triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 0, 9 * fromFaceNo + 0, rot(0, 0)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 0, 9 * fromFaceNo + 1, rot(1, 0)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 0, 9 * fromFaceNo + 2, rot(2, 0)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 0, 9 * toFaceNo + 0, -1));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 1, 9 * fromFaceNo + 3, rot(0, 0)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 1, 9 * fromFaceNo + 4, rot(1, 0)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 1, 9 * fromFaceNo + 5, rot(2, 0)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 1, 9 * toFaceNo + 3, -1));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 2, 9 * fromFaceNo + 6, rot(0, 0)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 2, 9 * fromFaceNo + 7, rot(1, 0)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 2, 9 * fromFaceNo + 8, rot(2, 0)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 2, 9 * toFaceNo + 6, -1));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 3, 9 * fromFaceNo + 0, rot(0, 1)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 3, 9 * fromFaceNo + 1, rot(1, 1)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 3, 9 * fromFaceNo + 2, rot(2, 1)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 3, 9 * toFaceNo + 1, -1));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 4, 9 * fromFaceNo + 3, rot(0, 1)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 4, 9 * fromFaceNo + 4, rot(1, 1)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 4, 9 * fromFaceNo + 5, rot(2, 1)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 4, 9 * toFaceNo + 4, -1));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 5, 9 * fromFaceNo + 6, rot(0, 1)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 5, 9 * fromFaceNo + 7, rot(1, 1)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 5, 9 * fromFaceNo + 8, rot(2, 1)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 5, 9 * toFaceNo + 7, -1));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 6, 9 * fromFaceNo + 0, rot(0, 2)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 6, 9 * fromFaceNo + 1, rot(1, 2)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 6, 9 * fromFaceNo + 2, rot(2, 2)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 6, 9 * toFaceNo + 2, -1));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 7, 9 * fromFaceNo + 3, rot(0, 2)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 7, 9 * fromFaceNo + 4, rot(1, 2)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 7, 9 * fromFaceNo + 5, rot(2, 2)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 7, 9 * toFaceNo + 5, -1));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 8, 9 * fromFaceNo + 6, rot(0, 2)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 8, 9 * fromFaceNo + 7, rot(1, 2)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 8, 9 * fromFaceNo + 8, rot(2, 2)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 8, 9 * toFaceNo + 8, -1));*/
				}
			}
//#pragma omp parallel for
			for (int i = 0; i < fixFaceIdx.size(); i++)
			{
				triple.push_back(Eigen::Triplet<double>(9 * nE + 9 * i + 0, 9 * fixFaceIdx[i] + 0, 1));
				triple.push_back(Eigen::Triplet<double>(9 * nE + 9 * i + 1, 9 * fixFaceIdx[i] + 1, 1));
				triple.push_back(Eigen::Triplet<double>(9 * nE + 9 * i + 2, 9 * fixFaceIdx[i] + 2, 1));
				triple.push_back(Eigen::Triplet<double>(9 * nE + 9 * i + 3, 9 * fixFaceIdx[i] + 3, 1));
				triple.push_back(Eigen::Triplet<double>(9 * nE + 9 * i + 4, 9 * fixFaceIdx[i] + 4, 1));
				triple.push_back(Eigen::Triplet<double>(9 * nE + 9 * i + 5, 9 * fixFaceIdx[i] + 5, 1));
				triple.push_back(Eigen::Triplet<double>(9 * nE + 9 * i + 6, 9 * fixFaceIdx[i] + 6, 1));
				triple.push_back(Eigen::Triplet<double>(9 * nE + 9 * i + 7, 9 * fixFaceIdx[i] + 7, 1));
				triple.push_back(Eigen::Triplet<double>(9 * nE + 9 * i + 8, 9 * fixFaceIdx[i] + 8, 1));
				y(9 * nE + 9 * i + 0) = fixFaceR[i](0, 0);
				y(9 * nE + 9 * i + 1) = fixFaceR[i](0, 1);
				y(9 * nE + 9 * i + 2) = fixFaceR[i](0, 2);
				y(9 * nE + 9 * i + 3) = fixFaceR[i](1, 0);
				y(9 * nE + 9 * i + 4) = fixFaceR[i](1, 1);
				y(9 * nE + 9 * i + 5) = fixFaceR[i](1, 2);
				y(9 * nE + 9 * i + 6) = fixFaceR[i](2, 0);
				y(9 * nE + 9 * i + 7) = fixFaceR[i](2, 1);
				y(9 * nE + 9 * i + 8) = fixFaceR[i](2, 2);
			}
			RCoef.setFromTriplets(triple.begin(), triple.end());
			Rsolver.compute(RCoef.transpose() * RCoef);
			Eigen::VectorXd x = Rsolver.solve(RCoef.transpose() * y);
			
			Eigen::VectorXd y_v(9 * nF + 3 * fixPointIdx.size());
//#pragma omp parallel for
			for (int faceNo = 0; faceNo < nF;  faceNo++)
			{
				Eigen::Matrix3d frameFace;
				frameFace << x(9 * faceNo + 0), x(9 * faceNo + 1), x(9 * faceNo + 2),
					x(9 * faceNo + 3), x(9 * faceNo + 4), x(9 * faceNo + 5),
					x(9 * faceNo + 6), x(9 * faceNo + 7), x(9 * faceNo + 8);
				Eigen::JacobiSVD<MatrixXd> svd(frameFace.transpose(), ComputeThinU | ComputeThinV);
				Eigen::MatrixXd u = svd.matrixU();
				Eigen::MatrixXd v = svd.matrixV();
				Eigen::Matrix3d U = (u * v.transpose()).transpose();
				Eigen::Matrix3d P = S_i[faceNo];
				Eigen::Matrix3d D = (U * P);

				y_v.segment(9 * faceNo + 0, 3) = D.row(0);
				y_v.segment(9 * faceNo + 3, 3) = D.row(1);
				y_v.segment(9 * faceNo + 6, 3) = D.row(2);
			}

			for (int i = 0; i < fixPointIdx.size(); i++)
			{
				y_v.segment(9 * nF + 3 * i + 0, 3) = fixPointPos[i];
			}
			Eigen::VectorXd vert = Vsolver.solve(DGCoef.transpose() * y_v);
			return vert;
		}

		SparsePCAForShapeGrad::SparsePCAForShapeGrad (PGMesh * _mesh, std::vector<std::vector<Eigen::VectorXd>> & sd)
		{
			mesh_ = _mesh;
			SpatialTemData = sd;
			nV = mesh_->n_vertices();
			nF = mesh_->n_faces();
			nE = mesh_->n_edges();
			DeformGrad.resize(nF);
			nSeq = (sd[0]).size();
			DFMatrix = Eigen::MatrixXd (nSeq, 9 * nF);
			Rs.resize(nF);
			Ss.resize(nF);
			DGs.resize(nF);
			V_inv.resize(nF);
			CMs.resize(nE);
		}

		SparsePCAForShapeGrad::SparsePCAForShapeGrad(PGMesh *_mesh)
		{
			mesh_ = _mesh;
			
			nV = mesh_->n_vertices();
			SpatialTemData.resize(nV);
			nF = mesh_->n_faces();
			nE = mesh_->n_edges();
			DeformGrad.resize(nF);
			nSeq = 2;
			DFMatrix = Eigen::MatrixXd(nSeq, 9 * nF);
			Rs.resize(nF);
			Ss.resize(nF);
			DGs.resize(nF);
			V_inv.resize(nF);
#pragma omp parallel for
			for (int i = 0; i < nV; i++)
			{
				SpatialTemData[i].resize(2);
				PGMesh::VertexHandle vh = mesh_->vertex_handle(i);
				PGMesh::Point &p = mesh_->point(vh);
				Eigen::VectorXd x(3);
				x << p[0], p[1], p[2];
				SpatialTemData[i][0] = x;
				SpatialTemData[i][1] = x;
			}
			CMs.resize(nE);
		}

		SparsePCAForShapeGrad::SparsePCAForShapeGrad(PGMesh *_mesh1, PGMesh *_mesh2)
		{
			mesh_ = _mesh1;

			nV = mesh_->n_vertices();
			SpatialTemData.resize(nV);
			nF = mesh_->n_faces();
			nE = mesh_->n_edges();
			DeformGrad.resize(nF);
			nSeq = 2;
			DFMatrix = Eigen::MatrixXd(nSeq, 9 * nF);
			Rs.resize(nF);
			Ss.resize(nF);
			DGs.resize(nF);
			V_inv.resize(nF);
			for (int i = 0; i < nV; i++)
			{
				SpatialTemData[i].resize(2);
				PGMesh::VertexHandle vh = _mesh1->vertex_handle(i);
				PGMesh::Point &p = _mesh1->point(vh);
				Eigen::VectorXd x(3);
				x << p[0], p[1], p[2];
				PGMesh::VertexHandle vh2 = _mesh2->vertex_handle(i);
				PGMesh::Point &p2 = _mesh2->point(vh);
				Eigen::VectorXd x2(3);
				x2 << p2[0], p2[1], p2[2];
				SpatialTemData[i][0] = x;
				SpatialTemData[i][1] = x2;
			}
			CMs.resize(nE);
		}

		Eigen::VectorXd SparsePCAForShapeGrad::reconstrtuction(PGMesh * _mesh, Eigen::VectorXd gradientComponent, int frameNumber){
			//const int nF = _mesh->n_faces();
			//const int nV = 9;			//每个向量中的元素个数，这里为9；
			const int vector_length = 9;	//每个向量中的元素个数，这里为9；
			

			//用来存储每个面上各个节点的编号
			struct Node_idx{
				int v0_idx;
				int v1_idx;
				int v2_idx;
			};

			//直接存V1-V0的值，2个为一组，代表一个面的V
			double *x = new double[ 2 * nF ];
			double *y = new double[ 2 * nF ];
			double *z = new double[ 2 * nF ];
			Node_idx* node_idx = new Node_idx[ nF ];		//每个面上对应三个点对应一个node_idx用于存储各个节点的编号

			int vCounter = 0;
			int tmp_faceCount = 0;

				
			//遍历每个面，并读取每个面上的所有点的位置信息
			int fNumber = 0;
			for( PGMesh::FaceIter f_it = _mesh->faces_begin(); f_it != _mesh->faces_end(); ++f_it,++fNumber ){
				double x_tmp, y_tmp, z_tmp;
				int i = 0;
				tmp_faceCount++;
				for( PGMesh::FaceVertexIter v_it = _mesh->fv_begin( *f_it ); v_it != _mesh->fv_end( *f_it ); ++v_it,++i ){
					auto point = _mesh->point( *v_it);
					if( i == 0 ){
					x_tmp = point.data()[0];
					y_tmp = point.data()[1];
					z_tmp = point.data()[2];
					node_idx[ fNumber ].v0_idx = v_it->idx();
					}else{
						x[vCounter] = point.data()[0] - x_tmp;
						y[vCounter] = point.data()[1] - y_tmp;
						z[vCounter] = point.data()[2] - z_tmp;
						if( i == 1 )
							node_idx[ fNumber ].v1_idx = v_it->idx();
						else node_idx[ fNumber ].v2_idx = v_it->idx();
						++vCounter;
					}
				}
			}
			
			tmp_faceCount;
			//int n = nV / 1000;		//fix the first vertex in every 1000 vertex
			//VectorXd Js( vector_length * nF + n * 3 );	//add 3n for the fixed vertexs
		//	VectorXd Js( vector_length * nF + 3 * 3 );	//add 3n for the fixed vertexs
			VectorXd Js( vector_length * nF + 9 * fixFacesNum );
			for( int i = 0; i < nF; ++i ){
				int startPos = i * vector_length;		//	start position of every component

				MatrixXd U_( 3, 3 );			// counterpart of U;
				U_.setZero();
				U_( 0, 1 ) = gradientComponent( startPos++ );
				U_( 0, 2 ) = gradientComponent( startPos++ );
				U_( 1, 2 ) = gradientComponent( startPos++ );
				U_( 1, 0 ) = 0 - U_( 0, 1 );
				U_( 2, 0 ) = 0 - U_( 0, 2 );
				U_( 2, 1 ) = 0 - U_( 1, 2 );

				MatrixXd U = U_.exp();
				MatrixXd P( 3, 3 );
				P( 0, 0 ) = gradientComponent( startPos++ ) + 1;
				P( 1, 1 ) = gradientComponent( startPos++ ) + 1;
				P( 2, 2 ) = gradientComponent( startPos++ ) + 1;
				P( 0, 1 ) = gradientComponent( startPos++ );
				P( 0, 2 ) = gradientComponent( startPos++ );
				P( 1, 2 ) = gradientComponent( startPos++ );
				P( 1, 0 ) = P( 0, 1 );
				P( 2, 0 ) = P( 0, 2 );
				P( 2, 1 ) = P( 1, 2 );
				MatrixXd J = U * P;

				//assign values to Js
				startPos = i * vector_length;
				Js.block( startPos, 0, 3, 1 ) = J.row( 0 ).transpose();
				startPos += 3;
				Js.block( startPos, 0, 3, 1 ) = J.row( 1 ).transpose();
				startPos += 3;
				Js.block( startPos, 0, 3, 1 ) = J.row( 2 ).transpose();
			}
			
			//Eigen::SparseMatrix<double> coefficient( nF * vector_length + 3 * n, 3 * nV );		//add 3n for  two fixed vertexs
		//	Eigen::SparseMatrix<double> coefficient( nF * vector_length + 3 * 3, 3 * nV );	
			Eigen::SparseMatrix<double> coefficient( nF * vector_length + 9 * fixFacesNum, 3 * nV );
			coefficient.setZero();
			Matrix< double, 3, 2 > V_input;
			HouseholderQR< Matrix< double, 3, 2 > > QR;
			std::vector<Eigen::Triplet<double> > triple;
			for( int i = 0; i < nF ; ++i){
				int v0_start = node_idx[i].v0_idx * 3;
				int v1_start = node_idx[i].v1_idx * 3;
				int v2_start = node_idx[i].v2_idx * 3;
				int startLine = i * vector_length;

				V_input( 0, 0 ) = x[ i * 2 ];
				V_input( 0, 1 ) = x[ i * 2 + 1];
				V_input( 1, 0 ) = y[ i * 2 ];
				V_input( 1, 1 ) = y[ i * 2 + 1];
				V_input( 2, 0 ) = z[ i * 2 ];
				V_input( 2, 1 ) = z[ i * 2 + 1];
				QR.compute( V_input );
				Matrix< double, 3, 2 > R_tmp = QR.matrixQR().triangularView< Upper>();
				Matrix3d Q_tmp = QR.householderQ();
				Matrix< double, 3, 2 > Q = Q_tmp.block( 0, 0, 3, 2 );
				Matrix2d R = R_tmp.block( 0, 0, 2, 2 );
				Matrix< double, 2, 3 > tmp = R.inverse() * Q.transpose();

				double a11 = tmp( 0, 0 ), a21 = tmp( 1, 0 ), a12 = tmp( 0, 1 ), 
					a22 = tmp( 1, 1 ), a13 = tmp( 0, 2 ), a23 = tmp( 1, 2 );
				double tmp_1 = 0.0 - a11 - a21;
				double tmp_2 = 0.0 - a12 - a22;
				double tmp_3 = 0.0 - a13 - a23;

				//v0
				triple.push_back( Eigen::Triplet<double>(startLine, v0_start, tmp_1));
				triple.push_back( Eigen::Triplet<double>(startLine + 1, v0_start, tmp_2));
				triple.push_back( Eigen::Triplet<double>(startLine + 2, v0_start++, tmp_3));

				triple.push_back( Eigen::Triplet<double>(startLine + 3, v0_start, tmp_1));
				triple.push_back( Eigen::Triplet<double>(startLine + 4, v0_start, tmp_2));
				triple.push_back( Eigen::Triplet<double>(startLine + 5, v0_start++, tmp_3));

				triple.push_back( Eigen::Triplet<double>(startLine + 6, v0_start, tmp_1));
				triple.push_back( Eigen::Triplet<double>(startLine + 7, v0_start, tmp_2));
				triple.push_back( Eigen::Triplet<double>(startLine + 8, v0_start, tmp_3));

				//v1
				triple.push_back( Eigen::Triplet<double>(startLine, v1_start, a11 ));
				triple.push_back( Eigen::Triplet<double>(startLine + 1, v1_start, a12));
				triple.push_back( Eigen::Triplet<double>(startLine + 2, v1_start++, a13));

				triple.push_back( Eigen::Triplet<double>(startLine + 3, v1_start, a11));
				triple.push_back( Eigen::Triplet<double>(startLine + 4, v1_start, a12));
				triple.push_back( Eigen::Triplet<double>(startLine + 5, v1_start++, a13));

				triple.push_back( Eigen::Triplet<double>(startLine + 6, v1_start, a11));
				triple.push_back( Eigen::Triplet<double>(startLine + 7, v1_start, a12));
				triple.push_back( Eigen::Triplet<double>(startLine + 8, v1_start, a13));

				//v2
				triple.push_back( Eigen::Triplet<double>(startLine, v2_start, a21 ));
				triple.push_back( Eigen::Triplet<double>(startLine + 1, v2_start, a22));
				triple.push_back( Eigen::Triplet<double>(startLine + 2, v2_start++, a23));

				triple.push_back( Eigen::Triplet<double>(startLine + 3, v2_start, a21));
				triple.push_back( Eigen::Triplet<double>(startLine + 4, v2_start, a22));
				triple.push_back( Eigen::Triplet<double>(startLine + 5, v2_start++, a23));

				triple.push_back( Eigen::Triplet<double>(startLine + 6, v2_start, a21));
				triple.push_back( Eigen::Triplet<double>(startLine + 7, v2_start, a22));
				triple.push_back( Eigen::Triplet<double>(startLine + 8, v2_start, a23));
			}
			coefficient.setFromTriplets( triple.begin(), triple.end() );
			//init the fix vertex by using the vertexes of the first face
			VectorXd fixPoints = fixFacePoints[frameNumber];
			for( int i = 0; i < fixFacesNum; ++i ){
				int idxPosition = i * 12;
				int beginLine = i * 9;
				int firstVertex = fixPoints( idxPosition );
				int secondVertex = fixPoints( idxPosition + 4 );
				int thirdVertex = fixPoints( idxPosition + 8 );
				++idxPosition;

				Js( vector_length * nF + beginLine ) = fixPoints( idxPosition++ );
				Js( vector_length * nF + beginLine + 1 ) = fixPoints( idxPosition++ );
				Js( vector_length * nF + beginLine + 2 ) = fixPoints( idxPosition++ );
				coefficient.insert( nF * vector_length + beginLine, firstVertex * 3 ) = 1;
				coefficient.insert( nF * vector_length + beginLine + 1, firstVertex * 3 + 1 ) = 1;
				coefficient.insert( nF * vector_length + beginLine + 2, firstVertex * 3 + 2 ) = 1;
				beginLine += 3;
				++idxPosition;

				Js( vector_length * nF + beginLine ) = fixPoints( idxPosition++ );
				Js( vector_length * nF + beginLine + 1 ) = fixPoints( idxPosition++ );
				Js( vector_length * nF + beginLine + 2 ) = fixPoints( idxPosition++ );
				coefficient.insert( nF * vector_length + beginLine, secondVertex * 3 ) = 1;
				coefficient.insert( nF * vector_length + beginLine + 1, secondVertex * 3 + 1 ) = 1;
				coefficient.insert( nF * vector_length + beginLine + 2, secondVertex * 3 + 2 ) = 1;
				beginLine += 3;
				++idxPosition;

				Js( vector_length * nF + beginLine ) = fixPoints( idxPosition++ );
				Js( vector_length * nF + beginLine + 1 ) = fixPoints( idxPosition++ );
				Js( vector_length * nF + beginLine + 2 ) = fixPoints( idxPosition++ );
				coefficient.insert( nF * vector_length + beginLine, thirdVertex * 3 ) = 1;
				coefficient.insert( nF * vector_length + beginLine + 1, thirdVertex * 3 + 1 ) = 1;
				coefficient.insert( nF * vector_length + beginLine + 2, thirdVertex * 3 + 2 ) = 1;

			}

			PGMesh::FaceIter tmp_iter = _mesh->faces_begin();
			PGMesh::FaceVertexIter tmp_iter2 = _mesh->fv_begin( *tmp_iter);
			int firstVertex = tmp_iter2->idx();
			++tmp_iter2;
			int secondVertex = tmp_iter2->idx();
			++tmp_iter2;
			int thirdVertex = tmp_iter2->idx();

			int tmp_idx = 0;
			VectorXd fixVertex = SpatialTemData[firstVertex][frameNumber];
			Js( vector_length * nF + tmp_idx ) = fixVertex( 0 );
			Js( vector_length * nF + tmp_idx + 1 ) = fixVertex( 1 );
			Js( vector_length * nF + tmp_idx + 2 ) = fixVertex( 2 );
			coefficient.insert( nF * vector_length + tmp_idx, firstVertex * 3 ) = 1;
			coefficient.insert( nF * vector_length + tmp_idx + 1, firstVertex * 3 + 1 ) = 1;
			coefficient.insert( nF * vector_length + tmp_idx + 2, firstVertex * 3 + 2 ) = 1;
			tmp_idx += 3;
			fixVertex = SpatialTemData[secondVertex][frameNumber];
			Js( vector_length * nF + tmp_idx ) = fixVertex( 0 );
			Js( vector_length * nF + tmp_idx + 1 ) = fixVertex( 1 );
			Js( vector_length * nF + tmp_idx + 2 ) = fixVertex( 2 );
			coefficient.insert( nF * vector_length + tmp_idx, secondVertex * 3 ) = 1;
			coefficient.insert( nF * vector_length + tmp_idx + 1, secondVertex * 3 + 1 ) = 1;
			coefficient.insert( nF * vector_length + tmp_idx + 2, secondVertex * 3 + 2 ) = 1;
			tmp_idx += 3;
			fixVertex = SpatialTemData[thirdVertex][frameNumber];
			Js( vector_length * nF + tmp_idx ) = fixVertex( 0 );
			Js( vector_length * nF + tmp_idx + 1 ) = fixVertex( 1 );
			Js( vector_length * nF + tmp_idx + 2 ) = fixVertex( 2 );
			coefficient.insert( nF * vector_length + tmp_idx, thirdVertex * 3 ) = 1;
			coefficient.insert( nF * vector_length + tmp_idx + 1, thirdVertex * 3 + 1 ) = 1;
			coefficient.insert( nF * vector_length + tmp_idx + 2, thirdVertex * 3 + 2 ) = 1;




			SimplicialCholesky<SparseMatrix<double>> solver;
			solver.compute( coefficient.transpose() * coefficient );
			VectorXd V;
			 V = solver.solve( coefficient.transpose() * Js);


			 //calculate the error for visualization
			 //VectorXd component = gradientComponent -  Xmean;
			 //double maxWeight = -100000.0;
			 //double * vertexWeight = new double[ nV ];
			 //memset( vertexWeight, 0.0, nV * sizeof( double ) );

			 //for( int i = 0; i < nF; ++i ){
				//	int tmp_idx = i * 9;
				//	VectorXd tmp = component.block( tmp_idx, 0, 9, 1 );
				//	/*
				//	double a = 0.0;
				//	for( int tmp_i = 0; tmp_i < 9; ++tmp_i )
				//		a += abs( tmp( tmp_i) );
				//	face_weight[i] = a;
				//	*/
				//	double weight = ( tmp.cwiseAbs() ).sum();
				//	if( weight > maxWeight )
				//		maxWeight = weight;

				//	PGMesh::FaceHandle fHandle = mesh_->face_handle(i);
				//	PGMesh::FaceVertexIter vertexIt = mesh_->fv_begin( fHandle );
				//	for( ; vertexIt != mesh_->fv_end( fHandle); ++ vertexIt )
				//	if( vertexWeight[vertexIt.handle().idx()] < weight )
				//		vertexWeight[vertexIt.handle().idx()] = weight;
				//}

			 // output a reconstructed mesh
			// char *pPath = new char[16];
			// sprintf( pPath,"%d_components%d_th.obj",k,frameNumber );
			// std::string tmp_str = pPath;
			// std::string str = "C:\\Users\\luhuina\\Desktop\\results\\reconstructModel\\" + tmp_str;
			// ofstream out(str);
			// for( int i = 0; i < nV; ++i ){
			//	 /*PGMesh::VertexHandle it = mesh_->vertex_handle(i);
			//	 PGMesh::Point point = mesh_->point( it );
			//	 int colorPersent = vertexWeight[i]  / maxWeight * 99;
			//	 unsigned char* colorPtr = colorMap + colorPersent * 3;*/
			//	 int tmp = i * 3;
			//	 out<<"v"<<" "<<V( tmp + 0 )<<" "<<V( tmp + 1 )<<" "<<V( tmp + 2)<<endl;
			////	 out<<(int)colorPtr[0]<<" "<<(int)colorPtr[1]<<" "<<(int)colorPtr[2]<<endl;
			// }
			// for( PGMesh::FaceIter fIter = _mesh->faces_begin(); fIter != _mesh->faces_end(); ++fIter ){
			//	 out<<"f"<<" ";
			//	 for( PGMesh::FaceVertexIter fvIter = _mesh->fv_begin( fIter.handle()); fvIter != _mesh->fv_end( fIter.handle()); ++fvIter)
			//		 out<<fvIter.handle().idx() + 1<<" ";
			//	 out<<endl;
			// }
			// out.close();
			
			 return V;
		}

		bool SparsePCAForShapeGrad::reconstrtuction(){
			if( W.size() <= 0 && KComp.size() <= 0 ){
				string str = "C:\\Users\\luhuina\\Desktop\\results\\component.txt";
				readKComponent( str );
				str = "C:\\Users\\luhuina\\Desktop\\results\\W.txt";
				readW( str );
				str = "C:\\Users\\luhuina\\Desktop\\results\\Xmean.txt";
				readXmean( str );
			}
			char* tmp = new char[130];
			for( int i = 0; i < nSeq/3; ++i )
				{
					Eigen::VectorXd cords = reconstrtuction( mesh_, getOneComponent( i * 3), i * 3);
					/*sprintf( tmp,"C:\\Users\\luhuina\\Desktop\\DGDATA\\DGDATA\\%d\\%d-th.obj", KComp.size(), 3*i );
					string str = tmp;
					outMesh(cords, str);*/
			}
			//reconstrtuction( mesh_, getComponent(), 0 );
			return true;
		}

		void SparsePCAForShapeGrad::subdivision(){
			mesh_->add_property( cogs );
			std::map < int, PGMesh::VertexHandle > iter_to_vHandle;			//将新mesh中的point与旧mesh中的vertex对应起来，避免重复创建相同的点。
			for( PGMesh::FaceIter f_iter = mesh_->faces_begin(); f_iter != mesh_->faces_end(); ++f_iter ){
				PGMesh::VertexHandle midVHandle[4];
				PGMesh::FaceHandle fHandle;
				PGMesh::FaceVertexIter fv_iter = mesh_->fv_begin( *f_iter);
				PGMesh::Point point1, point2, point3, point;

				if( iter_to_vHandle.find( fv_iter->idx() ) == iter_to_vHandle.end()){
					point1 = PGMesh::Point( mesh_->point( *fv_iter) );
					midVHandle[0] = subdivMesh.add_vertex( point1 );
					iter_to_vHandle[ fv_iter->idx() ] = midVHandle[0];
				}else{
					midVHandle[0] = iter_to_vHandle[ fv_iter->idx() ];
					point1 = subdivMesh.point( midVHandle[0] );
				}
				++fv_iter;

				if( iter_to_vHandle.find( fv_iter->idx() ) == iter_to_vHandle.end()){
					point2 = PGMesh::Point( mesh_->point( *fv_iter) );
					midVHandle[1] = subdivMesh.add_vertex( point2 );
					iter_to_vHandle[ fv_iter->idx() ] = midVHandle[1];
				}else{
					midVHandle[1] = iter_to_vHandle[ fv_iter->idx() ];
					point2 = subdivMesh.point( midVHandle[1] );
				}
				++fv_iter;

				if( iter_to_vHandle.find( fv_iter->idx() ) == iter_to_vHandle.end()){
					point3 = PGMesh::Point( mesh_->point( *fv_iter) );
					midVHandle[2] = subdivMesh.add_vertex( point3 );
					iter_to_vHandle[ fv_iter->idx() ] = midVHandle[2];
				}else{
					midVHandle[2] = iter_to_vHandle[ fv_iter->idx() ];
					point3 = subdivMesh.point( midVHandle[2] );
				}
				
				PGMesh::Point newPoint( ( point1 + point2 + point3 ) / 3.0 );
				midVHandle[3] = subdivMesh.add_vertex(( newPoint ));
				mesh_->property( cogs, *f_iter ) = midVHandle[3].idx();
				
				std::vector< PGMesh::VertexHandle > face_vHandle( 3 );
				face_vHandle.clear();
				face_vHandle.push_back( midVHandle[0]);
				face_vHandle.push_back( midVHandle[1]);
				face_vHandle.push_back( midVHandle[3]);
				subdivMesh.add_face( face_vHandle );

				face_vHandle.clear();
				face_vHandle.push_back( midVHandle[1]);
				face_vHandle.push_back( midVHandle[2]);
				face_vHandle.push_back( midVHandle[3]);
				subdivMesh.add_face( face_vHandle );

				face_vHandle.clear();
				face_vHandle.push_back( midVHandle[2]);
				//`face_vHandle.push_back( midVHandle[2]);
				face_vHandle.push_back( midVHandle[0]);
				face_vHandle.push_back( midVHandle[3]);
				subdivMesh.add_face( face_vHandle );
			}
			//OpenMesh::IO::write_mesh( subdivMesh, "C:\\Users\\luhuina\\Desktop\\outputnew_1.obj");
		}

		Eigen::VectorXd SparsePCAForShapeGrad::getOneComponent( int i ){
			/*
			VectorXd oneComponent( 9 * nF );
			if( i == -1)
				oneComponent.setZero();
			else
				for( int j = 0; j < nF; ++j ){
					int startPos = j * 9;
					oneComponent.block( startPos, 0, 9, 1 ) = DeformGrad[j][i];
				}
			return oneComponent;
			*/
			return  ( W[i].transpose() * KComp).transpose()  + Xmean ;
		}

		Eigen::VectorXd SparsePCAForShapeGrad::getComponent(){

			//VectorXd tmp = W[0] + W_single;
			return ( W_single.transpose() * KComp ).transpose()  + Xmean;
		}

		bool SparsePCAForShapeGrad::writeDeformGrad( string & path ){
			ofstream out( path );

			for( int i = 0; i < nSeq; ++i ){
				for( int j = 0; j < nF; ++j ){
					int k = 0;
				while( k < 9 ){
					double tmp = DeformGrad[j][i]( k );
					if( tmp >= -1.0e-10 && tmp <= 1.0e-10)
						out<<0.0<<" ";
					else
						out<<setprecision( 15 )<<tmp<<" ";
					++k;
					}
				}
				out<<endl;
			}
			out.close();
			return 1;
		}

		bool SparsePCAForShapeGrad::writeSubdivMesh_vertex( string & path ){
			ofstream out( path );
			PGMesh::VertexIter verticeBegin = subdivMesh.vertices_begin(), verticeEnd = subdivMesh.vertices_end();
			for( PGMesh::VertexIter verIte = verticeBegin; verIte != verticeEnd; ++verIte ){
				PGMesh::Point point = subdivMesh.point( *verIte);
				out<<point.data()[0]<<" "<<point.data()[1]<<" "<<point.data()[2]<<endl;
			}
			out.close();
			return 1;
		}

		bool SparsePCAForShapeGrad::writeSubdivMesh_face( string & path ){
			ofstream out( path );
			PGMesh::FaceIter faceBegin = subdivMesh.faces_begin(), faceEnd = subdivMesh.faces_end();
			for( PGMesh::FaceIter faceIte = faceBegin; faceIte != faceEnd; ++faceIte ){
				PGMesh::FaceVertexIter f_vIte = subdivMesh.fv_begin( *faceIte);
				out<<f_vIte->idx()<<" ";
				++f_vIte;
				out<<f_vIte->idx()<<" ";
				++f_vIte;
				out<<f_vIte->idx()<<endl;
			}
			out.close();
			return 1;
		}

		bool SparsePCAForShapeGrad::writeMeshProperty( string & filePath ){
			ofstream out( filePath );
			
			for( PGMesh::FaceIter fBegin = mesh_->faces_begin(); fBegin != mesh_->faces_end(); ++fBegin ){
				out<<mesh_->property( cogs, *fBegin )<<" ";
			}
			out.close();
			return true;
			
		}

		bool SparsePCAForShapeGrad::write_nSeq_nF( string & filePath ){
			ofstream out( filePath);
			out<<nSeq<<" "<<nF;
			out.close();
			return true;
		}

		bool SparsePCAForShapeGrad::readDeformGrad(){

			ifstream in( "C:\\Users\\luhuina\\Desktop\\results\\deformGrad.txt" );
			for( int i = 0; i < nSeq; ++i ){
				for( int j = 0; j < nF; ++j ){
					if( i == 0 )
						DeformGrad[j].resize( nSeq );
					VectorXd deformGrad_tmp( 9 );
					double tmp;
					int idx = 0;
					while( idx < 9 ){
						in>>tmp;
						deformGrad_tmp(idx) = tmp;
						++idx;
					}
					DeformGrad[j][i] = deformGrad_tmp;

				}
			}
			in.close();
			return true;
		}

		bool SparsePCAForShapeGrad::readKComponent( string & filePath ){
			ifstream in( filePath );
			KComp.resize( k, 9 * nE );
			int i = 0, j = 0;
			while( !in.eof() && i < k ){
				double tmp;
				in>>tmp;
				KComp( i , j ) = tmp;
				++j;
				if( j % ( 9 * nF ) == 0 ){
					++i;
					j = 0;
				}
			}
			
			in.close();
			return true;
		}

		bool SparsePCAForShapeGrad::readW( string & filePath ){
			W.clear();
			ifstream in( filePath );
			int j = 0;
			while( !in.eof() && j < nSeq ){
				VectorXd a( k );
				for( int i = 0; i < k; ++i )
					in>>a( i );
				W.push_back( a );
				++j;
			}
			in.close();
			return true;
		}
		
		bool SparsePCAForShapeGrad::readXmean( string & filePath ){
			ifstream in( filePath );
			Xmean.resize( 9 * nF );
			int i = 0;
			while( !in.eof() && i < 9 * nF )
				in>>Xmean( i++ );
			in.close();
			return true;
		}

		bool SparsePCAForShapeGrad::readFixFacePoints(){
			fixFacePoints.clear();
			ifstream in( "E:\\fixPoints.txt" );
			
			for( int i = 0; i < nSeq; ++i ){
				int idx;
				double x, y, z;
				int j = 0; 
				VectorXd tmp( 12 * fixFacesNum );
				while( j < fixFacesNum){
					for( int jj = 0; jj < 3; ++jj ){
						in>>idx>>x>>y>>z;
						VectorXd tt( 4 );
						tt<<idx, x, y, z;
						tmp.block( j * 12 + jj * 4, 0, 4, 1 ) = tt;
					}
					
					++j;
				}
				
				
				fixFacePoints.push_back( tmp );
			}
			
			in.close();
			return true;
		}

		bool SparsePCAForShapeGrad::writeVisualizeComponent_color( string & filePath ){
			ofstream out( filePath );
			double *face_weight = new double[nF];
			double * vertex_color = new double[ nV ];
			
			for( int seqNum = 0; seqNum < k; ++seqNum ){
				
				memset( vertex_color, 0.0, nV * sizeof( double ) );
				memset( face_weight, 0.0, nF * sizeof( double ) );

				VectorXd component = KComp.row( seqNum );
				double minWeight = 100000.0, maxWeight = -1000000.0;
				for( int i = 0; i < nF; ++i ){
					int tmp_idx = i * 9;
					VectorXd tmp = component.block( tmp_idx, 0, 9, 1 );
					/*
					double a = 0.0;
					for( int tmp_i = 0; tmp_i < 9; ++tmp_i )
						a += abs( tmp( tmp_i) );
					face_weight[i] = a;
					*/
					face_weight[i] = ( tmp.cwiseAbs() ).sum();
					if( face_weight[i] > maxWeight )
						maxWeight = face_weight[i];
					if( face_weight[i] < minWeight )
						minWeight = face_weight[i];

					PGMesh::FaceHandle fHandle = mesh_->face_handle(i);
					PGMesh::FaceVertexIter vertexIt = mesh_->fv_begin( fHandle );
					for( ; vertexIt != mesh_->fv_end( fHandle); ++ vertexIt )
					if( vertex_color[vertexIt->idx()] < face_weight[i] )
						vertex_color[vertexIt->idx()] = face_weight[i];
				}

				double range = maxWeight - minWeight;
				for( int i = 0; i < nF; ++i ){
					int colorPersent = (face_weight[ i ] - minWeight) / range * 99;
					unsigned char* colorPtr = colormap6 + colorPersent * 3;
					out<<(int)colorPtr[0]<<" "<<(int)colorPtr[1]<<" "<<(int)colorPtr[2]<<" ";
				}
				out<<endl;

				//output the components' visualization files
				char *pPath = new char[50];
				sprintf_s( pPath,50, "%d.obj",seqNum );
				std::string tmp_str = pPath;
				std::string str = "C:\\Users\\luhuina\\Desktop\\results\\visualize\\components\\" + tmp_str;
				delete []pPath;
				pPath = NULL;
				ofstream outputComponent(str);
				for( int i = 0; i < nV; ++i ){
					PGMesh::VertexHandle it = mesh_->vertex_handle(i);
					PGMesh::Point point = mesh_->point( it );
					int colorPersent = (vertex_color[i] - minWeight) / range * 99;
					unsigned char* colorPtr = colormap6 + colorPersent * 3;
					outputComponent<<"v "<<point.data()[0]<<" "<<point.data()[1]<<" "<<point.data()[2]<<" ";
					outputComponent<<(int) (colorPtr[0])  <<" "<< int(colorPtr[1])<< " " << int (colorPtr[2])<< endl;
				}
				for( PGMesh::FaceIter fIter = mesh_->faces_begin(); fIter != mesh_->faces_end(); ++fIter ){
					outputComponent<<"f"<<" ";
					for( PGMesh::FaceVertexIter fvIter = mesh_->fv_begin( *fIter); fvIter != mesh_->fv_end( *fIter); ++fvIter)
						outputComponent<<fvIter->idx() + 1<<" ";
					outputComponent<<endl;
				}
				outputComponent.close();

			}
			out.close();
			delete [] face_weight;
			face_weight = NULL;
			delete [] vertex_color;
			vertex_color = NULL;
			return true;
		}
		
		/** 根据点坐标导出网格*/
		void SparsePCAForShapeGrad::outMesh(Eigen::VectorXd cordComponent, string path)
		{
			ofstream offfile(path);
			offfile << "OFF"<<endl;
			offfile << nV << " " << nF<< " "<< 0 <<endl;
			
			for (int i = 0; i< nV; i++)
			{
				offfile << cordComponent(3 * i + 0) << " " <<cordComponent(3 * i + 1) << " " << cordComponent(3 * i + 2)  << endl;
			}
			for (PGMesh::FaceIter fit = mesh_->faces_begin(); fit != mesh_->faces_end(); ++ fit)
			{
				offfile << "3" << " ";
				for (PGMesh::FVIter fvit = mesh_->fv_begin(*fit); fvit != mesh_->fv_end(*fit); ++ fvit)
				{
					offfile << fvit->idx() << " ";
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

		void SparsePCAForShapeGrad::outInverFaceMesh()
		{
			for (int seqNo = 0; seqNo < nSeq; seqNo++)
			{
				char *pPath = new char[30];
				sprintf_s(pPath,30,"mesh%04d.off", seqNo);
				ofstream offfile(pPath);
				offfile << "OFF"<<endl;
				offfile << nV << " " << nF<< " "<< 0 <<endl;
				for (int i = 0; i< nV; i++)
				{
					Eigen::Vector3d cord = SpatialTemData[i][seqNo];
					offfile <<cord(0)<< " " <<cord(1) << " " << cord(2)  << endl;
				}
				for (PGMesh::FaceIter fit = mesh_->faces_begin(); fit != mesh_->faces_end(); ++ fit)
				{
					offfile << "3" << " ";
				
					PGMesh::FVIter fvit = mesh_->fv_begin(*fit);
					PGMesh::VertexHandle vv0 = *fvit;
					++fvit;
					PGMesh::VertexHandle vv1 = *fvit;
					++fvit;
					PGMesh::VertexHandle vv2 = *fvit;
					offfile << vv0.idx() << " "<< vv1.idx() << " " << vv2.idx() << " "  ;
					offfile << endl;
				}
				offfile.close();
			}
			
		}
		

		void SparsePCAForShapeGrad::writeFaceErrorOfComp(char * filePath, Eigen::VectorXd vertexCords, Eigen::VectorXd face_error)
		{
			double * vertex_color = new double[ nV ];
			memset( vertex_color, 0.0, nV * sizeof( double ) );
			double min = 1000000.0, max = -1000000.0;
			for( int i = 0; i < nF; ++i ){
				
				double sum = face_error(i) ;
				if( min > sum )	min = sum;
				if( max < sum )	max = sum;


				PGMesh::FaceHandle fHandle = mesh_->face_handle( i );
				PGMesh::FaceVertexIter vertexIt = mesh_->fv_begin( fHandle );
				for( ; vertexIt != mesh_->fv_end( fHandle); ++ vertexIt )
					if( vertex_color[vertexIt->idx()] < face_error(i) )
						vertex_color[vertexIt->idx()] = face_error(i);
			}
			double range = max - min;

			ofstream outputError(filePath);
			ofstream offfile1("DGcolor.txt");
			for( int i = 0; i < nV; ++i ){
				PGMesh::VertexHandle it = mesh_->vertex_handle(i);
				PGMesh::Point point = mesh_->point( it );
				double t = (vertex_color[i] - min) / range ;
				double rr, gg, bb;
				
				/*G = (0, 1,0), Y= (1,1,0), R = (1,0,0) 
                t在0~0.5用(1-2t)Y + 2tG， 在0.5~1用 (2t-1)R + (2-2t)Y*/
				if(t < 0.5)
				{
					rr = (1 - 2 * t) * 1 + 2 * t * 0;
					gg = (1 - 2 * t) * 1 + 2 * t * 1;
					bb = (1 - 2 * t) * 0 + 2 * t * 0;
				}
				else
				{
					rr = (2 * t - 1) * 1 + (2 - 2 * t) * 1;
					gg = (2 * t - 1) * 0 + (2 - 2 * t) * 1;
					bb = (2 * t - 1) * 0 + (2 - 2 * t) * 0;
				}
				int colorPersent = (vertex_color[i] - min) / range * 99;
				unsigned char* colorPtr = colormap6 + colorPersent * 3;
				//outputError<<"v "<<point.data()[0]<<" "<<point.data()[1]<<" "<<point.data()[2]<<" ";
				outputError<<"v "<<vertexCords(3 * i + 0)<<" "<<vertexCords(3 * i + 1)<<" "<<vertexCords(3 * i + 2)<<" ";
				outputError<<(int) (colorPtr[0])  <<" "<< int(colorPtr[1])<< " " << int (colorPtr[2])<< endl;
				offfile1<<(int) (colorPtr[0])  <<" "<< int(colorPtr[1])<< " " << int (colorPtr[2])<< endl;
				//outputError <<int (rr * 255)  <<" "<< int (gg * 255) << " " << int (bb * 255) << endl;
			}
			for( PGMesh::FaceIter fIter = mesh_->faces_begin(); fIter != mesh_->faces_end(); ++fIter ){
				outputError<<"f"<<" ";
				for( PGMesh::FaceVertexIter fvIter = mesh_->fv_begin( *fIter); fvIter != mesh_->fv_end( *fIter); ++fvIter)
				outputError<<fvIter->idx() + 1<<" ";
				/*PGMesh::FVIter fvit = mesh_->fv_begin(fIter.handle());
				PGMesh::VertexHandle vv0 = fvit.handle();
				++fvit;
				PGMesh::VertexHandle vv1 = fvit.handle();
				++fvit;
				PGMesh::VertexHandle vv2 = fvit.handle();
				outputError << vv2.idx() + 1 << " "<< vv1.idx() + 1 << " " << vv0.idx() + 1 << " "  ;*/

				outputError<<endl;
			}
			outputError.close();
			delete [] vertex_color;
			offfile1.close();
		}
		bool SparsePCAForShapeGrad::writeError( string & filePath ){
			ofstream out( filePath );
			double * face_error = new double[ nF ];
			double * vertex_color = new double[ nV ];
			memset( vertex_color, 0.0, nV * sizeof( double ) );
			double min = 1000000.0, max = -1000000.0;
			for( int i = 0; i < nF; ++i ){
				face_error[i] = 0.0;
				double sum = 0.0;
				VectorXd defGrad_Xmean = Xmean.block( i * 9, 0, 9, 1 );
				for( int j = 0; j < nSeq; ++j ){
					sum = sum + ( DeformGrad[i][j] - defGrad_Xmean ).cwiseAbs().sum();
				}
				face_error[i] = sum;
				if( min > sum )	min = sum;
				if( max < sum )	max = sum;


				PGMesh::FaceHandle fHandle = mesh_->face_handle( i );
				PGMesh::FaceVertexIter vertexIt = mesh_->fv_begin( fHandle );
				for( ; vertexIt != mesh_->fv_end( fHandle); ++ vertexIt )
					if( vertex_color[vertexIt->idx()] < face_error[i] )
						vertex_color[vertexIt->idx()] = face_error[i];

				

			}

			double range = max - min;
			for( int i = 0; i < nF; ++i ){
				int colorPersent = ( face_error[ i ] - min ) / range * 99;
				unsigned char* colorPtr = colormap6 + colorPersent * 3;
				out<<(int)colorPtr[0]<<" "<<(int)colorPtr[1]<<" "<<(int)colorPtr[2]<<" ";
			}
			out.close();

			//output the components' visualization files
			std::string str = "C:\\Users\\luhuina\\Desktop\\results\\visualize\\error.obj" ;
			ofstream outputError(str);
			for( int i = 0; i < nV; ++i ){
				PGMesh::VertexHandle it = mesh_->vertex_handle(i);
				PGMesh::Point point = mesh_->point( it );
				int colorPersent = (vertex_color[i] - min) / range * 99;
				unsigned char* colorPtr = colormap6 + colorPersent * 3;
				outputError<<"v "<<point.data()[0]<<" "<<point.data()[1]<<" "<<point.data()[2]<<" ";
				outputError<<(int) (colorPtr[0])  <<" "<< int(colorPtr[1])<< " " << int (colorPtr[2])<< endl;
			}
			for( PGMesh::FaceIter fIter = mesh_->faces_begin(); fIter != mesh_->faces_end(); ++fIter ){
				outputError<<"f"<<" ";
				for( PGMesh::FaceVertexIter fvIter = mesh_->fv_begin( *fIter); fvIter != mesh_->fv_end( *fIter); ++fvIter)
					outputError<<fvIter->idx() + 1<<" ";
				outputError<<endl;
			}
			outputError.close();

			delete [] face_error;
			delete [] vertex_color;

			return true;
		}

		bool SparsePCAForShapeGrad::visualizeComponent( int i_th ){
			VectorXd component = KComp.row( i_th );
			double * face_weight = new double[ nF ];
			double minWeight = 1000000.0, maxWeight = -100000.0;
		//	double *vertex_color = new double[ nV ];
		//	memset( vertex_color, 0, nV * sizeof( double ) );
			for( int i = 0; i < nF; ++i ){
				int tmp_idx = i * 9;
				VectorXd tmp = component.block( tmp_idx, 0, 9, 1 );
				double a = 0.0;
				for( int tmp_i = 0; tmp_i < 9; ++tmp_i )
					a += abs( tmp( tmp_i) );
				face_weight[i] = a;
				if( face_weight[i] > maxWeight )
					maxWeight = face_weight[i];
				if( face_weight[i] < minWeight )
					minWeight = face_weight[i];
			//	PGMesh::FaceHandle faceHandle( i );

			//	for( PGMesh::FaceVertexIter ite = mesh_->fv_begin( faceHandle ); ite != mesh_->fv_end( faceHandle); ++ite ){
			//		int idx = ite.handle().idx();
			//		if( vertex_color[idx] < face_weight[i] )
			//			vertex_color[idx] = face_weight[i];
			//	}

			}
			double range = maxWeight - minWeight;

			mesh_->request_face_colors();
			for( PGMesh::FaceIter it = mesh_->faces_begin(); it != mesh_->faces_end(); ++it ){
				PGMesh::Color c;
				int colorPersent = (face_weight[ it->idx() ] - minWeight) / range * 99;
				unsigned char* colorPtr = colormap6 + colorPersent * 3;
				c[0] = colorPtr[0];
				c[1] = colorPtr[1];
				c[2] = colorPtr[2];
				mesh_->set_color( *it, c);
			}

			/*
			char *pPath = new char[16];
			sprintf( pPath,"%d.obj",i_th );
			std::string tmp_str = pPath;
			std::string str = "C:\\Users\\luhuina\\Desktop\\results\\visualize\\" + tmp_str;
			ofstream out(str);
			for( int i = 0; i < nV; ++i ){
				PGMesh::VertexHandle it = mesh_->vertex_handle(i);
				PGMesh::Point point = mesh_->point( it );
				int colorPersent = (vertex_color[it.idx()] - minWeight) / range * 99;
				unsigned char* colorPtr = colorMap + colorPersent * 3;
				out<<"v "<<point.data()[0]<<" "<<point.data()[1]<<" "<<point.data()[2]<<" ";
				out<<(int) (colorPtr[0])  <<" "<< int(colorPtr[1])<< " " << int (colorPtr[2])<< endl;

			}
			for( PGMesh::FaceIter fIter = mesh_->faces_begin(); fIter != mesh_->faces_end(); ++fIter ){
				out<<"f"<<" ";
				for( PGMesh::FaceVertexIter fvIter = mesh_->fv_begin( fIter.handle()); fvIter != mesh_->fv_end( fIter.handle()); ++fvIter)
					out<<fvIter.handle().idx() + 1<<" ";
				out<<endl;
			}
			out.close();*/

			return true;
		}

bool SparsePCAForShapeGrad::doSPLOCS_py( string & deformGradPath, string & subdiviVertexPath, 
			string & subdiviFacePath , string & meshPropertyPath, int ndim ){

				return true;
}

void SparsePCAForShapeGrad::getColor( double color, double& r, double& g, double& b ){
			r = __min(__max(1.5 - 4 * abs(color - 0.75), 0.), 1.);
			g = __min(__max(1.5 - 4 * abs(color - 0.5), 0.), 1.);
			b = __min(__max(1.5 - 4 * abs(color - 0.25), 0.), 1.);
			r = r * 255;
			g = g * 255;
			b = b * 255;
}

bool SparsePCAForShapeGrad::setK( int k ){
		this->k = k;
		return true;
		/*W_single.resize( k );
		W_single.setZero();
		if( nSeq > 1 )
		return doSPLOCS_py( string( "C:\\Users\\luhuina\\Desktop\\results\\deformGrad.txt" ),
		string( "C:\\Users\\luhuina\\Desktop\\results\\subdiv_vertex.txt" ),
		string( "C:\\Users\\luhuina\\Desktop\\results\\subdiv_face.txt" ),
		string( "C:\\Users\\luhuina\\Desktop\\results\\MeshProperty.txt" ), 9 );*/
//		else return false;
}

bool SparsePCAForShapeGrad::set_nSeq( int nSeq ){
	this->nSeq = nSeq;
	readFixFacePoints();
	return true;
}

/** 提取第i帧的网格坐标向量*/
		Eigen::VectorXd SparsePCAForShapeGrad::getCords (int ith)
		{
			Eigen::VectorXd l(3 * nV);
			for(int i = 0; i < nV; i++)
			{
				l.segment(3 * i, 3) = SpatialTemData[i][ith];
			}
			return l;
		}

/** 可视化component(给定每条边的误差)*/
		void SparsePCAForShapeGrad::visualizeComp(Eigen::VectorXd cordComponent, Eigen::VectorXd edgeError, char* path)
		{
			ofstream offfile(path);
			ofstream offfile1("outColor.txt");
			offfile << "COFF"<<endl;
			offfile << nV << " " << nF<< " "<< 0 <<endl;
			std::vector<double> errors;
			errors.resize(nV);
			//遍历每个点
			for (int i = 0; i< nV; i++)
			{
				double pointerror = 0;
				PGMesh::VertexHandle vit = mesh_->vertex_handle(i);
				PGMesh::VertexOHalfedgeIter veit = mesh_->voh_begin(vit);
				int j = 0;
				errors[i] = 0;
				for (; veit != mesh_->voh_end(vit); ++veit)
				{
					int edgeidx = (*mesh_).edge_handle(*veit).idx();
					pointerror += edgeError(edgeidx);
					if (edgeError(edgeidx) > errors[i])
						errors[i] = edgeError(edgeidx);
					j++;
				}
				
			}
			double  maxElement = errors[0], minElement = errors[0];
			for (int i = 0; i < nV; ++i)
			{

				maxElement = maxElement > errors[i] ? maxElement : errors[i];
				minElement = minElement < errors[i] ? minElement : errors[i];
			}
			
			float range =maxElement - minElement;
			if (abs(range) < 1e-3)
			{
				range = 1;
			}

			//int ithElement = 0;
			for (int i = 0; i < nV; i++)
			{
				PGMesh::VertexHandle it = mesh_->vertex_handle(i);
				double v = errors[i];
// 				if (ithEigen == 1) //////////////////////////////////////201412/18
// 				{
// 					v = -pEigen[ithElement];
// 
// 				}				
				PGMesh::Color c;
				int colorPersent = (v - minElement) / range * 99;
				unsigned char* colorPtr = colormap6 + colorPersent * 3;
				//c.values_[0] = colorPtr[0];
				//c.values_[1] = colorPtr[1];
				//c.values_[2] = colorPtr[2];
				Eigen::Vector3d position = cordComponent.segment(3 * i, 3);
				offfile << "v" << " " << position( 0) << " " <<position(1) << " " << position(2)  <<" " ;
				offfile <<(int) (colorPtr[0])  <<" "<< int(colorPtr[1])<< " " << int (colorPtr[2])<< endl;
				offfile1<<(int) (colorPtr[0])  <<" "<< int(colorPtr[1])<< " " << int (colorPtr[2])<< endl;
			}
			for (PGMesh::FaceIter fit = mesh_->faces_begin(); fit != mesh_->faces_end(); ++ fit)
			{
				offfile << "f" << " ";
				for (PGMesh::FVIter fvit = mesh_->fv_begin(*fit); fvit != mesh_->fv_end(*fit); ++ fvit)
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
			offfile1.close();
		}

void SparsePCAForShapeGrad::setW_single(){
	ifstream in("E:\\GeometryProcessing6.0-20150803\\GeometryProcessing6.0\\src\\singleW.txt" );
	int idx;
	double values;
	W_single.setZero();
	while( in>>idx>>values )
		W_single( idx ) = values;
	in.close();

}
void SparsePCAForShapeGrad::setWs(std::vector<double>::iterator begin, std::vector< double >::iterator end ){
	int count = end - begin;
	settedW.resize( count );
	int i = 0;
	while( begin != end )
		settedW( i++ ) = *begin;

}

void SparsePCAForShapeGrad::reconstFormWsingle()
{
	Eigen::VectorXd component(3 * nF);
	component = W_single.transpose() * KComp;
	Eigen::VectorXd vertexCords(3 * nV );
		vertexCords = reconstrtuction(mesh_, 1 * component + Xmean, 0);
		//vertexCords = SpcaShapeGrad->reconstrtuction(SpcaShapeGrad->mesh_, x(index) * component + meanComp, 0);
		Eigen::VectorXd faceError (nF);
		for(int j = 0; j < nF; j++)
		{
			faceError(j) = component.segment(9 * j, 9).norm();
		}
		char *pPath = new char[30];
		//sprintf(pPath,"%6.2f _ %d-thDFComponent.obj",x(index), i);
		sprintf_s(pPath,30,"%reconstDGComponent.obj");
		writeFaceErrorOfComp(pPath, vertexCords, faceError);
}


void SparsePCAForShapeGrad::setAnchorPoints(std::vector<int> & _fixPointIdx, std::vector<Eigen::VectorXd> & _fixPointPos)
{
	fixPointIdx = _fixPointIdx;
	fixPointPos = _fixPointPos;
}
void SparsePCAForShapeGrad::setAnchorR(std::vector<int>& _fixFaceIdx, std::vector<Eigen::MatrixXd>& _fixFaceR)
{
	fixFaceIdx = _fixFaceIdx;
	fixFaceR = _fixFaceR;
}