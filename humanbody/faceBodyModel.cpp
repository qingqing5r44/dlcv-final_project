#include "faceBodyModel.h"

using namespace std;
using namespace Eigen;



FaceBodyModel::FaceBodyModel(int kS, int kP,  PGMesh *_mesh)
{
	spCM = new  SparseLocalizedForConnectMap(_mesh);
	spCM->calculateFacesNormals();
	spCM->calculateFacesFrame();
	spCM->computeDiEdge();
	spCM->presolve();
	k_Shapes = kS;
	k_Poses = kP;
	nV = _mesh->n_vertices();
	nE = _mesh->n_edges();
	nF = _mesh->n_faces();
	mean_shape = Eigen::VectorXd(3 * nV);
	shape_coeff = Eigen::VectorXd(k_Shapes);
	pose_coeff = Eigen::VectorXd(k_Poses);
	shape_coeff.setZero();
	pose_coeff.setZero();
	C_S = Eigen::MatrixXd(k_Shapes, 3 * nV);
	C_P = Eigen::MatrixXd(k_Poses, 2 * nE);
	ref_CM = Eigen::VectorXd(2 * nE);
	loadTensorflowModel();
}

Eigen::VectorXd FaceBodyModel::getVertexMatrix(PGMesh * _mesh)
{
	Eigen::VectorXd V(3 * nV);
	for (auto vit = _mesh->vertices_begin(); vit != _mesh->vertices_end(); ++vit)
	{
		PGMesh::Point p = _mesh->point(*vit);
		V(3 * vit->idx() + 0) = p[0];
		V(3 * vit->idx() + 1) = p[1];
		V(3 * vit->idx() + 2) = p[2];
	}
	return V;
}

Eigen::MatrixXd FaceBodyModel::getVertexMatrix1(PGMesh * _mesh)
{
	Eigen::MatrixXd V(_mesh->n_vertices(), 3);
	for (auto vit = _mesh->vertices_begin(); vit != _mesh->vertices_end(); ++vit)
	{
		PGMesh::Point p = _mesh->point(*vit);
		V(vit->idx() ,0) = p[0];
		V(vit->idx() , 1) = p[1];
		V(vit->idx() , 2) = p[2];
	}
	return V;
}

void FaceBodyModel::loadShapeBases(char* s_path)
{
	ifstream in(s_path);
	int i = 0, j = 0;
	while (!in.eof() && i < k_Shapes) {
		double tmp;
		in >> tmp;
		C_S(i, j) = tmp;
		//cout << C_S(i, j) << endl;
		++j;
		if (j % (3 * nV) == 0) {
			++i;
			j = 0;
		}
	}

	in.close();
}

void FaceBodyModel::loadPoseBases(char* p_path, int n)
{
	ifstream in(p_path);
	
	k_Poses = n;
	C_P = Eigen::MatrixXd(k_Poses, 2 * nE);
	int i = 0, j = 0;
	int iter = 0;
	while (!in.eof() && iter < k_Poses * nE * 2) {
		double tmp;
		in >> tmp;
		i = iter / (nE * 2);
		j = iter % (nE * 2);
		C_P(i, j) = tmp;
		iter++;
		//cout << i << " " <<j << endl;
	}

	in.close();
}

void FaceBodyModel::loadExpressionBases(char* p_path, int n)
{
	ifstream in(p_path);
	k_expressions = n;
	C_expressions = Eigen::MatrixXd(k_expressions, 2 * nE);
	int i = 0, j = 0;
	int iter = 0;
	while (!in.eof() && iter < k_expressions * nE * 2) {
		double tmp;
		in >> tmp;
		i = iter / (nE * 2);
		j = iter % (nE * 2);
		C_expressions(i, j) = tmp;
		iter++;
		//cout << i << " " <<j << endl;
	}
	in.close();
	Eigen::MatrixXd A_2 = C_expressions * (C_expressions.transpose());
	lltExpression = Eigen::LLT<MatrixXd>(A_2);
	return;
}

void FaceBodyModel::loadgestureBases(char* p_path, int n)
{
	ifstream in(p_path);
	k_gestures = n;
	C_gestures = Eigen::MatrixXd(k_gestures, 2 * nE);
	int i = 0, j = 0;
	int iter = 0;
	while (!in.eof() && iter < k_gestures * nE * 2) {
		double tmp;
		in >> tmp;
		i = iter / (nE * 2);
		j = iter % (nE * 2);
		C_gestures(i, j) = tmp;
		iter++;
		//cout << i << " " <<j << endl;
	}
	Eigen::MatrixXd A_2 = C_gestures * (C_gestures.transpose());
	lltGesture = Eigen::LLT<MatrixXd>(A_2);
	return;
}

void FaceBodyModel::optimizebodyfromde(Eigen::VectorXd &cm, Eigen::VectorXd &p_coef)
{
	Eigen::SparseMatrix<double> Coefficients(2 * nE + p_coef.size(), p_coef.size());
	Eigen::VectorXd y(2 * nE + p_coef.size());
	std::vector<Eigen::Triplet<double> > triple;
	for (int i = 0; i < nE; i++)
	{
		for (int j = 0; j < p_coef.size(); j++)
		{
			triple.push_back(Eigen::Triplet<double>(2 * i + 0, j, C_P(j, 2 * i + 0)));

			triple.push_back(Eigen::Triplet<double>(2 * i + 1, j, C_P(j, 2 * i + 1) * ref_CM(2 * i + 1)));
		}
		y(2 * i + 0) = cm(2 * i + 0) - ref_CM(2 * i + 0);
		y(2 * i + 1) = cm(2 * i + 1) - ref_CM(2 * i + 1);
	}
	for (int h = 0; h < p_coef.size(); h++)
	{
		triple.push_back(Eigen::Triplet<double>(2 * nE +  h , h, 3.0));
		y(2 * nE + h) = 3.0 * (p_coef(h));
	}
	Coefficients.setFromTriplets(triple.begin(), triple.end());
	SimplicialCholesky<SparseMatrix<double>> solver;
	solver.compute(Coefficients.transpose() * Coefficients);
	if (solver.info() != Success)
	{
		/// decomposit ion failed
		std::cout << "Decomposition failed" << std::endl;
		return;
	}
	p_coef = solver.solve(Coefficients.transpose() * y);
	cm = C_P.transpose() * p_coef;
	for (int i = 0; i < nE; i++)
	{
		cm(2 * i + 0) = ref_CM(2 * i + 0) + cm(2 * i + 0);
		cm(2 * i + 1) = ref_CM(2 * i + 1) * (1 + cm(2 * i + 1));
	}
	return;
}

void FaceBodyModel::optimizeexpressionfromde(Eigen::VectorXd &cm, Eigen::VectorXd &f_coef)
{
	Eigen::SparseMatrix<double> Coefficients(2 * nE, f_coef.size());
	Eigen::VectorXd y(2 * nE);
	std::vector<Eigen::Triplet<double> > triple;
	for (int i = 0; i < nE; i++)
	{
		for (int j = 0; j < f_coef.size(); j++)
		{
			triple.push_back(Eigen::Triplet<double>(2  * i + 0, j, C_expressions(j, 2 * i + 0)));

			triple.push_back(Eigen::Triplet<double>(2 * i + 1, j, C_expressions(j, 2 * i + 1) * ref_CM(2 * i + 1)));
		}
		y(2 * i + 0) = cm(2 * i + 0) - ref_CM(2 * i + 0);
		y(2 * i + 1) = cm(2 * i + 1) - ref_CM(2 * i + 1);
	}
	Coefficients.setFromTriplets(triple.begin(), triple.end());
	SimplicialCholesky<SparseMatrix<double>> solver;
	solver.compute(Coefficients.transpose() * Coefficients);
	if (solver.info() != Success)
	{
		/// decomposit ion failed
		std::cout << "Decomposition failed" << std::endl;
		return;
	}
	f_coef = solver.solve(Coefficients.transpose() * y);
	cm = C_expressions.transpose() * f_coef;
	for (int i = 0; i < nE; i++)
	{
		cm(2 * i + 0) = ref_CM(2 * i + 0) + cm(2 * i + 0);
		cm(2 * i + 1) = ref_CM(2 * i + 1) * (1 + cm(2 * i + 1));
	}
	return;
}

void FaceBodyModel::optimizegesturefromde(Eigen::VectorXd &cm, Eigen::VectorXd &g_coef)
{
	Eigen::SparseMatrix<double> Coefficients(2 * nE, g_coef.size());
	Eigen::VectorXd y(2 * nE);
	std::vector<Eigen::Triplet<double> > triple;
	for (int i = 0; i < nE; i++)
	{
		for (int j = 0; j < g_coef.size(); j++)
		{
			triple.push_back(Eigen::Triplet<double>(2 * i + 0, j, C_gestures(j, 2 * i + 0)));

			triple.push_back(Eigen::Triplet<double>(2 * i + 1, j, C_gestures(j, 2 * i + 1) * ref_CM(2 * i + 1)));
		}
		y(2 * i + 0) = cm(2 * i + 0) - ref_CM(2 * i + 0);
		y(2 * i + 1) = cm(2 * i + 1) - ref_CM(2 * i + 1);
	}
	Coefficients.setFromTriplets(triple.begin(), triple.end());
	SimplicialCholesky<SparseMatrix<double>> solver;
	solver.compute(Coefficients.transpose() * Coefficients);
	if (solver.info() != Success)
	{
		/// decomposit ion failed
		std::cout << "Decomposition failed" << std::endl;
		return;
	}
	g_coef = solver.solve(Coefficients.transpose() * y);
	cm = C_gestures.transpose() * g_coef;
	return;
}

void FaceBodyModel::loadMeanShape(char* m_path)
{
	ifstream in(m_path);

	int iter = 0;
	while (!in.eof() && iter < 3 * nV) {
		double tmp;
		in >> tmp;
		
		mean_shape(iter) = tmp;
		iter++;
		//cout << i << " " <<j << endl;
	}

	in.close();
}

void FaceBodyModel::outputBlendshapes()
{
	for (int i = 0; i < k_Shapes; i++)
	{
		char *compPath = new char[100];
		sprintf_s(compPath, 100, "shape%04d.obj", i);
		Eigen::VectorXd l = mean_shape +3 *  C_S.row(i).transpose();
		Eigen::VectorXd l1(3 * nV);
		spCM->reconstructionFromDiEdges(l, l1, 0, 0.5, 0);
		spCM->outMesh(l1, compPath);
	}
}


//读取body feature的idx和weights
void FaceBodyModel::readBodyFeatureidx(string p_path)
{
	body_feature_idx.resize(28);
	body_feature_weights.resize(28);
	ifstream file(p_path);
	int featureidx, vertexidx;
	double featureidx1, vertexidx1;
	double weight;
	file >> featureidx1;
	featureidx = int(featureidx1);
	while (featureidx != -1)
	{
		file >> vertexidx1;
		vertexidx = int(vertexidx1);
		file >> weight;
		body_feature_idx[featureidx].push_back(vertexidx);
		body_feature_weights[featureidx].push_back(weight);
		file >> featureidx1;
		featureidx = int(featureidx1);
	}
	return;
}

//读取body feature的idx和weights
void FaceBodyModel::readBodyFeatureidx1(string p_path)
{
	body_feature_idx.resize(33);
	body_feature_weights.resize(33);
	ifstream file(p_path);
	int featureidx, vertexidx;
	double featureidx1, vertexidx1;
	double weight;
	file >> featureidx1;
	featureidx = int(featureidx1);
	while (featureidx != -1)
	{
		file >> vertexidx1;
		vertexidx = int(vertexidx1);
		file >> weight;
		body_feature_idx[featureidx].push_back(vertexidx);
		body_feature_weights[featureidx].push_back(weight);
		file >> featureidx1;
		featureidx = int(featureidx1);
	}
	return;
}

void FaceBodyModel::outputPoseBases()
{
//	std::vector<int> fixPointIdx;
//	std::vector<Eigen::VectorXd> fixPointPos;
//	std::vector<int> fixFaceIdx;
//	std::vector<Eigen::MatrixXd> fixFaceR;
//	fixPointIdx.resize(1);
//	fixPointPos.resize(1);
//	fixFaceIdx.resize(1);
//	fixFaceR.resize(1);
//	for (int i = 0; i < 1; i++)
//	{
//		fixPointIdx[i] = i;
//		fixPointPos[i] = spGrad->SpatialTemData[i][0];
//	}
//	for (int i = 0; i < 1; i++)
//	{
//		fixFaceIdx[i] = i;
//		Eigen::MatrixXd R(3, 3);
//		R << 1, 0, 0,
//			0, 1, 0,
//			0, 0, 1;
//		fixFaceR[i] = R;
//	}
//	spGrad->setAnchorPoints(fixPointIdx, fixPointPos);
//	spGrad->setAnchorR(fixFaceIdx, fixFaceR);
//	spGrad->computeS_matrix();
//	spGrad->build_matrix();
////#pragma omp parallel for
//	for (int i = 0 ; i < k_Poses; i ++)
//	{
//		std::vector<Eigen::MatrixXd> R_ij;
//		std::vector<Eigen::MatrixXd> S_i;
//		R_ij.resize(spGrad->nE);
//		S_i.resize(spGrad->nF);
//		Eigen::VectorXd edge_comp = C_P.row(i);
//		spGrad->vecToMat(edge_comp, R_ij, S_i);
//		Eigen::VectorXd vert = spGrad->reconstructFromRS(spGrad->mesh_, R_ij, S_i);
//		char *meshpath1 = new char[100];
//		sprintf_s(meshpath1, 100, "Pose%04d.obj", i);
//		Eigen::VectorXd edgeError(nE);
//		double errors=0;
//		for (int i = 0; i < nE; i++)
//		{
//			//算第i帧与第0帧的差
//			double error = edge_comp.segment(9 * i, 9).norm();
//			edgeError(i) = error;
//			errors += error;
//		}
//		cout<< errors<<endl;
//		spGrad->visualizeComp(vert, edgeError, meshpath1);
//	}
}

Eigen::VectorXd FaceBodyModel::generateShapeWithoutPose(string shapepath)
{
	PGMesh *shapeMesh = new PGMesh();
	OpenMesh::IO::read_mesh(*shapeMesh, shapepath);
	Eigen::VectorXd shape = spCM->getVertexMatrix(shapeMesh);
	return shape;
}

Eigen::VectorXd FaceBodyModel::generateShapeWithoutPose(Eigen::VectorXd s_coeff)
{
	return mean_shape + C_S.transpose() * s_coeff;
}
		
void FaceBodyModel::generateShapeWithoutPose()
{
	shape_without_pose = mean_shape + C_S.transpose() * shape_coeff;
	return;
}

Eigen::MatrixXd FaceBodyModel::f_x(Eigen::VectorXd s_coeff, Eigen::VectorXd pose_parm, Eigen::VectorXd f_parm, Eigen::VectorXd g_parm)
{
	Eigen::VectorXd shape = generateShapeWithoutPose(s_coeff);
	spCM->outMesh1(shape, "refmesh.off");
	//std::ofstream fshape("shape.txt");
	//fshape << shape << endl;
	//fshape.close();
	updateRef(shape);
	std::vector<double> faceAnchor = getAnchors(s_coeff, pose_parm);
	//spCM->updateDeformed1(faceAnchor);
	getAnchorPos(s_coeff, pose_parm);
	//spCM->setAnchor(body_feature_idx, body_feature_weights, anchorpointpos);
	spCM->presolve();
	Eigen::VectorXd cm = getPoseDE(s_coeff, pose_parm);
	Eigen::VectorXd cm_face = C_expressions.transpose() * f_parm;
	Eigen::VectorXd cm_gesture = C_gestures.transpose() * g_parm;
	cm = cm + cm_face + cm_gesture;
	Eigen::VectorXd cm1(2 * nE);
	//std::cout << "cm" << cm(0) << cm(1) << cm(2) << std::endl;
	std::ofstream frefCM("ref_CM.txt");
	//ref_CM = spCM->DiEdgeDataMatrix.row(0);
	frefCM << ref_CM << endl;
	frefCM.close();
	for (int i = 0; i < nE; i++)
	{
		cm1(2 * i + 0) = ref_CM(2 * i + 0) + cm(2 * i + 0);
		cm1(2 * i + 1) = ref_CM(2 * i + 1) * (1 + cm(2 * i + 1));

	}
	std::ofstream fcm1("cm1.txt");
	fcm1 << cm1 << endl;
	fcm1.close();
	Eigen::VectorXd l1(3 * nV);

	spCM->reconstructionFromDiEdges(cm1, l1, 1, 0.5, 1);
	spCM->outMesh1(l1, "out_mesh.off");
	Eigen::MatrixXd L(nV, 3);
	for (int i = 0; i < nV; i++)
	{
		L(i, 0) = l1(3 * i + 0);
		L(i, 1) = l1(3 * i + 1);
		L(i, 2) = l1(3 * i + 2);
	}
	return L;
}
Eigen::VectorXd FaceBodyModel::getDE(string path)
{
	ifstream input(path);
	Eigen::VectorXd x(65380);
	//Eigen::VectorXd x(41328);
	for (int i = 0; i < x.size(); i++)
	{
		input >> x(i);
	}
	input.close();
	return x;
}

Eigen::MatrixXd FaceBodyModel::f_x(string shape_path, Eigen::VectorXd s_coeff, Eigen::VectorXd pose_parm)
{
	PGMesh *shapeMesh = new PGMesh();
	OpenMesh::IO::read_mesh(*shapeMesh, shape_path);
	Eigen::VectorXd shape = spCM->getVertexMatrix(shapeMesh);
	delete shapeMesh;
	spCM->outMesh1(shape, "refmesh.off");
	updateRef(shape);
	spCM->presolve();
	Eigen::VectorXd cm = getPoseDE(s_coeff, pose_parm);
	Eigen::VectorXd cm1(2 * nE);
	std::ofstream frefCM("ref_CM.txt");
	frefCM << ref_CM << endl;
	frefCM.close();
	for (int i = 0; i < nE; i++)
	{
		cm1(2 * i + 0) = ref_CM(2 * i + 0) + cm(2 * i + 0);
		cm1(2 * i + 1) = ref_CM(2 * i + 1) * (1 + cm(2 * i + 1));
	}
	std::ofstream fcm1("cm1.txt");
	fcm1 << cm1 << endl;
	fcm1.close();
	Eigen::VectorXd l1(3 * nV);

	spCM->reconstructionFromDiEdges(cm1, l1, 1, 0.5, 1);
	spCM->outMesh1(l1, "out_mesh.off");
	Eigen::MatrixXd L(nV, 3);
	for (int i = 0; i < nV; i++)
	{
		L(i, 0) = l1(3 * i + 0);
		L(i, 1) = l1(3 * i + 1);
		L(i, 2) = l1(3 * i + 2);
	}
	return L;
}

Eigen::MatrixXd FaceBodyModel::f_x(string shape_path, string de_path)
{
	PGMesh *shapeMesh = new PGMesh();
	OpenMesh::IO::read_mesh(*shapeMesh, shape_path);
	Eigen::VectorXd shape = spCM->getVertexMatrix(shapeMesh);
	delete shapeMesh;
	spCM->outMesh1(shape, "refmesh.off");
	updateRef(shape);
	spCM->presolve();
	Eigen::VectorXd cm = getDE(de_path);
	Eigen::VectorXd cm1(2 * nE);
	std::ofstream frefCM("ref_CM.txt");
	frefCM << ref_CM << endl;
	frefCM.close();
	for (int i = 0; i < nE; i++)
	{
		cm1(2 * i + 0) = ref_CM(2 * i + 0) + cm(2 * i + 0);
		cm1(2 * i + 1) = ref_CM(2 * i + 1) * (1 + cm(2 * i + 1));
	}
	std::ofstream fcm1("cm1.txt");
	fcm1 << cm1 << endl;
	fcm1.close();
	Eigen::VectorXd l1(3 * nV);

	spCM->reconstructionFromDiEdges(cm1, l1, 1, 0.5, 1);
	spCM->outMesh1(l1, "out_mesh.off");
	Eigen::MatrixXd L(nV, 3);
	for (int i = 0; i < nV; i++)
	{
		L(i, 0) = l1(3 * i + 0);
		L(i, 1) = l1(3 * i + 1);
		L(i, 2) = l1(3 * i + 2);
	}
	return L;
}

Eigen::MatrixXd FaceBodyModel::f_x(Eigen::VectorXd s_coeff, Eigen::VectorXd pose_parm)
{
	Eigen::VectorXd shape = generateShapeWithoutPose(s_coeff);
	spCM->outMesh1(shape, "refmesh.off");
	//std::ofstream fshape("shape.txt");
	//fshape << shape << endl;
	//fshape.close();
	updateRef(shape);
	//std::vector<double> faceAnchor = getAnchors(s_coeff, pose_parm);
	//spCM->updateDeformed1(faceAnchor);
	//getAnchorPos(s_coeff, pose_parm);
	//spCM->setAnchor(body_feature_idx, body_feature_weights, anchorpointpos);
	spCM->presolve();
	Eigen::VectorXd cm = getPoseDE(s_coeff, pose_parm);
	Eigen::VectorXd cm1(2 * nE);
	std::ofstream frefCM("ref_CM.txt");
	//ref_CM = spCM->DiEdgeDataMatrix.row(0);
	frefCM << ref_CM << endl;
	frefCM.close();
	for (int i = 0; i < nE; i++)
	{
		cm1(2 * i + 0) = ref_CM(2 * i + 0) + cm(2 * i + 0);
		cm1(2 * i + 1) = ref_CM(2 * i + 1) * (1 + cm(2 * i + 1));
	}
	std::ofstream fcm1("cm1.txt");
	fcm1 << cm1 << endl;
	fcm1.close();
	Eigen::VectorXd l1(3 * nV);

	spCM->reconstructionFromDiEdges(cm1, l1, 1, 0.5, 1);
	spCM->outMesh1(l1, "out_mesh.off");
	Eigen::MatrixXd L (nV, 3);
	for (int i = 0; i < nV; i++)
	{
		L(i, 0) = l1(3 * i + 0);
		L(i, 1) = l1(3 * i + 1);
		L(i, 2) = l1(3 * i + 2);
	}
	return L;
}

void FaceBodyModel::initFace(PGMesh *pg)
{
	spCM = new SparseLocalizedForConnectMap(pg);
	spCM->addAnchor(3050);
	spCM->addAnchor(3480);
	spCM->addAnchor(4811);
	spCM->addAnchor(4818);
	spCM->addAnchor(1103);
	spCM->addAnchor(877);
	spCM->addAnchor(5028);
	spCM->addAnchor(4952);
	spCM->addAnchor(1396);
	spCM->addAnchor(1578);
	spCM->addAnchor(1969);
	spCM->addAnchor(5642);
	spCM->addAnchor(8816);
	spCM->addAnchor(8771);
	spCM->addAnchor(5539);
	spCM->addAnchor(2243);
	spCM->calculateFacesNormals();
	spCM->calculateFacesFrame();
	spCM->computeDiEdge();
	spCM->presolve();

}

Eigen::MatrixXd FaceBodyModel::f_expression(Eigen::VectorXd f_parm)
{
	Eigen::VectorXd c1 = C_expressions.transpose() * f_parm;
	for (int i = 0; i < nE; i++)
	{
		c1(2 * i + 0) = c1(2 * i + 0) + ref_CM(2 * i + 0);
		c1(2 * i + 1) = (c1(2 * i + 1) + 1) * ref_CM(2 * i + 1);
	}
	Eigen::VectorXd l1(3 * spCM->nV);
	spCM->reconstructionFromDiEdges(c1, l1, 1, 0.5, 1);
	Eigen::MatrixXd V(nV, 3);
	for (int i = 0; i < nV; i++)
	{
		V.row(i) = l1.segment(3 * i + 0, 3);
	}
	return V;
}

Eigen::MatrixXd FaceBodyModel::f_gesture(Eigen::VectorXd g_parm)
{
	Eigen::VectorXd c1 = C_gestures.transpose() * g_parm;
	for (int i = 0; i < nE; i++)
	{
		c1(2 * i + 0) = c1(2 * i + 0) + ref_CM(2 * i + 0);
		c1(2 * i + 1) = (c1(2 * i + 1) + 1) * ref_CM(2 * i + 1);
	}
	Eigen::VectorXd l1(3 * spCM->nV);
	spCM->reconstructionFromDiEdges(c1, l1, 1, 0.5, 1);
	Eigen::MatrixXd V(nV, 3);
	for (int i = 0; i < nV; i++)
	{
		V.row(i) = l1.segment(3 * i + 0, 3);
	}
	return V;
}

Eigen::VectorXd FaceBodyModel::getFacePosfromCurrentV(Eigen::MatrixXd V)
{
	Eigen::VectorXd Vt(face_feature_idx.size() * 3);
	Vt.setZero();
	for (int i = 0; i < face_feature_idx.size(); i++)
	{
		int bf_idx = i;
		for (int j = 0; j < face_feature_idx[bf_idx].size(); j++)
		{
			Vt.segment(3 * i + 0, 3) += face_feature_weights[bf_idx][j] * V.row(face_feature_idx[bf_idx][j]);
		}
	}
	return Vt;
}

Eigen::VectorXd FaceBodyModel::getBodyPosfromCurrentV(Eigen::MatrixXd V, std::vector<std::vector<int>> idxset, std::vector<std::vector<double>> weightset)
{
	Eigen::VectorXd Vt(idxset.size() * 3);
	Vt.setZero();
	for (int i = 0; i < idxset.size(); i++)
	{
		
		for (int j = 0; j < idxset[i].size(); j++)
		{
			Vt.segment(3 * i + 0, 3) += weightset[i][j] * V.row(idxset[i][j]);
		}
	}
	return Vt;
}

Eigen::VectorXd FaceBodyModel::getBodyPosfromCurrentV(Eigen::MatrixXd V, std::vector<int> idxset)
{
	Eigen::VectorXd Vt(idxset.size() * 3);
	Vt.setZero();
	for (int i = 0; i < idxset.size(); i++)
	{
		int bf_idx = selectbodyfeature[ idxset[i]];
		for (int j = 0; j < body_feature_idx[bf_idx].size(); j++)
		{
			Vt.segment(3 * i + 0, 3) += body_feature_weights[bf_idx][j] * V.row(body_feature_idx[bf_idx][j]);
		}
	}
	return Vt;
}

Eigen::VectorXd FaceBodyModel::getBodyPosfromCurrentV1(Eigen::MatrixXd V, std::vector<int> idxset)
{
	Eigen::VectorXd Vt(idxset.size() * 3);
	for (int i = 0; i < idxset.size(); i++)
	{
		int bf_idx = idxset[i];
		
		Vt.segment(3 * i + 0, 3) =  V.row(bf_idx);
		
	}
	return Vt;
}

void FaceBodyModel::updateRef(Eigen::VectorXd x)
{
	//ref_CM = spCM->computeDiEdge(x);
	spCM->setMesh(x);
	spCM->calculateFacesNormals();
	spCM->calculateFacesFrame();
	spCM->computeDiEdge();
	ref_CM = spCM->DiEdgeDataMatrix.row(0);
}

std::vector<double> FaceBodyModel::getAnchors(Eigen::VectorXd s_coeff, Eigen::VectorXd motion_parm)
{
	std::vector<double> y(12 * 9);
	//定义入参和出参
	PyObject* pReturnValue;
	//两个入参分别为 x,y coordinate of each point
	PyObject* pyListX = PyList_New(s_coeff.size() + motion_parm.size());
	//#pragma omp parallel for
	for (int i = 0; i<s_coeff.size(); i++)
	{
		PyList_SetItem(pyListX, i, PyFloat_FromDouble(s_coeff(i)));
	}
	//#pragma omp parallel for
	for (int i = 0; i < motion_parm.size(); i++)
	{
		PyList_SetItem(pyListX, s_coeff.size() + i, PyFloat_FromDouble(motion_parm(i)));
	}

	PyTuple_SetItem(inputmotionfaceargs, 0, pyListX);
	pReturnValue = PyObject_CallObject(pFuncMotion2Face, inputmotionfaceargs);
	PyObject *it;
	//#pragma omp parallel for
	for (int i = 0; i < 12* 9; i++)
	{
		it = PyList_GetItem(pReturnValue, i);
		PyArg_Parse(it, "d", &(y[i]));
	}
	return y;
}

void FaceBodyModel::loadDenseCorrespondenceModule()
{
	pModule_coresspondence = NULL;
	pModule_coresspondence= PyImport_ImportModule("densecorresponding_face");
	pFuncLoadPointcloud = NULL;
	pFuncFindCoress = NULL;
	pointcloud = NULL;
	PyObject *pDict = PyModule_GetDict(pModule_coresspondence);
	pFuncLoadPointcloud  = PyDict_GetItemString(pDict, "load_pointcloud");
	pFuncFindCoress = PyDict_GetItemString(pDict, "find_coress_pairs");
}

void FaceBodyModel::load_pointcloud(int seqno)
{
	//定义入参和出参
	PyObject* pReturnValue;
	//两个入参分别为 seqlist
	PyObject* pyListX = PyList_New(1);
	PyList_SetItem(pyListX, 0, PyLong_FromLong(seqno));
	PyObject* pyList = PyTuple_New(1);
	PyTuple_SetItem(pyList, 0, pyListX);
	pReturnValue = PyObject_CallObject(pFuncLoadPointcloud, pyList);
	pointcloud = pReturnValue;
}

void FaceBodyModel::findCoressPairs(int seqno, int ith, int jth)
{
	//定义入参和出参
	PyObject* pReturnValue;
	//两个入参分别为 seqlist
	PyObject* pyListX = PyList_New(3);
	PyList_SetItem(pyListX, 0, PyLong_FromLong(seqno));
	PyList_SetItem(pyListX, 1, PyLong_FromLong(ith));
	PyList_SetItem(pyListX, 2, PyLong_FromLong(jth));
	PyObject* pyList = PyTuple_New(2);
	PyTuple_SetItem(pyList, 0, pyListX);
	PyTuple_SetItem(pyList, 1, pointcloud);
	pReturnValue = PyObject_CallObject(pFuncFindCoress, pyList);
}

void FaceBodyModel::loadTensorflowModel()
{
	Py_SetPythonHome(GetWC("F:/Anaconda3/envs/python37"));//D:/Program Files (x86)/Microsoft Visual Studio/Shared/Anaconda3_64/envs/tensorflow2
	Py_Initialize();//调用Py_Initialize()进行初始化
	PyRun_SimpleString("import sys");
	PyRun_SimpleString("sys.path.append('E:/Graphics/bachelor_thesis/Code/humanbody/humanbody/')");
	PyRun_SimpleString("sys.path.append('./')");
	pmaleModelMotion = NULL;
	pFuncMotion2DE = NULL;
	pModule2 = NULL;
	pFunc2 = NULL;
	sessPose = NULL;
	sessFace = NULL;
	sessHand = NULL;
	C_body = NULL;
	m2c_model = NULL;//the motion parameter to pose coefficient model
	ss_x = NULL;
	ss_y = NULL;
	bodyDEOutput = NULL;
	sessPoseMotion = NULL;
	tf_body_motion_x = NULL;
	inputmotionfaceargs = NULL;
	pFuncMotion2Face = NULL;
	//import math_test
	//pmaleModelMotion = PyImport_ImportModule("male_body_generate_DFAUST");
	pmaleModelMotion = PyImport_ImportModule("female_body_generate");
	//PyRun_SimpleString("import sys");
	//PyRun_SimpleString("import numpy as np");
	//PyRun_SimpleString("sys.path.append('./')");
	//PyRun_SimpleString("import tensorflow as tf");
	//PyRun_SimpleString("from sklearn.model_selection import train_test_split");
	//PyRun_SimpleString("from sklearn import preprocessing");
	//PyRun_SimpleString("hhh = np.load('DE_training_data.npz')");
	//PyRun_SimpleString("Input_ = hhh['Input_']");
	//PyRun_SimpleString("print(Input_)");
	PyObject *pDict = PyModule_GetDict(pmaleModelMotion);
	//对应math_test.py中的def add_func(a,b)函数
	pFucLoadbodymotion = PyDict_GetItemString(pDict, "load_body_motion");
	pFuncPose2DE = PyObject_GetAttrString(pmaleModelMotion, "recons");
	pFuncMotion2Face = PyObject_GetAttrString(pmaleModelMotion, "reconsFace");
	pFuncMotion2Pose = PyObject_GetAttrString(pmaleModelMotion, "outputPoseParm");
	pFuncPose2skeletonpos = PyObject_GetAttrString(pmaleModelMotion, "outputSkeletonPos");
	pFuncOptimizex = PyObject_GetAttrString(pmaleModelMotion, "optimizex");
	pFuncFindRigid = PyObject_GetAttrString(pmaleModelMotion, "findRigidtransform");
	pFuncgetJointsPos = PyObject_GetAttrString(pmaleModelMotion, "getJointsPos");
	PyObject * pReturn = NULL;
	//PyRun_SimpleString("from sklearn.externals import joblib");
	//PyRun_SimpleString("from sklearn import preprocessing");
	//PyRun_SimpleString("from sklearn.model_selection import train_test_split");
	//PyRun_SimpleString("import numpy as np");
	pReturn = PyObject_CallObject(pFucLoadbodymotion, NULL);
	if (!pReturn) {

		printf("import function result failed!!\n");

		return;
	}
	sessPose = PyTuple_GetItem(pReturn, 0);
	//m2c_model = PyTuple_GetItem(pReturn, 1);
	//ss_x = PyTuple_GetItem(pReturn, 2);
	//ss_y = PyTuple_GetItem(pReturn, 3);
	tf_body_motion_x = PyTuple_GetItem(pReturn, 1);
	//C_body = PyTuple_GetItem(pReturn, 5);
	bodyDEOutput = PyTuple_GetItem(pReturn, 2);
	
	inputmotionargs = PyTuple_New(4);
	PyTuple_SetItem(inputmotionargs, 1, tf_body_motion_x);
	PyTuple_SetItem(inputmotionargs, 2, bodyDEOutput);
	PyTuple_SetItem(inputmotionargs, 3, sessPose);

	inputmotionfaceargs = PyTuple_New(1);
	 
	inputmotionargs1 = PyTuple_New(5);
	PyTuple_SetItem(inputmotionargs1, 2, tf_body_motion_x);
	PyTuple_SetItem(inputmotionargs1, 3, bodyDEOutput);
	PyTuple_SetItem(inputmotionargs1, 4, sessPose);
}

Eigen::VectorXd FaceBodyModel::getPoseDE(Eigen::VectorXd s_coeff, Eigen::VectorXd pose_parm)
{

	Eigen::VectorXd y(2 * nE);
	//定义入参和出参
	PyObject* pReturnValue;
	//两个入参分别为 x,y coordinate of each point
	PyObject* pyListX = PyList_New(s_coeff.size() + pose_parm.size());
//#pragma omp parallel for
	for (int i = 0; i<s_coeff.size(); i++)
	{
		PyList_SetItem(pyListX, i, PyFloat_FromDouble(s_coeff(i)));
	}
//#pragma omp parallel for
	for (int i = 0; i < pose_parm.size(); i++)
	{
		PyList_SetItem(pyListX, s_coeff.size() + i, PyFloat_FromDouble(pose_parm(i)));
	}

	PyTuple_SetItem(inputmotionargs, 0, pyListX);
	pReturnValue = PyObject_CallObject(pFuncPose2DE, inputmotionargs);
	PyObject *it;
//#pragma omp parallel for
	for (int i = 0; i < 2 * nE; i++)
	{
		it = PyList_GetItem(pReturnValue, i);
		PyArg_Parse(it, "d", &(y(i)));
	}
	return y;
}

Eigen::VectorXd FaceBodyModel::getPoseParm(Eigen::VectorXd s_coeff, Eigen::VectorXd motion_parm)
{
	Eigen::VectorXd y(400);
	//定义入参和出参
	PyObject* pReturnValue;
	//两个入参分别为 x,y coordinate of each point
	PyObject* pyListX1 = PyList_New(s_coeff.size());
	PyObject* pyListX2 = PyList_New( motion_parm.size());
	//#pragma omp parallel for
	for (int i = 0; i<s_coeff.size(); i++)
	{
		PyList_SetItem(pyListX1, i, PyFloat_FromDouble(s_coeff(i)));
	}
	//#pragma omp parallel for
	for (int i = 0; i < motion_parm.size(); i++)
	{
		PyList_SetItem(pyListX2, i, PyFloat_FromDouble(motion_parm(i)));
	}
	PyObject* pyList = PyTuple_New(2);
	PyTuple_SetItem(pyList, 0, pyListX1);
	PyTuple_SetItem(pyList, 1, pyListX2);
	pReturnValue = PyObject_CallObject(pFuncMotion2Pose, pyList);
	PyObject *it;
	for (int i = 0; i < 400; i++)
	{
		it = PyList_GetItem(pReturnValue, i);
		PyArg_Parse(it, "d", &(y(i)));
	}
	return y;
}

std::vector<Eigen::Vector3d> FaceBodyModel::getJointsPos(Eigen::VectorXd s_coeff, Eigen::VectorXd motion_parm)
{
	Eigen::VectorXd y(72);
	std::vector<Eigen::Vector3d> out(24);
	//定义入参和出参
	PyObject* pReturnValue;
	//两个入参分别为 x,y coordinate of each point
	PyObject* pyListX1 = PyList_New(s_coeff.size());
	PyObject* pyListX2 = PyList_New(motion_parm.size());
	//#pragma omp parallel for
	for (int i = 0; i<s_coeff.size(); i++)
	{
		PyList_SetItem(pyListX1, i, PyFloat_FromDouble(s_coeff(i)));
	}
	//#pragma omp parallel for
	for (int i = 0; i < motion_parm.size(); i++)
	{
		PyList_SetItem(pyListX2, i, PyFloat_FromDouble(motion_parm(i)));
	}
	PyObject* pyList = PyTuple_New(2);
	PyTuple_SetItem(pyList, 0, pyListX1);
	PyTuple_SetItem(pyList, 1, pyListX2);
	pReturnValue = PyObject_CallObject(pFuncgetJointsPos, pyList);
	PyObject *it;
	for (int i = 0; i < 72; i++)
	{
		it = PyList_GetItem(pReturnValue, i);
		PyArg_Parse(it, "d", &(y(i)));
	}
	for (int i = 0; i < 24; i++)
	{
		out[i] = y.segment(3 * i, 3);
	}
	return out;
}


void FaceBodyModel::optimizexfromde(Eigen::VectorXd &cm, Eigen::VectorXd s_coef, Eigen::VectorXd &p_coef)
{
	//定义入参和出参
	PyObject* pReturnValue;
	//两个入参分别为 x,y coordinate of each point
	PyObject* pyListX1 = PyList_New(s_coef.size() + p_coef.size());
	PyObject* pyListX2 = PyList_New(2 * nE);

	for (int i = 0; i<s_coef.size(); i++)
	{
		PyList_SetItem(pyListX1, i, PyFloat_FromDouble(s_coef(i)));
	}
	//#pragma omp parallel for
	for (int i = 0; i < p_coef.size(); i++)
	{
		PyList_SetItem(pyListX1, s_coef.size() +  i, PyFloat_FromDouble(p_coef(i)));
	}
	for (int i = 0; i < 2 * nE; i++)
	{
		PyList_SetItem(pyListX2, i, PyFloat_FromDouble(cm(i)));
	}
	
	PyTuple_SetItem(inputmotionargs1, 0, pyListX1);
	PyTuple_SetItem(inputmotionargs1, 1, pyListX2);

	pReturnValue = PyObject_CallObject(pFuncOptimizex, inputmotionargs1);
	PyObject *out1, *out2;
	out1 = PyTuple_GetItem(pReturnValue, 0);
	out2 = PyTuple_GetItem(pReturnValue, 1);
	for (int i = 0; i < p_coef.size(); i++)
	{
		PyObject *it;
		it = PyList_GetItem(out1, s_coef.size() + i);
		PyArg_Parse(it, "d", &(p_coef(i)));
	}
}

void FaceBodyModel::findRigidTransform(Eigen::VectorXd Vt, Eigen::Matrix3d & R, Eigen::Vector3d &t, Eigen::MatrixXd layerpos)
{
	//定义入参和出参
	PyObject* pReturnValue;
	//两个入参分别为 x,y coordinate of each point
	PyObject* pyListX1 = PyList_New(Vt.size());
	PyObject* pyListX2 = PyList_New(9);
	PyObject* pyListX3 = PyList_New(3);
	PyObject* pyListX4 = PyList_New(Vt.size());
	PyList_SetItem(pyListX2, 0, PyFloat_FromDouble(R(0, 0)));
	PyList_SetItem(pyListX2, 1, PyFloat_FromDouble(R(0, 1)));
	PyList_SetItem(pyListX2, 2, PyFloat_FromDouble(R(0, 2)));
	PyList_SetItem(pyListX2, 3, PyFloat_FromDouble(R(1, 0)));
	PyList_SetItem(pyListX2, 4, PyFloat_FromDouble(R(1, 1)));
	PyList_SetItem(pyListX2, 5, PyFloat_FromDouble(R(1, 2)));
	PyList_SetItem(pyListX2, 6, PyFloat_FromDouble(R(2, 0)));
	PyList_SetItem(pyListX2, 7, PyFloat_FromDouble(R(2, 1)));
	PyList_SetItem(pyListX2, 8, PyFloat_FromDouble(R(2, 2)));
	PyList_SetItem(pyListX3, 0, PyFloat_FromDouble(t(0)));
	PyList_SetItem(pyListX3, 1, PyFloat_FromDouble(t(1)));
	PyList_SetItem(pyListX3, 2, PyFloat_FromDouble(t(2)));
	for (int i = 0; i < Vt.size()/3; i++)
	{
		PyList_SetItem(pyListX1, 3 * i + 0, PyFloat_FromDouble(Vt(3 * i + 0)));
		PyList_SetItem(pyListX1, 3 * i + 1, PyFloat_FromDouble(Vt(3 * i + 1)));
		PyList_SetItem(pyListX1, 3 * i + 2, PyFloat_FromDouble(Vt(3 * i + 2)));
		PyList_SetItem(pyListX4, 3 * i + 0, PyFloat_FromDouble(layerpos(i,0)));
		PyList_SetItem(pyListX4, 3 * i + 1, PyFloat_FromDouble(layerpos(i, 1)));
		PyList_SetItem(pyListX4, 3 * i + 2, PyFloat_FromDouble(layerpos(i, 0)));
	}
	PyObject* pyList = PyTuple_New(4);
	PyTuple_SetItem(pyList, 0, pyListX1);
	PyTuple_SetItem(pyList, 1, pyListX2);
	PyTuple_SetItem(pyList, 2, pyListX3);
	PyTuple_SetItem(pyList, 3, pyListX4);
	pReturnValue = PyObject_CallObject(pFuncFindRigid, pyList);
	PyObject *out1, *out2;
	out1 = PyTuple_GetItem(pReturnValue, 0);
	out2 = PyTuple_GetItem(pReturnValue, 1);
	PyObject *it;
	it = PyList_GetItem(out1, 0);
	PyArg_Parse(it, "d", &(R(0, 0)));
	it = PyList_GetItem(out1, 1);
	PyArg_Parse(it, "d", &(R(0, 1)));
	it = PyList_GetItem(out1, 2);
	PyArg_Parse(it, "d", &(R(0, 2)));
	it = PyList_GetItem(out1, 3);
	PyArg_Parse(it, "d", &(R(1, 0)));
	it = PyList_GetItem(out1, 4);
	PyArg_Parse(it, "d", &(R(1, 1)));
	it = PyList_GetItem(out1, 5);
	PyArg_Parse(it, "d", &(R(1, 2)));
	it = PyList_GetItem(out1, 6);
	PyArg_Parse(it, "d", &(R(2, 0)));
	it = PyList_GetItem(out1, 7);
	PyArg_Parse(it, "d", &(R(2, 1)));
	it = PyList_GetItem(out1, 8);
	PyArg_Parse(it, "d", &(R(2, 2)));
	it = PyList_GetItem(out2, 0);
	PyArg_Parse(it, "d", &(t(0)));
	it = PyList_GetItem(out2, 1);
	PyArg_Parse(it, "d", &(t(1)));
	it = PyList_GetItem(out2, 2);
	PyArg_Parse(it, "d", &(t(2)));
}

void FaceBodyModel::findRigidTransform(Eigen::VectorXd Vt, Eigen::Matrix3d & R, Eigen::Vector3d &t, std::vector<Eigen::Vector3d> layerpos)
{
	//定义入参和出参
	PyObject* pReturnValue;
	//两个入参分别为 x,y coordinate of each point
	PyObject* pyListX1 = PyList_New(Vt.size());
	PyObject* pyListX2 = PyList_New(9);
	PyObject* pyListX3 = PyList_New(3);
	PyObject* pyListX4 = PyList_New(Vt.size());
	PyList_SetItem(pyListX2, 0, PyFloat_FromDouble(R(0,0)));
	PyList_SetItem(pyListX2, 1, PyFloat_FromDouble(R(0, 1)));
	PyList_SetItem(pyListX2, 2, PyFloat_FromDouble(R(0, 2)));
	PyList_SetItem(pyListX2, 3, PyFloat_FromDouble(R(1, 0)));
	PyList_SetItem(pyListX2, 4, PyFloat_FromDouble(R(1, 1)));
	PyList_SetItem(pyListX2, 5, PyFloat_FromDouble(R(1, 2)));
	PyList_SetItem(pyListX2, 6, PyFloat_FromDouble(R(2, 0)));
	PyList_SetItem(pyListX2, 7, PyFloat_FromDouble(R(2, 1)));
	PyList_SetItem(pyListX2, 8, PyFloat_FromDouble(R(2, 2)));
	PyList_SetItem(pyListX3, 0, PyFloat_FromDouble(t(0)));
	PyList_SetItem(pyListX3, 1, PyFloat_FromDouble(t(1)));
	PyList_SetItem(pyListX3, 2, PyFloat_FromDouble(t(2)));
	for (int i = 0; i < layerpos.size(); i++)
	{
		PyList_SetItem(pyListX1, 3 * i + 0, PyFloat_FromDouble(Vt(3 * i + 0)));
		PyList_SetItem(pyListX1, 3 * i + 1, PyFloat_FromDouble(Vt(3 * i + 1)));
		PyList_SetItem(pyListX1, 3 * i + 2, PyFloat_FromDouble(Vt(3 * i + 2)));
		PyList_SetItem(pyListX4, 3 * i + 0, PyFloat_FromDouble(layerpos[i](0)));
		PyList_SetItem(pyListX4, 3 * i + 1, PyFloat_FromDouble(layerpos[i](1)));
		PyList_SetItem(pyListX4, 3 * i + 2, PyFloat_FromDouble(layerpos[i](2)));
	}
	PyObject* pyList = PyTuple_New(4);
	PyTuple_SetItem(pyList, 0, pyListX1);
	PyTuple_SetItem(pyList, 1, pyListX2);
	PyTuple_SetItem(pyList, 2, pyListX3);
	PyTuple_SetItem(pyList, 3, pyListX4);
	pReturnValue = PyObject_CallObject(pFuncFindRigid, pyList);
	PyObject *out1, *out2;
	out1 = PyTuple_GetItem(pReturnValue, 0);
	out2 = PyTuple_GetItem(pReturnValue, 1);
	PyObject *it;
	it = PyList_GetItem(out1, 0);
	PyArg_Parse(it, "d", &(R(0, 0)));
	it = PyList_GetItem(out1, 1);
	PyArg_Parse(it, "d", &(R(0, 1)));
	it = PyList_GetItem(out1, 2);
	PyArg_Parse(it, "d", &(R(0, 2)));
	it = PyList_GetItem(out1, 3);
	PyArg_Parse(it, "d", &(R(1, 0)));
	it = PyList_GetItem(out1, 4);
	PyArg_Parse(it, "d", &(R(1, 1)));
	it = PyList_GetItem(out1, 5);
	PyArg_Parse(it, "d", &(R(1, 2)));
	it = PyList_GetItem(out1, 6);
	PyArg_Parse(it, "d", &(R(2, 0)));
	it = PyList_GetItem(out1, 7);
	PyArg_Parse(it, "d", &(R(2, 1)));
	it = PyList_GetItem(out1, 8);
	PyArg_Parse(it, "d", &(R(2, 2)));
	it = PyList_GetItem(out2, 0);
	PyArg_Parse(it, "d", &(t(0)));
	it = PyList_GetItem(out2, 1);
	PyArg_Parse(it, "d", &(t(1)));
	it = PyList_GetItem(out2, 2);
	PyArg_Parse(it, "d", &(t(2)));
}


void FaceBodyModel::getAnchorPos(Eigen::VectorXd s_coeff, Eigen::VectorXd pose_parm)
{
	anchorpointpos.resize(27);
	PyObject* pReturnValue;
	//两个入参分别为 x,y coordinate of each point
	PyObject* pyListX1 = PyList_New(s_coeff.size());
	PyObject* pyListX2 = PyList_New(pose_parm.size());
	//#pragma omp parallel for
	for (int i = 0; i<s_coeff.size(); i++)
	{
		PyList_SetItem(pyListX1, i, PyFloat_FromDouble(s_coeff(i)));
	}
	//#pragma omp parallel for
	for (int i = 0; i < pose_parm.size(); i++)
	{
		PyList_SetItem(pyListX2,  i, PyFloat_FromDouble(pose_parm(i)));
	}
	PyObject* pyList = PyTuple_New(2);
	PyTuple_SetItem(pyList, 0, pyListX1);
	PyTuple_SetItem(pyList, 1, pyListX2);
	PyObject *it;
	pReturnValue = PyObject_CallObject(pFuncPose2skeletonpos, pyList);
	for (int i = 0; i < 27; i++)
	{
		Eigen::Vector3d content;
		it = PyList_GetItem(pReturnValue, 3 * i + 0);
		PyArg_Parse(it, "d", &(content(0)));
		it = PyList_GetItem(pReturnValue, 3 * i + 1);
		PyArg_Parse(it, "d", &(content(1)));
		it = PyList_GetItem(pReturnValue, 3 * i + 2);
		PyArg_Parse(it, "d", &(content(2)));
		anchorpointpos[i] = content;
	}
	return;
}