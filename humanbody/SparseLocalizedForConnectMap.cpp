#include "SparseLocalizedForConnectMap.h"


		unsigned char colormap7 [192] = 
		{
			0, 	0, 	143, 	
0, 	0, 	159, 	
0, 	0, 	175, 	
0, 	0, 	191, 	
0, 	0, 	207, 	
0, 	0, 	223, 	
0, 	0, 	239, 	
0, 	0, 	255, 	
0, 	16, 	255, 	
0, 	32, 	255, 	
0, 	48, 	255, 	
0, 	64, 	255, 	
0, 	80, 	255, 	
0, 	96, 	255, 	
0, 	112, 	255, 	
0, 	128, 	255, 	
0, 	143, 	255, 	
0, 	159, 	255, 	
0, 	175, 	255, 	
0, 	191, 	255, 	
0, 	207, 	255, 	
0, 	223, 	255, 	
0, 	239, 	255, 	
0, 	255, 	255, 	
16, 	255, 	239, 	
32, 	255, 	223, 	
48, 	255, 	207, 	
64, 	255, 	191, 	
80, 	255, 	175, 	
96, 	255, 	159, 	
112, 	255, 	143, 	
128, 	255, 	128, 	
143, 	255, 	112, 	
159, 	255, 	96, 	
175, 	255, 	80, 	
191, 	255, 	64, 	
207, 	255, 	48, 	
223, 	255, 	32, 	
239, 	255, 	16, 	
255, 	255, 	0, 	
255, 	239, 	0, 	
255, 	223, 	0, 	
255, 	207, 	0, 	
255, 	191, 	0, 	
255, 	175, 	0, 	
255, 	159, 	0, 	
255, 	143, 	0, 	
255, 	128, 	0, 	
255, 	112, 	0, 	
255, 	96, 	0, 	
255, 	80, 	0, 	
255, 	64, 	0, 	
255, 	48, 	0, 	
255, 	32, 	0, 	
255, 	16, 	0, 	
255, 	0, 	0, 	
239, 	0, 	0, 	
223, 	0, 	0, 	
207, 	0, 	0, 	
191, 	0, 	0, 	
175, 	0, 	0, 	
159, 	0, 	0, 	
143, 	0, 	0, 	
128, 	0, 	0, 	
		};
		unsigned char colormap5[300] = 
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
		unsigned char colorMap[300] = 
		{
			12 , 60 , 83
			, 12 , 62 , 83
			, 10 , 65 , 84
			,9 , 68 , 85
			,7 , 71 , 86
			,5 , 74 , 86
			,3 , 78 , 86
			,1 , 82 , 87
			,0 , 87 , 88
			,0 , 91 , 89
			,0 , 95 , 89
			,0 ,100 , 90
			,0 ,104 , 91
			,0 ,108 , 91
			,0 ,114 , 93
			,0 ,119 , 94
			,0 ,123 , 94
			,0 ,127 , 95
			,0 ,132 , 96
			,0 ,136 , 97
			,0 ,140 , 98
			,0 ,145 , 99
			,0 ,149 ,100
			,0 ,152 ,100
			,6 ,156 ,101
			, 11 ,160 ,102
			, 17 ,164 ,103
			, 24 ,168 ,105
			, 33 ,173 ,107
			, 42 ,177 ,108
			, 51 ,181 ,110
			, 60 ,186 ,112
			, 70 ,190 ,113
			, 81 ,194 ,116
			, 91 ,200 ,117
			,102 ,204 ,119
			,113 ,208 ,121
			,124 ,212 ,122
			,135 ,216 ,124
			,146 ,220 ,125
			,157 ,223 ,125
			,168 ,226 ,124
			,177 ,230 ,125
			,187 ,232 ,125
			,197 ,234 ,125
			,205 ,235 ,125
			,213 ,236 ,124
			,221 ,237 ,124
			,227 ,237 ,125
			,233 ,237 ,125
			,238 ,237 ,124
			,244 ,234 ,121
			,248 ,232 ,118
			,251 ,229 ,114
			,255 ,225 ,110
			,255 ,220 ,105
			,255 ,215 , 99
			,255 ,210 , 93
			,255 ,203 , 87
			,255 ,197 , 81
			,255 ,189 , 75
			,255 ,182 , 68
			,255 ,175 , 61
			,255 ,168 , 55
			,255 ,160 , 48
			,255 ,153 , 41
			,255 ,146 , 36
			,255 ,138 , 30
			,255 ,131 , 24
			,255 ,125 , 18
			,255 ,118 , 14
			,255 ,112 , 10
			,255 ,107  ,6
			,255 ,101  ,2
			,255 , 97  ,1
			,253 , 93  ,0
			,251 , 89  ,0
			,247 , 86  ,0
			,243 , 83  ,0
			,238 , 79  ,0
			,234 , 77  ,0
			,229 , 73  ,0
			,224 , 71  ,0
			,218 , 67  ,1
			,212 , 66  ,2
			,206 , 63  ,4
			,200 , 60  ,5
			,194 , 58  ,7
			,188 , 56  ,9
			,182 , 54 , 11
			,176 , 52 , 13
			,170 , 49 , 16
			,164 , 48 , 18
			,158 , 46 , 21
			,152 , 44 , 23
			,148 , 43 , 25
			,143 , 42 , 26
			,138 , 40 , 28
			,133 , 38 , 30
			,130 , 37 , 30
		};
		unsigned char colorMap2[300] ={
			    3,  86, 255
			,   5,  88, 255
			,   9,  90, 255
			,  12,  93, 255
			,  16,  95, 255
			,  20,  98, 255
			,  25, 101, 255
			,  29, 105, 255
			,  34, 108, 255
			,  38, 112, 255
			,  43, 116, 255
			,  49, 119, 255
			,  55, 124, 255
			,  60, 128, 255
			,  65, 132, 255
			,  71, 137, 255
			,  77, 141, 255
			,  83, 145, 255
			,  89, 149, 255
			,  95, 154, 255
			, 101, 159, 255
			, 107, 163, 255
			, 114, 168, 255
			, 121, 172, 255
			, 127, 178, 255
			, 133, 182, 255
			, 139, 187, 255
			, 146, 191, 255
			, 152, 195, 255
			, 158, 200, 255
			, 165, 204, 255
			, 171, 209, 255
			, 177, 213, 255
			, 182, 216, 255
			, 188, 220, 255
			, 194, 224, 255
			, 199, 228, 255
			, 205, 231, 255
			, 210, 234, 255
			, 215, 237, 255
			, 221, 240, 255
			, 225, 243, 255
			, 229, 245, 255
			, 234, 247, 255
			, 238, 249, 255
			, 242, 251, 255
			, 245, 252, 255
			, 248, 253, 255
			, 251, 254, 255
			, 255, 255, 255
			, 255, 255, 253
			, 255, 255, 250
			, 255, 255, 247
			, 255, 254, 243
			, 255, 253, 239
			, 255, 252, 235
			, 255, 251, 230
			, 255, 249, 225
			, 255, 247, 220
			, 255, 245, 215
			, 255, 242, 209
			, 255, 240, 203
			, 255, 237, 197
			, 255, 234, 192
			, 255, 231, 185
			, 255, 227, 178
			, 255, 224, 172
			, 255, 221, 165
			, 255, 217, 158
			, 255, 213, 152
			, 255, 210, 145
			, 255, 206, 138
			, 255, 202, 132
			, 255, 198, 124
			, 255, 194, 118
			, 255, 190, 111
			, 255, 186, 104
			, 255, 182,  97
			, 255, 178,  91
			, 255, 174,  84
			, 255, 171,  78
			, 255, 167,  72
			, 255, 162,  66
			, 255, 159,  59
			, 255, 156,  54
			, 255, 152,  48
			, 255, 149,  43
			, 255, 145,  37
			, 255, 142,  33
			, 255, 140,  29
			, 255, 137,  24
			, 255, 133,  20
			, 255, 131,  16
			, 255, 129,  13
			, 255, 127,  10
			, 255, 125,   7
			, 255, 124,   5
			, 255, 122,   3
			, 255, 122,   1
			, 255, 121,   1
		};

		void compress_index(
			const int *Ind,
			int nnz,
			int m,
			int *Ptr,
			int base)
		{
			int i;

			/* initialize everything to zero */
			for (i = 0; i<m + 1; i++) {
				Ptr[i] = 0;
			}
			/* count elements in every row */
			Ptr[0] = base;
			for (i = 0; i<nnz; i++) {
				Ptr[Ind[i] + (1 - base)]++;
			}
			/* add all the values */
			for (i = 0; i<m; i++) {
				Ptr[i + 1] += Ptr[i];
			}
		}

		void SparseLocalizedForConnectMap::setMesh(Eigen::VectorXd V)
		{
			for (int i = 0; i < nV; i++)
			{
				//SpatialTemData[i].resize(1);

				Eigen::VectorXd x = V.segment(3 * i, 3);
				SpatialTemData[i][0] = x;
			}
		}

		void SparseLocalizedForConnectMap::setMesh(Eigen::MatrixXd V)
		{
			for (int i = 0; i < nV; i++)
			{
				//SpatialTemData[i].resize(1);

				Eigen::VectorXd x = V.row(i).transpose();
				SpatialTemData[i][0] = x;
			}
		}

		void SparseLocalizedForConnectMap::setMesh(PGMesh * _mesh)
		{
#pragma omp parallel for
			for (int i = 0; i < nV; i++)
			{
				//SpatialTemData[i].resize(1);
				PGMesh::VertexHandle vh = _mesh->vertex_handle(i);
				PGMesh::Point &p = _mesh->point(vh);
				Eigen::VectorXd x(3);
				x << p[0], p[1], p[2];
				SpatialTemData[i][0] = x;
			}
		}

		void SparseLocalizedForConnectMap::updateMeshVertices(Eigen::VectorXd x, PGMesh *_mesh)
		{
			for (int i = 0; i < nV; i++)
			{
				//SpatialTemData[i].resize(2);
				PGMesh::VertexHandle vh = _mesh->vertex_handle(i);
				PGMesh::Point &p = _mesh->point(vh);
				p[0] = x(3 * i + 0);
				p[1] = x(3 * i + 1);
				p[2] = x(3 * i + 2);
			}
		}

		void SparseLocalizedForConnectMap::updateDeformed(PGMesh * _mesh)
		{
//#pragma omp parallel for
			for (int i = 0; i < nV; i++)
			{
				//SpatialTemData[i].resize(2);
				PGMesh::VertexHandle vh = _mesh->vertex_handle(i);
				PGMesh::Point &p = _mesh->point(vh);
				Eigen::VectorXd x(3);
				x << p[0], p[1], p[2];
				SpatialTemData[i][1] = x;
			}
				for (int j = 0; j < Faceanchors_.size(); j++)
				{
					faceAnchorFrames[1][j] = calculateFaceFrame(Faceanchors_[j], 1);
				}
				for (int j = 0; j < Vertexanchors_.size(); j++)
				{
					vertexAnchorCords[1][j].idx = Vertexanchors_[j];
					vertexAnchorCords[1][j].frameNo = 1;
					vertexAnchorCords[1][j].cord = SpatialTemData[Vertexanchors_[j]][1];
				}

		}

		void SparseLocalizedForConnectMap::updateDeformed1(PGMesh * _mesh)
		{
			/*Faceanchors_.resize(34);
			
			Faceanchors_[0] = 10817;
			Faceanchors_[1] = 10818;
			Faceanchors_[2] = 10833;
			Faceanchors_[3] = 10834;
			Faceanchors_[4] = 21688;
			Faceanchors_[5] = 21689;
			Faceanchors_[6] = 21704;
			Faceanchors_[7] = 21705;
			Faceanchors_[8] = 19939;
			Faceanchors_[9] = 19935;
			Faceanchors_[10] = 19111;
			Faceanchors_[11] = 19112;
			Faceanchors_[12] = 10657;
			Faceanchors_[13] = 5515;
			Faceanchors_[14] = 18687;
			Faceanchors_[15] = 18686;
			Faceanchors_[16] = 4629;
			Faceanchors_[17] = 4630;
			Faceanchors_[18] = 15305;
			Faceanchors_[19] = 20723;
			Faceanchors_[20] = 4432;
			Faceanchors_[21] = 9800;
			Faceanchors_[22] = 20293;
			Faceanchors_[23] = 20292;
			Faceanchors_[24] = 9449;
			Faceanchors_[25] = 9465;

			Faceanchors_[26] = 3742;
			Faceanchors_[27] = 14643;
			Faceanchors_[28] = 20827;
			Faceanchors_[29] = 3890;
			Faceanchors_[30] = 20337;
			Faceanchors_[31] = 9468;
			Faceanchors_[32] = 9794;
			Faceanchors_[33] = 20673;*/
			int a[1] = { 21623 };
			Faceanchors_.resize(1);
			copy(a, a + 1, Faceanchors_.begin());
			/*int a[49] = { 21623, 8440, 15834, 4962, 10857, 8711, 4174, 16313, 10741, 12069, 2194, 16361, 15539, 8252, 4794, 15542, 19119, 8246, 4715, 20802, 14756, 4486, 9846, 3847, 3730, 19577, 14750, 20320, 20368, 9465, 9501, 9214, 20328, 18721, 16879, 17090, 6149, 6008, 5887,  6223, 6169, 7297, 5841, 18211, 16883, 17019, 16686, 20573, 9297 };
			Faceanchors_.resize(49);
			copy(a, a + 49, Faceanchors_.begin());*/
			/*int a[12] = { 19522, 18827, 16651, 8026, 5785, 20689, 9815, 19705, 8834, 19903, 19270, 10080  };
			Faceanchors_.resize(12);
			copy(a, a + 12, Faceanchors_.begin());
			faceAnchorFrames[0].resize(Faceanchors_.size());
			faceAnchorFrames[1].resize(Faceanchors_.size());*/
			//#pragma omp parallel for
			for (int i = 0; i < nV; i++)
			{
				//SpatialTemData[i].resize(2);
				PGMesh::VertexHandle vh = _mesh->vertex_handle(i);
				PGMesh::Point &p = _mesh->point(vh);
				Eigen::VectorXd x(3);
				x << p[0], p[1], p[2];
				SpatialTemData[i][1] = x;
			}
			for (int j = 0; j < Faceanchors_.size(); j++)
			{
				faceAnchorFrames[1][j] = calculateFaceFrame(Faceanchors_[j], 1);
				faceAnchorFrames[0][j] = calculateFaceFrame(Faceanchors_[j], 0);
			}

		}
		void SparseLocalizedForConnectMap::updateDeformed1(std::vector<double> anchorfaces)
		{
			int a[12] = { 19522, 18827, 16651, 8026, 5785, 20689, 9815, 19705, 8834, 19903, 19270, 10080 };
			/*Faceanchors_.resize(12);
			copy(a, a + 12, Faceanchors_.begin());*/
			Faceanchors_.resize(1);
			copy(a, a + 1, Faceanchors_.begin());
			faceAnchorFrames[0].resize(Faceanchors_.size());
			faceAnchorFrames[1].resize(Faceanchors_.size());


			for (int i = 0; i < Faceanchors_.size(); i++)
			{

				faceanchor fanchor;
				fanchor.localframe << anchorfaces[9 * i + 0], anchorfaces[9 * i + 1], anchorfaces[9 * i + 2],
					anchorfaces[9 * i + 3], anchorfaces[9 * i + 4], anchorfaces[9 * i + 5], 
					anchorfaces[9 * i + 6], anchorfaces[9 * i + 7], anchorfaces[9 * i + 8];
				fanchor.frameNo = 1;
				fanchor.idx = Faceanchors_[i];
				faceAnchorFrames[1][i] = fanchor;
				faceAnchorFrames[0][i] = fanchor;
			}
		}
		void SparseLocalizedForConnectMap::updateDeformed1(string path_)
		{
			//int a[1] = { 21623 };
			/*int a[49] = { 21623, 8440, 15834, 4962, 10857, 8711, 4174, 16313, 10741, 12069, 2194, 16361, 15539, 8252, 4794, 15542, 19119, 8246, 4715, 20802, 14756, 4486, 9846, 3847, 3730, 19577, 14750, 20320, 20368, 9465, 9501, 9214, 20328, 18721, 16879, 17090, 6149, 6008, 5887,  6223, 6169, 7297, 5841, 18211, 16883, 17019, 16686, 20573, 9297 };
			Faceanchors_.resize(49);
			copy(a, a + 49, Faceanchors_.begin());*/
			int a[12] = { 19522, 18827, 16651, 8026, 5785, 20689, 9815, 19705, 8834, 19903, 19270, 10080 };
			Faceanchors_.resize(12);
			copy(a, a + 12, Faceanchors_.begin());
			/*Faceanchors_.resize(1);
			copy(a, a + 1, Faceanchors_.begin());*/
			faceAnchorFrames[0].resize(Faceanchors_.size());
			faceAnchorFrames[1].resize(Faceanchors_.size());
			std::ifstream input(path_);
			for (int i = 0; i < Faceanchors_.size(); i++)
			{
				double a0, a1, a2, a3, a4, a5, a6, a7, a8;
				input >> a0;
				input >> a1;
				input >> a2;
				input >> a3;
				input >> a4;
				input >> a5;
				input >> a6;
				input >> a7;
				input >> a8;
				/*if (i == 9 || i == 10)
					continue;*/
				faceanchor fanchor;
				fanchor.localframe << a0, a1, a2,
					a3, a4, a5,
					a6, a7, a8;
				fanchor.frameNo = 1;
				fanchor.idx = Faceanchors_[i];
				faceAnchorFrames[1][i] = fanchor;
				faceAnchorFrames[0][i] = calculateFaceFrame(Faceanchors_[i], 0);
			}
		}

		SparseLocalizedForConnectMap::SparseLocalizedForConnectMap(PGMesh *_mesh)
		{
			mesh_ = _mesh;
			//SpatialTemData = sd;
			mesh1 = new PGMesh();
			OpenMesh::IO::read_mesh(*mesh1, "F:\\scan\\smpl_H\\female_shape\\0000.obj");
			mesh1->request_face_normals();
			mesh1->update_face_normals();
			nV = mesh_->n_vertices();
			nF = mesh_->n_faces();
			nE = mesh_->n_edges();
			nSeq = 2;
			SpatialTemData.resize(nV);
//#pragma omp parallel for
			for (int i = 0; i < nV; i++)
			{
				SpatialTemData[i].resize(2);
				PGMesh::VertexHandle vh = _mesh->vertex_handle(i);
				PGMesh::Point &p = _mesh->point(vh);
				Eigen::VectorXd x(3);
				x << p[0], p[1], p[2];
				SpatialTemData[i][0] = x;
				SpatialTemData[i][1] = x;
			}
			FaceAreas.resize(nF);
			FaceNormals.resize(nF);
			for (int i = 0; i < nF; i++)
			{
				FaceAreas[i].resize(nSeq);
				FaceNormals[i].resize(nSeq);
			}
			EdgeLocalFrames.resize(nE);
			EdgeNormals.resize(nE);
			ConnectMapData.resize(nE);
			ConnectMapData1.resize(nE);
			faceAnchorFrames.resize(nSeq);
			vertexAnchorCords.resize(nSeq);
			for (int i = 0; i < 2; i++)
			{
				setAnchor(  i);
			}
			/*for (int i = 0; i < 40; i++)
			{
				setAnchor(250 * i);
			}
			setAnchor(nV - 1);*/
			for (int i = 0; i < nE; i++)
			{
				EdgeLocalFrames[i].resize(nSeq);
				EdgeNormals[i].resize(nSeq);
				ConnectMapData[i].resize(nSeq);
				ConnectMapData1[i].resize(nSeq);
			}
			for (int i = 0; i < nSeq; i++)
			{
				faceAnchorFrames[i].resize(Faceanchors_.size());
				vertexAnchorCords[i].resize(Vertexanchors_.size());
			}
			//setAnchor(0);
			/*calculateFaceArea();
			calculateEdgesNormals();
			calculateLocalFrame();
			calculateFacesNormals();*/
			for (int i = 0; i < nSeq; i++)
			{
				for (int j = 0; j < Faceanchors_.size(); j++)
				{
					faceAnchorFrames[i][j] = calculateFaceFrame(Faceanchors_[j], i);
				}
				for (int j = 0; j < Vertexanchors_.size(); j++)
				{
					vertexAnchorCords[i][j].idx = Vertexanchors_[j];
					vertexAnchorCords[i][j].frameNo = i;
					vertexAnchorCords[i][j].cord = SpatialTemData[Vertexanchors_[j]][i];
				}
			}
		}



		/** 初始化*/
		SparseLocalizedForConnectMap::SparseLocalizedForConnectMap(PGMesh * _mesh, std::vector<std::vector<Eigen::VectorXd>> & sd)
		{
			mesh_ = _mesh;
			SpatialTemData = sd;
			nV = mesh_->n_vertices();
			nF = mesh_->n_faces();
			nE = mesh_->n_edges();
			nSeq = SpatialTemData[0].size();
			FaceAreas.resize(nF);
			FaceNormals.resize(nF);
			for(int i = 0; i < nF; i++)
			{
				FaceAreas[i].resize(nSeq);
				FaceNormals[i].resize(nSeq);
			}
			EdgeLocalFrames.resize(nE);
			EdgeNormals.resize(nE);
			ConnectMapData.resize(nE);
			ConnectMapData1.resize(nE);
			faceAnchorFrames.resize(nSeq);
			vertexAnchorCords.resize(nSeq);
			for(int i = 0; i < 10; i++)
			{
				setAnchor(i);
			}
			for(int i = 0; i < nE; i++)
			{
				EdgeLocalFrames[i].resize(nSeq);
				EdgeNormals[i].resize(nSeq);
				ConnectMapData[i].resize(nSeq);
				ConnectMapData1[i].resize(nSeq);
			}
			for(int i = 0; i < nSeq; i++)
			{
				faceAnchorFrames[i].resize(Faceanchors_.size());
				vertexAnchorCords[i].resize(Vertexanchors_.size());
			}
			//setAnchor(0);
			/*calculateFaceArea();
			calculateEdgesNormals();
			calculateLocalFrame();
			calculateFacesNormals();*/
			for(int i = 0; i < nSeq; i++)
			{
				for(int j = 0; j < Faceanchors_.size(); j++ )
				{
					faceAnchorFrames[i][j] = calculateFaceFrame(Faceanchors_[j], i);
				}
				for(int j = 0; j < Vertexanchors_.size(); j ++)
				{
					vertexAnchorCords[i][j].idx = Vertexanchors_[j];
					vertexAnchorCords[i][j].frameNo = i;
					vertexAnchorCords[i][j].cord = SpatialTemData[Vertexanchors_[j]][i];
				}
			}
			
		//	for (int add = 0; add < 60; add++)
		//{
		//	/*vertexError.setZero();
		//	vertexErrorIncrement.setZero();
		//	for(int i = 0; i < ConnectMap->nSeq; i++)
		//	{
		//		ConnectMap->constructFrameFromComponentsError(i, k, str1, vertexErrorIncrement);
		//		vertexError = vertexError + vertexErrorIncrement;
		//	}
		//	int anchorno;
		//	double lllll = vertexError.maxCoeff(& anchorno);
		//	ConnectMap->addAnchor(anchorno);*/
		//	addAnchor(add * 166);
		//}
			//setAnchor(0);
		}
		/** 设置锚点*/
		void SparseLocalizedForConnectMap::setAnchor(unsigned int _idx)
		{
			OpenMesh::VertexHandle vh = mesh_->vertex_handle(_idx);
			anchorPointidx = vh.idx();
			anchorEidx = mesh_->ve_begin(vh)->idx();
			int faceNo;
			faceNo = mesh_->vf_begin(vh)->idx();
			//faceNo = 1;
			edgeanchors_.push_back(anchorEidx);
			Vertexanchors_.push_back(anchorPointidx);
			Faceanchors_.push_back(faceNo);
		}

		void SparseLocalizedForConnectMap::setAnchor(vector<vector<int>> idx_set, vector<vector<double>> weights, vector<Eigen::Vector3d> pos)
		{
			vertexAnchorCords.resize(nSeq);
			for (int l = 0; l < nSeq; l++)
			{
				vertexAnchorCords[l].resize(1);
				for (int i = 0; i < 1; i++)
				{
					vertexAnchorCords[l][i].weights = weights[27];
					vertexAnchorCords[l][i].idx_list = idx_set[27];
					vertexAnchorCords[l][i].cord = pos[27];
				}
			}
		}

		/** 添加锚点*/
		void SparseLocalizedForConnectMap::addAnchor(unsigned int _idx)
		{
			OpenMesh::VertexHandle vh = mesh_->vertex_handle(_idx);
			anchorPointidx = vh.idx();
			anchorEidx = mesh_->ve_begin(vh)->idx();
			int faceNo;
			faceNo = mesh_->vf_begin(vh)->idx();
			//faceNo = 1;
			edgeanchors_.push_back(anchorEidx);
			Vertexanchors_.push_back(anchorPointidx);
			Faceanchors_.push_back(faceNo);
			for(int i = 0; i < nSeq; i++)
			{
				faceAnchorFrames[i].push_back(calculateFaceFrame(faceNo, i));
				vertexanchor va;
				va.idx = anchorPointidx;
				va.frameNo = i;
				va.cord = SpatialTemData[anchorPointidx][i];
				vertexAnchorCords[i].push_back(va);
			}
		}

		/** 提取第i帧的网格坐标向量*/
		Eigen::VectorXd SparseLocalizedForConnectMap::getCords (int ith)
		{
			Eigen::VectorXd l(3 * nV);
			for(int i = 0; i < nV; i++)
			{
				l.segment(3 * i, 3) = SpatialTemData[i][ith];
			}
			return l;
		}

		/** 计算各帧每个面的面积以及法向量*/
		void SparseLocalizedForConnectMap::calculateFaceArea()
		{
			FaceAreas.resize(mesh_->n_faces());
			PGMesh::ConstFaceIter cfIt = mesh_->faces_begin();
			PGMesh::ConstFaceIter cfItEnd = mesh_->faces_end();
			for ( ; cfIt != cfItEnd; ++ cfIt )
			{
				FaceAreas[cfIt->idx()].resize(nSeq);
				for (int j = 0; j < nSeq; j++)
				{
					PGMesh::ConstFaceEdgeIter cfeIt = mesh_->cfe_iter(*cfIt);

					Eigen::VectorXd edgeVec1 = EDGE_POINT1( *cfeIt, j ) - EDGE_POINT0( *cfeIt, j);
					++cfeIt;
					Eigen::VectorXd edgeVec2 = EDGE_POINT1(*cfeIt, j) - EDGE_POINT0(*cfeIt, j);
					double area = (computeCross(edgeVec1, edgeVec2)).norm() * 0.5L;
					Eigen::Vector3d facenormal= computeCross(edgeVec1, edgeVec2) / (computeCross(edgeVec1, edgeVec2)).norm();
					FaceAreas[cfIt->idx()][j] = area;
					FaceNormals[cfIt->idx()][j] = facenormal;
				}
			}
		}

		/** 计算各帧每条边的法向量*/
		void SparseLocalizedForConnectMap::calculateEdgesNormals()
		{
			PGMesh::ConstEdgeIter ceIt = mesh_->edges_begin();
			PGMesh::ConstEdgeIter celtEnd = mesh_->edges_end();

			for ( int counter = 0; ceIt != celtEnd; ++ceIt, ++counter )
			{
				for(int j = 0; j < nSeq; j++)
				{
					Eigen::Vector3d edgeNormalSum(0.0,0.0,0.0);
					//double weightSum = .0L;
					double weightSum = 0.5;

					for ( int i = 0; i < 2; ++i )
					{
						OpenMesh::FaceHandle faceHandle = mesh_->face_handle(mesh_->halfedge_handle(*ceIt, i));
						/*if (ceIt.handle().idx() == 3771 && j == 1)
					{
						int  lll = faceHandle.idx();
						cout<< faceHandle<<lll;
					}*/
						// no idea about the underlying data type of 'normal' 
						if (faceHandle.idx() != -1)	// * Handle the case of an open mesh. * //
						{

							Eigen::Vector3d edgeNormal = FaceNormals[faceHandle.idx()][j];
							edgeNormal /= edgeNormal.norm();
							/*double weight = ptrMesh->property( fPropHandle, faceHandle )  ;
							edgeNormal *= weight;
							edgeNormalSum += edgeNormal;
							weightSum += weight;*/

							// * no weighted-average normal * //
							//edgeNormal *= weight;
							if((edgeNormalSum + edgeNormal).norm() == 0)
								edgeNormal << 0, 0, 0;
							edgeNormalSum += edgeNormal;
						}
					}

					/*if ( weightSum < 0.000000001L && weightSum >= 0.0L)
					{
						return -2;
					}

					edgeNormalSum /= weightSum;*/
					/*if (ceIt.handle().idx() == 3771 && j == 1)
					{
						Eigen::Vector3d  lll = edgeNormalSum;
						cout<< lll;
					}*/
					/*if (edgeNormalSum.norm()== 0)
					{
					}*/
					edgeNormalSum /= edgeNormalSum.norm();
					EdgeNormals[ceIt->idx()][j] = edgeNormalSum;	
				}
				
			}
		}

		/** 计算各帧各条边的局部标架*/
		void SparseLocalizedForConnectMap::calculateLocalFrame()
		{
			
			PGMesh::ConstEdgeIter ceIt = mesh_->edges_begin();
			PGMesh::ConstEdgeIter ceItEnd = mesh_->edges_end();
	
			//int i = 0;
			for ( ; ceIt != ceItEnd; ++ceIt )
			{
				for(int j = 0; j < nSeq; j++)
				{
					Eigen::Vector3d edgeVec = EDGE_POINT1(*ceIt, j) - EDGE_POINT0(*ceIt, j);
					Eigen::Vector3d zAxis = EdgeNormals[ceIt->idx()][j];
					/*if (ceIt.handle().idx() == 3771 && j == 1)
					{
						Eigen::Vector3d  lll = EdgeNormals [ceIt.handle().idx()][j];
						cout<< lll;
					}*/
					if (zAxis.norm() == 0)
					{
						cout<< "error"<<endl; 
					}
					zAxis /= zAxis.norm();		// The unit vector representing z axis of the local frame
	
					Eigen::Vector3d yAxis = computeCross (zAxis , edgeVec);
					if (yAxis.norm() == 0)
					{
						cout<< "error"<<endl; 
					}
					yAxis /= yAxis.norm();	// The unit vector representing y axis of the local frame
					
					Eigen::Vector3d xAxis = edgeVec / edgeVec.norm();	// The unit vector representing x axis of the local frame
					//Eigen::Vector3d xAxis =  computeCross (yAxis , zAxis);
					if (xAxis.norm() == 0)
					{
						cout<< "error"<<endl; 
					}
					Eigen::Matrix3d Frame;
					Frame<< xAxis(0), xAxis(1), xAxis(2),
						yAxis(0), yAxis(1), yAxis(2),
						zAxis(0), zAxis(1), zAxis(2);

					// * initialize p_prev_mat_X_ 
					/*p_prev_mat_X_->block(3 * i, 0, 1, 3) << xAxis[0], xAxis[1], xAxis[2];
					p_prev_mat_X_->block(3 * i + 1, 0, 1, 3) << yAxis[0], yAxis[1], yAxis[2];
					p_prev_mat_X_->block(3 * i + 2, 0, 1, 3) << zAxis[0], zAxis[1], zAxis[2];*/
					//++i;
					EdgeLocalFrames[ceIt->idx()][j] = Frame;
				}
				
			}

		}

		

		/** 给出旋转矩阵计算欧拉旋转角*/
		void SparseLocalizedForConnectMap::rotMat2EAngle(Eigen::Matrix3d  rotMat,Eigen::Vector3d & eulerAg)
		{
				double x_A = 0.0;
				double y_A = 0.0;
				double z_A = 0.0;

				//rotMat = rotMat.transpose();
				//Eigen::Matrix3d rotMat1 = rotMat.transpose();
				y_A = -asin(rotMat(2,0));

				assert(cos(y_A) > 1e-6 || cos(y_A)< -1e-6);

				x_A = atan2f(rotMat(2,1)/cos(y_A),rotMat(2,2)/cos(y_A));
				z_A = atan2f(rotMat(1,0)/cos(y_A),rotMat(0,0)/cos(y_A));

				eulerAg(0) = x_A;
				eulerAg(1) = y_A;
				eulerAg(2) = z_A;
				y_A = -asin(rotMat(0,2));

				//assert(cos(y_A) > 1e-6 || cos(y_A)< -1e-6);

				//x_A = atan2f(rotMat(1,2)/cos(y_A),rotMat(2,2)/cos(y_A));
				//z_A = atan2f(rotMat(0,1)/cos(y_A),rotMat(0,0)/cos(y_A));

				//eulerAg(0) = x_A;
				//eulerAg(1) = y_A;
				//eulerAg(2) = z_A;
			/*Eigen::Matrix3d log_U = rotMat.log();
			eulerAg(0) = log_U( 0, 1 );
			eulerAg(1) = log_U( 0, 2 );
			eulerAg(2) = log_U( 1, 2 );*/

		}

		/** 给出局部标价的变换矩阵计算旋转轴和旋转角*/
		void SparseLocalizedForConnectMap::rotMat2AxisAngle(Eigen::Matrix3d rotMat, Eigen::Vector3d & axis, double & angle)
		{
			
		}

		/** 通过欧拉角给出旋转矩阵*/
		void SparseLocalizedForConnectMap::ERAngle2rotMat(Eigen::Vector3d  eulerAg, Eigen::Matrix3d  &rotMat)
		{
			/*MatrixXd U_( 3, 3 );
			U_.setZero();
			U_( 0, 1 ) = eulerAg( 0 );
			U_( 0, 2 ) = eulerAg( 1 );
			U_( 1, 2 ) = eulerAg( 2 );
			U_( 1, 0 ) = - U_( 0, 1 );
			U_( 2, 0 ) = - U_( 0, 2 );
			U_( 2, 1 ) = - U_( 1, 2 );
			rotMat = U_.exp();*/
			double x = eulerAg(0);
			double y = eulerAg(1);
			double z = eulerAg(2);
			Eigen::Matrix3d x_rot, y_rot, z_rot;
			x_rot<< 1, 0, 0,
				0, cos(x), -sin(x),
				0, sin(x), cos(x);
			y_rot<< cos(y), 0, sin(y),
				0, 1, 0,
				-sin(y), 0, cos(y);
			z_rot<< cos(z), -sin(z), 0,
				sin(z), cos(z), 0,
				0,0,1;
			rotMat = (z_rot * y_rot * x_rot);

		}

		/** 建立各帧的连接映射*/
		void SparseLocalizedForConnectMap::buildEdgeGraph()
		{
			//ofstream file1("connectMapData.txt");
			edge_graph_.clear();
			PGMesh::EdgeIter ceit0 = mesh_->edges_begin();
			PGMesh::EdgeIter ceit_end0 = mesh_->edges_end();
			ConnectMapMatrix = Eigen::MatrixXd(nSeq,13 * nE);
			ConnectMapMatrix1 = Eigen::MatrixXd(nSeq,5 * nE);
			//ofstream onering("onering.txt");
			for (; ceit0 != ceit_end0; ++ceit0)
			{
				EdgeGraphNode egn;
				PGMesh::HalfedgeHandle heh0 = mesh_->halfedge_handle(*ceit0, 0);
				PGMesh::HalfedgeHandle heh1 = mesh_->halfedge_handle(*ceit0, 1);
				PGMesh::FaceHandle fh0 = mesh_->face_handle( heh0 );
				PGMesh::FaceHandle fh1 = mesh_->face_handle( heh1 );
				egn.center_node_id = ceit0->idx();
				if (fh0.is_valid())
				{
					PGMesh::HalfedgeHandle loop_heh = heh0;
			
					for (int i = 0; i < 2; ++i)
					{
						loop_heh = mesh_->next_halfedge_handle(loop_heh);
						PGMesh::EdgeHandle eh0 = mesh_->edge_handle(loop_heh);
						int idx = eh0.idx();
						egn.one_ring_idx.push_back(idx);
					}
				}
				if (fh1.is_valid())
				{
					PGMesh::HalfedgeHandle loop_heh = heh1;
			
					for (int i = 0; i < 2; ++i)
					{
						loop_heh = mesh_->next_halfedge_handle(loop_heh);
						PGMesh::EdgeHandle eh1 = mesh_->edge_handle(loop_heh);
						int idx = eh1.idx();
						egn.one_ring_idx.push_back(idx);
					}
				}
				if(egn.one_ring_idx.size() == 3)
				{
					egn.one_ring_idx.push_back(egn.one_ring_idx[0]);
				}
				if(egn.one_ring_idx.size() == 2)
				{
					egn.one_ring_idx.push_back(egn.one_ring_idx[0]);
					egn.one_ring_idx.push_back(egn.one_ring_idx[1]);
				}
				if(egn.one_ring_idx.size() == 1)
				{
					egn.one_ring_idx.push_back(egn.one_ring_idx[0]);
					egn.one_ring_idx.push_back(egn.one_ring_idx[0]);
					egn.one_ring_idx.push_back(egn.one_ring_idx[0]);
				}
				/*for (int i = 0; i < egn.one_ring_idx.size(); i ++)
				{
					int idx = egn.one_ring_idx[i];
					
				}*/
				//onering<< egn.one_ring_idx[0]<<" "<<egn.one_ring_idx[1]<<" "<<egn.one_ring_idx[2]<<" "<<egn.one_ring_idx[3]<<endl;
				int localEdgeIndex = ceit0->idx();
				egn.one_ring_local_rot.resize(nSeq);
				egn.rotationaxis.resize(nSeq);
				egn.one_ring_res_rot.resize(nSeq);
				for (int j = 0; j < nSeq; j++)
				{
					egn.one_ring_local_rot[j].resize(4);
					egn.rotationaxis[j].resize(4);
					egn.one_ring_res_rot[j].resize(4);
					Eigen::VectorXd  connectMapVector(13);
					//Eigen::VectorXd connectVector(5);
					// * Compute local rotation relative to the center's. * //
					for (int i = 0; i < 4; i++)
					{
						int comparedIndex = egn.one_ring_idx[i];
						Eigen::Matrix3d rot;
						Eigen::Matrix3d rot1 = EdgeLocalFrames[comparedIndex][j];
						Eigen::Matrix3d rot2 = EdgeLocalFrames[localEdgeIndex][j];
						rot = rot2 * rot1.inverse();
						egn.one_ring_local_rot[j][i] = rot;
						//Eigen::Matrix3d rottrans = rot.transpose();
						/*Eigen::Matrix3d rottrans = rot.transpose();*/
						//Matrix3d log_U = rottrans.log();
						//Eigen::Vector3d  eulerAg;
						////rotMat2EAngle(rottrans, eulerAg);
						//eulerAg<< log_U( 0, 1 ), log_U( 0, 2 ), log_U( 1, 2 );
						//connectMapVector.segment(3 * i, 3) = eulerAg;
						//egn.one_ring_local_rot[j][i] = rot;
						Eigen::Matrix3d rot0 = egn.one_ring_local_rot[0][i];
						Eigen::Matrix3d rot3 = rot * rot0.inverse();
						egn.one_ring_res_rot[j][i] = rot3;
						Eigen::Vector3d  eulerAg;
						rotMat2EAngle(rot3.transpose(), eulerAg);
						/*Matrix3d log_U = (rot3).log();
						eulerAg<< log_U( 0, 1 ), log_U( 0, 2 ), log_U( 1, 2 );*/
						connectMapVector.segment(3 * i, 3) = eulerAg;
						/*GeometryProcess::Math::Geometry::axisAndAngle(quat, axis_and_angle);
						angle = axis_and_angle.norm();
						if(axis_and_angle.norm() == 0)
							angle = 0;
						else
						{
							angle = axis_and_angle.norm();
							axis_and_angle = axis_and_angle/angle;
						}
						egn.rotationaxis[j][i] = axis_and_angle;
						connectVector(i) = angle;*/
					}
					double elength = (EDGE_POINT1(*ceit0, j) - EDGE_POINT0(*ceit0, j)).norm();
					egn.currentLength = elength;
					connectMapVector(12) = elength;
					//connectVector(4) = elength;
					ConnectMapData[ceit0->idx()][j] = connectMapVector;
					//ConnectMapData1[ceit0.handle().idx()][j] = connectVector;
					ConnectMapMatrix.block<1, 13>(j, 13 * ceit0->idx() ) = connectMapVector.transpose();
					//ConnectMapMatrix1.block<1, 5>(j, 5 * ceit0.handle().idx() ) = connectVector.transpose();
					//file1<<connectMapVector<<endl;
					
				}
				edge_graph_.push_back(egn);
			}
			
 		}

		/** 可视化第i帧连接映射*/
		void SparseLocalizedForConnectMap::visualizeConnectMap(int ith)
		{
			Eigen::VectorXd edgeError(nE);
			for(int i = 0; i < nE; i++)
			{
				//算第i帧与第0帧的差
				double error = (IICsMatrix.block<1,3>(ith,3 * i)- IICsMatrix.block<1,3>(0,3 * i)).norm();
				edgeError(i) = error;
			}
			char *pPath = new char[50];
			sprintf_s(pPath,50, "%d-thConMapVisual.obj", ith);
			string offname = pPath;
			ofstream offfile(offname);
			delete []pPath;
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
				if(veit->is_valid())
				{
					pointerror=0;
				for (; veit != mesh_->voh_end(vit); ++veit)
				{
					int edgeidx = mesh_->edge_handle(*veit).idx();
					
					pointerror += edgeError(edgeidx);
					
					//pointerror = 0;
					j++;
				}
				
				//求该点对应相邻边的误差的平均值
				errors[i] = pointerror/j;
				}
				else
					errors[i] = 0;
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
				unsigned char* colorPtr = colorMap + colorPersent * 3;
				//c.values_[0] = colorPtr[0];
				//c.values_[1] = colorPtr[1];
				//c.values_[2] = colorPtr[2];
				Eigen::Vector3d position = SpatialTemData[i][ith];
				offfile << "v" << " " << position( 0) << " " <<position(1) << " " << position(2)  <<" " ;
				offfile <<(int) (colorPtr[0])  <<" "<< int(colorPtr[1])<< " " << int (colorPtr[2])<< endl;
			}


			for (PGMesh::FaceIter fit = mesh_->faces_begin(); fit !=		 mesh_->faces_end(); ++ fit)
			{
				offfile << "f" << " ";
				for (PGMesh::FVIter fvit = mesh_->fv_begin(*fit); fvit != mesh_->fv_end(*fit); ++ fvit)
				{
					offfile << fvit->idx() + 1 << " ";
				}
				offfile << endl;
			}
			offfile.close();
		}

		/** 细分网格*/
		/* 每条边增加一个顶点，每个三角形变为4个三角形*/
		void SparseLocalizedForConnectMap::subdivision()
		{
			subdivMesh = new PGMesh();
			
			mesh_->add_property(vhandle, "The new handle on each vertex of old mesh");

			for (PGMesh::VertexIter vit = mesh_->vertices_begin(); vit != mesh_->vertices_end(); ++ vit)
			{
				PGMesh::Point p = mesh_->point(*vit);
				mesh_->property(vhandle, *vit) = subdivMesh->add_vertex(p);
			}
			
			mesh_->add_property(ehandle, "the new vertex handle on each edge");
			for (PGMesh::EdgeIter eit = mesh_->edges_begin(); eit != mesh_->edges_end(); ++ eit)
			{
				PGMesh::HalfedgeHandle heit = mesh_->halfedge_handle(*eit, 0);
				PGMesh::VertexHandle vit0 = mesh_->from_vertex_handle(heit);
				PGMesh::VertexHandle vit1 = mesh_->to_vertex_handle(heit);
				PGMesh::Point p0 = mesh_->point(vit0);
				PGMesh::Point p1 = mesh_->point(vit1);
				PGMesh::Point p((p0 + p1) / 2.0);
				mesh_->property(ehandle, *eit) = subdivMesh->add_vertex(p);
			}
			for (PGMesh::FaceIter fit = mesh_->faces_begin(); fit != mesh_->faces_end(); ++ fit)
			{
				PGMesh::FaceHalfedgeIter fhit = mesh_->fh_begin(*fit);
				PGMesh::HalfedgeHandle h0 = *fhit;
				++fhit;
				PGMesh::HalfedgeHandle h1 = *fhit;
				++fhit;
				PGMesh::HalfedgeHandle h2 = *fhit;
				PGMesh::VertexHandle v0 = mesh_->from_vertex_handle(h0);
				PGMesh::VertexHandle v1 = mesh_->to_vertex_handle(h0);
				PGMesh::VertexHandle v2 = mesh_->to_vertex_handle(h1);
				PGMesh::VertexHandle p0 = mesh_->property(vhandle, v0);
				PGMesh::VertexHandle p1 = mesh_->property(vhandle, v1);
				PGMesh::VertexHandle p2 = mesh_->property(vhandle, v2);
				PGMesh::VertexHandle p3 = mesh_->property(ehandle, mesh_->edge_handle(h0));
				PGMesh::VertexHandle p4 = mesh_->property(ehandle, mesh_->edge_handle(h1));
				PGMesh::VertexHandle p5 = mesh_->property(ehandle, mesh_->edge_handle(h2));
				std::vector< PGMesh::VertexHandle > face_vHandle( 3 );
				face_vHandle.clear();
				face_vHandle.push_back( p0);
				face_vHandle.push_back( p3);
				face_vHandle.push_back( p5);
				subdivMesh->add_face( face_vHandle );

				face_vHandle.clear();
				face_vHandle.push_back( p3);
				face_vHandle.push_back( p1);
				face_vHandle.push_back( p4);
				subdivMesh->add_face( face_vHandle );

				face_vHandle.clear();
				face_vHandle.push_back( p4);
				face_vHandle.push_back( p2);
				face_vHandle.push_back( p5);
				subdivMesh->add_face( face_vHandle );

				face_vHandle.clear();
				face_vHandle.push_back( p3);
				face_vHandle.push_back( p4);
				face_vHandle.push_back( p5);
				subdivMesh->add_face( face_vHandle );
			}
			nSubV = subdivMesh->n_vertices();
			nSubF = subdivMesh->n_faces();
			//OpenMesh::IO::write_mesh( *subdivMesh, "outputnew_1.obj");
		}

		/** 载入连接映射文件*/
		void SparseLocalizedForConnectMap::loadLA()
		{
			IICsMatrix.resize(nSeq, 2 * nE);
			ifstream ConMaps("E:\\GeometryProcessing6.0-20150803\\GeometryProcessing6.0\\src\\verts.txt");
			for(int i = 0; i < nSeq; i++)
			{
				for (int j = 0; j < nE; j++)
				{
					for (int k = 0; k < 2; k++)
					{
						// PyList_SetItem(verts, 13 * nE * i + 13 * j + k, Py_BuildValue("d",IICsData[j][i](k)));
						ConMaps >> IICsMatrix(i, 2 * j + k) ;
					}
				}
			}
		}

		/** 载入锚点数据*/
		void SparseLocalizedForConnectMap::loadAnchorData()
		{
		}

		void SparseLocalizedForConnectMap::loadAnchorData(string str1)
		{
			
			int  int_tmp;
			double dec_tmp;
			std::set<int> idx_set;
			std::vector<int> idx_vec;
			Vertexanchors_.resize(0);
			Vertexanchors_.clear();
			vertexAnchorCords.clear();
			vertexAnchorCords.resize(2);
			vertexAnchorCords[0].resize(0);
			vertexAnchorCords[1].resize(0);
			vertexAnchorCords[0].clear();
			vertexAnchorCords[1].clear();
			std::ifstream ifile(str1, std::ios::in);
			while (ifile >> int_tmp)
			{
				if (int_tmp < 0)
				{
					ifile.close();
					return ;
				}
				else
				{
					Vertexanchors_.push_back(int_tmp);
					vertexanchor v1;
					v1.idx = int_tmp;
					Eigen::VectorXd ll(3);
					ifile >> ll(0);
					ifile >> ll(1);
					ifile >> ll(2);
					v1.cord = ll;
					v1.frameNo = 1;
					vertexAnchorCords[0].push_back(v1);
					vertexAnchorCords[1].push_back(v1);
				}
			}
		}

		/** 初始化回归*/
		void SparseLocalizedForConnectMap::initialCongression()
		{
			edge_regressions.clear();
			PGMesh::EdgeIter ceit0 = mesh_->edges_begin();
			PGMesh::EdgeIter ceit_end0 = mesh_->edges_end();
			for (; ceit0 != ceit_end0; ++ceit0)
			{
				regreCoef egn;
				PGMesh::HalfedgeHandle heh0 = mesh_->halfedge_handle( *ceit0, 0 );
				PGMesh::HalfedgeHandle heh1 = mesh_->halfedge_handle( *ceit0, 1 );
				PGMesh::FaceHandle fh0 = mesh_->face_handle( heh0 );
				PGMesh::FaceHandle fh1 = mesh_->face_handle( heh1 );
				egn.center_node_id = (*ceit0).idx();
				if (fh0.is_valid())
				{
					PGMesh::HalfedgeHandle loop_heh = heh0;
			
					for (int i = 0; i < 2; ++i)
					{
						loop_heh = mesh_->next_halfedge_handle(loop_heh);
						PGMesh::EdgeHandle eh0 = mesh_->edge_handle(loop_heh);
						int idx = eh0.idx();
						egn.one_ring_idx.push_back(idx);
					}
				}
				if (fh1.is_valid())
				{
					PGMesh::HalfedgeHandle loop_heh = heh1;
			
					for (int i = 0; i < 2; ++i)
					{
						loop_heh = mesh_->next_halfedge_handle(loop_heh);
						PGMesh::EdgeHandle eh1 = mesh_->edge_handle(loop_heh);
						int idx = eh1.idx();
						egn.one_ring_idx.push_back(idx);
					}
				}
				if(egn.one_ring_idx.size() == 3)
				{
					egn.one_ring_idx.push_back(egn.one_ring_idx[0]);
				}
				if(egn.one_ring_idx.size() == 2)
				{
					egn.one_ring_idx.push_back(egn.one_ring_idx[0]);
					egn.one_ring_idx.push_back(egn.one_ring_idx[1]);
				}
				if(egn.one_ring_idx.size() == 1)
				{
					egn.one_ring_idx.push_back(egn.one_ring_idx[0]);
					egn.one_ring_idx.push_back(egn.one_ring_idx[0]);
					egn.one_ring_idx.push_back(egn.one_ring_idx[0]);
				}
				
				edge_regressions.push_back(egn);
			}
		}

		/** 估计回归系数*/
		void SparseLocalizedForConnectMap::estimateCoeff()
		{
			ofstream conCoef("E:\\GeometryProcessing6.0-20150803\\GeometryProcessing6.0\\src\\coeffs.txt");
			Eigen::VectorXd y(nSeq);
			Eigen::MatrixXd X(nSeq, 6);
			Eigen::VectorXd unit(nSeq);
			Eigen::VectorXd beta(6);
			Eigen::MatrixXd XX(nSeq - 1, 6);
			Eigen::VectorXd yy(nSeq-1);
			for(int i = 0; i < nSeq; i ++)
			{
				unit(i) = 1;
			}
			for(int edgeNo = 0; edgeNo < nE; edgeNo++)
			{
				int idx0 = edge_regressions[edgeNo].one_ring_idx[0];
				int idx1 = edge_regressions[edgeNo].one_ring_idx[1];
				int idx2 = edge_regressions[edgeNo].one_ring_idx[2];
				int idx3 = edge_regressions[edgeNo].one_ring_idx[3];
				X.col(1) = IICsMatrix.col(3 * idx0 + 1);
				X.col(2) = IICsMatrix.col(3 * idx1 + 1);
				X.col(3) = IICsMatrix.col(3 * idx2 + 1);
				X.col(4) = IICsMatrix.col(3 * idx3 + 1);
				X.col(5) = IICsMatrix.col(3 * edgeNo + 1);
				X.col(0) = unit;
				y = IICsMatrix.col(3 * edgeNo + 2);
				for(int i = 1; i < nSeq; i++)
				{
					X.row(i) = X.row(i) - X.row(0);
					y(i) = y(i) - y(0);
				}
				X.row(0) = X.row(0) - X.row(0);
				y(0) = 0;
				X.col(0) = unit;
				XX = X.block(1, 0, nSeq - 1, 6);
				yy = y.segment(1, nSeq - 1);
				beta = (XX.transpose() * XX).jacobiSvd(ComputeThinU | ComputeThinV).solve(XX.transpose() * yy);
				edge_regressions[edgeNo].one_ring_Coef = beta;
				conCoef<< beta <<endl;
			}
		}

		/** 提取形状基*/
		void SparseLocalizedForConnectMap::getBasis(string str1, string str2)
		{
			ifstream CBase(str1);
			for(int i = 0; i < kCom; i++)
			{
				for (int j = 0; j < 2 * nE; j++)
				{
					CBase >> shapeBasis(i, j);
				}
			}

			ifstream WCof(str2);
			for (int i = 0; i < nSeq; i ++)
			{
				for (int j = 0; j < kCom; j++)
				{
					WCof >> W(i, j);
				}
			}
		}

		/** 载入数据到文件中*/
		void SparseLocalizedForConnectMap::loadData()
		{
			
			//传细分网格面的坐标参数
			ofstream facestxt("faces1.txt");
			for (PGMesh::FaceIter fit = subdivMesh->faces_begin(); fit != subdivMesh->faces_end(); ++ fit)
			{
			    PGMesh::FVIter fvit = subdivMesh->fv_begin(*fit);
				facestxt<< fvit->idx()<<" ";
			    ++ fvit;
				facestxt<< fvit->idx()<<" ";
				++ fvit;
				facestxt<< fvit->idx()<< endl;
			}
			//传第一帧细分网格坐标
			ofstream verttxt("vert1.txt");
			for (int i = 0; i < nSubV; i++)
			{
				PGMesh::Point p = subdivMesh->point(subdivMesh->vertex_handle(i));
				for (int j = 0; j < 3; j++)
				{
					verttxt << p[j] <<" ";
				}
				verttxt<<endl;
			}
			//传原来网格的边在细分网格中点的index
			ofstream indextxt("index1.txt");
			for(PGMesh::EdgeIter eit = mesh_->edges_begin(); eit != mesh_->edges_end(); ++eit)
			{
				PGMesh::VertexHandle v = mesh_->property(ehandle, *eit);
				indextxt << v.idx()<<endl;
			}
		}

		/** 提取细节*/
		void SparseLocalizedForConnectMap::loadDetail()
		{

		}

		/** 局部稀疏分解*/
		void SparseLocalizedForConnectMap::sparseLocalizedDecomp(int kComp)
		{
			// 初始化Python
			//在使用Python系统前，必须使用Py_Initialize对其
			//进行初始化。它会载入Python的内建模块并添加系统路
			//径到模块搜索路径中。这个函数没有返回值，检查系统
			//是否初始化成功需要使用Py_IsInitialized。

			
			//Py_Initialize();

			//// 检查初始化是否成功
			//if ( !Py_IsInitialized() )
			//{
			//	//state = 0;
			//	return ;
			//}

			//// 添加当前路径
			//PyRun_SimpleString("import sys");
			//PyRun_SimpleString("sys.path.append('./')");
			//PyObject *pName,*pModule,*pDict,*pFunc,*verts, *faces, *pnSeq, *pnV, *pnF, *pndim, *firstFrame, *index, *arglist, * results;

			// 载入名为pytest的脚本
			//pName = PyString_FromString("sparse");
			//pModule = PyImport_Import(pName);
			//if ( !pModule )
			//{
			//	//state = 0;
			//	printf("can't find pytest.py");
			//	getchar();
			//	return;
			//}
			// pDict = PyModule_GetDict(pModule);
			//if ( !pDict )
			//{
			//	//state = 0;
			//	return ;
			//}
			// 找出函数名为main的函数
			//pFunc = PyDict_GetItemString(pDict, "transmain");
			//if ( !pFunc || !PyCallable_Check(pFunc) )
			//{
			//printf("can't find function [transmain]");
			////state = 0;
			//getchar();
			//return;
			//}
			//传connection map参数
			//verts = PyList_New(13 * nE * nSeq); // new reference
			
			//传connection map数据, 细分面坐标, 帧数，细分点个数，细分面个数，connection map维数, 边数, 第一帧坐标
			
			/*arglist = PyTuple_New(10);
			PyTuple_SetItem(arglist, 0, verts);
			PyTuple_SetItem(arglist, 1, faces);
			PyTuple_SetItem(arglist, 2, Py_BuildValue("i", nSeq));
			PyTuple_SetItem(arglist, 3, Py_BuildValue("i", nSubV));
			PyTuple_SetItem(arglist, 4, Py_BuildValue("i", nSubF));
			PyTuple_SetItem(arglist, 5, Py_BuildValue("i", 13));
			PyTuple_SetItem(arglist, 6, Py_BuildValue("i", kComp));
			PyTuple_SetItem(arglist, 7, Py_BuildValue("i", nE));
			PyTuple_SetItem(arglist, 8, firstFrame);
			PyTuple_SetItem(arglist, 9, index);*/

			////调用函数
			////results = PyEval_CallObject(pFunc, arglist);
			////if(results == NULL)
			////{
			////	printf("call failed!");
			////	//state = 0;
			////	return;
			////}
			shapeBasis = Eigen::MatrixXd(kComp, 2 * nE);
			W = Eigen::MatrixXd(nSeq, kComp);
			kCom = kComp;
			W_single = Eigen::VectorXd (kComp);
		}

		/** 设置单独的W*/
		void SparseLocalizedForConnectMap::setW_single()
		{
			W_single.setZero();
			ifstream in("edit.txt");
			int idx;
			double values;
			while( in>>idx>>values )
				W_single( idx ) = values;
			in.close();
		}

		void SparseLocalizedForConnectMap::setW_single(string str)
		{
			W_single.setZero();
			ifstream in(str);
			
			for (int i = 0; i < this->kCom; i++)
			{
				in >> W_single(i);
			}
			in.close();
		}
		
		/** 设置multiple的W*/
		void SparseLocalizedForConnectMap::setW_multiple()
		{
			W_single.setZero();
			ifstream in("E:\\GeometryProcessing6.0-20150803\\GeometryProcessing6.0\\src\\edit1.txt" );
			int edit_num = 0;
			in>>edit_num;
			for (int i = 0; i < edit_num; i++)
			{				
				int idx;
				double values;
				in>>idx>>values;
				W_single( idx ) = values;
			}
			in.close();
		}

		/** 根据W_single重构形状*/
		void SparseLocalizedForConnectMap::reconstFormWsingle()
		{
			Eigen::VectorXd l = W_single.transpose() * shapeBasis;
			Eigen::VectorXd edgeError (nE);
			for(int j = 0; j < nE; j++)
			{
				edgeError(j) = l.segment(2 * j, 2).norm();
			}
			char *pPath = new char[30];
			sprintf_s(pPath,30,"editreconst.obj");
			Eigen::VectorXd l0 = DiEdgeDataMatrix.row(0);
			Eigen::VectorXd l1(3 * nV);
			l += l0;
			reconstructionFromDiEdges(l, l1, 0, 0.5, 0);
			visualizeComponent(l1, edgeError, pPath);
			//outMesh(l1, pPath);
		}

		/** 可视化形状基的局部范围*/
		void SparseLocalizedForConnectMap::visualizeBasisLocal (int kth)
		{
			Eigen::VectorXd detaiBasis = shapeBasis.row(kth);
			Eigen::VectorXd edgeError (nE);
			for(int i = 0; i < nE; i++)
			{
				//算第i帧与第0帧的差
				double error = detaiBasis.segment(13 * i, 13).norm();
				edgeError(i) = error;
			}
			char *pPath = new char[50];
			sprintf_s(pPath,50,"%d-thBasisVisual.obj",kth);
			string offname = pPath;
			ofstream offfile(offname);
			delete []pPath;
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
				for (; veit != mesh_->voh_end(vit); ++veit)
				{
					int edgeidx = mesh_->edge_handle(*veit).idx();
					pointerror += edgeError(edgeidx);
					j++;
				}
				//求该点对应相邻边的误差的平均值
				errors[i] = pointerror/j;
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
				unsigned char* colorPtr = colorMap + colorPersent * 3;
				//c.values_[0] = colorPtr[0];
				//c.values_[1] = colorPtr[1];
				//c.values_[2] = colorPtr[2];
				Eigen::Vector3d position = SpatialTemData[i][0];
				offfile << "v" << " " << position( 0) << " " <<position(1) << " " << position(2)  <<" " ;
				offfile <<(int) (colorPtr[0])  <<" "<< int(colorPtr[1])<< " " << int (colorPtr[2])<< endl;
			}


			for (PGMesh::FaceIter fit = mesh_->faces_begin(); fit != mesh_->faces_end(); ++ fit)
			{
				offfile << "f" << " ";
				for (PGMesh::FVIter fvit = mesh_->fv_begin(*fit); fvit != mesh_->fv_end(*fit); ++ fvit)
				{
					offfile << fvit->idx() + 1 << " ";
				}
				offfile << endl;
			}
			offfile.close();
		}
		
		void SparseLocalizedForConnectMap::reconstruction(Eigen::VectorXd edgeComponent, Eigen::VectorXd & cordComponent, int anchorIthSeq, double radius, int ithComp)
		{
			//if( edgeComponent.size() != 13  * nE)
			if( edgeComponent.size() != 13  * nE)
				return;
			/** 首先重构各个标架,先采用求解线性方程组*/
			//Eigen::VectorXd frameCor(9 * nE);
			//Eigen::SparseMatrix<double> Coefficients(36 * nE + 9 ,9 * nE);
			Eigen::SparseMatrix<double> Coefficients(36 * nE + 5 * 9 ,9 * nE);
			Coefficients.setZero();
			Eigen::VectorXd frameCor(36 * nE + 5 * 9);
			//Eigen::VectorXd frameCor(36 * nE + 9);
			frameCor.setZero();
			std::vector<Eigen::Triplet<double> > triple;  
			for (int i = 0; i < nE; i++)
			{
				/*先把与四条相邻边的标架的变换矩阵得到*/
				edge_graph_[i].current_one_ring_local_rot.resize(4);
				for (int ringidx = 0 ; ringidx < 4; ringidx ++)
				{
					Eigen::Vector3d Angaxis; 
					Eigen::Matrix3d rot,trform, rot2Trans;
					Eigen::Matrix3d rot0 = edge_graph_[i].one_ring_local_rot[0][ringidx];
					//Eigen::Matrix3d rotNew;
					//double angle = edgeComponent(5 * i + ringidx);
					//Eigen::Vector3d axis = edge_graph_[i].rotationaxis[anchorIthSeq][ringidx];
					//axis = angle * axis;
					Eigen::Vector3d eulerAg;
					Eigen::Matrix3d U;
					eulerAg = edgeComponent.segment(13 * i + 3 * ringidx, 3);
					/*U.setZero();
					U(0,1) = eulerAg(0);
					U(0,2) = eulerAg(1);
					U(1,2) = eulerAg(2);
					U(1,0) = - U(0,1);
					U(2,0) = - U(0,2);
					U(2,1) = - U(1,2);
					trform = (U.exp());*/
					ERAngle2rotMat(eulerAg, trform);
					//ERAngle2rotMat(eulerAg, rotTrans);
					rot = trform.transpose() * rot0;
					//rot = rotTrans.transpose();
					//rot = GeometryProcess::Math::Geometry::rotationMatrix<Eigen::Matrix3d, Eigen::Vector3d>(axis, angle);
					//rot = rotNew.transpose();
					//EuAng << edgeComponent(13 * i + 3 * ringidx + 0), edgeComponent(13 * i + 3 * ringidx + 1), edgeComponent(13 * i + 3 * ringidx + 2);
					//ERAngle2rotMat(EuAng, rot);
					//rot = edge_graph_[i].one_ring_local_rot[0][ringidx];
					//edge_graph_[i].current_one_ring_local_rot[ringidx] = rot;
					int comparedring =  edge_graph_[i].one_ring_idx[ringidx];
					int local = i;
					triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 0, 9 * comparedring +0, rot(0,0)));
					triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 0, 9 * comparedring +3, rot(0,1)));
					triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 0, 9 * comparedring +6,rot(0,2)));
					triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 0, 9 * local +0, -1));
					triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 1, 9 * comparedring +1,rot(0,0)));
					triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 1, 9 * comparedring +4,rot(0,1)));
					triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 1, 9 * comparedring +7,rot(0,2)));
					triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 1, 9 * local +1, -1));
					triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 2, 9 * comparedring +2,rot(0,0)));
					triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 2, 9 * comparedring +5,rot(0,1)));
					triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 2, 9 * comparedring +8,rot(0,2)));
					triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 2, 9 * local +2,-1));
					triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 3, 9 * comparedring +0,rot(1,0)));
					triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 3, 9 * comparedring +3,rot(1,1)));
					triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 3, 9 * comparedring +6,rot(1,2)));
					triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 3, 9 * local +3, -1));
					triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 4, 9 * comparedring +1, rot(1,0)));
					triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 4, 9 * comparedring +4, rot(1,1)));
					triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 4, 9 * comparedring +7, rot(1,2)));
					triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 4, 9 * local +4, -1));
					triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 5, 9 * comparedring +2, rot(1,0)));
					triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 5, 9 * comparedring +5, rot(1,1)));
					triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 5, 9 * comparedring +8, rot(1,2)));
					triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 5, 9 * local +5, -1));
					triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 6, 9 * comparedring +0, rot(2,0)));
					triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 6, 9 * comparedring +3, rot(2,1)));
					triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 6, 9 * comparedring +6, rot(2,2)));
					triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 6, 9 * local +6, -1));
				    triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 7, 9 * comparedring +1, rot(2,0)));
					triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 7, 9 * comparedring +4, rot(2,1)));
					triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 7, 9 * comparedring +7, rot(2,2)));
					triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 7, 9 * local +7, -1));
					triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 8, 9 * comparedring +2, rot(2,0)));
					triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 8, 9 * comparedring +5, rot(2,1)));
					triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 8, 9 * comparedring +8, rot(2,2)));
					triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 8, 9 * local +8, -1));
				}
			}
			//设置锚点
			Eigen::Matrix3d anchorFrame =  EdgeLocalFrames[anchorEidx][anchorIthSeq];
			int anchorEidx0 = edge_graph_[anchorEidx].one_ring_idx[0];
			int anchorEidx1 = edge_graph_[anchorEidx].one_ring_idx[1];
			int anchorEidx2 = edge_graph_[anchorEidx].one_ring_idx[2];
			int anchorEidx3 = edge_graph_[anchorEidx].one_ring_idx[3];
			Eigen::MatrixXd anchorFrame0 =  EdgeLocalFrames[anchorEidx0][anchorIthSeq];
			Eigen::MatrixXd anchorFrame1 =  EdgeLocalFrames[anchorEidx1][anchorIthSeq];
			Eigen::MatrixXd anchorFrame2 =  EdgeLocalFrames[anchorEidx2][anchorIthSeq];
			Eigen::MatrixXd anchorFrame3 =  EdgeLocalFrames[anchorEidx3][anchorIthSeq];
			int l = 0;
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 0, 9 * anchorEidx + 0, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 1, 9 * anchorEidx + 1, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 2, 9 * anchorEidx + 2, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 3, 9 * anchorEidx + 3, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 4, 9 * anchorEidx + 4, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 5, 9 * anchorEidx + 5, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 6, 9 * anchorEidx + 6, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 7, 9 * anchorEidx + 7, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 8, 9 * anchorEidx + 8, 1));
			frameCor(36 * nE + 0) = anchorFrame(0,0);
			frameCor(36 * nE + 1) = anchorFrame(0,1);
			frameCor(36 * nE + 2) = anchorFrame(0,2);
			frameCor(36 * nE + 3) = anchorFrame(1,0);
			frameCor(36 * nE + 4) = anchorFrame(1,1);
			frameCor(36 * nE + 5) = anchorFrame(1,2);
			frameCor(36 * nE + 6) = anchorFrame(2,0);
			frameCor(36 * nE + 7) = anchorFrame(2,1);
			frameCor(36 * nE + 8) = anchorFrame(2,2);
			l++;
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 0, 9 *  anchorEidx0 + 0, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 1, 9 * anchorEidx0 + 1, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 2, 9 * anchorEidx0 + 2, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 3, 9 * anchorEidx0 + 3, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 4, 9 * anchorEidx0 + 4, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 5, 9 * anchorEidx0 + 5, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 6, 9 * anchorEidx0 + 6, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 7, 9 * anchorEidx0 + 7, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 8, 9 * anchorEidx0 + 8, 1));
			frameCor(36 * nE + l * 9 + 0) = anchorFrame0(0,0);
			frameCor(36 * nE + l * 9 + 1) = anchorFrame0(0,1);
			frameCor(36 * nE + l * 9 + 2) = anchorFrame0(0,2);
			frameCor(36 * nE + l * 9 + 3) = anchorFrame0(1,0);
			frameCor(36 * nE + l * 9 + 4) = anchorFrame0(1,1);
			frameCor(36 * nE + l * 9 + 5) = anchorFrame0(1,2);
			frameCor(36 * nE + l * 9 + 6) = anchorFrame0(2,0);
			frameCor(36 * nE + l * 9 + 7) = anchorFrame0(2,1);
			frameCor(36 * nE + l * 9 + 8) = anchorFrame0(2,2);
			l++;
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 0, 9 *  anchorEidx1 + 0, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 1, 9 * anchorEidx1 + 1, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 2, 9 * anchorEidx1 + 2, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 3, 9 * anchorEidx1 + 3, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 4, 9 * anchorEidx1 + 4, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 5, 9 * anchorEidx1 + 5, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 6, 9 * anchorEidx1 + 6, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 7, 9 * anchorEidx1 + 7, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 8, 9 * anchorEidx1 + 8, 1));
			frameCor(36 * nE + l * 9 + 0) = anchorFrame1(0,0);
			frameCor(36 * nE + l * 9 + 1) = anchorFrame1(0,1);
			frameCor(36 * nE + l * 9 + 2) = anchorFrame1(0,2);
			frameCor(36 * nE + l * 9 + 3) = anchorFrame1(1,0);
			frameCor(36 * nE + l * 9 + 4) = anchorFrame1(1,1);
			frameCor(36 * nE + l * 9 + 5) = anchorFrame1(1,2);
			frameCor(36 * nE + l * 9 + 6) = anchorFrame1(2,0);
			frameCor(36 * nE + l * 9 + 7) = anchorFrame1(2,1);
			frameCor(36 * nE + l * 9 + 8) = anchorFrame1(2,2);
			l++;
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 0, 9 *  anchorEidx2 + 0, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 1, 9 * anchorEidx2 + 1, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 2, 9 * anchorEidx2 + 2, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 3, 9 * anchorEidx2 + 3, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 4, 9 * anchorEidx2 + 4, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 5, 9 * anchorEidx2 + 5, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 6, 9 * anchorEidx2 + 6, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 7, 9 * anchorEidx2 + 7, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 8, 9 * anchorEidx2 + 8, 1));
			frameCor(36 * nE + l * 9 + 0) = anchorFrame2(0,0);
			frameCor(36 * nE + l * 9 + 1) = anchorFrame2(0,1);
			frameCor(36 * nE + l * 9 + 2) = anchorFrame2(0,2);
			frameCor(36 * nE + l * 9 + 3) = anchorFrame2(1,0);
			frameCor(36 * nE + l * 9 + 4) = anchorFrame2(1,1);
			frameCor(36 * nE + l * 9 + 5) = anchorFrame2(1,2);
			frameCor(36 * nE + l * 9 + 6) = anchorFrame2(2,0);
			frameCor(36 * nE + l * 9 + 7) = anchorFrame2(2,1);
			frameCor(36 * nE + l * 9 + 8) = anchorFrame2(2,2);
			l++;
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 0, 9 *  anchorEidx3 + 0, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 1, 9 * anchorEidx3 + 1, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 2, 9 * anchorEidx3 + 2, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 3, 9 * anchorEidx3 + 3, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 4, 9 * anchorEidx3 + 4, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 5, 9 * anchorEidx3 + 5, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 6, 9 * anchorEidx3 + 6, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 7, 9 * anchorEidx3 + 7, 1));
			triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 8, 9 * anchorEidx3 + 8, 1));
			frameCor(36 * nE + l * 9 + 0) = anchorFrame3(0,0);
			frameCor(36 * nE + l * 9 + 1) = anchorFrame3(0,1);
			frameCor(36 * nE + l * 9 + 2) = anchorFrame3(0,2);
			frameCor(36 * nE + l * 9 + 3) = anchorFrame3(1,0);
			frameCor(36 * nE + l * 9 + 4) = anchorFrame3(1,1);
			frameCor(36 * nE + l * 9 + 5) = anchorFrame3(1,2);
			frameCor(36 * nE + l * 9 + 6) = anchorFrame3(2,0);
			frameCor(36 * nE + l * 9 + 7) = anchorFrame3(2,1);
			frameCor(36 * nE + l * 9 + 8) = anchorFrame3(2,2);
			Coefficients.setFromTriplets(triple.begin(), triple.end());  
			Eigen::VectorXd test(9 * nE);
			/*for(int i = 0; i < nE; i++)
			{
				for(int j = 0; j < 3; j++)
				{
					for( int k = 0; k < 3; k++)
					{
						test(9 * i + 3 * j + k) = EdgeLocalFrames[i][0](j, k);
					}
				}
			}*/
			//Eigen::VectorXd testresult = Coefficients * test - frameCor;
			//ofstream tr("testresult.txt");
			//tr << testresult;
			SimplicialCholesky<SparseMatrix<double>> solver;
			//ofstream solvertxt("solve1.txt");
			//solvertxt<<Coefficients.transpose() * Coefficients;
			solver.compute(Coefficients.transpose() * Coefficients);
			if(solver.info()!=Success)
			{
				/// decomposit ion failed
				std::cout<<"Decomposition failed" << std::endl;
				return;
			}
			Eigen::VectorXd x_up; 
			x_up = solver.solve(Coefficients.transpose() * frameCor);
			//ofstream solvert("solve.txt");
			//solvert<< x_up;
			//solvert.close();
			//delete &Coefficients;
			//delete &frameCor;
			/* 试采用牛顿法求解各个标架*/
			/** 根据标架和边长得到每条边的向量，通过边向量重构得到顶点位置*/
			//Eigen::SparseMatrix<double> Coeffs(3 * nE + 3 * Vertexanchors_.size() , 3 * nV);
			 Eigen::SparseMatrix<double> Coeffs(3 * nE + 3 *3 , 3 * nV);
			Coeffs.setZero();
			//Eigen::VectorXd b(3 * nE + 3 * Vertexanchors_.size());
			Eigen::VectorXd b(3 * nE + 3 * 3);
			b.setZero();
			PGMesh::EdgeIter ceit0 = mesh_->edges_begin();
			PGMesh::EdgeIter ceit_end0 = mesh_->edges_end();
			triple.clear();
			for (; ceit0 != ceit_end0; ++ceit0)
			{
				int i = ceit0->idx();
				Eigen::VectorXd Edgeframe(3);
				Edgeframe = EdgeLocalFrames[i][0].row(0);
				Edgeframe << x_up(9 * ceit0->idx() + 0), x_up(9 * ceit0->idx() + 1), x_up(9 * i + 2);
				double edgelength = edgeComponent(13 * ceit0->idx() + 12);
				//double edgelength = ConnectMapData1[ceit0.handle().idx()][0](4);
				Eigen::Vector3d vectorEdge = (edgelength * Edgeframe) / (Edgeframe.norm());
				int from = POINT0_HANDLE(*ceit0).idx();
				int to = POINT1_HANDLE(*ceit0).idx();
				//Eigen::Vector3d vectorEdge = SpatialTemData[to][0] - SpatialTemData[from][0];
				triple.push_back(Eigen::Triplet<double>(3 * i + 0, 3 * to + 0, 1));
				triple.push_back(Eigen::Triplet<double>(3 * i + 0, 3 * from + 0, -1));
				b(3 * i + 0) = vectorEdge(0);
				triple.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * to + 1, 1));
				triple.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * from + 1, -1));
				b(3 * i + 1) = vectorEdge(1);
				triple.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * to + 2, 1));
				triple.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * from + 2, -1));
				b(3 * i + 2) = vectorEdge(2);
			}
			//for (int i = 0; i < Vertexanchors_.size(); i ++)
				for (int i = 0; i <3; i ++)
			{
				triple.push_back(Eigen::Triplet<double>(3 * nV + 3 * i + 0, 3 * Vertexanchors_[i] + 0, 1));
				b(3 * nV + 3 * i + 0) = SpatialTemData[0][anchorIthSeq](0);
				triple.push_back(Eigen::Triplet<double>(3 * nV + 3 * i + 1, 3 * Vertexanchors_[i] + 1, 1));
				b(3 * nV + 3 * i + 1) = SpatialTemData[0][anchorIthSeq](1);
				triple.push_back(Eigen::Triplet<double>(3 * nV + 3 * i + 2, 3 * Vertexanchors_[i] + 2, 1));
				b(3 * nV + 3 * i + 2) = SpatialTemData[0][anchorIthSeq](2);
			}
			/*for (int i = 0; i < nV; i++)
			{
				test3( 3 * i + 0 ) = SpatialTemData[i][anchorIthSeq](0);
				test3( 3 * i + 1 ) = SpatialTemData[i][anchorIthSeq](1);
				test3( 3 * i + 2 ) = SpatialTemData[i][anchorIthSeq](2);
			}*/
			/*ofstream test3out("testnV.txt");
			test3out << Coeffs * test3 - b;
			test3out.close();*/
			Coeffs.setFromTriplets(triple.begin(), triple.end()); 
			SimplicialCholesky<SparseMatrix<double>> solver1;
			solver1.compute(Coeffs.transpose() * Coeffs);
			if(solver1.info()!=Success)
			{
				/// decomposit ion failed
				std::cout<<"Decomposition failed" << std::endl;
				return;
			}
			Eigen::VectorXd y_up; 
			y_up = solver1.solve(Coeffs.transpose() * b);
			//Eigen::VectorXd test1 = Coeffs * y_up - b;
			//ofstream outy("outy_up.txt");
			//outy<< test1;
			/*output the aligned mesh to a file: seqNo.obj*/
			char *pPath = new char[30];
			sprintf_s(pPath,30,"%d-thConMapFrame.obj", ithComp);
			string offname = pPath;
			ofstream offfile(offname);
			delete []pPath;
			offfile << "COFF"<<endl;
				offfile << nV << " " << nF<< " "<< 0 <<endl;
			Eigen::VectorXd errorVec(nV);
			Eigen::Vector3d error1;
			Eigen::Vector3d error2;
			double maxerror = 0;
			double minerror = 1000;
			for(int i = 0; i < nV; i++)
			{
				 error1 = y_up.segment(3*i, 3) ;
				 error2 = SpatialTemData[i][anchorIthSeq];
				 errorVec(i) = (( y_up.segment(3*i, 3) - SpatialTemData[i][anchorIthSeq]).norm());
				if(errorVec(i) > maxerror)
					maxerror = errorVec(i);
				if(errorVec(i) < minerror)
					minerror = errorVec(i);
			}
			double meanerror = errorVec.mean();
			//double maxerror = errorVec.
			//ofstream offfile("connectOut.obj");
			for (int i = 0; i< nV; i++)
			{
				int color =(int) ((errorVec(i))  * 255/radius);
				offfile << "v" << " " << y_up(3 * i + 0) << " " <<y_up(3 * i + 1) << " " << y_up(3 * i + 2)  <<" " << 255 <<" "<< 255<< " " << 255-color<< endl;
			}
			for (PGMesh::FaceIter fit = mesh_->faces_begin(); fit != mesh_->faces_end(); ++ fit)
			{
				offfile << "f" << " ";
				for (PGMesh::FVIter fvit = mesh_->fv_begin(*fit); fvit != mesh_->fv_end(*fit); ++ fvit)
				{
					offfile << fvit->idx() + 1 << " ";
				}
				offfile << endl;
			}
			offfile << "# maxerror " << maxerror/radius << endl;
			offfile << "# meanerror "<< meanerror/radius<<endl;
			offfile.close();
		}

		/** 利用简化线性方程求解重构网格*/
		void SparseLocalizedForConnectMap::reconstruction1(Eigen::VectorXd edgeComponent, Eigen::VectorXd & cordComponent, int anchorIthSeq, double radius, int ithComp)
		{
			//if( edgeComponent.size() != 13  * nE)
			if( edgeComponent.size() != 17  * nE)
				return;
			/** 首先重构各个标架,先采用求解线性方程组*/
			Eigen::VectorXd frameCor1(12 * nE+ 3 * edgeanchors_.size());
			Eigen::VectorXd frameCor2(12 * nE+ 3 * edgeanchors_.size());
			Eigen::VectorXd frameCor3(12 * nE+ 3 * edgeanchors_.size());
			Eigen::SparseMatrix<double> Coefficients(12 * nE + 3 * edgeanchors_.size() ,3 * nE);
			frameCor1.setZero();
			frameCor2.setZero();
			frameCor3.setZero();

			/** 根据标架和边长得到每条边的向量，通过边向量重构得到顶点位置*/
			Eigen::SparseMatrix<double> Coeffs(3 * nE + 3 * Vertexanchors_.size(), 3 * nV);

			std::vector<Eigen::Triplet<double> > triple;  
			for (int i = 0; i < nE; i++)
			{
				/*先把与四条相邻边的标架的变换矩阵得到*/
				edge_graph_[i].current_one_ring_local_rot.resize(4);
				for (int ringidx = 0 ; ringidx < 4; ringidx ++)
				{
					Eigen::Vector3d Angaxis; 
					Eigen::Matrix3d rot, rotTrans;
					//Eigen::Matrix3d rotNew;
					//double angle = edgeComponent(5 * i + ringidx);
					//Eigen::Vector3d axis = edge_graph_[i].rotationaxis[anchorIthSeq][ringidx];
					//axis = angle * axis;
					Eigen::Vector4d quat;
					quat = edgeComponent.segment(17 * i + 4 * ringidx, 4);
					Math::Geometry::axisAndAngle(quat, Angaxis);
					rotTrans = Math::Geometry::rotationMatrix<Eigen::Matrix3d,Eigen::Vector3d>(Angaxis);
					Eigen::Matrix3d rotFirst = edge_graph_[i].one_ring_local_rot[0][ringidx];
					rot = rotTrans * rotFirst;
					triple.push_back(Eigen::Triplet<double>(12 * i + 3 * ringidx + 0, 3 * i +0, rot(0,0)));
					triple.push_back(Eigen::Triplet<double>(12 * i + 3 * ringidx + 0, 3 * i +1, rot(0,1)));
					triple.push_back(Eigen::Triplet<double>(12 * i + 3 * ringidx + 0, 3 * i +2,rot(0,2)));
					triple.push_back(Eigen::Triplet<double>(12 * i + 3 * ringidx + 0, 3 * edge_graph_[i].one_ring_idx[ringidx] +0, -1));
					triple.push_back(Eigen::Triplet<double>(12 * i + 3 * ringidx + 1, 3 * i +0,rot(1,0)));
					triple.push_back(Eigen::Triplet<double>(12 * i + 3 * ringidx + 1, 3 * i +1,rot(1,1)));
					triple.push_back(Eigen::Triplet<double>(12 * i + 3 * ringidx + 1, 3 * i +2,rot(1,2)));
					triple.push_back(Eigen::Triplet<double>(12 * i + 3 * ringidx + 1, 3 * edge_graph_[i].one_ring_idx[ringidx] +1, -1));
					triple.push_back(Eigen::Triplet<double>(12 * i + 3 * ringidx + 2, 3 * i +0, rot(2,0)));
					triple.push_back(Eigen::Triplet<double>(12 * i + 3 * ringidx + 2, 3 * i +1, rot(2,1)));
					triple.push_back(Eigen::Triplet<double>(12 * i + 3 * ringidx + 2, 3 * i +2, rot(2,2)));
					triple.push_back(Eigen::Triplet<double>(12 * i + 3 * ringidx + 2, 3 * edge_graph_[i].one_ring_idx[ringidx] +2, -1));
				}
			}
			//设置锚点
			for(int l = 0; l < edgeanchors_.size(); l++)
			{
				Eigen::Matrix3d anchorFrame =  EdgeLocalFrames[edgeanchors_[l]][anchorIthSeq];
				int anchorEidx0 = edge_graph_[edgeanchors_[l]].one_ring_idx[0];
				int anchorEidx1 = edge_graph_[edgeanchors_[l]].one_ring_idx[1];
				int anchorEidx2 = edge_graph_[edgeanchors_[l]].one_ring_idx[2];
				int anchorEidx3 = edge_graph_[edgeanchors_[l]].one_ring_idx[3];
				triple.push_back(Eigen::Triplet<double>(12 * nE + l * 3 + 0, 3 * anchorEidx + 0, 1));
				triple.push_back(Eigen::Triplet<double>(12 * nE + l * 3 + 1, 3 * anchorEidx + 1, 1));
				triple.push_back(Eigen::Triplet<double>(12 * nE + l * 3 + 2, 3 * anchorEidx + 2, 1));
				frameCor1(12 * nE + 3 * l + 0) = anchorFrame(0,0);
				frameCor1(12 * nE + 3 * l + 1) = anchorFrame(1,0);
				frameCor1(12 * nE + 3 * l + 2) = anchorFrame(2,0);
				frameCor2(12 * nE + 3 * l + 0) = anchorFrame(0,1);
				frameCor2(12 * nE + 3 * l + 1) = anchorFrame(1,1);
				frameCor2(12 * nE + 3 * l + 2) = anchorFrame(2,1);
				frameCor3(12 * nE + 3 * l + 0) = anchorFrame(0,2);
				frameCor3(12 * nE + 3 * l + 1) = anchorFrame(1,2);
				frameCor3(12 * nE + 3 * l + 2) = anchorFrame(2,2);
			}
			
			Coefficients.setFromTriplets(triple.begin(), triple.end());  
			SimplicialCholesky<SparseMatrix<double>> solver;
			//ofstream solvertxt("solve1.txt");
			//solvertxt<<Coefficients.transpose() * Coefficients;
			solver.compute(Coefficients.transpose() * Coefficients);
			if(solver.info()!=Success)
			{
				/// decomposit ion failed
				std::cout<<"Decomposition failed" << std::endl;
				return;
			}
			Eigen::VectorXd x_up1; 
			x_up1 = solver.solve(Coefficients.transpose() * frameCor1);
			Eigen::VectorXd x_up2; 
			x_up2 = solver.solve(Coefficients.transpose() * frameCor2);
			Eigen::VectorXd x_up3; 
			x_up3 = solver.solve(Coefficients.transpose() * frameCor3);
			Eigen::VectorXd x_up(9 * nE);
			for(int i = 0; i < nE; i++)
			{
				x_up(9 * i + 0) = x_up1(3 * i + 0);
				x_up(9 * i + 3) = x_up1(3 * i + 1);
				x_up(9 * i + 6) = x_up1(3 * i + 2);
				x_up(9 * i + 1) = x_up2(3 * i + 0);
				x_up(9 * i + 4) = x_up2(3 * i + 1);
				x_up(9 * i + 7) = x_up2(3 * i + 2);
				x_up(9 * i + 2) = x_up3(3 * i + 0);
				x_up(9 * i + 5) = x_up3(3 * i + 1);
				x_up(9 * i + 8) = x_up3(3 * i + 2);
			}
			Coeffs.setZero();
			Eigen::VectorXd b(3 * nE + 3 * Vertexanchors_.size());
			b.setZero();
			PGMesh::EdgeIter ceit0 = mesh_->edges_begin();
			PGMesh::EdgeIter ceit_end0 = mesh_->edges_end();
			triple.clear();
			for (; ceit0 != ceit_end0; ++ceit0)
			{
				int i = ceit0->idx();
				Eigen::VectorXd Edgeframe(3);
				Edgeframe = EdgeLocalFrames[i][0].row(0);
				Edgeframe << x_up(9 * ceit0->idx() + 0), x_up(9 * ceit0->idx() + 1), x_up(9 * i + 2);
				double edgelength = edgeComponent(17 * ceit0->idx() + 16);
				//double edgelength = ConnectMapData1[ceit0.handle().idx()][0](4);
				Eigen::Vector3d vectorEdge = (edgelength * Edgeframe) / (Edgeframe.norm());
				int from = POINT0_HANDLE(*ceit0).idx();
				int to = POINT1_HANDLE(*ceit0).idx();
				//Eigen::Vector3d vectorEdge = SpatialTemData[to][0] - SpatialTemData[from][0];
				triple.push_back(Eigen::Triplet<double>(3 * i + 0, 3 * to + 0, 1));
				triple.push_back(Eigen::Triplet<double>(3 * i + 0, 3 * from + 0, -1));
				b(3 * i + 0) = vectorEdge(0);
				triple.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * to + 1, 1));
				triple.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * from + 1, -1));
				b(3 * i + 1) = vectorEdge(1);
				triple.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * to + 2, 1));
				triple.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * from + 2, -1));
				b(3 * i + 2) = vectorEdge(2);
			}
			for (int i = 0; i < Vertexanchors_.size(); i ++)
			{
				triple.push_back(Eigen::Triplet<double>(3 * nE + 3 * i + 0, 3 * Vertexanchors_[i] + 0, 1));
				b(3 * nE + 3 * i + 0) = SpatialTemData[0][anchorIthSeq](0);
				triple.push_back(Eigen::Triplet<double>(3 * nE + 3 * i + 1, 3 * Vertexanchors_[i] + 1, 1));
				b(3 * nE + 3 * i + 1) = SpatialTemData[0][anchorIthSeq](1);
				triple.push_back(Eigen::Triplet<double>(3 * nE + 3 * i + 2, 3 * Vertexanchors_[i] + 2, 1));
				b(3 * nE + 3 * i + 2) = SpatialTemData[0][anchorIthSeq](2);
			}

			//Eigen::VectorXd test3(3 * nV);
			/*for (int i = 0; i < nV; i++)
			{
				test3( 3 * i + 0 ) = SpatialTemData[i][anchorIthSeq](0);
				test3( 3 * i + 1 ) = SpatialTemData[i][anchorIthSeq](1);
				test3( 3 * i + 2 ) = SpatialTemData[i][anchorIthSeq](2);
			}*/
			/*ofstream test3out("testnV.txt");
			test3out << Coeffs * test3 - b;
			test3out.close();*/
			Coeffs.setFromTriplets(triple.begin(), triple.end()); 
			SimplicialCholesky<SparseMatrix<double>> solver1;
			solver1.compute(Coeffs.transpose() * Coeffs);
			if(solver1.info()!=Success)
			{
				/// decomposit ion failed
				std::cout<<"Decomposition failed" << std::endl;
				return;
			}
			Eigen::VectorXd y_up; 
			y_up = solver1.solve(Coeffs.transpose() * b);
			//Eigen::VectorXd test1 = Coeffs * y_up - b;
			//ofstream outy("outy_up.txt");
			//outy<< test1;
			/*output the aligned mesh to a file: seqNo.obj*/
			char *pPath = new char[30];
			sprintf_s(pPath,30,"%d-thConMapFrame.obj", ithComp);
			string offname = pPath;
			ofstream offfile(offname);
			delete []pPath;
			offfile << "COFF"<<endl;
				offfile << nV << " " << nF<< " "<< 0 <<endl;
			Eigen::VectorXd errorVec(nV);
			Eigen::Vector3d error1;
			Eigen::Vector3d error2;
			double maxerror = 0;
			double minerror = 1000;
			for(int i = 0; i < nV; i++)
			{
				 error1 = y_up.segment(3*i, 3) ;
				 error2 = SpatialTemData[i][anchorIthSeq];
				 errorVec(i) = (( y_up.segment(3*i, 3) - SpatialTemData[i][anchorIthSeq]).norm());
				if(errorVec(i) > maxerror)
					maxerror = errorVec(i);
				if(errorVec(i) < minerror)
					minerror = errorVec(i);
			}
			double meanerror = errorVec.mean();
			//double maxerror = errorVec.
			//ofstream offfile("connectOut.obj");
			for (int i = 0; i< nV; i++)
			{
				int color =(int) ((errorVec(i))  * 255/radius);
				offfile << "v" << " " << y_up(3 * i + 0) << " " <<y_up(3 * i + 1) << " " << y_up(3 * i + 2)  <<" " << 255 <<" "<< 255<< " " << 255-color<< endl;
			}
			for (PGMesh::FaceIter fit = mesh_->faces_begin(); fit != mesh_->faces_end(); ++ fit)
			{
				offfile << "f" << " ";
				for (PGMesh::FVIter fvit = mesh_->fv_begin(*fit); fvit != mesh_->fv_end(*fit); ++ fvit)
				{
					offfile << fvit->idx() + 1 << " ";
				}
				offfile << endl;
			}
			offfile << "# maxerror " << maxerror/radius << endl;
			offfile << "# meanerror "<< meanerror/radius<<endl;
			offfile.close();
		}

		/** 利用牛顿法重构网格*/
		void SparseLocalizedForConnectMap::reconstructionNewton(Eigen::VectorXd edgeComponent, Eigen::VectorXd & cordComponent, int anchorIthSeq, double radius, int ithComp)
		{
			//if( edgeComponent.size() != 13  * nE)
			if( edgeComponent.size() != 17  * nE)
				return;
			/** 首先重构各个标架,先采用求解线性方程组*/
			//Eigen::VectorXd frameCor(9 * nE);
			//Eigen::SparseMatrix<double> Coefficients(36 * nE + 9 ,9 * nE);
			Eigen::SparseMatrix<double> Coefficients(36 * nE + 9 * edgeanchors_.size() ,9 * nE);
			Eigen::SparseMatrix<double> A_rot(6 * nE, 9 * nE);
			Eigen::VectorXd f_rot( 6 * nE);
			Coefficients.setZero();
			Eigen::VectorXd frameCor(36 * nE + 9 * edgeanchors_.size());
			//Eigen::VectorXd frameCor(36 * nE + 9);
			frameCor.setZero();
			f_rot.setZero();
			std::vector<Eigen::Triplet<double> > triple;  
			Eigen::VectorXd x_up(9 *  nE);
			/*initialization of the x_up*/
			for(int i = 0; i < nE; i++)
			{
				//the rotation matrix is the identical matrix
				x_up(i * 9 + 0) = 1;
				x_up(i * 9 + 1) = 0;
				x_up(i * 9 + 2) = 0;
				x_up(i * 9 + 3) = 0;
				x_up(i * 9 + 4) = 1;
				x_up(i * 9 + 5) = 0;
				x_up(i * 9 + 6) = 0;
				x_up(i * 9 + 7) = 0;
				x_up(i * 9 + 8) = 1;
			}

			int currentIter = 0;
			int totalIter = 6;

			while(1)
			{
				if (currentIter >= totalIter)
					break;
				currentIter++;
				/*Set the two sparse matrices zero*/
				A_rot.setZero();
				Coefficients.setZero();
				/*calculate the value of f_rot*/
				for(int i = 0; i < nE; i++)
				{
					f_rot(i * 6 + 0) = x_up(i * 9 + 0) * x_up(i * 9 + 3) + x_up(i * 9 + 1) * x_up(i * 9 + 4) + x_up(i * 9 + 2) * x_up(i * 9 + 5);//
					f_rot(i * 6 + 1) = x_up(i * 9 + 0) * x_up(i * 9 + 6) + x_up(i * 9 + 1) * x_up(i * 9 + 7) + x_up(i * 9 + 2) * x_up(i * 9 + 8);//
					f_rot(i * 6 + 2) = x_up(i * 9 + 3) * x_up(i * 9 + 6) + x_up(i * 9 + 4) * x_up(i * 9 + 7) + x_up(i * 9 + 5) * x_up(i * 9 + 8);//
					f_rot(i * 6 + 3) = x_up(i * 9 + 0) * x_up(i * 9 + 0) + x_up(i * 9 + 1) * x_up(i * 9 + 1) + x_up(i * 9 + 2) * x_up(i * 9 + 2) - 1;
					f_rot(i * 6 + 4) = x_up(i * 9 + 3) * x_up(i * 9 + 3) + x_up(i * 9 + 4) * x_up(i * 9 + 4) + x_up(i * 9 + 5) * x_up(i * 9 + 5) - 1;
					f_rot(i * 6 + 5) = x_up(i * 9 + 6) * x_up(i * 9 + 6) + x_up(i * 9 + 7) * x_up(i * 9 + 7) + x_up(i * 9 + 8) * x_up(i * 9 + 8) - 1;
				}

				/*calculate the Jacobian matrtix of f_rot*/
				for(int i = 0; i < nE; i++)
				{
					triple.push_back(Eigen::Triplet<double>(i * 6 + 0, i * 9 + 0, x_up(i * 9 + 3)));
					triple.push_back(Eigen::Triplet<double>(i * 6 + 0, i * 9 + 3, x_up(i * 9 + 0)));
					triple.push_back(Eigen::Triplet<double>(i * 6 + 0, i * 9 + 1, x_up(i * 9 + 4)));
					triple.push_back(Eigen::Triplet<double>(i * 6 + 0, i * 9 + 4, x_up(i * 9 + 1)));
					triple.push_back(Eigen::Triplet<double>(i * 6 + 0, i * 9 + 2, x_up(i * 9 + 5)));
					triple.push_back(Eigen::Triplet<double>(i * 6 + 0, i * 9 + 5, x_up(i * 9 + 2)));
					triple.push_back(Eigen::Triplet<double>(i * 6 + 1, i * 9 + 0, x_up(i * 9 + 6)));
					triple.push_back(Eigen::Triplet<double>(i * 6 + 1, i * 9 + 6, x_up(i * 9 + 0)));
					triple.push_back(Eigen::Triplet<double>(i * 6 + 1, i * 9 + 1, x_up(i * 9 + 7)));
					triple.push_back(Eigen::Triplet<double>(i * 6 + 1, i * 9 + 7, x_up(i * 9 + 1)));
					triple.push_back(Eigen::Triplet<double>(i * 6 + 1, i * 9 + 2, x_up(i * 9 + 8)));
					triple.push_back(Eigen::Triplet<double>(i * 6 + 1, i * 9 + 8, x_up(i * 9 + 2)));
					triple.push_back(Eigen::Triplet<double>(i * 6 + 2, i * 9 + 3, x_up(i * 9 + 6)));
					triple.push_back(Eigen::Triplet<double>(i * 6 + 2, i * 9 + 6, x_up(i * 9 + 3)));
					triple.push_back(Eigen::Triplet<double>(i * 6 + 2, i * 9 + 4, x_up(i * 9 + 7)));
					triple.push_back(Eigen::Triplet<double>(i * 6 + 2, i * 9 + 7, x_up(i * 9 + 4)));
					triple.push_back(Eigen::Triplet<double>(i * 6 + 2, i * 9 + 5, x_up(i * 9 + 8)));
					triple.push_back(Eigen::Triplet<double>(i * 6 + 2, i * 9 + 8, x_up(i * 9 + 5)));
					triple.push_back(Eigen::Triplet<double>(i * 6 + 3, i * 9 + 0, 2 * x_up(i * 9 + 0)));
					triple.push_back(Eigen::Triplet<double>(i * 6 + 3, i * 9 + 1, 2 * x_up(i * 9 + 1)));
					triple.push_back(Eigen::Triplet<double>(i * 6 + 3, i * 9 + 2, 2 * x_up(i * 9 + 2)));
					triple.push_back(Eigen::Triplet<double>(i * 6 + 4, i * 9 + 3, 2 * x_up(i * 9 + 3)));
					triple.push_back(Eigen::Triplet<double>(i * 6 + 4, i * 9 + 4, 2 * x_up(i * 9 + 4)));
					triple.push_back(Eigen::Triplet<double>(i * 6 + 4, i * 9 + 5, 2 * x_up(i * 9 + 5)));
					triple.push_back(Eigen::Triplet<double>(i * 6 + 5, i * 9 + 6, 2 * x_up(i * 9 + 6)));
					triple.push_back(Eigen::Triplet<double>(i * 6 + 5, i * 9 + 7, 2 * x_up(i * 9 + 7)));
					triple.push_back(Eigen::Triplet<double>(i * 6 + 5, i * 9 + 8, 2 * x_up(i * 9 + 8)));
				}
				A_rot.setFromTriplets(triple.begin(), triple.end());  
				triple.clear();

				for (int i = 0; i < nE; i++)
				{
					/*先把与四条相邻边的标架的变换矩阵得到*/
					edge_graph_[i].current_one_ring_local_rot.resize(4);
					for (int ringidx = 0 ; ringidx < 4; ringidx ++)
					{
						Eigen::Vector3d Angaxis; 
						Eigen::Matrix3d rot, rotTrans;
						//Eigen::Matrix3d rotNew;
						//double angle = edgeComponent(5 * i + ringidx);
						//Eigen::Vector3d axis = edge_graph_[i].rotationaxis[anchorIthSeq][ringidx];
						//axis = angle * axis;
						Eigen::Vector4d quat;
						quat = edgeComponent.segment(17 * i + 4 * ringidx, 4);
						Math::Geometry::axisAndAngle(quat, Angaxis);
						rotTrans = Math::Geometry::rotationMatrix<Eigen::Matrix3d,Eigen::Vector3d>(Angaxis);
						Eigen::Matrix3d rotFirst = edge_graph_[i].one_ring_local_rot[0][ringidx];
						rot = rotTrans * rotFirst;
						Eigen::Matrix3d frame, frame1, error;
						frame << x_up(9 * i +0) , x_up(9 * i +1) , x_up(9 * i +2) ,
							x_up(9 * i +3) , x_up(9 * i +4) , x_up(9 * i +5) ,
							x_up(9 * i +6) , x_up(9 * i +7) , x_up(9 * i +8) ;
						frame = rot * frame;
						int new1 = 9 * edge_graph_[i].one_ring_idx[ringidx] ;
						frame1 << x_up(new1 +0) , x_up(new1 +1) , x_up(new1 +2) ,
							x_up(new1 +3) , x_up(new1 +4) , x_up(new1 +5) ,
							x_up(new1 +6) , x_up(new1 +7) , x_up(new1 +8) ;
						error = frame - frame1;
						frameCor(36 * i + 9 * ringidx + 0) = error(0,0);
						frameCor(36 * i + 9 * ringidx + 1) = error(0,1);
						frameCor(36 * i + 9 * ringidx + 2) = error(0,2);
						frameCor(36 * i + 9 * ringidx + 3) = error(1,0);
						frameCor(36 * i + 9 * ringidx + 4) = error(1,1);
						frameCor(36 * i + 9 * ringidx + 5) = error(1,2);
						frameCor(36 * i + 9 * ringidx + 6) = error(2,0);
						frameCor(36 * i + 9 * ringidx + 7) = error(2,1);
						frameCor(36 * i + 9 * ringidx + 8) = error(2,2);
						//rot = GeometryProcess::Math::Geometry::rotationMatrix<Eigen::Matrix3d, Eigen::Vector3d>(axis, angle);
						//rot = rotNew.transpose();
						//EuAng << edgeComponent(13 * i + 3 * ringidx + 0), edgeComponent(13 * i + 3 * ringidx + 1), edgeComponent(13 * i + 3 * ringidx + 2);
						//ERAngle2rotMat(EuAng, rot);
						//rot = edge_graph_[i].one_ring_local_rot[0][ringidx];
						//edge_graph_[i].current_one_ring_local_rot[ringidx] = rot;
						triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 0, 9 * i +0, rot(0,0)));
						triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 0, 9 * i +3, rot(0,1)));
						triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 0, 9 * i +6,rot(0,2)));
						triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 0, 9 * edge_graph_[i].one_ring_idx[ringidx] +0, -1));
						triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 1, 9 * i +1,rot(0,0)));
						triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 1, 9 * i +4,rot(0,1)));
						triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 1, 9 * i +7,rot(0,2)));
						triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 1, 9 * edge_graph_[i].one_ring_idx[ringidx] +1, -1));
						triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 2, 9 * i +2,rot(0,0)));
						triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 2, 9 * i +5,rot(0,1)));
						triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 2, 9 * i +8,rot(0,2)));
						triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 2, 9 * edge_graph_[i].one_ring_idx[ringidx] +2,-1));
						triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 3, 9 * i +0,rot(1,0)));
						triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 3, 9 * i +3,rot(1,1)));
						triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 3, 9 * i +6,rot(1,2)));
						triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 3, 9 * edge_graph_[i].one_ring_idx[ringidx] +3, -1));
						triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 4, 9 * i +1, rot(1,0)));
						triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 4, 9 * i +4, rot(1,1)));
						triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 4, 9 * i +7, rot(1,2)));
						triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 4, 9 * edge_graph_[i].one_ring_idx[ringidx] +4, -1));
						triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 5, 9 * i +2, rot(1,0)));
						triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 5, 9 * i +5, rot(1,1)));
						triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 5, 9 * i +8, rot(1,2)));
						triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 5, 9 * edge_graph_[i].one_ring_idx[ringidx] +5, -1));
						triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 6, 9 * i +0, rot(2,0)));
						triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 6, 9 * i +3, rot(2,1)));
						triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 6, 9 * i +6, rot(2,2)));
						triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 6, 9 * edge_graph_[i].one_ring_idx[ringidx] +6, -1));
						triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 7, 9 * i +1, rot(2,0)));
						triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 7, 9 * i +4, rot(2,1)));
						triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 7, 9 * i +7, rot(2,2)));
						triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 7, 9 * edge_graph_[i].one_ring_idx[ringidx] +7, -1));
						triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 8, 9 * i +2, rot(2,0)));
						triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 8, 9 * i +5, rot(2,1)));
						triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 8, 9 * i +8, rot(2,2)));
						triple.push_back(Eigen::Triplet<double>(36 * i + 9 * ringidx + 8, 9 * edge_graph_[i].one_ring_idx[ringidx] +8, -1));
					}

				}
				//设置锚点
				for (int l = 0; l < edgeanchors_.size(); l++)
				{
					Eigen::Matrix3d anchorFrame =  EdgeLocalFrames[edgeanchors_[l]][anchorIthSeq];
					triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 0, 9 * edgeanchors_[l] + 0, 1));
					triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 1, 9 * edgeanchors_[l] + 1, 1));
					triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 2, 9 * edgeanchors_[l] + 2, 1));
					triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 3, 9 * edgeanchors_[l] + 3, 1));
					triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 4, 9 * edgeanchors_[l] + 4, 1));
					triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 5, 9 * edgeanchors_[l] + 5, 1));
					triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 6, 9 * edgeanchors_[l] + 6, 1));
					triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 7, 9 * edgeanchors_[l] + 7, 1));
					triple.push_back(Eigen::Triplet<double>(36 * nE + l * 9 + 8, 9 * edgeanchors_[l] + 8, 1));
					frameCor(36 * nE + l * 9 + 0) = x_up( 9 * edgeanchors_[l] + 0) - anchorFrame(0,0);
					frameCor(36 * nE + l * 9 + 1) = x_up( 9 * edgeanchors_[l] + 1) - anchorFrame(0,1);
					frameCor(36 * nE + l * 9 + 2) = x_up( 9 * edgeanchors_[l] + 2) - anchorFrame(0,2);
					frameCor(36 * nE + l * 9 + 3) = x_up( 9 * edgeanchors_[l] + 3) - anchorFrame(1,0);
					frameCor(36 * nE + l * 9 + 4) = x_up( 9 * edgeanchors_[l] + 4) - anchorFrame(1,1);
					frameCor(36 * nE + l * 9 + 5) = x_up( 9 * edgeanchors_[l] + 5) - anchorFrame(1,2);
					frameCor(36 * nE + l * 9 + 6) = x_up( 9 * edgeanchors_[l] + 6) - anchorFrame(2,0);
					frameCor(36 * nE + l * 9 + 7) = x_up( 9 * edgeanchors_[l] + 7) - anchorFrame(2,1);
					frameCor(36 * nE + l * 9 + 8) = x_up( 9 * edgeanchors_[l] + 8) - anchorFrame(2,2);
				}

				Coefficients.setFromTriplets(triple.begin(), triple.end());  
				triple.clear();
				double w_rot = 1.0;
				double w_reg = 10.0;
				
				SimplicialCholesky<SparseMatrix<double>> solver;
				//ofstream solvertxt("solve1.txt");
				//solvertxt<<Coefficients.transpose() * Coefficients;
				solver.compute(w_rot * A_rot.transpose() * A_rot + w_reg * Coefficients.transpose() * Coefficients);
				if(solver.info()!=Success)
				{
					/// decomposit ion failed
					std::cout<<"Decomposition failed" << std::endl;
					return;
				}
				VectorXd x_update;                                                                                                                 
				VectorXd Axb ;

				 
				 Axb = -1.0f * (w_rot  * A_rot.transpose() * f_rot + w_reg * Coefficients.transpose() * frameCor );
				 //Axb = -1.0f * (w_rot  * A_rot.transpose() * f_rot + w_con * A_con.transpose() * f_con);

				 x_update = solver.solve(Axb);
				 //b_reg = b_reg + x_update;
			
				x_up = x_up + x_update;
			}
			
			//ofstream solvert("solve.txt");
			//solvert<< x_up;
			//solvert.close();
			//delete &Coefficients;
			//delete &frameCor;
			/* 试采用牛顿法求解各个标架*/
			/** 根据标架和边长得到每条边的向量，通过边向量重构得到顶点位置*/
			Eigen::SparseMatrix<double> Coeffs(3 * nE + 3 * Vertexanchors_.size() , 3 * nV);
			 
			Coeffs.setZero();
			Eigen::VectorXd b(3 * nE + 3 * Vertexanchors_.size());
			b.setZero();
			PGMesh::EdgeIter ceit0 = mesh_->edges_begin();
			PGMesh::EdgeIter ceit_end0 = mesh_->edges_end();
			triple.clear();
			for (; ceit0 != ceit_end0; ++ceit0)
			{
				int i = ceit0->idx();
				Eigen::VectorXd Edgeframe(3);
				Edgeframe = EdgeLocalFrames[i][0].row(0);
				Edgeframe << x_up(9 * ceit0->idx() + 0), x_up(9 * ceit0->idx() + 1), x_up(9 * i + 2);
				double edgelength = edgeComponent(17 * ceit0->idx() + 16);
				//double edgelength = ConnectMapData1[ceit0.handle().idx()][0](4);
				Eigen::Vector3d vectorEdge = (edgelength * Edgeframe) / (Edgeframe.norm());
				int from = POINT0_HANDLE(*ceit0).idx();
				int to = POINT1_HANDLE(*ceit0).idx();
				//Eigen::Vector3d vectorEdge = SpatialTemData[to][0] - SpatialTemData[from][0];
				triple.push_back(Eigen::Triplet<double>(3 * i + 0, 3 * to + 0, 1));
				triple.push_back(Eigen::Triplet<double>(3 * i + 0, 3 * from + 0, -1));
				b(3 * i + 0) = vectorEdge(0);
				triple.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * to + 1, 1));
				triple.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * from + 1, -1));
				b(3 * i + 1) = vectorEdge(1);
				triple.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * to + 2, 1));
				triple.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * from + 2, -1));
				b(3 * i + 2) = vectorEdge(2);
			}
			for (int i = 0; i < Vertexanchors_.size(); i ++)
			{
				triple.push_back(Eigen::Triplet<double>(3 * nV + 3 * i + 0, 3 * Vertexanchors_[i] + 0, 1));
				b(3 * nV + 3 * i + 0) = SpatialTemData[0][anchorIthSeq](0);
				triple.push_back(Eigen::Triplet<double>(3 * nV + 3 * i + 1, 3 * Vertexanchors_[i] + 1, 1));
				b(3 * nV + 3 * i + 1) = SpatialTemData[0][anchorIthSeq](1);
				triple.push_back(Eigen::Triplet<double>(3 * nV + 3 * i + 2, 3 * Vertexanchors_[i] + 2, 1));
				b(3 * nV + 3 * i + 2) = SpatialTemData[0][anchorIthSeq](2);
			}
			/*for (int i = 0; i < nV; i++)
			{
				test3( 3 * i + 0 ) = SpatialTemData[i][anchorIthSeq](0);
				test3( 3 * i + 1 ) = SpatialTemData[i][anchorIthSeq](1);
				test3( 3 * i + 2 ) = SpatialTemData[i][anchorIthSeq](2);
			}*/
			/*ofstream test3out("testnV.txt");
			test3out << Coeffs * test3 - b;
			test3out.close();*/
			Coeffs.setFromTriplets(triple.begin(), triple.end()); 
			SimplicialCholesky<SparseMatrix<double>> solver1;
			solver1.compute(Coeffs.transpose() * Coeffs);
			if(solver1.info()!=Success)
			{
				/// decomposit ion failed
				std::cout<<"Decomposition failed" << std::endl;
				return;
			}
			Eigen::VectorXd y_up; 
			y_up = solver1.solve(Coeffs.transpose() * b);
			//Eigen::VectorXd test1 = Coeffs * y_up - b;
			//ofstream outy("outy_up.txt");
			//outy<< test1;
			/*output the aligned mesh to a file: seqNo.obj*/
			char *pPath = new char[30];
			sprintf_s(pPath,30, "%d-thConMapFrame.obj", ithComp);
			string offname = pPath;
			ofstream offfile(offname);
			delete []pPath;
			offfile << "COFF"<<endl;
				offfile << nV << " " << nF<< " "<< 0 <<endl;
			Eigen::VectorXd errorVec(nV);
			Eigen::Vector3d error1;
			Eigen::Vector3d error2;
			double maxerror = 0;
			double minerror = 1000;
			for(int i = 0; i < nV; i++)
			{
				 error1 = y_up.segment(3*i, 3) ;
				 error2 = SpatialTemData[i][anchorIthSeq];
				 errorVec(i) = (( y_up.segment(3*i, 3) - SpatialTemData[i][anchorIthSeq]).norm());
				if(errorVec(i) > maxerror)
					maxerror = errorVec(i);
				if(errorVec(i) < minerror)
					minerror = errorVec(i);
			}
			double meanerror = errorVec.mean();
			//double maxerror = errorVec.
			//ofstream offfile("connectOut.obj");
			for (int i = 0; i< nV; i++)
			{
				int color =(int) ((errorVec(i))  * 255/radius);
				offfile << "v" << " " << y_up(3 * i + 0) << " " <<y_up(3 * i + 1) << " " << y_up(3 * i + 2)  <<" " << 255 <<" "<< 255<< " " << 255-color<< endl;
			}
			for (PGMesh::FaceIter fit = mesh_->faces_begin(); fit != mesh_->faces_end(); ++ fit)
			{
				offfile << "f" << " ";
				for (PGMesh::FVIter fvit = mesh_->fv_begin(*fit); fvit != mesh_->fv_end(*fit); ++ fvit)
				{
					offfile << fvit->idx() + 1 << " ";
				}
				offfile << endl;
			}
			offfile << "# maxerror " << maxerror/radius << endl;
			offfile << "# meanerror "<< meanerror/radius<<endl;
			offfile.close();
		}

		/** 计算各帧每个面的法向量*/
		void SparseLocalizedForConnectMap::calculateFacesNormals()
		{
			//nSeq = 1;
			/*初始化FaceNormals*/
			FaceNormals.resize(nF);
			for(int i = 0; i < nF; i++)
			{
				FaceNormals[i].resize(nSeq);
			}
			/*遍历每个面的每一帧，计算出其面的法向量*/
			for(PGMesh::FaceIter fit = mesh_->faces_begin(); fit != mesh_->faces_end(); ++ fit)
			{
				int faceNo = fit->idx();
				PGMesh::VertexHandle v0;
				PGMesh::VertexHandle v1;
				PGMesh::VertexHandle v2;
				PGMesh::FVIter fvit = mesh_->fv_begin(*fit);
				v0 = *fvit;
				++fvit;
				v1 = *fvit;
				++fvit;
				v2 = *fvit;
				for(int frameNo = 0; frameNo < nSeq; frameNo ++)
				{
					Eigen::Vector3d e0 = POINT(v1, frameNo) - POINT(v0, frameNo);
					Eigen::Vector3d e1 = POINT(v2, frameNo) - POINT(v1, frameNo);
					Eigen::Vector3d Normal = computeCross(e0, e1);
					Normal = Normal / Normal.norm();
					FaceNormals[faceNo][frameNo] = Normal;
				}
			}
			nSeq = SpatialTemData[0].size();
		}

		Eigen::MatrixXd SparseLocalizedForConnectMap::calculateFacesNormals(Eigen::VectorXd x)
		{
			Eigen::MatrixXd FN(nF, 3);
			ofstream normal3("normal3.txt");
			/*遍历每个面的每一帧，计算出其面的法向量*/
//#pragma omp parallel for 
			for (int faceNo = 0 ; faceNo < nF; ++faceNo)
			{
				auto fit = mesh_->face_handle(faceNo);
				//int faceNo = fit->idx();
				PGMesh::VertexHandle v0;
				PGMesh::VertexHandle v1;
				PGMesh::VertexHandle v2;
				PGMesh::FVIter fvit = mesh_->fv_begin(fit);
				v0 = *fvit;
				++fvit;
				v1 = *fvit;
				++fvit;
				v2 = *fvit;

				Eigen::Vector3d e0 = POINT(v1, x) - POINT(v0, x);
				Eigen::Vector3d e1 = POINT(v2, x) - POINT(v1, x);
				Eigen::Vector3d Normal = computeCross(e0, e1);
				Normal = Normal / Normal.norm();
				FN.row(faceNo) = Normal;
				normal3 << Normal(0) << " " << Normal(1) << " " << Normal(2) << endl;
			}
			normal3.close();
			return FN;
		}

		/** 计算单帧单个面的法向量*/
		Eigen::Vector3d SparseLocalizedForConnectMap::calculateFaceNormal(int faceNo, int frameNo)
		{
			PGMesh::FaceHandle fit = mesh_->face_handle(faceNo);
			PGMesh::VertexHandle v0;
				PGMesh::VertexHandle v1;
				PGMesh::VertexHandle v2;
				PGMesh::FVIter fvit = mesh_->fv_begin(fit);
				v0 = *fvit;
				++fvit;
				v1 = *fvit;
				++fvit;
				v2 = *fvit;
				Eigen::Vector3d e0 = POINT(v1, frameNo) - POINT(v0, frameNo);
					Eigen::Vector3d e1 = POINT(v2, frameNo) - POINT(v1, frameNo);
					Eigen::Vector3d Normal = computeCross(e0, e1);
					Normal = Normal / Normal.norm();
					return Normal;
		}

		/** 计算各帧各个面的局部标架*/
		/*以列向量存储*/
		void SparseLocalizedForConnectMap::calculateFacesFrame()
		{
			//nSeq = 1;
			/*初始化FaceLocalFrames*/
			FaceLocalFrames.resize(nF);
			for(int i = 0; i < nF; i++)
			{
				FaceLocalFrames[i].resize(nSeq);
			}
			/*遍历每个面的每一帧，得到面的局部标架*/
			for(PGMesh::FaceIter fit = mesh_->faces_begin(); fit != mesh_->faces_end(); ++fit)
			{
				int faceNo = fit->idx();
				PGMesh::VertexHandle v0;
				PGMesh::VertexHandle v1;
				PGMesh::VertexHandle v2;
				PGMesh::FVIter fvit = mesh_->fv_begin(*fit);
				/*v0 = fvit.handle();
				++fvit;
				v1 = fvit.handle();
				++fvit;
				v2 = fvit.handle();*/
				PGMesh::HalfedgeHandle he0;
				PGMesh::FHIter fheit = mesh_->fh_begin(*fit);
				he0 = *fheit;
				v0 = mesh_->from_vertex_handle(he0);
				v1 = mesh_->to_vertex_handle(he0);
				for(int frameNo = 0; frameNo < nSeq; frameNo++)
				{
					Eigen::Matrix3d faceFrame;
					Eigen::Vector3d x_axis, y_axis, z_axis;
					x_axis = POINT(v1, frameNo) - POINT(v0, frameNo);
					x_axis = x_axis / x_axis.norm();
					z_axis = FaceNormals[faceNo][frameNo];
					y_axis = computeCross(z_axis, x_axis);
					y_axis = y_axis / y_axis.norm();
					faceFrame.col(0) = x_axis;
					faceFrame.col(1) = y_axis;
					faceFrame.col(2) = z_axis;
					FaceLocalFrames[faceNo][frameNo] = faceFrame;
				}
			}
			nSeq = SpatialTemData[0].size();
		}

		/** 计算单帧单个面的标架*/
		SparseLocalizedForConnectMap::faceanchor SparseLocalizedForConnectMap::calculateFaceFrame(int faceNo, int frameNo)
		{
			faceanchor fanchor;
			PGMesh::FaceHandle fit = mesh_->face_handle(faceNo);
			PGMesh::VertexHandle v0;
				PGMesh::VertexHandle v1;
				PGMesh::VertexHandle v2;
				PGMesh::FVIter fvit = mesh_->fv_begin(fit);
				/*v0 = fvit.handle();
				++fvit;
				v1 = fvit.handle();
				++fvit;
				v2 = fvit.handle();*/
				PGMesh::HalfedgeHandle he0;
				PGMesh::FHIter fheit = mesh_->fh_begin(fit);
				he0 = *fheit;
				v0 = mesh_->from_vertex_handle(he0);
				v1 = mesh_->to_vertex_handle(he0);
			//计算第一条边
				Eigen::Vector3d x_axis, y_axis, z_axis;
				x_axis = POINT(v1, frameNo) - POINT(v0, frameNo);
					x_axis = x_axis / x_axis.norm();
					z_axis = calculateFaceNormal(faceNo, frameNo);
					y_axis = computeCross(z_axis, x_axis);
					y_axis = y_axis / y_axis.norm();
					fanchor.localframe.col(0) = x_axis;
					fanchor.localframe.col(1) = y_axis;
					fanchor.localframe.col(2) = z_axis;
					fanchor.frameNo = frameNo;
					fanchor.idx = faceNo;
					return fanchor;
		}

		/** 计算IICs数据*/
		void SparseLocalizedForConnectMap::computeIICs()
		{
			/*初始化IICsData*/
			//IICsData.resize(nE);
			//LocalTransforms.resize(nE);
			//for(int i = 0; i < nE; i++)
			//{
			//	IICsData[i].resize(nSeq);
			//	//LocalTransforms[i].resize(nSeq);
			//}
			//ofstream file1("project.txt");
			/*初始化IICsMatrix*/
			//nSeq = 1;
			IICsMatrix = Eigen::MatrixXd(nSeq, 3 * nE);
			/*遍历每条边，计算IIC*/
			for(PGMesh::EdgeIter eit = mesh_->edges_begin(); eit != mesh_->edges_end(); ++eit)
			{
				int edgeNo = eit->idx();
				PGMesh::HalfedgeHandle heh0 = mesh_->halfedge_handle( *eit, 0 );
				PGMesh::HalfedgeHandle heh1 = mesh_->halfedge_handle( *eit, 1 );
				PGMesh::VertexHandle vj = mesh_->from_vertex_handle(heh0);
				PGMesh::VertexHandle vi = mesh_->to_vertex_handle(heh0);
				PGMesh::FaceHandle fh0 = mesh_->face_handle( heh0 );
				PGMesh::FaceHandle fh1 = mesh_->face_handle( heh1 );
				
				/*找出fh0中的三条半边的handle*/
				PGMesh::HalfedgeHandle he0;
				PGMesh::HalfedgeHandle he1;
				PGMesh::HalfedgeHandle he2;
				PGMesh::EdgeHandle e0;
				PGMesh::EdgeHandle e1;
				PGMesh::EdgeHandle e2;
				PGMesh::FHIter fheit = mesh_->fh_begin(fh0);
				he0 = *fheit;
				e0 = mesh_->edge_handle(he0);
				++fheit;
				he1 = *fheit;
				e1 = mesh_->edge_handle(he1);
				++fheit;
				he2 = *fheit;
				e2 = mesh_->edge_handle(he2);
				///*找出fh1中的三个点的handle*/
				//PGMesh::VertexHandle v0;
				//PGMesh::VertexHandle v1;
				//PGMesh::VertexHandle v2;
				//PGMesh::FVIter fvit = mesh_->fv_begin(fh1);
				//v0 = fvit.handle();
				//++fvit;
				//v1 = fvit.handle();
				//++fvit;
				//v2 = fvit.handle();
				
				for(int frameNo = 0; frameNo < nSeq; frameNo ++)
				{
					double length = VECTOR(*eit, frameNo).norm();
					/*若为边缘边，则标架无变化，因此两角度分别为0*/
					if(!(fh0.is_valid() && fh1.is_valid()))
					{
						//double length = VECTOR(eit.handle(), frameNo).norm();
						Eigen::Matrix3d trans;
						trans<< 1, 0, 0,
							0, 1, 0,
							0, 0, 1;
						Eigen::Vector3d Iic (0,0,length);
						//IICsData[edgeNo][frameNo] = Iic;
						IICsMatrix.block<1, 3>(frameNo, 3 * edgeNo) = Iic.transpose();
						//LocalTransforms[edgeNo][frameNo] = trans;
					}
					else
					{
						Eigen::Matrix3d frame1 = FaceLocalFrames[fh0.idx()][frameNo];
						Eigen::Matrix3d frame2 = FaceLocalFrames[fh1.idx()][frameNo];
						Eigen::Matrix3d frame3;
						/*先求出二面角，即phi*/
						double phi;
						Eigen::Vector3d normal1 = frame1.col(2);
						Eigen::Vector3d normal2 = frame2.col(2);
						Eigen::Matrix3d trans = frame1.transpose() * frame2;
						//LocalTransforms[edgeNo][frameNo] = trans;
						Eigen::Vector3d e_ij = POINT(vj, frameNo) - POINT(vi, frameNo);
						e_ij = e_ij / (e_ij.norm());
						Eigen::Vector3d crossnormal = computeCross(normal1, normal2);
						double dotcross = e_ij.dot(crossnormal);
						if(dotcross >= 0)
						{
							double dotNorm1Norm2 = normal1.dot(normal2);
							if(dotNorm1Norm2 < 1 - 1e-6 && dotNorm1Norm2 > -1 + 1e-6)
							{
								phi = acos(dotNorm1Norm2);
							}
							else if(dotNorm1Norm2 >= 1 - 1e-6 )
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
							if(dotNorm1Norm2 < 1 - 1e-6 && dotNorm1Norm2 > -1 + 1e-6)
							{
								phi = -acos(dotNorm1Norm2);
							}   
							else if(dotNorm1Norm2 >= 1 - 1e-6 )
							{
								phi = 0;
							}
							else
							{
								phi = -3.1415926;
							}
						}
						/*将第二个三角形的第一条边旋转到与第一个三角形共面*/
						Eigen::Vector3d edgevector = frame2.col(0);
						OpenMesh::Vec3d edgevector0 = eigenTransToVec3d(edgevector);
						
						Eigen::Vector3d rotateaxis = -e_ij;
						OpenMesh::Vec3d rotateaxis0 = eigenTransToVec3d(rotateaxis);
						rotateaxis0 = rotateaxis0 * phi;
						OpenMesh::Vec3d edgevectorproj0;
						edgevectorproj0 = Math::rotate<OpenMesh::Vec3d>(rotateaxis0, edgevector0);
						Eigen::Vector3d edgevectorproj = Vec3dToEigen(edgevectorproj0);
						Eigen::Vector3d projectedge = frame1.transpose() * edgevectorproj; 
						//file1<<projectedge(0) <<" "<< projectedge(1)<< " "<< projectedge(2)<<endl;
						Eigen::Vector3d projectRotationaxistolocal = frame1.transpose() * e_ij;
						
						int judge;
						if(heh0 == he0)
							judge = 1;
						else if(heh0 == he1)
							judge = 2;
						else
							judge = 3;
						/*求theta*/
						//file1<<projectRotationaxistolocal(0) <<" "<< projectRotationaxistolocal(1)<< " "<< projectRotationaxistolocal(2)<<" "<< judge <<endl;
						double theta;
						double x = projectedge(0);
						double y = projectedge(1);
						//file1<<projectedge(0) <<" "<< projectedge(1)<< " "<< projectedge(2)<<" "<< judge <<endl;
						if (y >= 0)
						{
							if(x >= -1 + 1e-6)
							{
								theta = acos(x);
							}
							else
							{
								theta = 3.1415926;
							}
						}
						else if ( x < 0)
						{
							if(-x < 1 - 1e-6)
							{
								theta = 3.1415926 + acos(-x);
							}
							else
								theta = 3.1415926;
						}
						else
						{
							if (x < 1 - 1e-6)
							{
								theta = 2 * 3.1415926 - acos(x);
							}
							else
								theta = 2 * 3.1415926;
						}
						Eigen::Vector3d Iic(theta, phi, length);
						//file1<<Iic(0) <<" "<< Iic(1)<< " "<< Iic(2)<<endl;
						//IICsData[edgeNo][frameNo] = Iic;
						IICsMatrix.block<1, 3>(frameNo, 3 * edgeNo) = Iic.transpose();
						/*Eigen::Matrix3d frame2Inframe1;
						Eigen::Matrix3d matrixNormal;
						matrixNormal << 0, -1, 0,
							1, 0, 0,
							0, 0, 0;
						Eigen::Matrix3d expNormal = (theta * matrixNormal).exp();
						Eigen::Matrix3d matrixedge;
						matrixedge << 0, -projectRotationaxistolocal(2), projectRotationaxistolocal(1),
							projectRotationaxistolocal(2), 0, -projectRotationaxistolocal(0),
							-projectRotationaxistolocal(1), projectRotationaxistolocal(0), 0;
						Eigen::Matrix3d expedge = (phi * matrixedge).exp();
						frame2Inframe1 = expedge * expNormal;
						frame3 = frame1 * frame2Inframe1;*/
						//file1<<frame2Inframe1 <<endl<<LocalTransforms[edgeNo][frameNo]<<endl<<endl;
					}
				}
				
			}
			nSeq = SpatialTemData[0].size();
			//file1 .close();
		}
		 
		/** 计算某网格的二面角和边长数据*/
		Eigen::VectorXd SparseLocalizedForConnectMap::computeDiEdge(Eigen::VectorXd x)
		{
			Eigen::VectorXd DEdata = Eigen::VectorXd(2 * nE);
			std::vector<double> DEdata1(2 * nE);
			Eigen::MatrixXd FN = calculateFacesNormals(x);
			/*updateMeshVertices(x, mesh1);
			OpenMesh::IO::write_mesh(*mesh1, "update.obj");
			mesh1->update_face_normals();*/
			/*遍历每条边，计算IIC*/

			for (int edgeNo = 0; edgeNo < mesh_->n_edges(); edgeNo++)
			{
				PGMesh::EdgeHandle eit = mesh_->edge_handle(edgeNo);
				//int edgeNo = eit->idx();
				PGMesh::HalfedgeHandle heh0 = mesh_->halfedge_handle(eit, 0);
				PGMesh::HalfedgeHandle heh1 = mesh_->halfedge_handle(eit, 1);
				PGMesh::VertexHandle vj = mesh_->from_vertex_handle(heh0);
				PGMesh::VertexHandle vi = mesh_->to_vertex_handle(heh0);
				PGMesh::FaceHandle fh0 = mesh_->face_handle(heh0);
				PGMesh::FaceHandle fh1 = mesh_->face_handle(heh1);

				/*找出fh0中的三条半边的handle*/
				PGMesh::HalfedgeHandle he0;
				PGMesh::HalfedgeHandle he1;
				PGMesh::HalfedgeHandle he2;
				PGMesh::EdgeHandle e0;
				PGMesh::EdgeHandle e1;
				PGMesh::EdgeHandle e2;
				PGMesh::FHIter fheit = mesh_->fh_begin(fh0);
				he0 = *fheit;
				e0 = mesh_->edge_handle(he0);
				++fheit;
				he1 = *fheit;
				e1 = mesh_->edge_handle(he1);
				++fheit;
				he2 = *fheit;
				e2 = mesh_->edge_handle(he2);



				double length = sqrt((x(3 * vi.idx() + 0) - x(3 * vj.idx() + 0))* (x(3 * vi.idx() + 0) - x(3 * vj.idx() + 0)) + (x(3 * vi.idx() + 1) - x(3 * vj.idx() + 1))*  (x(3 * vi.idx() + 1) - x(3 * vj.idx() + 1)) + (x(3 * vi.idx() + 2) - x(3 * vj.idx() + 2))* (x(3 * vi.idx() + 2) - x(3 * vj.idx() + 2)));
				//cout << length << endl;
				/*若为边缘边，则标架无变化，因此两角度分别为0*/
				{
					if (!(fh0.is_valid() && fh1.is_valid()))
					{
						Eigen::Vector2d Iic(0, length);
						//IICsData[edgeNo][frameNo] = Iic;
						DEdata1[2 * edgeNo + 0] = 0;
						DEdata1[2 * edgeNo + 1] = length;
						//DEdata.segment( 2 * edgeNo, 2) = Iic.transpose();
						//LocalTransforms[edgeNo][frameNo] = trans;
					}
					else
					{

						Eigen::Matrix3d frame3;
						/*先求出二面角，即phi*/
						double phi;
						Eigen::VectorXd normal1 = FN.row(fh0.idx());
						Eigen::VectorXd normal2 = FN.row(fh1.idx());
						Eigen::Vector3d e_ij = POINT(vj, x) - POINT(vi, x);
						e_ij = e_ij / (e_ij.norm());
						Eigen::Vector3d crossnormal = computeCross(normal1, normal2);
						double dotcross = e_ij.dot(crossnormal);
						//e_ij(0) * crossnormal[0] + e_ij(1) *crossnormal[1] +e_ij(2) * crossnormal[2];
						double dotNorm1Norm2 = normal1(0) * normal2(0) + normal1(1) * normal2(1) + normal1(2) * normal2(2);
						//cout << normal1 << endl;
						if (dotcross >= 0)
						{

							if (dotNorm1Norm2 < 1 - 1e-6 && dotNorm1Norm2 > -1 + 1e-6)
							{
								phi = acos(dotNorm1Norm2);
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
							if (dotNorm1Norm2 < 1 - 1e-6 && dotNorm1Norm2 > -1 + 1e-6)
							{
								phi = -acos(dotNorm1Norm2);
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
						//DEdata1[2 * edgeNo + 0] = phi;
						//DEdata1[2 * edgeNo + 1] = length;
						Eigen::Vector2d Iic(phi, length);
						DEdata.segment(2 * edgeNo, 2) = Iic;
					}
				}
			}
			return DEdata;
		}

		
		Eigen::VectorXd SparseLocalizedForConnectMap::getVertexMatrix(PGMesh * _mesh)
		{
			Eigen::VectorXd V(3 * nV);
			for (auto vit = _mesh->vertices_begin(); vit != _mesh->vertices_end(); ++vit)
			{
				PGMesh::Point p = _mesh->point(*vit);
				V(3 * vit->idx()+ 0) = p[0];
				V(3 * vit->idx()+ 1) = p[1];
				V(3 * vit->idx()+ 2) = p[2];
			}
			return V;
		}


		/** 计算二面角和边长数据*/
		void SparseLocalizedForConnectMap::computeDiEdge()
		{
			
			/*初始化IICsMatrix*/
			//nSeq = 1;
			DiEdgeDataMatrix = Eigen::MatrixXd(nSeq, 2 * nE);
			/*遍历每条边，计算IIC*/
			for(PGMesh::EdgeIter eit = mesh_->edges_begin(); eit != mesh_->edges_end(); ++eit)
			{
				int edgeNo = eit->idx();
				PGMesh::HalfedgeHandle heh0 = mesh_->halfedge_handle( *eit, 0 );
				PGMesh::HalfedgeHandle heh1 = mesh_->halfedge_handle( *eit, 1 );
				PGMesh::VertexHandle vj = mesh_->from_vertex_handle(heh0);
				PGMesh::VertexHandle vi = mesh_->to_vertex_handle(heh0);
				PGMesh::FaceHandle fh0 = mesh_->face_handle( heh0 );
				PGMesh::FaceHandle fh1 = mesh_->face_handle( heh1 );
				
				/*找出fh0中的三条半边的handle*/
				PGMesh::HalfedgeHandle he0;
				PGMesh::HalfedgeHandle he1;
				PGMesh::HalfedgeHandle he2;
				PGMesh::EdgeHandle e0;
				PGMesh::EdgeHandle e1;
				PGMesh::EdgeHandle e2;
				PGMesh::FHIter fheit = mesh_->fh_begin(fh0);
				he0 = *fheit;
				e0 = mesh_->edge_handle(he0);
				++fheit;
				he1 = *fheit;
				e1 = mesh_->edge_handle(he1);
				++fheit;
				he2 = *fheit;
				e2 = mesh_->edge_handle(he2);
				
				
				for(int frameNo = 0; frameNo < nSeq; frameNo ++)
				{
					double length = VECTOR(*eit, frameNo).norm();
					/*若为边缘边，则标架无变化，因此两角度分别为0*/
					if(!(fh0.is_valid() && fh1.is_valid()))
					{
						//double length = VECTOR(eit.handle(), frameNo).norm();
						Eigen::Matrix3d trans;
						trans<< 1, 0, 0,
							0, 1, 0,
							0, 0, 1;
						Eigen::Vector2d Iic (0,length);
						//IICsData[edgeNo][frameNo] = Iic;
						DiEdgeDataMatrix.block<1, 2>(frameNo, 2 * edgeNo) = Iic.transpose();
						//LocalTransforms[edgeNo][frameNo] = trans;
					}
					else
					{
						Eigen::Matrix3d frame1 = FaceLocalFrames[fh0.idx()][frameNo];
						Eigen::Matrix3d frame2 = FaceLocalFrames[fh1.idx()][frameNo];
						Eigen::Matrix3d frame3;
						/*先求出二面角，即phi*/
						double phi;
						Eigen::Vector3d normal1 = frame1.col(2);
						Eigen::Vector3d normal2 = frame2.col(2);
						Eigen::Matrix3d trans = frame1.transpose() * frame2;
						//LocalTransforms[edgeNo][frameNo] = trans;
						Eigen::Vector3d e_ij = POINT(vj, frameNo) - POINT(vi, frameNo);
						e_ij = e_ij / (e_ij.norm());
						Eigen::Vector3d crossnormal = computeCross(normal1, normal2);
						double dotcross = e_ij.dot(crossnormal);
						if(dotcross >= 0)
						{
							double dotNorm1Norm2 = normal1.dot(normal2);
							if(dotNorm1Norm2 < 1 - 1e-6 && dotNorm1Norm2 > -1 + 1e-6)
							{
								phi = acos(dotNorm1Norm2);
							}
							else if(dotNorm1Norm2 >= 1 - 1e-6 )
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
							if(dotNorm1Norm2 < 1 - 1e-6 && dotNorm1Norm2 > -1 + 1e-6)
							{
								phi = -acos(dotNorm1Norm2);
							}   
							else if(dotNorm1Norm2 >= 1 - 1e-6 )
							{
								phi = 0;
							}
							else
							{
								phi = -3.1415926;
							}
						}
						
						Eigen::Vector2d Iic( phi, length);
						DiEdgeDataMatrix.block<1, 2>(frameNo, 2 * edgeNo) = Iic.transpose();
						
					}
				}
				
			}
			nSeq = SpatialTemData[0].size();
			//file1 .close();
		}

		/** 根据IIC重构网格*/
		/*采用论文Computer Graphics Forum 2012*/
		/*"Linear surface reconstruction from discrete fundamental forms on triangle meshes."*/
		/**/
		void SparseLocalizedForConnectMap::reconstructionFromIICs(Eigen::VectorXd edgeComponent, Eigen::VectorXd & cordComponent, int anchorIthSeq, double radius, int ithComp)
		{
			/** 通过角度利用线性回归估计边长*/
			/*std::vector<double> edgelengths(nE);
			for(int edgeNo = 0; edgeNo < nE; edgeNo++)
			{
				Eigen::VectorXd x(6);
				x(0) = 1;
				for(int i = 0; i < 4; i++)
				{
					int index = edge_regressions[edgeNo].one_ring_idx[i];
					x(i+1) = edgeComponent(3 * index + 1) - IICsMatrix(0, 3 * index + 1);
				}
				x(5) = edgeComponent(3 * edgeNo + 1) - IICsMatrix(0, 3 * edgeNo + 1);
				Eigen::VectorXd coef = edge_regressions[edgeNo].one_ring_Coef;
				edgelengths[edgeNo] = coef.transpose() * x; 
				edgelengths[edgeNo]+= IICsMatrix(0, 3 * edgeNo + 2);
			}*/
			//ofstream file1("cossin.txt");
			if( edgeComponent.size() != 3  * nE)
				return;
			Eigen::SparseMatrix<double> Coefficients(9 * nE + Faceanchors_.size() * 9 ,9 * nF);
			Coefficients.setZero();
			Eigen::VectorXd frameCor(9 * nE + Faceanchors_.size() * 9);
			frameCor.setZero();
			std::vector<Eigen::Triplet<double> > triple;  
			triple.clear();
			//ofstream file2("inversetrans.txt");
			//ofstream file3("reconstrans.txt");
			//ofstream file4("axis.txt");
			for(PGMesh::EdgeIter eit = mesh_->edges_begin(); eit != mesh_->edges_end(); ++ eit)
			{ 
				int edgeNo = eit->idx();
				PGMesh::HalfedgeHandle heh0 = mesh_->halfedge_handle( *eit, 0 );
				PGMesh::HalfedgeHandle heh1 = mesh_->halfedge_handle( *eit, 1 );
				PGMesh::VertexHandle vj = mesh_->from_vertex_handle(heh0);
				PGMesh::VertexHandle vi = mesh_->to_vertex_handle(heh0);
				PGMesh::FaceHandle fh0 = mesh_->face_handle( heh0 );
				PGMesh::FaceHandle fh1 = mesh_->face_handle( heh1 );
				int fromFaceNo = fh0.idx();
				int toFaceNo = fh1.idx();
				double phi = edgeComponent(3 * edgeNo + 0);
				double theta = edgeComponent(3 * edgeNo + 1);
				Eigen::Matrix3d matrixNormal;
				matrixNormal << 0, -1, 0,
					1, 0, 0,
					0, 0, 0;
				Eigen::Matrix3d expNormal;
				int loop = phi/(2 * 3.1415926);
				phi = phi - loop * 2 * 3.1415926;
				if(phi < 3.141592653)
				{
					expNormal = (phi * matrixNormal).exp();
				}
				else
				{
					expNormal = (-(2 * 3.141592653 - phi) * matrixNormal).exp();
				}
				//file4 << expNormal<<endl<<phi<<endl;
				/*找出fh0中的三条半边的handle*/
				PGMesh::HalfedgeHandle he0;
				PGMesh::HalfedgeHandle he1;
				PGMesh::HalfedgeHandle he2;
				PGMesh::EdgeHandle e0;
				PGMesh::EdgeHandle e1;
				PGMesh::EdgeHandle e2;
				PGMesh::FHIter fheit = mesh_->fh_begin(fh0);
				he0 = *fheit;
				e0 = mesh_->edge_handle(he0);
				++fheit;
				he1 = *fheit;
				e1 = mesh_->edge_handle(he1);
				++fheit;
				he2 = *fheit;
				e2 = mesh_->edge_handle(he2);
				Eigen::Vector3d axis;
				int judge;
				if(heh0 == he0)
				{
					axis << -1, 0, 0;
					judge = 1;
				}
				else if(heh0 == he1)
				{

					double a = edgeComponent(3 * e0.idx() + 2);
					double b = edgeComponent(3 * e1.idx() + 2);
					double c = edgeComponent(3 * e2.idx() + 2);
					/*double a = IICsMatrix.row(0) (3 * e0.idx() + 2);
					double b = IICsMatrix.row(0)(3 * e1.idx() + 2);
					double c = IICsMatrix.row(0)(3 * e2.idx() + 2);*/
					/*double a = edgelengths [e0.idx()];
					double b = edgelengths [e1.idx()];
					double c = edgelengths [e2.idx()];*/

					double sqcos = (a * a + b * b - c * c)/ (2 * a * b);
					double squaresin = 1 - sqcos * sqcos;
					axis << sqcos, -sqrt(squaresin), 0;
					if(squaresin < 0)
						axis(1) = 0;
					if(sqcos > 1)
						axis(0) = 1;
					if(sqcos < -1)
						axis(0) = -1;
					//file1 << axis(0)<<" "<< axis(1)<<" "<< axis(2)<<endl;
					judge = 2;
				}
				else
				{
					double a = edgeComponent(3 * e0.idx() + 2);
					double b = edgeComponent(3 * e2.idx() + 2);
					double c = edgeComponent(3 * e1.idx() + 2);
					/*double a = IICsMatrix.row(0) (3 * e0.idx() + 2);
					double b = IICsMatrix.row(0)(3 * e1.idx() + 2);
					double c = IICsMatrix.row(0)(3 * e2.idx() + 2);*/
					/*double a = edgelengths [e0.idx()];
					double b = edgelengths [e1.idx()];
					double c = edgelengths [e2.idx()];*/

					double sqcos = (a * a + b * b - c * c)/ (2 * a * b);
					double squaresin = 1 - sqcos * sqcos;

					axis << sqcos, sqrt(squaresin), 0;
					if(squaresin < 0)
						axis(1) = 0;
					if(sqcos > 1)
						axis(0) = 1;
					if(sqcos < -1)
						axis(0) = -1;
					//file1 << axis(0)<<" "<< axis(1)<<" "<< axis(2)<<endl;
					judge = 3;
				}
				//file4<< axis(0)<<" "<<axis(1)<<" "<<axis(2)<< " "<<judge<<endl;
				Eigen::Matrix3d matrixedge;
				matrixedge << 0, -axis(2), axis(1),
					axis(2), 0, -axis(0),
					-axis(1), axis(0), 0;
				Eigen::Matrix3d expedge = (theta * matrixedge).exp();
				Eigen::Matrix3d rot = expedge * expNormal;
				//file1<< rot << endl << endl;
				//file3<< LocalTransforms[edgeNo][ithComp]<<endl<<endl;
				if((fh0.is_valid() && fh1.is_valid()))
				{
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 0, 9 * fromFaceNo +0, rot(0,0)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 0, 9 * fromFaceNo +1, rot(1,0)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 0, 9 * fromFaceNo +2,rot(2,0)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 0, 9 * toFaceNo +0, -1));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 1, 9 * fromFaceNo +3,rot(0,0)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 1, 9 * fromFaceNo +4,rot(1,0)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 1, 9 * fromFaceNo +5,rot(2,0)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 1, 9 * toFaceNo +3, -1));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 2, 9 * fromFaceNo +6,rot(0,0)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 2, 9 * fromFaceNo +7,rot(1,0)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 2, 9 * fromFaceNo +8,rot(2,0)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 2, 9 * toFaceNo +6,-1));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 3, 9 * fromFaceNo +0,rot(0,1)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 3, 9 * fromFaceNo +1,rot(1,1)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 3, 9 * fromFaceNo +2,rot(2,1)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 3, 9 * toFaceNo +1, -1));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 4, 9 * fromFaceNo +3, rot(0,1)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 4, 9 * fromFaceNo +4, rot(1,1)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 4, 9 * fromFaceNo +5, rot(2,1)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 4, 9 * toFaceNo +4, -1));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 5, 9 * fromFaceNo +6, rot(0,1)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 5, 9 * fromFaceNo +7, rot(1,1)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 5, 9 * fromFaceNo +8, rot(2,1)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 5, 9 * toFaceNo +7, -1));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 6, 9 * fromFaceNo +0, rot(0,2)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 6, 9 * fromFaceNo +1, rot(1,2)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 6, 9 * fromFaceNo +2, rot(2,2)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 6, 9 * toFaceNo +2, -1));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 7, 9 * fromFaceNo +3, rot(0,2)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 7, 9 * fromFaceNo +4, rot(1,2)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 7, 9 * fromFaceNo +5, rot(2,2)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 7, 9 * toFaceNo +5, -1));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 8, 9 * fromFaceNo +6, rot(0,2)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 8, 9 * fromFaceNo +7, rot(1,2)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 8, 9 * fromFaceNo +8, rot(2,2)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 8, 9 * toFaceNo +8, -1));
				}
			}
			//设置锚点
			for(int l = 0; l < Faceanchors_.size(); l++)
			{
				Eigen::Matrix3d anchorFrame = faceAnchorFrames[anchorIthSeq][l].localframe; 
				triple.push_back(Eigen::Triplet<double>(9 * nE + l * 9 + 0, 9 * Faceanchors_[l] + 0, 1));
				triple.push_back(Eigen::Triplet<double>(9 * nE + l * 9 + 1, 9 * Faceanchors_[l] + 1, 1));
				triple.push_back(Eigen::Triplet<double>(9 * nE + l * 9 + 2, 9 * Faceanchors_[l] + 2, 1));
				triple.push_back(Eigen::Triplet<double>(9 * nE + l * 9 + 3, 9 * Faceanchors_[l] + 3, 1));
				triple.push_back(Eigen::Triplet<double>(9 * nE + l * 9 + 4, 9 * Faceanchors_[l] + 4, 1));
				triple.push_back(Eigen::Triplet<double>(9 * nE + l * 9 + 5, 9 * Faceanchors_[l] + 5, 1));
				triple.push_back(Eigen::Triplet<double>(9 * nE + l * 9 + 6, 9 * Faceanchors_[l] + 6, 1));
				triple.push_back(Eigen::Triplet<double>(9 * nE + l * 9 + 7, 9 * Faceanchors_[l] + 7, 1));
				triple.push_back(Eigen::Triplet<double>(9 * nE + l * 9 + 8, 9 * Faceanchors_[l] + 8, 1));
				frameCor(9 * nE + l * 9 + 0) = anchorFrame(0,0);
				frameCor(9 * nE + l * 9 + 1) = anchorFrame(0,1);
				frameCor(9 * nE + l * 9 + 2) = anchorFrame(0,2);
				frameCor(9 * nE + l * 9 + 3) = anchorFrame(1,0);
				frameCor(9 * nE + l * 9 + 4) = anchorFrame(1,1);
				frameCor(9 * nE + l * 9 + 5) = anchorFrame(1,2);
				frameCor(9 * nE + l * 9 + 6) = anchorFrame(2,0);
				frameCor(9 * nE + l * 9 + 7) = anchorFrame(2,1);
				frameCor(9 * nE + l * 9 + 8) = anchorFrame(2,2);
			}
			Coefficients.setFromTriplets(triple.begin(), triple.end());  
			triple.clear();
			SimplicialCholesky<SparseMatrix<double>> solver;
			//ofstream solvertxt("solve1.txt");
			//solvertxt<<Coefficients.transpose() * Coefficients;
			solver.compute(Coefficients.transpose() * Coefficients);
			if(solver.info()!=Success)
			{
				/// decomposit ion failed
				std::cout<<"Decomposition failed" << std::endl;
				return;
			}
			Eigen::VectorXd x_up; 
			x_up = solver.solve(Coefficients.transpose() * frameCor);
			/*solvertxt<< x_up;
			solvertxt.close();*/
			//solvertxt.close();
			/*Given the frames of all faces, reconstruct the vertex coordinates*/
			//Eigen::SparseMatrix<double> vertexCoff(9 * nF + 3 , 3 * nV);
			Eigen::VectorXd fvertex(9 * nF + 3);
			Eigen::VectorXd y_up(3 * nV);
			//ofstream axisview("axis.txt");
			//ofstream edgeview("edgevector.txt");
			for(PGMesh::FaceIter fit = mesh_->faces_begin(); fit != mesh_->faces_end(); ++fit)
			{
				PGMesh::FaceHandle fh = *fit;
				int faceNo = fh.idx();
				/*找出fh0中的三条半边的handle*/
				PGMesh::HalfedgeHandle he0;
				PGMesh::HalfedgeHandle he1;
				PGMesh::HalfedgeHandle he2;
				PGMesh::EdgeHandle e0;
				PGMesh::EdgeHandle e1;
				PGMesh::EdgeHandle e2;
				PGMesh::FHIter fheit = mesh_->fh_begin(fh);
				he0 = *fheit;
				e0 = mesh_->edge_handle(he0);
				++fheit;
				he1 = *fheit;
				e1 = mesh_->edge_handle(he1);
				++fheit;
				he2 = *fheit;
				e2 = mesh_->edge_handle(he2);
				/*找出三条边的长度*/
				double a, b, c;
				a = edgeComponent(3 * e0.idx() + 2);
				b = edgeComponent(3 * e1.idx() + 2);
				c = edgeComponent(3 * e2.idx() + 2);
				 /*a = IICsMatrix.row(0) (3 * e0.idx() + 2);
				 b = IICsMatrix.row(0)(3 * e1.idx() + 2);
				 c = IICsMatrix.row(0)(3 * e2.idx() + 2);*/
				 /*a = edgelengths [e0.idx()];
				 b = edgelengths [e1.idx()];
				 c = edgelengths [e2.idx()];*/
				/*找出fh1中的三个点的handle*/
				PGMesh::VertexHandle v0;
				PGMesh::VertexHandle v1;
				PGMesh::VertexHandle v2;
				v0 = mesh_->from_vertex_handle(he0);
				v1 = mesh_->to_vertex_handle(he0);
				v2 = mesh_->to_vertex_handle(he1);
				
				/*提取出该面的标架中的a，b轴*/
				Eigen::Matrix3d frameFace;
				frameFace << x_up(9 * faceNo + 0), x_up(9 * faceNo + 1), x_up(9 * faceNo + 2),
					x_up(9 * faceNo + 3), x_up(9 * faceNo + 4), x_up(9 * faceNo + 5),
					x_up(9 * faceNo + 6), x_up(9 * faceNo + 7), x_up(9 * faceNo + 8);
				Eigen::JacobiSVD<MatrixXd> svd( frameFace, ComputeThinU | ComputeThinV);
				Eigen::MatrixXd u = svd.matrixU();
				Eigen::MatrixXd v = svd.matrixV();
				Eigen::Matrix3d U = u * v.transpose();
				Eigen::Vector3d a_axis, b_axis;
				a_axis = U.col(0);
				b_axis = U.col(1);
				a_axis = a_axis/ a_axis.norm();
				b_axis = b_axis /b_axis.norm();
				//axisview<< a_axis(0)<<" "<<a_axis(1)<<" "<<a_axis(2)<<endl;
				//axisview<< b_axis(0)<<" "<<b_axis(1)<<" "<<b_axis(2)<<endl;
				/*提取出角度 alpha, beta*/
				double cosalpha, sinalpha, cosbeta, sinbeta;
				cosalpha = (a * a + c * c - b * b)/ (2 * a * c);
				sinalpha = sqrt(1 - cosalpha * cosalpha);
				cosbeta = (a * a + b * b - c * c)/ (2 * a * b);
				sinbeta = sqrt(1 - cosbeta * cosbeta);
				
				if(cosalpha > 1)
				{
					cosalpha = 1;
					sinalpha = 0;
				}
				if(cosalpha < -1)
				{
					cosalpha = -1;
					sinalpha = 0;
				}
					
				if(cosbeta > 1)
					{cosbeta = 1;
				sinbeta = 0;}
				if(cosbeta < -1)
					{cosbeta = -1;
				sinbeta = 0;}
				//axisview<< cosalpha<<" "<<sinalpha<< " " << cosbeta<<" "<< sinbeta<<endl;
				/*得到各边向量*/
				Eigen::Vector3d v0v1, v1v2, v2v0;
				v0v1 =  a * a_axis;
				v1v2 = -b * (cosbeta * a_axis - sinbeta * b_axis);
				v2v0 = -c * (cosalpha * a_axis + sinalpha * b_axis);
				//edgeview<< v0v1(0)<<" "<<v0v1(1)<<" "<<v0v1(2)<<endl;
				//edgeview<< v1v2(0)<<" "<<v1v2(1)<<" "<<v1v2(2)<<endl;
				//edgeview<< v2v0(0)<<" "<<v2v0(1)<<" "<<v2v0(2)<<endl;
				/*填入矩阵和向量元素*/
				//triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 0, 3 * v1.idx() + 0, 1));
				//triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 0, 3 * v0.idx() + 0, -1));
				fvertex(faceNo * 9 + 0) = v0v1(0);
				//triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 1, 3 * v1.idx() + 1, 1));
				//triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 1, 3 * v0.idx() + 1, -1));
				fvertex(faceNo * 9 + 1) = v0v1(1);
				//triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 2, 3 * v1.idx() + 2, 1));
				//triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 2, 3 * v0.idx() + 2, -1));
				fvertex(faceNo * 9 + 2) = v0v1(2);
				//triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 3, 3 * v2.idx() + 0, 1));
				//triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 3, 3 * v1.idx() + 0, -1));
				fvertex(faceNo * 9 + 3) = v1v2(0);
				//triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 4, 3 * v2.idx() + 1, 1));
				//triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 4, 3 * v1.idx() + 1, -1));
				fvertex(faceNo * 9 + 4) = v1v2(1);
				//triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 5, 3 * v2.idx() + 2, 1));
				//triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 5, 3 * v1.idx() + 2, -1));
				fvertex(faceNo * 9 + 5) = v1v2(2);
				//triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 6, 3 * v0.idx() + 0, 1));
				//triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 6, 3 * v2.idx() + 0, -1));
				fvertex(faceNo * 9 + 6) = v2v0(0);
				//triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 7, 3 * v0.idx() + 1, 1));
				//triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 7, 3 * v2.idx() + 1, -1));
				fvertex(faceNo * 9 + 7) = v2v0(1);
				//triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 8, 3 * v0.idx() + 2, 1));
				//triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 8, 3 * v2.idx() + 2, -1));
				fvertex(faceNo * 9 + 8) = v2v0(2);
			}
			//axisview.close();
			//edgeview.close();
			/*添加锚点*/
			
				//triple.push_back(Eigen::Triplet<double>(9 * nF + 0, 3 * 0 + 0, 1));
			fvertex(9 * nF + 0) = vertexAnchorCords[anchorIthSeq][0].cord(0);
					
				//triple.push_back(Eigen::Triplet<double>(9 * nF  + 1, 3 * 0+ 1, 1));
				fvertex(9 * nF  + 1) = vertexAnchorCords[anchorIthSeq][0].cord(1);
				//triple.push_back(Eigen::Triplet<double>(9 * nF + 2, 3 * 0 + 2, 1));
				fvertex(9 * nF + 2) = vertexAnchorCords[anchorIthSeq][0].cord(2);
			
			//vertexCoff.setFromTriplets(triple.begin(), triple.end());  
			//triple.clear();
			//ofstream solvertxt("solve1.txt");
			//solvertxt<<Coefficients.transpose() * Coefficients;
			//SimplicialCholesky<SparseMatrix<double>> solver1;
			//solver1.compute(vertexCoff.transpose() * vertexCoff);
			/*if(solver1.info()!=Success)
			{
				/// decomposit ion failed
				std::cout<<"Decomposition failed" << std::endl;
				return;
			}*/
			y_up = presolver.solve(vertexCoff.transpose() * fvertex);
			cordComponent = y_up;
			
		}

		/** 根据二面角和边长重构网格*/
		void SparseLocalizedForConnectMap::reconstructionFromDiEdges(Eigen::VectorXd edgeComponent, Eigen::VectorXd & cordComponent, int anchorIthSeq, double radius, int ithComp)
		{
			
			if( edgeComponent.size() != 2  * nE)
				return;
			Eigen::SparseMatrix<double> Coefficients(9 * nE + Faceanchors_.size() * 9 ,9 * nF);
			Coefficients.setZero();
			Eigen::VectorXd frameCor(9 * nE + Faceanchors_.size() * 9);
			frameCor.setZero();
			std::vector<Eigen::Triplet<double> > triple;  
			triple.clear();


			//int *h_cooRowIndA = NULL;
			//int *h_coo_csrColIndA = NULL;
			//double *h_coo_csrValA = NULL;
			//h_coo_csrValA = (double *)malloc(sizeof(double)* 1171620);
			//h_cooRowIndA = (int *)malloc(sizeof(int)* 1171620);
			//h_coo_csrColIndA = (int *)malloc(sizeof(int)* 1171620);
			//int * tmp_row = h_cooRowIndA;
			//int *temp_col = h_coo_csrColIndA;
			//double * tmp_val = h_coo_csrValA;
			//ofstream file2("inversetrans.txt");
			//ofstream file3("reconstrans.txt");
			//ofstream file4("axis.txt");
			//ofstream efile("efile.txt");
			//ofstream efile1("efile1.txt");
			for(PGMesh::EdgeIter eit = mesh_->edges_begin(); eit != mesh_->edges_end(); ++ eit)
			{  
				if(mesh_->is_boundary(*eit))
					continue;
				int edgeNo = eit->idx();
				PGMesh::HalfedgeHandle heh0 = mesh_->halfedge_handle( *eit, 0 );
				PGMesh::HalfedgeHandle heh1 = mesh_->halfedge_handle( *eit, 1 );
				PGMesh::VertexHandle vj = mesh_->from_vertex_handle(heh0);
				PGMesh::VertexHandle vi = mesh_->to_vertex_handle(heh0);
				PGMesh::FaceHandle fh0 = mesh_->face_handle( heh0 );
				PGMesh::FaceHandle fh1 = mesh_->face_handle( heh1 );
				int fromFaceNo = fh0.idx();
				int toFaceNo = fh1.idx();
				//double phi = edgeComponent(3 * edgeNo + 0);
				double theta = edgeComponent(2 * edgeNo + 0); 
				
				/*int loop = phi/(2 * 3.1415926);
				phi = phi - loop * 2 * 3.1415926;
				if(phi < 3.141592653)
				{
					expNormal = (phi * matrixNormal).exp();
				}
				else
				{
					expNormal = (-(2 * 3.141592653 - phi) * matrixNormal).exp();
				}*/
				//file4 << expNormal<<endl<<phi<<endl;
				/*找出fh0中的三条半边的handle*/
				PGMesh::HalfedgeHandle he0;
				PGMesh::HalfedgeHandle he1;
				PGMesh::HalfedgeHandle he2;
				PGMesh::EdgeHandle e0;
				PGMesh::EdgeHandle e1;
				PGMesh::EdgeHandle e2;
				PGMesh::FHIter fheit = mesh_->fh_begin(fh0);
				he0 = *fheit;
				e0 = mesh_->edge_handle(he0);
				++fheit;
				he1 = *fheit;
				e1 = mesh_->edge_handle(he1);
				++fheit;
				he2 = *fheit;
				e2 = mesh_->edge_handle(he2);
				Eigen::Vector3d axis;
				int judge;
				if(heh0 == he0)
				{
					axis << -1, 0, 0;
					judge = 1;
					//efile << -1 << " " << -1<< " " <<-1 << endl;
				}
				else if(heh0 == he1)
				{

					double a = edgeComponent(2 * e0.idx() + 1);
					double b = edgeComponent(2 * e1.idx() + 1);
					double c = edgeComponent(2 * e2.idx() + 1);
					//efile << 2 * e0.idx()+1 << " " << 2 * e1.idx()+1 << " " << 2 * e2.idx()+1 << endl;
					/*double a = IICsMatrix.row(0) (3 * e0.idx() + 2);
					double b = IICsMatrix.row(0)(3 * e1.idx() + 2);
					double c = IICsMatrix.row(0)(3 * e2.idx() + 2);*/
					/*double a = edgelengths [e0.idx()];
					double b = edgelengths [e1.idx()];
					double c = edgelengths [e2.idx()];*/

					double sqcos = (a * a + b * b - c * c)/ (2 * a * b);
					double squaresin = 1 - sqcos * sqcos;
					axis << sqcos, -sqrt(squaresin), 0;
					if(squaresin < 0)
						axis(1) = 0;
					if(sqcos > 1)
						axis(0) = 1;
					if(sqcos < -1)
						axis(0) = -1;
					//file1 << axis(0)<<" "<< axis(1)<<" "<< axis(2)<<endl;
					judge = 2;
				}
				else
				{
					double a = edgeComponent(2 * e0.idx() + 1);
					double b = edgeComponent(2 * e2.idx() + 1);
					double c = edgeComponent(2 * e1.idx() + 1);
					/*double a = IICsMatrix.row(0) (3 * e0.idx() + 2);
					double b = IICsMatrix.row(0)(3 * e1.idx() + 2);
					double c = IICsMatrix.row(0)(3 * e2.idx() + 2);*/
					/*double a = edgelengths [e0.idx()];
					double b = edgelengths [e1.idx()];
					double c = edgelengths [e2.idx()];*/
					//efile << 2 * e0.idx()+1 << " " << 2 * e2.idx()+1 << " " << 2 * e1.idx()+1 << endl;
					double sqcos = (a * a + b * b - c * c)/ (2 * a * b);
					double squaresin = 1 - sqcos * sqcos;

					axis << sqcos, sqrt(squaresin), 0;
					if(squaresin < 0)
						axis(1) = 0;
					if(sqcos > 1)
						axis(0) = 1;
					if(sqcos < -1)
						axis(0) = -1;
					//file1 << axis(0)<<" "<< axis(1)<<" "<< axis(2)<<endl;
					judge = 3;
				}
				//file4<< axis(0)<<" "<<axis(1)<<" "<<axis(2)<< " "<<judge<<endl;
				Eigen::Matrix3d matrixedge;
				matrixedge << 0, -axis(2), axis(1),
					axis(2), 0, -axis(0),
					-axis(1), axis(0), 0;
				Eigen::Matrix3d expedge = (theta * matrixedge).exp();
				/*求axis对应的第二个坐标轴*/
				Eigen::Vector3d axisY = computeCross (Eigen::Vector3d(0,0,1), axis);
				Eigen::Matrix3d matrixNormal;
				matrixNormal<< axis(0), axisY(0), 0,
				axis(1), axisY(1), 0,
				axis(2), axisY(2), 1;
				/*以axis作为坐标轴，求出第二个三角形的第一条边在该坐标轴下的坐标*/

				/*找出fh1中的三条半边的handle*/
				PGMesh::HalfedgeHandle he3;
				PGMesh::HalfedgeHandle he4;
				PGMesh::HalfedgeHandle he5;
				PGMesh::EdgeHandle e3;
				PGMesh::EdgeHandle e4;
				PGMesh::EdgeHandle e5;
				fheit = mesh_->fh_begin(fh1);
				he3 = *fheit;
				e3 = mesh_->edge_handle(he3);
				++fheit;
				he4 = *fheit;
				e4 = mesh_->edge_handle(he4);
				++fheit;
				he5 = *fheit;
				e5 = mesh_->edge_handle(he5);
				Eigen::Vector3d axis1;
				if(heh1 == he3)
				{
					axis1 << 1, 0, 0;
					//efile1 << -1 << " " << -1 << " " << -1 << endl;

				}
				else if(heh1 == he4)
				{
					double a = edgeComponent(2 * e3.idx() + 1);
					double b = edgeComponent(2 * e4.idx() + 1);
					double c = edgeComponent(2 * e5.idx() + 1);
					//efile1 << 2 * e3.idx() + 1 << " " << 2 * e4.idx() + 1 << " " << 2 * e5.idx() + 1 << endl;
					double sqcos = (a * a + b * b - c * c)/ (2 * a * b);
					double squaresin = 1 - sqcos * sqcos;
					axis1 << -sqcos, -sqrt(squaresin), 0;
					if(squaresin < 0)
						axis1(1) = 0;
					if(sqcos > 1)
						axis1(0) = -1;
					if(sqcos < -1)
						axis1(0) = 1;
				}
				else
				{
					double a = edgeComponent(2 * e3.idx() + 1);
					double b = edgeComponent(2 * e5.idx() + 1);
					double c = edgeComponent(2 * e4.idx() + 1);
					//efile1 << 2 * e3.idx() + 1 << " " << 2 * e5.idx() + 1 << " " << 2 * e4.idx() + 1 << endl;
					double sqcos = (a * a + b * b - c * c)/ (2 * a * b);
					double squaresin = 1 - sqcos * sqcos;
					axis1 << -sqcos, sqrt(squaresin), 0;
					if(squaresin < 0)
						axis1(1) = 0;
					if(sqcos > 1)
						axis1(0) = -1;
					if(sqcos < -1)
						axis1(0) = 1;
				}
				Eigen::Vector3d axisY1 = computeCross (Eigen::Vector3d(0,0,1), axis1);
				Eigen::Matrix3d matrixNormal1;
				matrixNormal1<< axis1(0), axisY1(0), 0,
				axis1(1), axisY1(1), 0,
				axis1(2), axisY1(2), 1;
				Eigen::Matrix3d rot = expedge * matrixNormal* matrixNormal1;
				//file1<< rot << endl << endl;
				//file3<< LocalTransforms[edgeNo][ithComp]<<endl<<endl;
				
				if((fh0.is_valid() && fh1.is_valid()))
				{
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 0, 9 * fromFaceNo +0, rot(0,0)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 0, 9 * fromFaceNo +1, rot(1,0)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 0, 9 * fromFaceNo +2,rot(2,0)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 0, 9 * toFaceNo +0, -1));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 1, 9 * fromFaceNo +3,rot(0,0)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 1, 9 * fromFaceNo +4,rot(1,0)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 1, 9 * fromFaceNo +5,rot(2,0)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 1, 9 * toFaceNo +3, -1));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 2, 9 * fromFaceNo +6,rot(0,0)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 2, 9 * fromFaceNo +7,rot(1,0)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 2, 9 * fromFaceNo +8,rot(2,0)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 2, 9 * toFaceNo +6,-1));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 3, 9 * fromFaceNo +0,rot(0,1)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 3, 9 * fromFaceNo +1,rot(1,1)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 3, 9 * fromFaceNo +2,rot(2,1)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 3, 9 * toFaceNo +1, -1));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 4, 9 * fromFaceNo +3, rot(0,1)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 4, 9 * fromFaceNo +4, rot(1,1)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 4, 9 * fromFaceNo +5, rot(2,1)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 4, 9 * toFaceNo +4, -1));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 5, 9 * fromFaceNo +6, rot(0,1)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 5, 9 * fromFaceNo +7, rot(1,1)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 5, 9 * fromFaceNo +8, rot(2,1)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 5, 9 * toFaceNo +7, -1));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 6, 9 * fromFaceNo +0, rot(0,2)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 6, 9 * fromFaceNo +1, rot(1,2)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 6, 9 * fromFaceNo +2, rot(2,2)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 6, 9 * toFaceNo +2, -1));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 7, 9 * fromFaceNo +3, rot(0,2)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 7, 9 * fromFaceNo +4, rot(1,2)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 7, 9 * fromFaceNo +5, rot(2,2)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 7, 9 * toFaceNo +5, -1));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 8, 9 * fromFaceNo +6, rot(0,2)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 8, 9 * fromFaceNo +7, rot(1,2)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 8, 9 * fromFaceNo +8, rot(2,2)));
					triple.push_back(Eigen::Triplet<double>(9 * edgeNo + 8, 9 * toFaceNo +8, -1));
					//*tmp_val++ = rot(0, 0);
					//*tmp_val++ = rot(1, 0);
					//*tmp_val++ = rot(2, 0);
					//*tmp_val++ = -1;
					//*tmp_val++ = rot(0, 0);
					//*tmp_val++ = rot(1, 0);
					//*tmp_val++ = rot(2, 0);
					//*tmp_val++ = -1;
					//*tmp_val++ = rot(0, 0);
					//*tmp_val++ = rot(1, 0);
					//*tmp_val++ = rot(2, 0);
					//*tmp_val++ = -1;
					//*tmp_val++ = rot(0, 1);
					//*tmp_val++ = rot(1, 1);
					//*tmp_val++ = rot(2, 1);
					//*tmp_val++ = -1;
					//*tmp_val++ = rot(0, 1);
					//*tmp_val++ = rot(1, 1);
					//*tmp_val++ = rot(2, 1);
					//*tmp_val++ = -1;
					//*tmp_val++ = rot(0, 1);
					//*tmp_val++ = rot(1, 1);
					//*tmp_val++ = rot(2, 1);
					//*tmp_val++ = -1;
					//*tmp_val++ = rot(0, 2);
					//*tmp_val++ = rot(1, 2);
					//*tmp_val++ = rot(2, 2);
					//*tmp_val++ = -1;
					//*tmp_val++ = rot(0, 2);
					//*tmp_val++ = rot(1, 2);
					//*tmp_val++ = rot(2, 2);
					//*tmp_val++ = -1;
					//*tmp_val++ = rot(0, 2);
					//*tmp_val++ = rot(1, 2);
					//*tmp_val++ = rot(2, 2);
					//*tmp_val++ = -1;

					//*tmp_row++ = 9 * edgeNo + 0;
					//*tmp_row++ = 9 * edgeNo + 0;
					//*tmp_row++ = 9 * edgeNo + 0;
					//*tmp_row++ = 9 * edgeNo + 0;
					//*tmp_row++ = 9 * edgeNo + 1;
					//*tmp_row++ = 9 * edgeNo + 1;
					//*tmp_row++ = 9 * edgeNo + 1;
					//*tmp_row++ = 9 * edgeNo + 1;
					//*tmp_row++ = 9 * edgeNo + 2;
					//*tmp_row++ = 9 * edgeNo + 2;
					//*tmp_row++ = 9 * edgeNo + 2;
					//*tmp_row++ = 9 * edgeNo + 2;
					//*tmp_row++ = 9 * edgeNo + 3;
					//*tmp_row++ = 9 * edgeNo + 3;
					//*tmp_row++ = 9 * edgeNo + 3;
					//*tmp_row++ = 9 * edgeNo + 3;
					//*tmp_row++ = 9 * edgeNo + 4;
					//*tmp_row++ = 9 * edgeNo + 4;
					//*tmp_row++ = 9 * edgeNo + 4;
					//*tmp_row++ = 9 * edgeNo + 4;
					//*tmp_row++ = 9 * edgeNo + 5;
					//*tmp_row++ = 9 * edgeNo + 5;
					//*tmp_row++ = 9 * edgeNo + 5;
					//*tmp_row++ = 9 * edgeNo + 5;
					//*tmp_row++ = 9 * edgeNo + 6;
					//*tmp_row++ = 9 * edgeNo + 6;
					//*tmp_row++ = 9 * edgeNo + 6;
					//*tmp_row++ = 9 * edgeNo + 6;
					//*tmp_row++ = 9 * edgeNo + 7;
					//*tmp_row++ = 9 * edgeNo + 7;
					//*tmp_row++ = 9 * edgeNo + 7;
					//*tmp_row++ = 9 * edgeNo + 7;
					//*tmp_row++ = 9 * edgeNo + 8;
					//*tmp_row++ = 9 * edgeNo + 8;
					//*tmp_row++ = 9 * edgeNo + 8;
					//*tmp_row++ = 9 * edgeNo + 8;

					//*temp_col++ = 9 * fromFaceNo + 0;
					//*temp_col++ = 9 * fromFaceNo + 1;
					//*temp_col++ = 9 * fromFaceNo + 2;
					//*temp_col++ = 9 * toFaceNo + 0;
					//*temp_col++ = 9 * fromFaceNo + 3;
					//*temp_col++ = 9 * fromFaceNo + 4;
					//*temp_col++ = 9 * fromFaceNo + 5;
					//*temp_col++ = 9 * toFaceNo + 3;
					//*temp_col++ = 9 * fromFaceNo + 6;
					//*temp_col++ = 9 * fromFaceNo + 7;
					//*temp_col++ = 9 * fromFaceNo + 8;
					//*temp_col++ = 9 * toFaceNo + 6;
					//*temp_col++ = 9 * fromFaceNo + 0;
					//*temp_col++ = 9 * fromFaceNo + 1;
					//*temp_col++ = 9 * fromFaceNo + 2;
					//*temp_col++ = 9 * toFaceNo + 1;
					//*temp_col++ = 9 * fromFaceNo + 3;
					//*temp_col++ = 9 * fromFaceNo + 4;
					//*temp_col++ = 9 * fromFaceNo + 5;
					//*temp_col++ = 9 * toFaceNo + 4;
					//*temp_col++ = 9 * fromFaceNo + 6;
					//*temp_col++ = 9 * fromFaceNo + 7;
					//*temp_col++ = 9 * fromFaceNo + 7;
					//*temp_col++ = 9 * toFaceNo + 7;
					//*temp_col++ = 9 * fromFaceNo + 0;
					//*temp_col++ = 9 * fromFaceNo + 1;
					//*temp_col++ = 9 * fromFaceNo + 2;
					//*temp_col++ = 9 * toFaceNo + 2;
					//*temp_col++ = 9 * fromFaceNo + 3;
					//*temp_col++ = 9 * fromFaceNo + 4;
					//*temp_col++ = 9 * fromFaceNo + 5;
					//*temp_col++ = 9 * toFaceNo + 5;
					//*temp_col++ = 9 * fromFaceNo + 6;
					//*temp_col++ = 9 * fromFaceNo + 7;
					//*temp_col++ = 9 * fromFaceNo + 8;
					//*temp_col++ = 9 * toFaceNo + 8;
				}
			}
			//设置锚点
			for(int l = 0; l < Faceanchors_.size(); l++)
			{
				Eigen::Matrix3d anchorFrame = faceAnchorFrames[anchorIthSeq][l].localframe; 
				triple.push_back(Eigen::Triplet<double>(9 * nE + l * 9 + 0, 9 * Faceanchors_[l] + 0, (1.0 / std::sqrt(vertexAnchorCords[0].size())))); 
				triple.push_back(Eigen::Triplet<double>(9 * nE + l * 9 + 1, 9 * Faceanchors_[l] + 1, (1.0 / std::sqrt(vertexAnchorCords[0].size()))));
				triple.push_back(Eigen::Triplet<double>(9 * nE + l * 9 + 2, 9 * Faceanchors_[l] + 2, (1.0 / std::sqrt(vertexAnchorCords[0].size()))));
				triple.push_back(Eigen::Triplet<double>(9 * nE + l * 9 + 3, 9 * Faceanchors_[l] + 3, (1.0 / std::sqrt(vertexAnchorCords[0].size()))));
				triple.push_back(Eigen::Triplet<double>(9 * nE + l * 9 + 4, 9 * Faceanchors_[l] + 4, (1.0 / std::sqrt(vertexAnchorCords[0].size()))));
				triple.push_back(Eigen::Triplet<double>(9 * nE + l * 9 + 5, 9 * Faceanchors_[l] + 5, (1.0 / std::sqrt(vertexAnchorCords[0].size()))));
				triple.push_back(Eigen::Triplet<double>(9 * nE + l * 9 + 6, 9 * Faceanchors_[l] + 6, (1.0 / std::sqrt(vertexAnchorCords[0].size()))));
				triple.push_back(Eigen::Triplet<double>(9 * nE + l * 9 + 7, 9 * Faceanchors_[l] + 7, (1.0 / std::sqrt(vertexAnchorCords[0].size()))));
				triple.push_back(Eigen::Triplet<double>(9 * nE + l * 9 + 8, 9 * Faceanchors_[l] + 8, (1.0 / std::sqrt(vertexAnchorCords[0].size()))));
				frameCor(9 * nE + l * 9 + 0) = (1.0 / std::sqrt(vertexAnchorCords[0].size())) * anchorFrame(0,0);
				frameCor(9 * nE + l * 9 + 1) = (1.0 / std::sqrt(vertexAnchorCords[0].size()))* anchorFrame(0,1);
				frameCor(9 * nE + l * 9 + 2) = (1.0 / std::sqrt(vertexAnchorCords[0].size()))*anchorFrame(0,2);
				frameCor(9 * nE + l * 9 + 3) = (1.0 / std::sqrt(vertexAnchorCords[0].size()))*anchorFrame(1,0);
				frameCor(9 * nE + l * 9 + 4) = (1.0 / std::sqrt(vertexAnchorCords[0].size()))* anchorFrame(1,1);
				frameCor(9 * nE + l * 9 + 5) = (1.0 / std::sqrt(vertexAnchorCords[0].size()))* anchorFrame(1,2);
				frameCor(9 * nE + l * 9 + 6) = (1.0 / std::sqrt(vertexAnchorCords[0].size()))* anchorFrame(2,0);
				frameCor(9 * nE + l * 9 + 7) = (1.0 / std::sqrt(vertexAnchorCords[0].size()))* anchorFrame(2,1);
				frameCor(9 * nE + l * 9 + 8) = (1.0 / std::sqrt(vertexAnchorCords[0].size()))* anchorFrame(2,2);

				//*tmp_val++ = (1.0 / std::sqrt(vertexAnchorCords[0].size()));
				//*tmp_val++ = (1.0 / std::sqrt(vertexAnchorCords[0].size()));
				//*tmp_val++ = (1.0 / std::sqrt(vertexAnchorCords[0].size()));
				//*tmp_val++ = (1.0 / std::sqrt(vertexAnchorCords[0].size()));
				//*tmp_val++ = (1.0 / std::sqrt(vertexAnchorCords[0].size()));
				//*tmp_val++ = (1.0 / std::sqrt(vertexAnchorCords[0].size()));
				//*tmp_val++ = (1.0 / std::sqrt(vertexAnchorCords[0].size()));
				//*tmp_val++ = (1.0 / std::sqrt(vertexAnchorCords[0].size()));
				//*tmp_val++ = (1.0 / std::sqrt(vertexAnchorCords[0].size()));

				//*tmp_row++ = 9 * nE + l * 9 + 0;
				//*tmp_row++ = 9 * nE + l * 9 + 1;
				//*tmp_row++ = 9 * nE + l * 9 + 2;
				//*tmp_row++ = 9 * nE + l * 9 + 3;
				//*tmp_row++ = 9 * nE + l * 9 + 4;
				//*tmp_row++ = 9 * nE + l * 9 + 5;
				//*tmp_row++ = 9 * nE + l * 9 + 6;
				//*tmp_row++ = 9 * nE + l * 9 + 7;
				//*tmp_row++ = 9 * nE + l * 9 + 8;

				//*temp_col++ = 9 * Faceanchors_[l] + 0;
				//*temp_col++ = 9 * Faceanchors_[l] + 1;
				//*temp_col++ = 9 * Faceanchors_[l] + 2;
				//*temp_col++ = 9 * Faceanchors_[l] + 3;
				//*temp_col++ = 9 * Faceanchors_[l] + 4;
				//*temp_col++ = 9 * Faceanchors_[l] + 5;
				//*temp_col++ = 9 * Faceanchors_[l] + 6;
				//*temp_col++ = 9 * Faceanchors_[l] + 7;
				//*temp_col++ = 9 * Faceanchors_[l] + 8;

			}
			
			Coefficients.setFromTriplets(triple.begin(), triple.end());  
			triple.clear();
			SimplicialCholesky<SparseMatrix<double>> solver;
			//ofstream solvertxt("solve1.txt");
			//solvertxt<<Coefficients.transpose() * Coefficients;
			solver.compute(Coefficients.transpose() * Coefficients);
			if(solver.info()!=Success)
			{
				/// decomposit ion failed
				std::cout<<"Decomposition failed" << std::endl;
				return;
			}
			Eigen::VectorXd x_up; 
			x_up = solver.solve(Coefficients.transpose() * frameCor);
			/*solvertxt<< x_up;
			solvertxt.close();*/
			//solvertxt.close();
			/*Given the frames of all faces, reconstruct the vertex coordinates*/
			//Eigen::SparseMatrix<double> vertexCoff(9 * nF + 3 , 3 * nV);
			Eigen::VectorXd fvertex(9 * nF + 3 * vertexAnchorCords[0].size());
			Eigen::VectorXd y_up(3 * nV);
			//ofstream axisview("axis.txt");
			//ofstream edgeview("edgevector.txt");
			for(PGMesh::FaceIter fit = mesh_->faces_begin(); fit != mesh_->faces_end(); ++fit)
			{
				PGMesh::FaceHandle fh = *fit;
				int faceNo = fh.idx();
				/*找出fh0中的三条半边的handle*/
				PGMesh::HalfedgeHandle he0;
				PGMesh::HalfedgeHandle he1;
				PGMesh::HalfedgeHandle he2;
				PGMesh::EdgeHandle e0;
				PGMesh::EdgeHandle e1;
				PGMesh::EdgeHandle e2;
				PGMesh::FHIter fheit = mesh_->fh_begin(fh);
				he0 = *fheit;
				e0 = mesh_->edge_handle(he0);
				++fheit;
				he1 = *fheit;
				e1 = mesh_->edge_handle(he1);
				++fheit;
				he2 = *fheit;
				e2 = mesh_->edge_handle(he2);
				/*找出三条边的长度*/
				double a, b, c;
				a = edgeComponent(2 * e0.idx() + 1);
				b = edgeComponent(2 * e1.idx() + 1);
				c = edgeComponent(2 * e2.idx() + 1);
				 /*a = IICsMatrix.row(0) (3 * e0.idx() + 2);
				 b = IICsMatrix.row(0)(3 * e1.idx() + 2);
				 c = IICsMatrix.row(0)(3 * e2.idx() + 2);*/
				 /*a = edgelengths [e0.idx()];
				 b = edgelengths [e1.idx()];
				 c = edgelengths [e2.idx()];*/
				/*找出fh1中的三个点的handle*/
				PGMesh::VertexHandle v0;
				PGMesh::VertexHandle v1;
				PGMesh::VertexHandle v2;
				v0 = mesh_->from_vertex_handle(he0);
				v1 = mesh_->to_vertex_handle(he0);
				v2 = mesh_->to_vertex_handle(he1);
				
				/*提取出该面的标架中的a，b轴*/
				Eigen::Matrix3d frameFace;
				frameFace << x_up(9 * faceNo + 0), x_up(9 * faceNo + 1), x_up(9 * faceNo + 2),
					x_up(9 * faceNo + 3), x_up(9 * faceNo + 4), x_up(9 * faceNo + 5),
					x_up(9 * faceNo + 6), x_up(9 * faceNo + 7), x_up(9 * faceNo + 8);
				Eigen::JacobiSVD<MatrixXd> svd( frameFace, ComputeThinU | ComputeThinV);
				Eigen::MatrixXd u = svd.matrixU();
				Eigen::MatrixXd v = svd.matrixV();
				Eigen::Matrix3d U = u * v.transpose();
				Eigen::Vector3d a_axis, b_axis;
				a_axis = U.col(0);
				b_axis = U.col(1);
				a_axis = a_axis/ a_axis.norm();
				b_axis = b_axis /b_axis.norm();
				//axisview<< a_axis(0)<<" "<<a_axis(1)<<" "<<a_axis(2)<<endl;
				//axisview<< b_axis(0)<<" "<<b_axis(1)<<" "<<b_axis(2)<<endl;
				/*提取出角度 alpha, beta*/
				double cosalpha, sinalpha, cosbeta, sinbeta;
				cosalpha = (a * a + c * c - b * b)/ (2 * a * c);
				sinalpha = sqrt(1 - cosalpha * cosalpha);
				cosbeta = (a * a + b * b - c * c)/ (2 * a * b);
				sinbeta = sqrt(1 - cosbeta * cosbeta);
				
				if(cosalpha > 1)
				{
					cosalpha = 1;
					sinalpha = 0;
				}
				if(cosalpha < -1)
				{
					cosalpha = -1;
					sinalpha = 0;
				}
					
				if(cosbeta > 1)
					{cosbeta = 1;
				sinbeta = 0;}
				if(cosbeta < -1)
					{cosbeta = -1;
				sinbeta = 0;}
				//axisview<< cosalpha<<" "<<sinalpha<< " " << cosbeta<<" "<< sinbeta<<endl;
				/*得到各边向量*/
				Eigen::Vector3d v0v1, v1v2, v2v0;
				v0v1 =  a * a_axis;
				v1v2 = -b * (cosbeta * a_axis - sinbeta * b_axis);
				v2v0 = -c * (cosalpha * a_axis + sinalpha * b_axis);
				//edgeview<< v0v1(0)<<" "<<v0v1(1)<<" "<<v0v1(2)<<endl;
				//edgeview<< v1v2(0)<<" "<<v1v2(1)<<" "<<v1v2(2)<<endl;
				//edgeview<< v2v0(0)<<" "<<v2v0(1)<<" "<<v2v0(2)<<endl;
				/*填入矩阵和向量元素*/
				//triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 0, 3 * v1.idx() + 0, 1));
				//triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 0, 3 * v0.idx() + 0, -1));
				fvertex(faceNo * 9 + 0) = v0v1(0);
				//triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 1, 3 * v1.idx() + 1, 1));
				//triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 1, 3 * v0.idx() + 1, -1));
				fvertex(faceNo * 9 + 1) = v0v1(1);
				//triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 2, 3 * v1.idx() + 2, 1));
				//triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 2, 3 * v0.idx() + 2, -1));
				fvertex(faceNo * 9 + 2) = v0v1(2);
				//triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 3, 3 * v2.idx() + 0, 1));
				//triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 3, 3 * v1.idx() + 0, -1));
				fvertex(faceNo * 9 + 3) = v1v2(0);
				//triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 4, 3 * v2.idx() + 1, 1));
				//triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 4, 3 * v1.idx() + 1, -1));
				fvertex(faceNo * 9 + 4) = v1v2(1);
				//triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 5, 3 * v2.idx() + 2, 1));
				//triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 5, 3 * v1.idx() + 2, -1));
				fvertex(faceNo * 9 + 5) = v1v2(2);
				//triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 6, 3 * v0.idx() + 0, 1));
				//triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 6, 3 * v2.idx() + 0, -1));
				fvertex(faceNo * 9 + 6) = v2v0(0);
				//triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 7, 3 * v0.idx() + 1, 1));
				//triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 7, 3 * v2.idx() + 1, -1));
				fvertex(faceNo * 9 + 7) = v2v0(1);
				//triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 8, 3 * v0.idx() + 2, 1));
				//triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 8, 3 * v2.idx() + 2, -1));
				fvertex(faceNo * 9 + 8) = v2v0(2);
			}
			//axisview.close();
			//edgeview.close();
			/*添加锚点*/
			for(int j = 0; j < vertexAnchorCords[0].size(); j++)
			{
				//triple.push_back(Eigen::Triplet<double>(9 * nF + 0, 3 * 0 + 0, 1));
			fvertex(9 * nF + 3 * j+ 0) = vertexAnchorCords[anchorIthSeq][j].cord(0);
					
				//triple.push_back(Eigen::Triplet<double>(9 * nF  + 1, 3 * 0+ 1, 1));
				fvertex(9 * nF  + 3 * j +1) = vertexAnchorCords[anchorIthSeq][j].cord(1);
				//triple.push_back(Eigen::Triplet<double>(9 * nF + 2, 3 * 0 + 2, 1));
				fvertex(9 * nF + 3 * j + 2) = vertexAnchorCords[anchorIthSeq][j].cord(2);
			}
			//vertexCoff.setFromTriplets(triple.begin(), triple.end());  
			//triple.clear();
			//ofstream solvertxt("solve1.txt");
			//solvertxt<<Coefficients.transpose() * Coefficients;
			//SimplicialCholesky<SparseMatrix<double>> solver1;
			//solver1.compute(vertexCoff.transpose() * vertexCoff);
			/*if(solver1.info()!=Success)
			{
				/// decomposit ion failed
				std::cout<<"Decomposition failed" << std::endl;
				return;
			}*/
			y_up = presolver.solve(vertexCoff.transpose() * fvertex);
			cordComponent = y_up;
		}
		/** 先解presolver(考虑锚点由点的线性混合)*/
		void SparseLocalizedForConnectMap::presolve1()
		{
			vertexCoff = Eigen::SparseMatrix<double>(9 * nF + 3 * vertexAnchorCords[0].size(), 3 * nV);
			//presolver = SimplicialCholesky<SparseMatrix<double>>();
			std::vector<Eigen::Triplet<double> > triple;
			triple.clear();
			for (PGMesh::FaceIter fit = mesh_->faces_begin(); fit != mesh_->faces_end(); ++fit)
			{
				PGMesh::FaceHandle fh = *fit;
				int faceNo = fh.idx();
				/*找出fh0中的三条半边的handle*/
				PGMesh::HalfedgeHandle he0;
				PGMesh::HalfedgeHandle he1;
				PGMesh::HalfedgeHandle he2;
				PGMesh::EdgeHandle e0;
				PGMesh::EdgeHandle e1;
				PGMesh::EdgeHandle e2;
				PGMesh::FHIter fheit = mesh_->fh_begin(fh);
				he0 = *fheit;
				e0 = mesh_->edge_handle(he0);
				++fheit;
				he1 = *fheit;
				e1 = mesh_->edge_handle(he1);
				++fheit;
				he2 = *fheit;
				e2 = mesh_->edge_handle(he2);
				PGMesh::VertexHandle v0;
				PGMesh::VertexHandle v1;
				PGMesh::VertexHandle v2;
				v0 = mesh_->from_vertex_handle(he0);
				v1 = mesh_->to_vertex_handle(he0);
				v2 = mesh_->to_vertex_handle(he1);


				/*填入矩阵和向量元素*/
				triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 0, 3 * v1.idx() + 0, 1));
				triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 0, 3 * v0.idx() + 0, -1));
				//fvertex(faceNo * 9 + 0) = v0v1(0);
				triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 1, 3 * v1.idx() + 1, 1));
				triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 1, 3 * v0.idx() + 1, -1));
				//fvertex(faceNo * 9 + 1) = v0v1(1);
				triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 2, 3 * v1.idx() + 2, 1));
				triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 2, 3 * v0.idx() + 2, -1));
				//fvertex(faceNo * 9 + 2) = v0v1(2);
				triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 3, 3 * v2.idx() + 0, 1));
				triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 3, 3 * v1.idx() + 0, -1));
				//fvertex(faceNo * 9 + 3) = v1v2(0);
				triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 4, 3 * v2.idx() + 1, 1));
				triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 4, 3 * v1.idx() + 1, -1));
				//fvertex(faceNo * 9 + 4) = v1v2(1);
				triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 5, 3 * v2.idx() + 2, 1));
				triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 5, 3 * v1.idx() + 2, -1));
				//fvertex(faceNo * 9 + 5) = v1v2(2);
				triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 6, 3 * v0.idx() + 0, 1));
				triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 6, 3 * v2.idx() + 0, -1));
				//fvertex(faceNo * 9 + 6) = v2v0(0);
				triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 7, 3 * v0.idx() + 1, 1));
				triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 7, 3 * v2.idx() + 1, -1));
				//fvertex(faceNo * 9 + 7) = v2v0(1);
				triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 8, 3 * v0.idx() + 2, 1));
				triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 8, 3 * v2.idx() + 2, -1));
				//fvertex(faceNo * 9 + 8) = v2v0(2);
			}
			/*添加锚点*/
			for (int j = 0; j < vertexAnchorCords[0].size(); j++)
			{
				for (int l = 0; l < vertexAnchorCords[0][j].idx_list.size(); l++)
				{
					triple.push_back(Eigen::Triplet<double>(9 * nF + 3 * j + 0, 3 * vertexAnchorCords[0][j].idx_list[l] + 0, 1.0 * vertexAnchorCords[0][j].weights[l]));
					//fvertex(9 * nF + 0) = SpatialTemData[0][anchorIthSeq](0);
					triple.push_back(Eigen::Triplet<double>(9 * nF + 3 * j + 1, 3 * vertexAnchorCords[0][j].idx_list[l] + 1, 1.0 * vertexAnchorCords[0][j].weights[l]));
					//fvertex(9 * nF  + 1) = SpatialTemData[0][anchorIthSeq](1);
					triple.push_back(Eigen::Triplet<double>(9 * nF + 3 * j + 2, 3 * vertexAnchorCords[0][j].idx_list[l] + 2, 1.0 * vertexAnchorCords[0][j].weights[l]));
					//fvertex(9 * nF + 2) = SpatialTemData[0][anchorIthSeq](2);
					vertexCoff.setFromTriplets(triple.begin(), triple.end());
				}
				
			}

			triple.clear();
			//ofstream solvertxt("solve1.txt");
			//solvertxt<<Coefficients.transpose() * Coefficients;
			presolver.compute(vertexCoff.transpose() * vertexCoff);
			if (presolver.info() != Success)
			{
				/// decomposit ion failed
				std::cout << "Decomposition failed" << std::endl;
				return;
			}
		}
		/** 先解presolver*/
		void SparseLocalizedForConnectMap::presolve()
		{
			vertexCoff = Eigen::SparseMatrix<double> (9 * nF + 3 * vertexAnchorCords[0].size() , 3 * nV);
			//presolver = SimplicialCholesky<SparseMatrix<double>>();
			std::vector<Eigen::Triplet<double> > triple;  
			triple.clear();
			for(PGMesh::FaceIter fit = mesh_->faces_begin(); fit != mesh_->faces_end(); ++fit)
			{
				PGMesh::FaceHandle fh = *fit;
				int faceNo = fh.idx();
				/*找出fh0中的三条半边的handle*/
				PGMesh::HalfedgeHandle he0;
				PGMesh::HalfedgeHandle he1;
				PGMesh::HalfedgeHandle he2;
				PGMesh::EdgeHandle e0;
				PGMesh::EdgeHandle e1;
				PGMesh::EdgeHandle e2;
				PGMesh::FHIter fheit = mesh_->fh_begin(fh);
				he0 = *fheit;
				e0 = mesh_->edge_handle(he0);
				++fheit;
				he1 = *fheit;
				e1 = mesh_->edge_handle(he1);
				++fheit;
				he2 = *fheit;
				e2 = mesh_->edge_handle(he2);
				PGMesh::VertexHandle v0;
				PGMesh::VertexHandle v1;
				PGMesh::VertexHandle v2;
				v0 = mesh_->from_vertex_handle(he0);
				v1 = mesh_->to_vertex_handle(he0);
				v2 = mesh_->to_vertex_handle(he1);
				
				
				/*填入矩阵和向量元素*/
				triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 0, 3 * v1.idx() + 0, 1));
				triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 0, 3 * v0.idx() + 0, -1));
				//fvertex(faceNo * 9 + 0) = v0v1(0);
				triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 1, 3 * v1.idx() + 1, 1));
				triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 1, 3 * v0.idx() + 1, -1));
				//fvertex(faceNo * 9 + 1) = v0v1(1);
				triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 2, 3 * v1.idx() + 2, 1));
				triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 2, 3 * v0.idx() + 2, -1));
				//fvertex(faceNo * 9 + 2) = v0v1(2);
				triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 3, 3 * v2.idx() + 0, 1));
				triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 3, 3 * v1.idx() + 0, -1));
				//fvertex(faceNo * 9 + 3) = v1v2(0);
				triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 4, 3 * v2.idx() + 1, 1));
				triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 4, 3 * v1.idx() + 1, -1));
				//fvertex(faceNo * 9 + 4) = v1v2(1);
				triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 5, 3 * v2.idx() + 2, 1));
				triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 5, 3 * v1.idx() + 2, -1));
				//fvertex(faceNo * 9 + 5) = v1v2(2);
				triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 6, 3 * v0.idx() + 0, 1));
				triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 6, 3 * v2.idx() + 0, -1));
				//fvertex(faceNo * 9 + 6) = v2v0(0);
				triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 7, 3 * v0.idx() + 1, 1));
				triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 7, 3 * v2.idx() + 1, -1));
				//fvertex(faceNo * 9 + 7) = v2v0(1);
				triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 8, 3 * v0.idx() + 2, 1));
				triple.push_back(Eigen::Triplet<double>(faceNo * 9 + 8, 3 * v2.idx() + 2, -1));
				//fvertex(faceNo * 9 + 8) = v2v0(2);
			}
			/*添加锚点*/
			for (int j = 0; j < vertexAnchorCords[0].size(); j++)
			{
				triple.push_back(Eigen::Triplet<double>(9 * nF + 3 * j + 0, 3 *vertexAnchorCords[0][j].idx  + 0, 1.0));
				//fvertex(9 * nF + 0) = SpatialTemData[0][anchorIthSeq](0);
				triple.push_back(Eigen::Triplet<double>(9 * nF + 3 * j + 1, 3 * vertexAnchorCords[0][j].idx + 1, 1.0));
				//fvertex(9 * nF  + 1) = SpatialTemData[0][anchorIthSeq](1);
				triple.push_back(Eigen::Triplet<double>(9 * nF + 3 * j + 2, 3 * vertexAnchorCords[0][j].idx  + 2, 1.0));
				//fvertex(9 * nF + 2) = SpatialTemData[0][anchorIthSeq](2);
				vertexCoff.setFromTriplets(triple.begin(), triple.end());  
			}
			
			triple.clear();
			//ofstream solvertxt("solve1.txt");
			//solvertxt<<Coefficients.transpose() * Coefficients;
			presolver.compute(vertexCoff.transpose() * vertexCoff);
			if(presolver.info()!=Success)
			{
				/// decomposit ion failed
				std::cout<<"Decomposition failed" << std::endl;
				return;
			}
		}
		void SparseLocalizedForConnectMap::outMesh2(Eigen::VectorXd cordComponent, string path)
		{
			ofstream offfile(path);
			//offfile << "COFF"<<endl;
			//offfile << nV << " " << nF<< " "<< 0 <<endl;
			
			for (int i = 0; i< nV; i++)
			{
				offfile << "v" << " " << cordComponent(3 * i + 0) << " " <<cordComponent(3 * i + 1) << " " << cordComponent(3 * i + 2)  << endl;
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
		}
		/** 根据点坐标导出网格*/
		void SparseLocalizedForConnectMap::outMesh1(Eigen::VectorXd cordComponent, string path)
		{
			ofstream offfile(path);
			offfile << "OFF"<<endl;
			offfile << nV << " " << nF<< " "<< 0 <<endl;
			
			for (int i = 0; i< nV; i++)
			{
				//std::cout << cordComponent(3 * i + 0) << " " << cordComponent(3 * i + 1) << " " << cordComponent(3 * i + 2) << endl;
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

		void SparseLocalizedForConnectMap::outMesh1(Eigen::MatrixXd cordComponent, string path)
		{
			ofstream offfile(path);
			offfile << "OFF" << endl;
			offfile << nV << " " << nF << " " << 0 << endl;

			for (int i = 0; i< nV; i++)
			{
				offfile << cordComponent(i , 0) << " " << cordComponent(i , 1) << " " << cordComponent(i , 2) << endl;
			}
			for (PGMesh::FaceIter fit = mesh_->faces_begin(); fit != mesh_->faces_end(); ++fit)
			{
				offfile << "3" << " ";
				for (PGMesh::FVIter fvit = mesh_->fv_begin(*fit); fvit != mesh_->fv_end(*fit); ++fvit)
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

		/** 根据点坐标导出网格*/
		void SparseLocalizedForConnectMap::outMesh(Eigen::VectorXd cordComponent, char* path)
		{
			ofstream offfile(path);
			//offfile << "COFF"<<endl;
			//offfile << nV << " " << nF<< " "<< 0 <<endl;
			
			for (int i = 0; i< nV; i++)
			{
				offfile << "v" << " " << cordComponent(3 * i + 0) << " " <<cordComponent(3 * i + 1) << " " << cordComponent(3 * i + 2)  << endl;
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
		}
		void SparseLocalizedForConnectMap::outMesh(Eigen::MatrixXd cordComponent, char* path)
		{
			ofstream offfile(path);
			//offfile << "COFF"<<endl;
			//offfile << nV << " " << nF<< " "<< 0 <<endl;

			for (int i = 0; i< nV; i++)
			{
				offfile << "v" << " " << cordComponent( i , 0) << " " << cordComponent(i,  1) << " " << cordComponent(i , 2) << endl;
			}
			for (PGMesh::FaceIter fit = mesh_->faces_begin(); fit != mesh_->faces_end(); ++fit)
			{
				offfile << "f" << " ";
				for (PGMesh::FVIter fvit = mesh_->fv_begin(*fit); fvit != mesh_->fv_end(*fit); ++fvit)
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

		/** 可视化component(给定每条边的误差)*/
		void SparseLocalizedForConnectMap::visualizeComponent(Eigen::VectorXd cordComponent, Eigen::VectorXd edgeError, char* path)
		{
			
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
					int edgeidx = mesh_->edge_handle(*veit).idx();
					pointerror += edgeError(edgeidx);
					if (edgeError(edgeidx) > errors[i])
						errors[i] = edgeError(edgeidx);
					j++;
				}
				//求该点对应相邻边的误差的平均值
				//errors[i] = pointerror/j;
			}
			colorize(errors, cordComponent, path);
		}

		/** 给出顶点权重，给模型上色*/
		void SparseLocalizedForConnectMap::colorize(std::vector<double> &errorss, Eigen::VectorXd cordComponent, char* path)
		{
			
			std::vector<double> errors = errorss;
			ofstream offfile(path);
			ofstream offfile1("colorout.txt");
			offfile << "COFF"<<endl;
			offfile << nV << " " << nF<< " "<< 0 <<endl;
			double  maxElement = errors[0], minElement = errors[0];
			for (int i = 0; i < nV; ++i)
			{

				maxElement = maxElement > errors[i] ? maxElement : errors[i];
				minElement = minElement < errors[i] ? minElement : errors[i];
			}
			
			float range =(maxElement - minElement);
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
				double t = (v - minElement) / range ;
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
					if (t > 1)
						t = 1;
					rr = (2 * t - 1) * 1 + (2 - 2 * t) * 1;
					gg = (2 * t - 1) * 0 + (2 - 2 * t) * 1;
					bb = (2 * t - 1) * 0 + (2 - 2 * t) * 0;
				}
				//int color1 = (v - minElement) * 192 / range;
				//int color2 = (v - minElement) * 63 / range;
				int colorPersent1 = (v - minElement) / range * 99;
				unsigned char* colorPtr = colormap5 + colorPersent1 * 3;
				//c.values_[0] = colorPtr[0];
				//c.values_[1] = colorPtr[1];
				//c.values_[2] = colorPtr[2];
				Eigen::Vector3d position = cordComponent.segment(3 * i, 3);
				offfile << "v" << " " << position( 0) << " " <<position(1) << " " << position(2)  <<" " ;
				offfile <<(int) (colorPtr[0])  <<" "<< int(colorPtr[1])<< " " << int (colorPtr[2])<< endl;
				offfile1 <<(int) (colorPtr[0])  <<" "<< int(colorPtr[1])<< " " << int (colorPtr[2])<< endl;
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

		/** 可视化网格()*/
		void SparseLocalizedForConnectMap::visualizeMesh(Eigen::VectorXd cordComponent, int ith, char* path, double radius)
		{
			ofstream offfile(path);
			offfile << "COFF"<<endl;
				offfile << nV << " " << nF<< " "<< 0 <<endl;
			Eigen::VectorXd errorVec(nV);
			Eigen::Vector3d error1;
			Eigen::Vector3d error2;
			double maxerror = 0;
			double minerror = 1000;
			for(int i = 0; i < nV; i++)
			{
				 error1 = cordComponent.segment(3*i, 3) ;
				 error2 = SpatialTemData[i][ith];
				 errorVec(i) = (( cordComponent.segment(3*i, 3) - SpatialTemData[i][ith]).norm());
				if(errorVec(i) > maxerror)
					maxerror = errorVec(i);
				if(errorVec(i) < minerror)
					minerror = errorVec(i);
			}
			double meanerror = errorVec.mean();
			//double maxerror = errorVec.
			//ofstream offfile("connectOut.obj");
			double r;
			if(maxerror > radius)
				r = maxerror;
			else
				r = radius;
			for (int i = 0; i< nV; i++)
			{
				//int color =(int) ((errorVec(i))  * 255/radius);
				int colorPersent1 = (errorVec(i) - 0) / r * 64;
				unsigned char* colorPtr = colormap7 + colorPersent1 * 3;
				offfile << "v" << " " << cordComponent(3 * i + 0) << " " <<cordComponent(3 * i + 1) << " " << cordComponent(3 * i + 2)  <<" ";
				offfile <<(int) (colorPtr[0])  <<" "<< int(colorPtr[1])<< " " << int (colorPtr[2])<< endl;
			}
			for (PGMesh::FaceIter fit = mesh_->faces_begin(); fit != mesh_->faces_end(); ++ fit)
			{
				offfile << "f" << " ";
				/*for (PGMesh::FVIter fvit = mesh_->fv_begin(fit.handle()); fvit != mesh_->fv_end(fit.handle()); ++ fvit)
				{
					offfile << fvit.handle().idx() + 1 << " ";
				} 
				
				offfile << endl;*/
				PGMesh::FVIter fvit = mesh_->fv_begin(*fit);
				PGMesh::VertexHandle vv0 = *fvit;
				++fvit;
				PGMesh::VertexHandle vv1 = *fvit;
				++fvit;
				PGMesh::VertexHandle vv2 = *fvit;
				offfile << vv2.idx() + 1 << " "<< vv1.idx() + 1 << " " << vv0.idx() + 1 << " "  ;
				offfile << endl;
			}
			offfile << "# maxerror " << maxerror/radius << endl;
			offfile << "# meanerror "<< meanerror/radius<<endl;
			offfile.close();
		}

		/** 利用IIC components重构第ith帧*/
		void SparseLocalizedForConnectMap::constructFrameFromComponents(int ith, int compNum, string str)
		{
			char *tmp = new char[120];
			sprintf_s(tmp,120,  "E:\\new_STED\\original\\samba1\\mesh_%04d.off", ith);//原始序列的第i帧
			PGMesh *pg = new PGMesh();
			OpenMesh::IO::read_mesh(*pg, tmp);
			updateDeformed(pg);
			Eigen::VectorXd reconst(2 * nE);
			reconst = W.row(ith) * shapeBasis;
			char *pPath = new char[30];
			Eigen::VectorXd l0 = DiEdgeDataMatrix.row(0);
			Eigen::VectorXd l1(3 * nV);
			reconst += l0;
			reconstructionFromDiEdges(reconst, l1, 1, 0.5, 1);
			outMesh1(l1, str);
		}
		/** 判断误差IIC components重构第ith帧*/
		void SparseLocalizedForConnectMap::constructFrameFromComponentsError(int ith, int compNum, string str, Eigen::VectorXd & errorIncrement)
		{
			Eigen::VectorXd reconst(2 * nE);
			reconst = W.row(ith) * shapeBasis;
			char *pPath = new char[30];
			Eigen::VectorXd l0 = DiEdgeDataMatrix.row(0);
			Eigen::VectorXd l1(3 * nV);
			reconst += l0;
			reconstructionFromDiEdges(reconst, l1, ith, 0.5, ith);
			for(int i = 0; i < nV; i++)
			{
				double error = (SpatialTemData[i][ith] - l1.segment(3 * i, 3)).norm();
				errorIncrement(i) = error;
			}
		}

		/** 输入w和e_ij，计算出能量*/
		double SparseLocalizedForConnectMap::computeEnergy(Eigen::VectorXd Weights, Eigen::VectorXd v_ij)
		{
			/*计算出相应的边长和二面角*/
			Eigen::VectorXd shapeVec = shapeBasis.transpose() * Weights;
			/*求出所有面的标架*/
			//for()
			/*计算出所有的Q矩阵*/
			return 1;
		}

