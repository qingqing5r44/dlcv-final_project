//stl
#include <iostream> 
#include <fstream> 
#include <string>
#include <vector>
#include <memory>
#include <cmath>
#include <algorithm>
#include <map>
#include <set>
#include <stdio.h> 

#include <time.h> 
#include <omp.h>


//IGL

//#include <igl/unproject_onto_mesh.h>
////#include <igl/viewer/Viewer.h>
//#include <igl/embree/unproject_onto_mesh.h>
//#include <igl/opengl2/draw_rectangular_marquee.h>
//#include <igl/opengl2/unproject.h>
//#include <igl/arap.h>
//#include <igl/biharmonic_coordinates.h>
//#include <igl/cat.h>
//#include <igl/cotmatrix.h>
//#include <igl/massmatrix.h>
//#include <igl/matrix_to_list.h>
//#include <igl/parula.h>
//#include <igl/point_mesh_squared_distance.h>
//#include <igl/remove_unreferenced.h>
//#include <igl/slice.h>

//openmesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include "OpenMesh/Core/Mesh/Handles.hh"

//cuBLAS
//#include "cuda_runtime.h"
//#include <cusparse.h>
//#include "cublas_v2.h"
//#include <helper_functions.h>  // helper for shared functions common to CUDA Samples
//#include <helper_cuda.h>       // helper for CUDA error checking

//eigen
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/SparseCholesky>
#include <Eigen/SVD>

//glut
//#include <GL/glut.h>
//#include <Windows.h>

//MATLAB
//#include "mclmcr.h"
//#include "matrix.h"
//#include "mclcppclass.h"
//#include"lars.h"
//
//#pragma comment(lib,"lars.lib")

typedef double decimal;
typedef int integer;

typedef std::vector<integer> int_vector;
typedef std::vector<decimal> dec_vector;
typedef std::map<integer, integer> int_map;
//

typedef OpenMesh::TriMesh_ArrayKernelT<> MyMesh;
//typedef OpenMesh::PolyMesh_ArrayKernelT<> MyMesh;
typedef std::vector<MyMesh::Point> pt_vector;

typedef Eigen::SparseMatrix<decimal> SMatXd;
typedef Eigen::Triplet<decimal> T;

#define M_PI       3.14159265358979323846
