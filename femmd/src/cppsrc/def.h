#ifndef __DEF_H__
#define __DEF_H__

#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#define EIGEN_RUNTIME_NO_MALLOC
#include "Eigen/Dense"
#include "Eigen/StdVector"

//using namespace Eigen;

#define USE_DOUBLE
//#define USE_3D
#define USE_2D

#ifdef USE_DOUBLE
typedef double real;
#endif

#ifdef USE_2D
#define NDIM 2
typedef Eigen::Vector2d VectorR;
typedef Eigen::Vector2i VectorI;
typedef Eigen::Vector3i elementIndexes;
typedef Eigen::Matrix2d MatrixDxD;

#else
#define NDIM 3
typedef Eigen::Vector3d VectorR;
typedef Eigen::Vector3i VectorI;
typedef Eigen::Vector4i elementIndexes;
typedef Eigen::Matrix3d MatrixDxD;

#endif

extern real PI;
extern real one_over_ndim_factorial;

#define RANDOM01 ((real) rand() / (RAND_MAX))

enum PBCTYPE {NOPBC, XPBC, XYPBC, XYZPBC};

#define Max(x1, x2)  (((x1) > (x2)) ? (x1) : (x2))
#define Min(x1, x2)  (((x1) < (x2)) ? (x1) : (x2))
	

typedef VectorR Box;

#define hvector std::vector

#define minImage(dr, box, x)                           \
	if (dr[x] >= 0.5*box[x]) dr[x] -= box[x];              \
		else if (dr[x] < -0.5*box[x]) dr[x] += box[x];

#ifdef __CUDACC__
__host__ __device__ __forceinline__
#endif 
static void nearestImage(VectorR &dr, VectorR &box, PBCTYPE pbcType){
	switch (pbcType){
	case XYZPBC:
		minImage(dr, box, 0);
		minImage(dr, box, 1);
		minImage(dr, box, 2);
		break;
	case XYPBC:
		minImage(dr, box, 0);
		minImage(dr, box, 1);
		break;
	case XPBC:
		minImage(dr, box, 0);
		break;
	case NOPBC:
		break;
	}
};

#define pbcCalc(r, box, x) r[x] = r[x] - floor(r[x] / box[x])*box[x];

static void applyBoundaryCondition(VectorR &r, VectorR &box, PBCTYPE pbcType){
	switch (pbcType){
	case XYZPBC:
		pbcCalc(r, box, 0);
		pbcCalc(r, box, 1);
		pbcCalc(r, box, 2);
		break;
	case XYPBC:
		pbcCalc(r, box, 0);
		pbcCalc(r, box, 1);
		break;
	case XPBC:
		pbcCalc(r, box, 0);
		break;
	case NOPBC:
		break;
	}
}

//class ForceOnWalls { public: real x0, x1, y0, y1; };

inline real nint(real a){
	if (a >= 0.0) return floor(a + 0.5);
	else return floor(a - 0.5);
}

struct matNxN {
	real *data;
	real **m;
	void init(int n1) {
		n = n1;
		data = new real[n*n];
		m = new real*[n];
		m[0] = data;
		for (int i = 1; i < n; i++) m[i] = m[i - 1] + n;//m[i] = &data[i*n];
	}

	void clear() {
		delete[] data;
		delete[] m;
		n = 0;
	}
	real at(int i, int j) { return data[i*n + j]; }
	void set(int i, int j, real v) { data[i*n + j] = v; }
	real* operator[] (int i) { return m[i]; }
	int n;
};


//struct mat3x3 {
//	real m[3][3];
//};

//struct tetraIndexes{
//struct elementIndexes{
//	int ind[3];
//	inline int& operator[](int idx) { return ind[idx]; }
//};
/*


inline void a_dot_b(mat3x3 a, mat3x3 b, mat3x3 *c) {
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			c->m[i][j] = 0;
			for (int k = 0; k < 3; k++) {
				c->m[i][j] = c->m[i][j] + a.m[i][k] * b.m[k][j];
			}
		}
	}
}

//bool invertMatrix(const double m[16], double invOut[16])

static bool invertMatrix(real *m, real *invOut)
{
	real inv[9], det;
	int i;

	det = -m[2]*m[4]*m[6] + m[1]*m[5]*m[6] + m[2]*m[3]*m[7] - m[0]*m[5]*m[7] - m[1]*m[3]*m[8] + m[0]*m[4]*m[8];
	if (det == 0)
		return false;

	inv[0] = -m[5] * m[7] + m[4] * m[8];
	inv[1] =  m[2] * m[7] - m[1] * m[8];
	inv[2] = -m[2] * m[4] + m[1] * m[5];
	inv[3] =  m[5] * m[6] - m[3] * m[8];
	inv[4] = -m[2] * m[6] + m[0] * m[8];
	inv[5] =  m[2] * m[3] - m[0] * m[5];
	inv[6] = -m[4] * m[6] + m[3] * m[7];
	inv[7] =  m[1] * m[6] - m[0] * m[7];
	inv[8] = -m[1] * m[3] + m[0] * m[4];

	det = 1.0f / det;

	for (i = 0; i < 9; i++)
		invOut[i] = inv[i] * det;

	return true;
}
*/

#endif