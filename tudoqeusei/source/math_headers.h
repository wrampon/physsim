#ifndef _MATH_HEADERS_H_
#define _MATH_HEADERS_H_

//used for conversion
#include "limits.h"
#include "global_headers.h"

// glm
#include "glm.hpp"
#include "gtc\matrix_transform.hpp"
#include "gtc\quaternion.hpp"

//Eigen
#include <Core>
#include <Dense>
#include <Sparse>

#define minimun(a, b) a>b? b : a
#define maximun(a, b) a>b? a : b

//print-me
#define printme(e) {std::cout<<#e<<"="<<(e)<<std::endl;}

typedef int IndexType;
typedef Eigen::Matrix<ScalarType, 3, 3, 0, 3, 3> EigenMatrix3;
typedef Eigen::Matrix<ScalarType, 3, 1, 0, 3, 1> EigenVector3;
typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> VectorX;
typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> Matrix;
typedef Eigen::SparseMatrix<ScalarType> SparseMatrix;
typedef Eigen::Triplet<ScalarType, IndexType> SparseMatrixTriplet;

// eigen vector accessor
#define block_vector(a) block<3,1>(3*(a), 0)
#define block_3vector(a) block<1,3>(a, 0)

// eigen 2 glm, glm 2 eigen
glm::vec3 Eigen2GLM(const EigenVector3& eigen_vector);
EigenVector3 GLM2Eigen(const glm::vec3& glm_vector);

#endif