#ifndef HELPER_H
#define HELPER_H


#include <iostream>
#include <matrix.h>
#include <vector.h>

using namespace ASC_bla;

template<typename T>
Matrix<T> ToMatrix(VectorView<T> vec) {
  Matrix<T> res(3, 3);
  res.Row(0) = vec.segment(0, 3);  // First row
  res.Row(1) = vec.segment(3, 3);  // Second row
  res.Row(2) = vec.segment(6, 3);  // Third row

  return res;
}

template<typename T>
Vector<T> ToVector(MatrixView<T> m) {
  Vector<T> res(9);
  res.setConstant(0);

  res.segment(0, 3) = m.Row(0);
  res.segment(3, 3) = m.Row(1);
  res.segment(6, 3) = m.Row(2);

  return res;
}

template<typename T, typename S>
auto ApplyTransformation(const VectorView<T> q, const VectorView<S> pos) {
  Matrix<decltype(q(0)* pos(0))> R = ToMatrix(q.segment(3, 9));
  Vector<decltype(q(0)* pos(0))> res = q.segment(0, 3) + R * pos;
  return res;
}
// Function to convert a 3D vector to a skew-symmetric matrix (hat operator)
template<typename T>
Matrix<T> vectorToSkewSymmetric(const VectorView<T> v) {
    Matrix<T> skewSymmetric(3, 3);
    skewSymmetric(0, 0) = 0;
    skewSymmetric(0, 1) = -v(2);
    skewSymmetric(0, 2) = v(1);
    skewSymmetric(1, 0) = v(2);
    skewSymmetric(1, 1) = 0;
    skewSymmetric(1, 2) = -v(0);
    skewSymmetric(2, 0) = -v(1);
    skewSymmetric(2, 1) = v(0);
    skewSymmetric(2, 2) = 0;
    
    return skewSymmetric;
}

// Function to convert a skew-symmetric matrix to a 3D vector (inverse hat, vee operator)
template<typename T>
auto skewSymmetricToVector(const MatrixExpr<T>& skewMatrix) -> Vector<decltype(skewMatrix(0, 0)*1.0)> {
    Vector<decltype(1.0*skewMatrix(0, 0))> res = {skewMatrix(0, 0)*1.0, skewMatrix(0, 0)*1.0, skewMatrix(0, 0)*1.0};

    res(0) = (skewMatrix(2, 1) - skewMatrix(1, 2))/2;
    res(1) = (skewMatrix(0, 2) - skewMatrix(2, 0))/2;
    res(2) = (skewMatrix(1, 0) - skewMatrix(0, 1))/2;

    
    return res;
}

size_t dim_per_body();

size_t dim_per_motion();



#endif