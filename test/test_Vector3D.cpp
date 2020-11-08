// File to test Vector3D
#include <iostream>

#include "Vector3D.hpp"

void test_Vector3D()
{
  std::cout << "Testing Vector3D ..." << std::endl;

  // Construct a Vector3D
  const int N=10;
  Vector3D vec(N);

  // Test Vector3D methods
  assert(vec.size() == N);
  int dim=vec.dim();
  assert(dim == 3);

  // Check vector<double> access for each dimension
  for (int d=0; d < dim; ++d)
  {
    std::vector<double>& dv = vec[d];
    assert(dv.size() == N);
    dv[0] = 1;
    dv[N-1] = 1;
  }

  // Check const 
  const Vector3D& cvec = vec;
  for (int d=0; d < dim; ++d)
  {
    const std::vector<double>& cdv = cvec[d];
    assert(cdv.size() == N);
    double val = cdv[0];
    val = cdv[N-1];
  }

  std::cout << ".. Vector3D tests passed!" << std::endl;
}
