// File to test SimpleParticleList functionality
#include <iostream>

#include "SimpleParticleList.hpp"

void test_SimpleParticleList()
{
  std::cout << "Testing SimpleParticleList ..." << std::endl;
 
  // Test SimpleParticleList constructor
  const int N=10;
  const double m=.1;
  SimpleParticleList spl(N, m);

  // Test SimpleParticleList methods
  assert(spl.get_pos().size() == N);
  int dim=spl.get_pos().dim();
  assert(dim == 3);
  assert(spl.get_mass() == m);

  // Check vector<double> access for each dimension
  for (int d=0; d < dim; ++d)
  {
    std::vector<double>& dv = spl.get_pos()[d];
    assert(dv.size() == N);
  }

  // Check const versions
  const SimpleParticleList& cspl = spl;
  assert(cspl.get_mass() == m);
  for (int d=0; d < dim; ++d)
  {
    const std::vector<double>& cdv = cspl.get_pos()[d];
    assert(cdv.size() == N);
  }

  std::cout << ".. SimpleParticleList tests passed!" << std::endl;
}
