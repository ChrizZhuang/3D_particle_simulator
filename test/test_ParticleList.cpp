// File to test ParticleList functionality
#include <iostream>

#include "SimpleParticleList.hpp"

void test_ParticleList()
{
  std::cout << "Testing ParticleList ..." << std::endl;
 
  const int N=10;
  const double m=.1;
  // Construct a concrete SimpleParticleList to test ParticleList
  SimpleParticleList spl(N, m);

  // Test ParticleList methods (not SimpleParticleList!)
  assert(spl.ParticleList::get_pos().size() == N);
  int dim=spl.ParticleList::get_pos().dim();
  assert(dim == 3);

  // Check vector<double> access for each dimension
  for (int d=0; d < dim; ++d)
  {
    std::vector<double>& dv = spl.ParticleList::get_pos()[d];
    assert(dv.size() == N);
  }

  // Check const 
  const ParticleList& pl = spl;
  for (int d=0; d < dim; ++d)
  {
    const std::vector<double>& cdv = pl.ParticleList::get_pos()[d];
    assert(cdv.size() == N);
  }

  std::cout << ".. ParticleList tests passed!" << std::endl;
}
