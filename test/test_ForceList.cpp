// File to test LatticeParticleForce functionality
#include <iostream>
#include <memory>

#include "GravityForce.hpp"
#include "ForceList.hpp"
#include "SimpleParticleList.hpp"

void test_ForceList()
{
  std::cout << "Testing ForceList ..." << std::endl;
 
  const int N=27; // Number of particles in each direction
  double m=.1;
  SimpleParticleList spl(N, m);

  // Construct a concrete GravityForce
  double g1[3]={0., 0., -9.8};
  std::shared_ptr<GravityForce> gf1 = std::make_shared<GravityForce>(N, g1);
  double g2[3]={9.8, 0., 0.}; // create a different gravity vector
  std::shared_ptr<GravityForce> gf2 = std::make_shared<GravityForce>(N, g2);

  ForceList flist(N);
  flist.push_back(gf1);
  flist.push_back(gf2);
  const Vector3D& vec3d = flist.calc_force(spl);

  // Check that the sum of forces F=mg
  for (int d=0; d < vec3d.dim(); ++d)
    for (int i=0; i < vec3d.size(); ++i)
      assert(vec3d[d][i] == m*(g1[d] + g2[d]));

  std::cout << ".. ForceList tests passed!" << std::endl;
}
