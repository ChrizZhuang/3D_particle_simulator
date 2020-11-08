// File to test GravityForce
#include <iostream>

#include "GravityForce.hpp"
#include "SimpleParticleList.hpp"
#include "LatticeParticleList.hpp"

void test_GravityForce()
{
  std::cout << "Testing GravityForce ..." << std::endl;
 
  const int N=27; // Number of particles, 3 in each direction

  // Check for SimpleParticleList
  {
    std::cout << ".. GravityForce for SimpleParticleList" << std::endl;
    // Construct a concrete GravityForce
    double g1[3]={0., 0., -9.8};
    GravityForce gf1(N, g1);

    double m=.1;
    SimpleParticleList spl(N, m);
    const Vector3D& vec3d = gf1.calc_force(spl);

    // Check that the force F=mg
    for (int d=0; d < vec3d.dim(); ++d)
      for (int i=0; i < vec3d.size(); ++i)
        assert(vec3d[d][i] == m*g1[d]);
  }

  // Check for LatticeParticleList
  {
    std::cout << ".. GravityForce for LatticeParticleList" << std::endl;
    double g2[3]={9.8, 0., 0.}; // create a different gravity vector
    GravityForce gf2(N, g2);
    double m2=.2;
    std::vector<double> vecm(N,m2);
    LatticeParticleList lpl(N, vecm);
    const Vector3D& vec3d = gf2.calc_force(lpl);

    // Check that the force F=mg
    for (int d=0; d < vec3d.dim(); ++d)
      for (int i=0; i < vec3d.size(); ++i)
        assert(vec3d[d][i] == vecm[i]*g2[d]);
  }

  std::cout << ".. GravityForce tests passed!" << std::endl;
}
