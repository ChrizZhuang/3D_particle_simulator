// File to test FrictionForce functionality
#include <iostream>

#include "FrictionForce.hpp"
#include "SimpleParticleList.hpp"

void test_FrictionForce()
{
  std::cout << "Testing FrictionForce ..." << std::endl;
 
  const int N=27; // Number of particles
  double m=.1;
  SimpleParticleList spl(N, m);
  // Initialize the velocity to 1
  Vector3D& vel = spl.get_vel();
  for (int d=0; d < vel.dim(); ++d)
    for (int i=0; i < vel.size(); ++i)
      vel[d][i]=1.;

  // Construct a FrictionForce
  double mu=0.1;
  FrictionForce f(N,mu);

  // Check that the force F= -mu * vel
  const Vector3D& force = f.calc_force(spl);
  for (int d=0; d < force.dim(); ++d)
    for (int i=0; i < force.size(); ++i)
      assert(force[d][i] == -mu*vel[d][i]);

  std::cout << ".. FrictionForce tests passed!" << std::endl;
}
