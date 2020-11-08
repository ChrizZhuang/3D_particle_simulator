// File to test LatticeParticleForce functionality
#include <iostream>

#include "LatticeParticleList.hpp"
#include "LatticeParticleForce.hpp"

void test_LatticeParticleForce()
{
  std::cout << "Testing LatticeParticleForce ..." << std::endl;
 
  const int N=3; // Number of particles in each direction
  const int Np=N*N*N; // Note, must be a perfect cube

  // Construct a concrete LatticeParticleForce
  double c=1;
  double a=1./(N-1);
  LatticeParticleForce lpf(Np, c, a);

  // Calculate the equillibrium force on a LatticeParticleList
  double m2=.2;
  std::vector<double> vecm(Np,m2);
  LatticeParticleList lpl(Np, vecm);
  // Set the positions to be on the lattice, spacing=a on [0,1]^3
  std::cout << ".. testing spacing = spring equillibrium!" << std::endl;
  for (int k=0; k < N; ++k)
  for (int j=0; j < N; ++j)
  for (int i=0; i < N; ++i)
  {
    lpl.pos(0,i,j,k) = a*i;
    lpl.pos(1,i,j,k) = a*j;
    lpl.pos(2,i,j,k) = a*k;
  }

  // Check that the resulting forces are round-off
  const Vector3D& f=lpf.calc_force(lpl);
  double tol=1e-15;
  for (int k=0; k < N; ++k)
  for (int j=0; j < N; ++j)
  for (int i=0; i < N; ++i)
  {
    int ix=i+N*(j + N*k);
    // std::cout << "f[0][" << ix << "]=" << f[0][ix] << std::endl;
    assert(fabs(f[0][ix]) < tol);
  }

  // Set the positions to be on the lattice, spacing=a/2 on [0,1/2]^3
  std::cout << ".. testing spacing = 1/2 spring equillibrium!" << std::endl;
  for (int k=0; k < N; ++k)
  for (int j=0; j < N; ++j)
  for (int i=0; i < N; ++i)
  {
    lpl.pos(0,i,j,k) = a*i/2.;
    lpl.pos(1,i,j,k) = a*j/2.;
    lpl.pos(2,i,j,k) = a*k/2.;
  }

  // Check that the resulting forces are -/+a/2 on face centers
  const Vector3D& f2=lpf.calc_force(lpl);
  for (int k=0; k < N; ++k)
  for (int j=0; j < N; ++j)
  for (int i=0; i < N; ++i)
  {
    int ix=i+N*(j + N*k);
    /*
    std::cout << "f2[0][" << ix << "]=" << f2[0][ix] << std::endl;
    std::cout << "f2[1][" << ix << "]=" << f2[1][ix] << std::endl;
    std::cout << "f2[2][" << ix << "]=" << f2[2][ix] << std::endl;
    */
    if (ix==12) // i-1
      assert(fabs(f2[0][ix] + a/2) < tol);
    else if (ix==14) // i+1
      assert(fabs(f2[0][ix] - a/2) < tol);
    else if (ix==10) // j-1
      assert(fabs(f2[1][ix] + a/2) < tol);
    else if (ix==16) // j+1
      assert(fabs(f2[1][ix] - a/2) < tol);
    else if (ix==4) // k-1
      assert(fabs(f2[2][ix] + a/2) < tol);
    else if (ix==22) // k+1
      assert(fabs(f2[2][ix] - a/2) < tol);
    else
      assert(fabs(f2[0][ix]) < tol);
  }

  std::cout << ".. LatticeParticleForce tests passed!" << std::endl;
}
