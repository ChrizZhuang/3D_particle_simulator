// File to test main functions
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <cmath>
#include <assert.h>
#include <string>

#include "Vector3D.hpp"
#include "SimpleParticleList.hpp"
#include "Spring1DForce.hpp"
#include "LatticeParticleList.hpp"
#include "LatticeParticleForce.hpp"
#include "GravityForce.hpp"
#include "FrictionForce.hpp"
#include "verlet_integrator.hpp"
#include "test_main.hpp"

void test_func()
{
  struct Args
  {
    int N = 64; // -n <N> (int), number of interior particles
    int Nstep = 1e4; // -nstep <Nstep> (int), number of time steps
    std::string test; // -test <test> (string), the type of test
    double L = 2; // -l <L> (double), length of the cubic simulation box
    double t_end = 2; // -time (double), end of simulation time.
  };

  Args args;
  args.test = "equil";

  // args.test = "shift";

  // args.test = "moving";

  /*
  * Struct to hold constants needed for the particle calculation
  * NOTE: we declare it here "static const" to avoid updating values
  */
  struct ParticleConst
  {
    double gamma = 5;
    double c=0.5; // spring constant N/m
    double mass = 0.1; // mass per particle, kg
  };
  static const ParticleConst pc;

  // Defensive programming: check that the input N is an integer to the third power
  int index = 0; // initiate a index to mark when the input N is not an integer to the third power
  for (int i=0; i<=args.N; i++)
  {
    if (args.N == pow(i, 3))
    {
      index = 1; // change the index if the input N is an integer to the third power
    }
  }

  if (index == 0)
  {
    std::cout << "error: The input N is not an integer to the third power!" << std::endl; // print out message to indicate the possible mistakes
  }
  assert(index == 1);

  // Convert the mass value to a vector
  std::vector<double> mass_vec;
  mass_vec.assign(args.N, pc.mass); 
  // Instantiate the LatticeParticleList
  LatticeParticleList lpl(args.N, mass_vec);

  // Instantiate verlet_integrator
  verlet_integrator vi(args.N, args.Nstep, args.test, args.L, args.t_end, pc.gamma, pc.c, pc.mass, lpl);
  
  // Initialize the particle positions, velocities and accelerations regarding the condition 
  vi.init_particles(lpl);
  
  // Test init
  test_main tm(args.N, args.Nstep, args.test, args.L, args.t_end, pc.gamma, pc.c, pc.mass);
  tm.test_init(lpl);
  
  // Instantiate some objects
  double equil_distance = args.L/(std::cbrt(args.N)-1);
  LatticeParticleForce lpf(args.N, pc.c, equil_distance);
  double drag = 0.1;
  FrictionForce ff(args.N, drag);
  double g[3] = {0, 0, -9.8}; // gravity in z direction
  GravityForce gf(args.N, g);

  // Calculate the particle final position, velocity and acceleration
  vi.do_time_integration(lpl, lpf, ff, gf);

  // Test integration
  tm.test_results(lpl);
  
  std::cout << ".. Main tests passed!" << std::endl;
}
