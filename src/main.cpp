#ifndef TEST

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

struct InputArgs
{
  int N; // -n <N> (int), number of interior particles
  int Nstep; // -nstep <Nstep> (int), number of time steps
  std::string filename; // -f <filename> (string), output file name
  std::string test; // -test <test> (string), the type of test
  double L; // -l <L> (double), length of the cubic simulation box
  double t_end; // -time (double), end of simulation time.
};

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

/* 
 * Implement a command line argument parser 
 */
InputArgs parse_args(int argc, char* argv[])
{
  InputArgs args={.N=1, .Nstep=1, .test = "equil", .L=1.0, .t_end=1.0}; // set default values

  // Convert argv into vector of strings if there is only the exacutable's name
  if (argc == 1)
    return args;

  // Otherwise parse for arguments we know about
  std::vector<std::string> all_args;
  all_args.assign(argv+1, argv+argc);

  for (auto a=all_args.begin(); a != all_args.end(); ++a)
  {
    auto arg=*a;
    if (arg=="-n")
    {
      ++a; // read the next argument
      auto arg=*a;
      args.N=stoi(arg);
    }
    else if (arg=="-nstep")
    {
      ++a; // read the next argument
      auto arg=*a;
      args.Nstep=stoi(arg);
    }
    else if (arg=="-f")
    {
      ++a; // read the next argument
      auto arg=*a;
      args.filename=arg;
    }
    else if (arg=="-test")
    {
      ++a; // read the next argument
      auto arg=*a;
      args.test=arg;
    }
    else if (arg=="-l")
    {
      ++a; // read the next argument
      auto arg=*a;
      args.L=stoi(arg);
    }
    else if (arg=="-time")
    {
      ++a; // read the next argument
      auto arg=*a;
      args.t_end=stoi(arg);
    }
  }
  return args;
}
void print_results(const InputArgs& args,
                  LatticeParticleList& lpl)
{
  int cbrtN = std::cbrt(args.N); // calculate the number of particles each side (include the fixed ones)
  std::cout << std::endl << "Positions:" << std::endl;
  for (int k=0; k<cbrtN; ++k) // z level iteration 
  {
    std::cout << "The " << k+1 <<"th layer: " << std::endl;
    for (int j=0; j<cbrtN; ++j) // y level iteration
    {
      for (int i=0; i<cbrtN; ++i) // x level iteration
      {
        std::cout<< lpl.get_pos()[0][i + cbrtN*(j + cbrtN*k)] << ' ' << lpl.get_pos()[1][i + cbrtN*(j + cbrtN*k)] << ' ' << lpl.get_pos()[2][i + cbrtN*(j + cbrtN*k)] << std::endl;
      }
    }
    std::cout << std::endl;
  }
  
  std::cout << std::endl << "Velocities:" << std::endl;
  for (int k=0; k<cbrtN; ++k) // z level iteration 
  {
    std::cout << "The " << k+1 <<"th layer: " << std::endl;
    for (int j=0; j<cbrtN; ++j) // y level iteration
    {
      for (int i=0; i<cbrtN; ++i) // x level iteration
      {
        std::cout<< lpl.get_vel()[0][i + cbrtN*(j + cbrtN*k)] << ' ' << lpl.get_vel()[1][i + cbrtN*(j + cbrtN*k)] << ' ' << lpl.get_vel()[2][i + cbrtN*(j + cbrtN*k)] << std::endl;
      }
    }
    std::cout << std::endl;
  }

  std::cout << std::endl << "Accelerations:" << std::endl;
  for (int k=0; k<cbrtN; ++k) // z level iteration 
  {
    std::cout << "The " << k+1 <<"th layer: " << std::endl;
    for (int j=0; j<cbrtN; ++j) // y level iteration
    {
      for (int i=0; i<cbrtN; ++i) // x level iteration
      {
        std::cout<< lpl.get_acc()[0][i + cbrtN*(j + cbrtN*k)] << ' ' << lpl.get_acc()[1][i + cbrtN*(j + cbrtN*k)] << ' ' << lpl.get_acc()[2][i + cbrtN*(j + cbrtN*k)] << std::endl;
      }
    }
    std::cout << std::endl;
  }
}

void output_final(const InputArgs& args, 
                  LatticeParticleList& lpl)
{
  std::ofstream myfile;
  myfile.open(args.filename);
  if (myfile.is_open())
  {
    std::cout << "Opened: " << args.filename << std::endl;

    // print out the title
    myfile << "Results of condition \"" << args.test << "\" at " << args.t_end << " s." << std::endl << std::endl << std::endl; 

    // print out the basic information about the system
    myfile << "Information about this system: " << std::endl << std::endl;
    myfile << "   Number of total particles: " << args.N << ';' << std::endl;
    myfile << "   Number of interior particles: " << std::cbrt(args.N)-2 << ';' << std::endl;
    myfile << "   Spring constant: 0.5 N/m;" << std::endl << std::endl;

    // print out the details of simulation
    myfile << "Information about the simulation: " << std::endl;
    myfile << "   Method: Verlet integration; " << std::endl; 
    myfile << "   Type of test: " << args.test << ';' << std::endl;
    myfile << "   Number of steps: " << args.Nstep << " steps;" << std::endl << std::endl << std::endl;

    // print out the simlation results
    int cbrtN = std::cbrt(lpl.get_N()); // calculate the number of particles each side (include the fixed ones)
    myfile << "Simulation result: " << std::endl << std::endl;
    myfile << "Positions: " << std::endl << std::endl;
    myfile << std::setprecision(8); // set the precision to be 8 digits
    for (int k=0; k<cbrtN; ++k) // z level iteration 
    {
      myfile << "The " << k+1 <<"th layer: " << std::endl;
      for (int j=0; j<cbrtN; ++j) // y level iteration
      {
        for (int i=0; i<cbrtN; ++i) // x level iteration
        {
          myfile << "   " << lpl.get_pos()[0][i + cbrtN*(j + cbrtN*k)] << ' ' << lpl.get_pos()[1][i + cbrtN*(j + cbrtN*k)] << ' ' << lpl.get_pos()[2][i + cbrtN*(j + cbrtN*k)] << std::endl;
        }
      }
      myfile << std::endl;
    }
    
    myfile << std::endl << "Velocities:" << std::endl << std::endl;
    for (int k=0; k<cbrtN; ++k) // z level iteration 
    {
      myfile << "The " << k+1 <<"th layer: " << std::endl;
      for (int j=0; j<cbrtN; ++j) // y level iteration
      {
        for (int i=0; i<cbrtN; ++i) // x level iteration
        {
          myfile << "   " << lpl.get_vel()[0][i + cbrtN*(j + cbrtN*k)] << ' ' << lpl.get_vel()[1][i + cbrtN*(j + cbrtN*k)] << ' ' << lpl.get_vel()[2][i + cbrtN*(j + cbrtN*k)] << std::endl;
        }
      }
      myfile << std::endl;
    }

    myfile << std::endl << "Accelerations:" << std::endl << std::endl;
    for (int k=0; k<cbrtN; ++k) // z level iteration 
    {
      myfile << "The " << k+1 <<"th layer: " << std::endl;
      for (int j=0; j<cbrtN; ++j) // y level iteration
      {
        for (int i=0; i<cbrtN; ++i) // x level iteration
        {
          myfile << "   " << lpl.get_acc()[0][i + cbrtN*(j + cbrtN*k)] << ' ' << lpl.get_acc()[1][i + cbrtN*(j + cbrtN*k)] << ' ' << lpl.get_acc()[2][i + cbrtN*(j + cbrtN*k)] << std::endl;
        }
      }
      myfile << std::endl;
    }
    myfile << std::endl;
    myfile.close(); // close the file to prevent any further change
  }
}

int main(int argc, char* argv[])
{
  // Parse the input arguments
  auto args = parse_args(argc, argv); 

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

  // Print out the command line arguments we parsed
  std::cout << argv[0] << " Input args:" << std::endl;
  std::cout << " N = " << args.N << std::endl; // total number of particles (include the fixed particles)
  std::cout << " Nstep = " << args.Nstep << std::endl; // number of iteration you wish to take
  std::cout << " Filename: " << args.filename << std::endl; // output filename
  std::cout << " Run: " << args.test << std::endl; // wich kind of condition you would like to run
  std::cout << " Length: " << args.L << std::endl; // length of the simulation box
  std::cout << " End time: " <<args.t_end << std::endl << std::endl; // length of time for simulation

  // Convert the mass value to a vector
  
  std::vector<double> mass_vec;
  mass_vec.assign(args.N, pc.mass); 
  // Instantiate the LatticeParticleList
  LatticeParticleList lpl(args.N, mass_vec);
  // The number passed inside the instance lpl include the fixed particles on the boundaries, so for free particle, it is necessary to -2
  std::cout << " Number of total particles in instance lpl: " << lpl.get_N() << std::endl;
  std::cout << " Mass of each particle in instance lpl: " << lpl.get_vec_mass()[0] << std::endl << std::endl;

  // Instantiate verlet_integrator
  verlet_integrator vi(args.N, args.Nstep, args.test, args.L, args.t_end, pc.gamma, pc.c, pc.mass, lpl);
  
  // Initialize the particle positions, velocities and accelerations regarding the condition 
  vi.init_particles(lpl);

  // Print out the initialization results
  // std::cout << "Initialization results: " << std::endl;
  // print_results(args, lpl);
  
  // Instantiate some objects
  double equil_distance = args.L/(std::cbrt(args.N)-1);
  LatticeParticleForce lpf(args.N, pc.c, equil_distance);
  double drag = 0.1;
  FrictionForce ff(args.N, drag);
  double g[3] = {0, 0, -9.8}; // gravity in z direction
  GravityForce gf(args.N, g);

  // Calculate the particle final position, velocity and acceleration
  vi.do_time_integration(lpl, lpf, ff, gf);
  
  // Print out the integration results
  // std::cout << "Integration results: " << std::endl;
  // print_results(args, lpl);

  // Write output files
  output_final(args, lpl);
  
  return 0;
}

#endif // ndef TEST
