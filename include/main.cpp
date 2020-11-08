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

/*
 * Main routine for our simulation
 */
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
  InputArgs args={.N=1, .Nstep=1, .L=1.0, .t_end=1.0}; // set default values

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

/* 
 * Initilize the positions and velocities for each particle with respect to time
*/
void init_particles(const InputArgs& args, 
                    const ParticleConst& pc,
                    LatticeParticleList& lpl)
{ 
  int cbrtN = std::cbrt(lpl.get_N()); // calculate the number of particles each side (include the fixed ones)
  double unit_l = args.L/(cbrtN - 1); // calculate the equil distance with respect to box length L
  int dim_total = 3; // get the total dimensions

  // for equil condition
  if (args.test == "equil")
  {
    for (int k=0; k<cbrtN; ++k) // z level iteration 
    {
      for (int j=0; j<cbrtN; ++j) // y level iteration
      {
        for (int i=0; i<cbrtN; ++i) // x level iteration
        {
          // Initialize the positions
          lpl.get_pos()[0][i + cbrtN*(j + cbrtN*k)] = i*unit_l; // x coordinates
          lpl.get_pos()[1][i + cbrtN*(j + cbrtN*k)] = j*unit_l; // y coordinates
          lpl.get_pos()[2][i + cbrtN*(j + cbrtN*k)] = k*unit_l; // z coordinates

          // Initialize the velocities and accelerations
          for (int dim=0; dim<dim_total; ++dim)
          {
            lpl.get_vel()[dim][i + cbrtN*(j + cbrtN*k)] = 0; // get velocities for 3 dimensions
            lpl.get_acc()[dim][i + cbrtN*(j + cbrtN*k)] = 0; // get accelerations for 3 dimensions
          }
          
        }
      }
    }
  }

  // for shift condition
  else if (args.test == "shift")
  {
    double shift = 0.5 * unit_l;
    // assign the inital positions, velocities and accelerations
    for (int k=0; k<cbrtN; ++k) // z level iteration 
    {
      for (int j=0; j<cbrtN; ++j) // y level iteration
      {
        for (int i=0; i<cbrtN; ++i) // x level iteration
        {
          // Initialize the positions
          lpl.get_pos()[0][i + cbrtN*(j + cbrtN*k)] = i*unit_l; // x coordinates
          lpl.get_pos()[1][i + cbrtN*(j + cbrtN*k)] = j*unit_l; // y coordinates
          lpl.get_pos()[2][i + cbrtN*(j + cbrtN*k)] = k*unit_l; // z coordinates

          // Initialize velocities and accelerations
          for (int dim=0; dim < dim_total; ++dim) // for each dimension
          {
            // Initalize the velocities
            lpl.get_vel()[dim][i + cbrtN*(j + cbrtN*k)] = 0; 
            // Initalize the acclerations
            lpl.get_acc()[dim][i + cbrtN*(j + cbrtN*k)] = 0; 
          }
        }
      }
    }
    // Add extra shift to the x coordinates of interior particles
    if (cbrtN>=3) // if cbrt<3 (0, 1 or 2), no interior particles occur.
    {
      for (int k=1; k<cbrtN-1; ++k) // z level iteration 
      {
        for (int j=1; j<cbrtN-1; ++j) // y level iteration
        {
          for (int i=1; i<cbrtN-1; ++i) // x level iteration
          {
            // add extra shift to the x coordinates of interior particles
            lpl.get_pos()[0][i + cbrtN*(j + cbrtN*k)] += shift; 
          }
        }
      }
    }
  }

  // for moving condition
  else if (args.test == "moving")
  {
    // Calculate the velocity
    double lambda_1 = -2 + 2 * cos(M_PI * args.L / (args.N+1));
    double omega_1 = sqrt(-pc.gamma * lambda_1);
    double delta_x = args.L / (args.N+1); 
    double velocity = 0.5 * delta_x * omega_1;

    // assign the inital positions, velocities and accelerations
    for (int k=0; k<cbrtN; ++k) // z level iteration 
    {
      for (int j=0; j<cbrtN; ++j) // y level iteration
      {
        for (int i=0; i<cbrtN; ++i) // x level iteration
        {
          // Initialize positions
          lpl.get_pos()[0][i + cbrtN*(j + cbrtN*k)] = i*unit_l; // x coordinates
          lpl.get_pos()[1][i + cbrtN*(j + cbrtN*k)] = j*unit_l; // y coordinates
          lpl.get_pos()[2][i + cbrtN*(j + cbrtN*k)] = k*unit_l; // z coordinates

          // Initialize velocities and accelerations
          for (int dim=0; dim < dim_total; ++dim) // for each dimension
          {
            // Initalize the velocities
            lpl.get_vel()[dim][i + cbrtN*(j + cbrtN*k)] = 0; 
            // Initalize the acclerations
            lpl.get_acc()[dim][i + cbrtN*(j + cbrtN*k)] = 0; 
          }
        }
      }
    }
    // Add extra velocities to the interior particles in x direction
    if (cbrtN>=3) // if cbrt<3 (0, 1 or 2), no interior particles occur.
    {
      for (int k=1; k<cbrtN-1; ++k) // z level iteration 
      {
        for (int j=1; j<cbrtN-1; ++j) // y level iteration
        {
          for (int i=1; i<cbrtN-1; ++i) // x level iteration
          {
            // add extra velocities to the interior particles in x direction
            lpl.get_vel()[0][i + cbrtN*(j + cbrtN*k)] += velocity; 
          }
        }
      }
    }

  }
  
  else // print out some useful information to point out the possible error
  {
    std::cout << "Test should be among equil, shift or moving. " << std::endl;
  }
}

void test_init(const InputArgs& args,
               LatticeParticleList& lpl)
{
  // Test the size of initialized positions, velocities and accelerations
  assert(lpl.get_pos().size() == args.N);
  assert(lpl.get_vel().size() == args.N);
  assert(lpl.get_acc().size() == args.N);

  // Test values of positions, velocities and accelerations 
  int cbrtN = std::cbrt(lpl.get_N()); // calculate the number of particles each side (include the fixed ones)
  double unit_l = args.L/(cbrtN - 1); // calculate the equil distance with respect to box length L
  int dim_total = 3; // get the total dimensions

  // for equil condition
  if (args.test == "equil")
  {
    for (int k=0; k<cbrtN; ++k) // z level iteration 
    {
      for (int j=0; j<cbrtN; ++j) // y level iteration
      {
        for (int i=0; i<cbrtN; ++i) // x level iteration
        {
          // Test the positions
          assert(lpl.get_pos()[0][i + cbrtN*(j + cbrtN*k)] == i*unit_l); // x coordinates
          assert(lpl.get_pos()[1][i + cbrtN*(j + cbrtN*k)] == j*unit_l); // y coordinates
          assert(lpl.get_pos()[2][i + cbrtN*(j + cbrtN*k)] == k*unit_l); // z coordinates

          // Test the velocities and accelerations
          for (int dim=0; dim<dim_total; ++dim)
          {
            assert(lpl.get_vel()[dim][i + cbrtN*(j + cbrtN*k)] == 0); 
            assert(lpl.get_acc()[dim][i + cbrtN*(j + cbrtN*k)] == 0); 
          }
          
        }
      }
    }
  }

  // for shift condition
  else if (args.test == "shift")
  {
    double shift = 0.5 * unit_l;
    // Test the inital positions, velocities and accelerations
    for (int k=0; k<cbrtN; ++k) // z level iteration 
    {
      for (int j=0; j<cbrtN; ++j) // y level iteration
      {
        for (int i=0; i<cbrtN; ++i) // x level iteration
        {
          // Test the positions
          // for x coordinates
          if(k == 0 || k == cbrtN-1 || j == 0 || j == cbrtN-1 || i == 0 || i == cbrtN-1) // for outer particles
          {
            assert(lpl.get_pos()[0][i + cbrtN*(j + cbrtN*k)] == i*unit_l); 
          }
          else // for interior particles
          {
            assert(lpl.get_pos()[0][i + cbrtN*(j + cbrtN*k)] == i*unit_l + shift); 
          }
          
          assert(lpl.get_pos()[1][i + cbrtN*(j + cbrtN*k)] == j*unit_l); // y coordinates
          assert(lpl.get_pos()[2][i + cbrtN*(j + cbrtN*k)] == k*unit_l); // z coordinates

          // Test velocities and accelerations
          for (int dim=0; dim < dim_total; ++dim) // for each dimension
          {
            // Test the velocities
            assert(lpl.get_vel()[dim][i + cbrtN*(j + cbrtN*k)] == 0); 
            // Test the acclerations
            assert(lpl.get_acc()[dim][i + cbrtN*(j + cbrtN*k)] == 0); 
          }
        }
      }
    }
  }

  // for moving condition
  else if (args.test == "moving")
  {
    // Calculate the velocity
    double lambda_1 = -2 + 2 * cos(M_PI * args.L / (args.N+1));
    double omega_1 = sqrt(-pc.gamma * lambda_1);
    double delta_x = args.L / (args.N+1); 
    double velocity = 0.5 * delta_x * omega_1;

    // Test the inital positions, velocities and accelerations
    for (int k=0; k<cbrtN; ++k) // z level iteration 
    {
      for (int j=0; j<cbrtN; ++j) // y level iteration
      {
        for (int i=0; i<cbrtN; ++i) // x level iteration
        {
          // Test positions
          assert(lpl.get_pos()[0][i + cbrtN*(j + cbrtN*k)] == i*unit_l); // x coordinates
          assert(lpl.get_pos()[1][i + cbrtN*(j + cbrtN*k)] == j*unit_l); // y coordinates
          assert(lpl.get_pos()[2][i + cbrtN*(j + cbrtN*k)] == k*unit_l); // z coordinates

          // Test velocities and accelerations
          for (int dim=0; dim < dim_total; ++dim) // for each dimension
          {
            // Test the velocities
            if (dim == 0) // for x dimension
            {
              if(k == 0 || k == cbrtN-1 || j == 0 || j == cbrtN-1 || i == 0 || i == cbrtN-1) // for outer particles
              {
                assert(lpl.get_vel()[dim][i + cbrtN*(j + cbrtN*k)] == 0); 
              }
              else // for interior particles
              {
                assert(fabs(lpl.get_vel()[dim][i + cbrtN*(j + cbrtN*k)] == velocity)); 
              }
            }
            else // for other dimensions
            {
              assert(lpl.get_vel()[dim][i + cbrtN*(j + cbrtN*k)] == 0); 
            }
            // Test the acclerations
            assert(lpl.get_acc()[dim][i + cbrtN*(j + cbrtN*k)] == 0); 
          }
        }
      }
    }
  }

}
/*
 * Function to add force with Vector3D type
 */
Vector3D add_force(const Vector3D f1,
                   const Vector3D f2,
                   const Vector3D f3,
                   const InputArgs& args)
{ 
  int dim_total = 3; // total dimensions
  int cbrtN = std::cbrt(args.N); // calculate the number of particles each side (include the fixed ones)
  Vector3D f_total = f1;

  for (int k=0; k<cbrtN; ++k) // z level iteration 
  {
    for (int j=0; j<cbrtN; ++j) // y level iteration
    {
      for (int i=0; i<cbrtN; ++i) // x level iteration
      {
        for (int dim=0; dim<dim_total; ++dim)
        {
          f_total[dim][i + cbrtN*(j + cbrtN*k)] += f2[dim][i + cbrtN*(j + cbrtN*k)];
          f_total[dim][i + cbrtN*(j + cbrtN*k)] += f3[dim][i + cbrtN*(j + cbrtN*k)];
        }
      }
    }
  }
  return f_total;
}
/*
 * Function to calculate acceleration
 */
void calculate_acceleration(const InputArgs& args,
                            const ParticleConst& pc,
                            LatticeParticleList& lpl,
                            LatticeParticleForce& lpf,
                            FrictionForce& ff,
                            GravityForce& gf)
{
  Vector3D Spring_force = lpf.calc_force(lpl); // calculate the spring force
  Vector3D Friction_force = ff.calc_force(lpl); // calculate the friction
  Vector3D Gravity_force = gf.calc_force(lpl); // calculate the gravity
  Vector3D force = add_force(Spring_force, Friction_force, Gravity_force, args); // add all force together ('+' is operator overloading)
  lpl.calc_accel(force); // update the acceleration
}

/* 
 * Update x, v, a for every step.
 * partial update for velocity, then update the position
 */
void do_verlet_step(const InputArgs& args, 
                    const ParticleConst& pc,
                    double dt, 
                    LatticeParticleList& lpl,
                    LatticeParticleForce& lpf,
                    FrictionForce& ff,
                    GravityForce& gf)
{
  int dim_total = 3;

  for (int i=0; i < args.N; ++i) // update just interior points
  {
    for (int dim=0; dim<dim_total; ++dim)
    {
      lpl.get_pos()[dim][i] = lpl.get_pos()[dim][i] + dt*(lpl.get_vel()[dim][i] + 0.5*dt*lpl.get_acc()[dim][i]); // add new position
      lpl.get_vel()[dim][i] = lpl.get_vel()[dim][i] + 0.5*dt*lpl.get_acc()[dim][i]; // update the partial velocity
    }
  }

  // update the acceleration of interior particles
  calculate_acceleration(args, pc, lpl, lpf, ff, gf);

  // final update of the velocity with the new acceleration
  for (int i=0; i < args.N; ++i)
  {
    for (int dim=0; dim<dim_total; ++dim)
    {
      lpl.get_vel()[dim][i] = lpl.get_vel()[dim][i] + 0.5*dt*lpl.get_acc()[dim][i]; // update the velocity
    }
  }

}


// Do the requested number of time steps of the Verlet integrator
void do_time_integration(const InputArgs& args, 
                         const ParticleConst& pc,
                         LatticeParticleList& lpl,
                         LatticeParticleForce& lpf,
                         FrictionForce& ff,
                         GravityForce& gf)
{
  double dt=args.t_end/args.Nstep; // end time divided by total time steps
  int cbrtN = std::cbrt(lpl.get_N()); // calculate the number of particles each side (include the fixed ones)
  double unit_l = args.L/(cbrtN - 1); // calculate the equil distance with respect to box length L
  int dim_total = 3;

  // Calculate initial accelerations
  calculate_acceleration(args, pc, lpl, lpf, ff, gf);
  
  // Print out some useful information
  std::cout << "Beginning time stepping ..." << std::endl;
  std::cout << "  dt=" << dt << std::endl;
  for (int n=0; n < args.Nstep; ++n)
  {
    do_verlet_step(args, pc, dt, lpl, lpf, ff, gf);// update x, v, a for each step
  }
  std::cout << "... completed time step: " << args.Nstep << std::endl;

  // Fixed the boundary particles
  for (int k=0; k<cbrtN; ++k) // z level iteration 
  {
    for (int j=0; j<cbrtN; ++j) // y level iteration
    {
      for (int i=0; i<cbrtN; ++i) // x level iteration
      {
        if (k<1 || j<1 || i<1)
        {
          // Initialize the positions
          lpl.get_pos()[0][i + cbrtN*(j + cbrtN*k)] = i*unit_l; // x coordinates
          lpl.get_pos()[1][i + cbrtN*(j + cbrtN*k)] = j*unit_l; // y coordinates
          lpl.get_pos()[2][i + cbrtN*(j + cbrtN*k)] = k*unit_l; // z coordinates

          // Initialize the velocities and accelerations
          for (int dim=0; dim<dim_total; ++dim)
          {
            lpl.get_vel()[dim][i + cbrtN*(j + cbrtN*k)] = 0; // get velocities for 3 dimensions
            lpl.get_acc()[dim][i + cbrtN*(j + cbrtN*k)] = 0; // get accelerations for 3 dimensions
          }
        }
        else if (k==cbrtN-1 || j==cbrtN-1 || i==cbrtN-1)
        {
          // Initialize the positions
          lpl.get_pos()[0][i + cbrtN*(j + cbrtN*k)] = i*unit_l; // x coordinates
          lpl.get_pos()[1][i + cbrtN*(j + cbrtN*k)] = j*unit_l; // y coordinates
          lpl.get_pos()[2][i + cbrtN*(j + cbrtN*k)] = k*unit_l; // z coordinates

          // Initialize the velocities and accelerations
          for (int dim=0; dim<dim_total; ++dim)
          {
            lpl.get_vel()[dim][i + cbrtN*(j + cbrtN*k)] = 0; // get velocities for 3 dimensions
            lpl.get_acc()[dim][i + cbrtN*(j + cbrtN*k)] = 0; // get accelerations for 3 dimensions
          }
        }

      }
    }
  }
  // Check if the calculated interior points have z coordinates lower than 0
  // If it is true, assign the corresponding particles' positions, velocities and accelerations to 0
  for (int k=1; k<cbrtN-1; ++k) // z level iteration 
  {
    for (int j=1; j<cbrtN-1; ++j) // y level iteration
    {
      for (int i=1; i<cbrtN-1; ++i) // x level iteration
      {
        if (lpl.get_pos()[2][i + cbrtN*(j + cbrtN*k)] < 0)
        {
          lpl.get_pos()[2][i + cbrtN*(j + cbrtN*k)] = 0; // position
          lpl.get_vel()[2][i + cbrtN*(j + cbrtN*k)] = 0; // velocity
          lpl.get_acc()[2][i + cbrtN*(j + cbrtN*k)] = 0; // acceleration
        }
      }
    }
  }
}

void test_results(const InputArgs& args,
                  LatticeParticleList& lpl)
{
  // Test the size of initialized positions, velocities and accelerations
  assert(lpl.get_pos().size() == args.N);
  assert(lpl.get_vel().size() == args.N);
  assert(lpl.get_acc().size() == args.N);

  // Test values of positions, velocities and accelerations 
  int cbrtN = std::cbrt(lpl.get_N()); // calculate the number of particles each side (include the fixed ones)
  double unit_l = args.L/(cbrtN - 1); // calculate the equil distance with respect to box length L
  int dim_total = 3; // get the total dimensions

  // for equil condition
  if (args.test == "equil")
  {
    for (int k=0; k<cbrtN; ++k) // z level iteration 
    {
      for (int j=0; j<cbrtN; ++j) // y level iteration
      {
        for (int i=0; i<cbrtN; ++i) // x level iteration
        {
          // Test the positions
          if(k == 0 || k == cbrtN-1 || j == 0 || j == cbrtN-1 || i == 0 || i == cbrtN-1) // for outer particles
          {
            assert(fabs(lpl.get_pos()[0][i + cbrtN*(j + cbrtN*k)] - i*unit_l) < 1e-14); // x coordinates
            assert(fabs(lpl.get_pos()[1][i + cbrtN*(j + cbrtN*k)] - j*unit_l) < 1e-14); // y coordinates
            assert(fabs(lpl.get_pos()[2][i + cbrtN*(j + cbrtN*k)] - k*unit_l) < 1e-14); // z coordinates
            
            // Test the velocities and accelerations
            for (int dim=0; dim<dim_total; ++dim)
            {
              assert(fabs(lpl.get_vel()[dim][i + cbrtN*(j + cbrtN*k)]) < 1e-14); // get velocities for 3 dimensions
              assert(fabs(lpl.get_acc()[dim][i + cbrtN*(j + cbrtN*k)]) < 1e-14); // get accelerations for 3 dimensions
            }
          }
          
        }
      }
    }
  }

  // for shift condition
  else if (args.test == "shift")
  {
    double shift = 0.5 * unit_l;
    // Test positions, velocities and accelerations
    for (int k=0; k<cbrtN; ++k) // z level iteration 
    {
      for (int j=0; j<cbrtN; ++j) // y level iteration
      {
        for (int i=0; i<cbrtN; ++i) // x level iteration
        {
          // Test the positions
          if(k == 0 || k == cbrtN-1 || j == 0 || j == cbrtN-1 || i == 0 || i == cbrtN-1) // for outer particles
          {
            assert(lpl.get_pos()[0][i + cbrtN*(j + cbrtN*k)] == i*unit_l); // x coordinates
            assert(lpl.get_pos()[1][i + cbrtN*(j + cbrtN*k)] == j*unit_l); // y coordinates
            assert(lpl.get_pos()[2][i + cbrtN*(j + cbrtN*k)] == k*unit_l); // z coordinates
          }

          // Test velocities and accelerations in y direction
          assert(fabs(lpl.get_vel()[1][i + cbrtN*(j + cbrtN*k)]) < 1e-14); // y velocities
          assert(fabs(lpl.get_acc()[1][i + cbrtN*(j + cbrtN*k)]) < 1e-14); // y accelerations
        }
      }
    }
  }

  // for moving condition
  else if (args.test == "moving")
  {
    // Calculate the velocity
    double lambda_1 = -2 + 2 * cos(M_PI * args.L / (args.N+1));
    double omega_1 = sqrt(-pc.gamma * lambda_1);
    double delta_x = args.L / (args.N+1); 
    double velocity = 0.5 * delta_x * omega_1;

    // Test the inital positions, velocities and accelerations
    for (int k=0; k<cbrtN; ++k) // z level iteration 
    {
      for (int j=0; j<cbrtN; ++j) // y level iteration
      {
        for (int i=0; i<cbrtN; ++i) // x level iteration
        {
          // Test the positions
          if(k == 0 || k == cbrtN-1 || j == 0 || j == cbrtN-1 || i == 0 || i == cbrtN-1) // for outer particles
          {
            assert(fabs(lpl.get_pos()[0][i + cbrtN*(j + cbrtN*k)] - i*unit_l) < 1e-14); // x coordinates
            assert(fabs(lpl.get_pos()[1][i + cbrtN*(j + cbrtN*k)] - j*unit_l) < 1e-14); // x coordinates
            assert(fabs(lpl.get_pos()[2][i + cbrtN*(j + cbrtN*k)] - k*unit_l) < 1e-14); // x coordinates

            for (int dim=0; dim<dim_total; ++dim)
            {
              assert(lpl.get_vel()[dim][i + cbrtN*(j + cbrtN*k)] == 0); // velocities
              assert(lpl.get_acc()[dim][i + cbrtN*(j + cbrtN*k)] == 0); // accelerations
            }
          }
        }
      }
    }
  }

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
    myfile << "The results of condition \"" << args.test << "\" at " << args.t_end << " s." << std::endl << std::endl << std::endl; 

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
    myfile << "Positions: " << std::endl << "   ";
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
    
    myfile << std::endl << "Velocities:" << std::endl << "   ";
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

    myfile << std::endl << "Accelerations:" << std::endl << "   ";
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

  // Initialize the particle positions, velocities and accelerations regarding the condition 
  init_particles(args, pc, lpl);

  // Test init
  std::cout << " Testing init..." << std::endl;
  test_init(args, lpl);
  std::cout << " Initialization tests passed!" << std::endl << std::endl;

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
  do_time_integration(args, pc, lpl, lpf, ff, gf);

  // Test integration
  std::cout << std::endl << " Test integration results..." << std::endl;
  test_results(args, lpl);
  std::cout << " Integration tests passed!" << std::endl << std::endl;

  // Print out the integration results
  // std::cout << "Integration results: " << std::endl;
  // print_results(args, lpl);

  // Write output files
  output_final(args, lpl);
  
  return 0;
}


#endif // ndef TEST
