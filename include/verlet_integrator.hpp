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
 * Initilize the positions and velocities for each particle with respect to time
*/
class verlet_integrator
{
  public:
 
  struct Args
  {
    int m_N; //  number of interior particles
    int m_Nstep; // number of time steps
    std::string m_test; // the type of test
    double m_L; // length of the cubic simulation box
    double m_t_end; // end of simulation time
  };

  /*
  * Struct to hold constants needed for the particle calculation
  * NOTE: we declare it here "static const" to avoid updating values
  */
  struct ParticleConst
  {  
    double m_gamma;
    double m_c; // spring constant N/m
    double m_mass; // mass per particle, kg
  };

  verlet_integrator(const int& N, const int& Nstep, 
              const std::string& test, const double& L,
              const double& t_end, const double& gamma,
              const double& c, const double& mass,
              LatticeParticleList& lpl)
  {
    m_args.m_L = L;
    m_args.m_N = N; 
    m_args.m_Nstep = Nstep;
    m_args.m_t_end = t_end;
    m_args.m_test = test;
    m_pc.m_c = c;
    m_pc.m_gamma = gamma;
    m_pc.m_mass = mass;
    mass_vec.assign(N, mass);
  }
  
  void init_particles(LatticeParticleList& lpl)
  { 
    int cbrtN = std::cbrt(lpl.get_N()); // calculate the number of particles each side (include the fixed ones)
    double unit_l = m_args.m_L/(cbrtN - 1); // calculate the equil distance with respect to box length L
    int dim_total = 3; // get the total dimensions

    // for equil condition
    if (m_args.m_test == "equil")
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
    else if (m_args.m_test == "shift")
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
    else if (m_args.m_test == "moving")
    {
      // Calculate the velocity
      double lambda_1 = -2 + 2 * cos(M_PI * m_args.m_L / (m_args.m_N+1));
      double omega_1 = sqrt(-m_pc.m_gamma * lambda_1);
      double delta_x = m_args.m_L / (m_args.m_N+1); 
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
  
  /*
 * Function to add force with Vector3D type
 */
Vector3D add_force(const Vector3D f1,
                   const Vector3D f2,
                   const Vector3D f3)
{ 
  int dim_total = 3; // total dimensions
  int cbrtN = std::cbrt(m_args.m_N); // calculate the number of particles each side (include the fixed ones)
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
void calculate_acceleration(LatticeParticleList& lpl,
                            LatticeParticleForce& lpf,
                            FrictionForce& ff,
                            GravityForce& gf)
{
  Vector3D Spring_force = lpf.calc_force(lpl); // calculate the spring force
  Vector3D Friction_force = ff.calc_force(lpl); // calculate the friction
  Vector3D Gravity_force = gf.calc_force(lpl); // calculate the gravity
  Vector3D force = add_force(Spring_force, Friction_force, Gravity_force);// add all force together ('+' is operator overloading)
  lpl.calc_accel(force); // update the acceleration
}

/* 
 * Update x, v, a for every step.
 * partial update for velocity, then update the position
 */
void do_verlet_step(double dt, 
                    LatticeParticleList& lpl,
                    LatticeParticleForce& lpf,
                    FrictionForce& ff,
                    GravityForce& gf)
{
  int dim_total = 3;

  for (int i=0; i < m_args.m_N; ++i) // update just interior points
  {
    for (int dim=0; dim<dim_total; ++dim)
    {
      lpl.get_pos()[dim][i] = lpl.get_pos()[dim][i] + dt*(lpl.get_vel()[dim][i] + 0.5*dt*lpl.get_acc()[dim][i]); // add new position
      //if (i == 4)
      //{
        //  std::cout<<"delta: " << lpl.get_pos()[0][i] << ' ' << lpl.get_pos()[1][i] << ' ' << lpl.get_pos()[2][i] << std::endl;
      //}
      lpl.get_vel()[dim][i] = lpl.get_vel()[dim][i] + 0.5*dt*lpl.get_acc()[dim][i]; // update the partial velocity
    }
  }

  // update the acceleration of interior particles
  calculate_acceleration(lpl, lpf, ff, gf);

  // final update of the velocity with the new acceleration
  for (int i=0; i < m_args.m_N; ++i)
  {
    for (int dim=0; dim<dim_total; ++dim)
    {
      lpl.get_vel()[dim][i] = lpl.get_vel()[dim][i] + 0.5*dt*lpl.get_acc()[dim][i]; // update the velocity
    }
  }

}


// Do the requested number of time steps of the Verlet integrator
void do_time_integration(LatticeParticleList& lpl,
                         LatticeParticleForce& lpf,
                         FrictionForce& ff,
                         GravityForce& gf)
{
  double dt=m_args.m_t_end/m_args.m_Nstep; // end time divided by total time steps
  int cbrtN = std::cbrt(lpl.get_N()); // calculate the number of particles each side (include the fixed ones)
  double unit_l = m_args.m_L/(cbrtN - 1); // calculate the equil distance with respect to box length L
  int dim_total = 3;

  // Calculate initial accelerations
  calculate_acceleration(lpl, lpf, ff, gf);
  
  // Print out some useful information
  std::cout << "Beginning time stepping ..." << std::endl;
  std::cout << "  dt=" << dt << std::endl;
  for (int n=0; n < m_args.m_Nstep; ++n)
  {
    do_verlet_step(dt, lpl, lpf, ff, gf);// update x, v, a for each step
  }
  std::cout << "... completed time step: " << m_args.m_Nstep << std::endl;

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

  private:
    Args m_args;
    ParticleConst m_pc;
    std::vector<double> mass_vec;

};