#pragma once

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


class test_main
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

  test_main(const int& N, const int& Nstep, 
            const std::string& test, const double& L,
            const double& t_end, const double& gamma,
            const double& c, const double& mass)
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

  // test_results
  void test_results(LatticeParticleList& lpl)
  {
    // Test the size of initialized positions, velocities and accelerations
    assert(lpl.get_pos().size() == m_args.m_N);
    assert(lpl.get_vel().size() == m_args.m_N);
    assert(lpl.get_acc().size() == m_args.m_N);

    // Test values of positions, velocities and accelerations 
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
    else if (m_args.m_test == "shift")
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
            assert(fabs(lpl.get_acc()[1][i + cbrtN*(j + cbrtN*k)]) < 1e-13); // y accelerations
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


  void test_init(LatticeParticleList& lpl)
  {
    // Test the size of initialized positions, velocities and accelerations
    assert(lpl.get_pos().size() == m_args.m_N);
    assert(lpl.get_vel().size() == m_args.m_N);
    assert(lpl.get_acc().size() == m_args.m_N);

    // Test values of positions, velocities and accelerations 
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
    else if (m_args.m_test == "shift")
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
    else if (m_args.m_test == "moving")
    {
      // Calculate the velocity
      double lambda_1 = -2 + 2 * cos(M_PI * m_args.m_L / (m_args.m_N+1));
      double omega_1 = sqrt(-m_pc.m_gamma * lambda_1);
      double delta_x = m_args.m_L / (m_args.m_N+1); 
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
  private:
    Args m_args;
    ParticleConst m_pc;
    std::vector<double> mass_vec;

};