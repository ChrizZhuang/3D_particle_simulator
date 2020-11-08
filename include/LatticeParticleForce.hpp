#pragma once

#include <cmath>
#include <cassert>

#include "Force.hpp"
// #include "LatticeParticleList.hpp"

/**
 * Class to represent forces between particles on a cubic lattice.
 * Forces are calculated only between neighbors, end point forces are zero.
 */
class LatticeParticleForce : public Force
{
public:
  /**
   * Constructor: initialize the lattice particle force model
   * @param N - the number of particles (should be a perfect cube)
   * @param c - spring constant, units "N/m"
   * @param a - spring equilibrium distance, units "m"
   */
  LatticeParticleForce(const int& N, const double& c, const double& a)
    : Force(N)
  {
    // Check that N is a perfect integer cube
    m_cbrtN = std::cbrt(N);
    assert(m_cbrtN*m_cbrtN*m_cbrtN==N);
    assert(c >= 0);
    assert(a >= 0);
    m_c = c;
    m_a = a;
  }
  
  /**
   * Return the calculated forces based on a LatticeParticleList
   * @overload
   * @param pl - ref to ParticleList, but must be a LatticeParticleList
   */
  virtual const Vector3D& calc_force(const ParticleList& pl)
  {
    // Make sure this particle list is the same size
    int N = pl.get_pos().size();
    assert(N == m_force.size());

    // Do a dynamic cast to LatticeParticleList, since that's the 
    // only kind of particle list we can calculate a force for
    const LatticeParticleList* ptr_lpl = 
      dynamic_cast<const LatticeParticleList*>(&pl);
    assert(ptr_lpl != nullptr);
    const LatticeParticleList& lpl=*ptr_lpl;

    // Initialize all force values to 0 so we can accumulate them
    for (int d=0; d < m_force.dim(); ++d)
    for (int ix=0; ix < N; ++ix)
      m_force[d][ix]=0;

    // Loop over interior particles, calculate the forces 
    // for each particle with its neighbors in the lattice
    for (int k=1; k < m_cbrtN-1; ++k)
    for (int j=1; j < m_cbrtN-1; ++j)
    for (int i=1; i < m_cbrtN-1; ++i)
    {
      int ix=i + m_cbrtN*(j + m_cbrtN*k); // index using ijk ordering
      // Position of current particle
      double x1[3]={lpl.pos(0,i,j,k), lpl.pos(1,i,j,k), lpl.pos(2,i,j,k)}; 
      // std::cout << "x1[0][" << ix << "]=" << x1[0] << std::endl;

      // Neighbor force for i+/-1, j+/-1, k+/-1 neighbors
      const int ixoff[6][3]={ {-1,0,0},{1,0,0}, // i offsets
                              {0,-1,0},{0,1,0}, // j offsets
                              {0,0,-1},{0,0,1}}; // k offsets
      // Loop over neighbor offsets
      for (int ixf=0; ixf<6; ++ixf)
      {
        int i2=i+ixoff[ixf][0];
        int j2=j+ixoff[ixf][1];
        int k2=k+ixoff[ixf][2];
        int ix2 = i2 + m_cbrtN*(j2 + m_cbrtN*k2); // index of p2
        // Position of particle 2
        double x2[3]={lpl.pos(0,i2,j2,k2),
                      lpl.pos(1,i2,j2,k2),
                      lpl.pos(2,i2,j2,k2)};
        // std::cout << "x2[0][" << ix2 << "]=" << x2[0] << std::endl;
        // Difference x2 - x1
        double dx[3]={x2[0]-x1[0], x2[1]-x1[1], x2[2]-x1[2]};
        // std::cout << "dx[0][" << ix2 << "]=" << dx[0] << std::endl;
        // std::cout << "dx[1][" << ix2 << "]=" << dx[1] << std::endl;
        // std::cout << "dx[2][" << ix2 << "]=" << dx[2] << std::endl;
        // Magnitude |x2 - x1|
        double magDx=std::sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
        // std::cout << "magDx=" << magDx << std::endl;
        // Finally calculate the force on particle 2 for d = x,y,z component
        for (int d=0; d<3; ++d)
        {
          double f12=-m_c*(dx[d] - m_a*dx[d] / magDx);
          // NOTE: only add the force to particle 2 if exterior,
          // otherwise we double-add the force to interior particles!
          if ((i2==0) || (i2==m_cbrtN-1) ||
              (j2==0) || (j2==m_cbrtN-1) ||
              (k2==0) || (k2==m_cbrtN-1))
          {
            m_force[d][ix2] += f12;
          }
          m_force[d][ix] += -f12; // opposite force on this particle
        }
      }
    }

    return m_force;
  }

protected:
  int m_cbrtN; /// number of particles in each direction, cube root of N
  double m_c; /// spring constant, units "N/m"
  double m_a; /// spring equilibrium distance, units "m"
};

