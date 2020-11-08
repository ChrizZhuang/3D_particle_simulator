#pragma once

#include <cmath>
#include <vector>

#include "ParticleList.hpp"

/**
 * A derived class for particles on a lattice with different masses
 */
class LatticeParticleList : public ParticleList
{
  public:
    /**
     * Constructor: sets particle masses and calls ParticleList constructor
     * @param N - number of particles, must be a perfect cube, N=n^3
     * @param vecm - vector of masses of each particle
     */
    LatticeParticleList(const int& N, const std::vector<double>& vecm)
      : ParticleList(N)
    {
      int n = std::cbrt(N);
      assert(n*n*n==N); // Make sure N is a perfect cube
      m_cbrtN = n;
      assert(vecm.size() == N); // Make sure N masses were passed in
      for (int i=0; i < N; ++i)
        assert(vecm[i] > 0); // Make sure each mass is ok
      m_vec_mass=vecm;
    }

    /**
     * Vector of masses of each particle
     */
    const std::vector<double>& get_vec_mass() const
    {
      return m_vec_mass;
    }

    /**
     * Non-const indexing into pos using the lattice indices
     */
    inline double& pos(const int& d, const int& i, const int& j, const int& k)
    {
      int ix = i + m_cbrtN*(j + m_cbrtN*k); // using ijk ordering
      return m_pos[d][ix];
    }

    /**
     * Non-const indexing into vel using the lattice indices
     */
    inline double& vel(const int& d, const int& i, const int& j, const int& k)
    {
      int ix = i + m_cbrtN*(j + m_cbrtN*k); // using ijk ordering
      return m_vel[d][ix];
    }

    /**
     * Non-const indexing into acc using the lattice indices
     */
    inline double& acc(const int& d, const int& i, const int& j, const int& k)
    {
      int ix = i + m_cbrtN*(j + m_cbrtN*k); // using ijk ordering
      return m_acc[d][ix];
    }


    /**
     * Const indexing into pos using the lattice indices
     */
    inline const double& pos(const int& d, const int& i, const int& j, 
                             const int& k) const
    {
      int ix = i + m_cbrtN*(j + m_cbrtN*k); // using ijk ordering
      return m_pos[d][ix];
    }

    /**
     * Const indexing into vel using the lattice indices
     */
    inline const double& vel(const int& d, const int& i, const int& j, 
                             const int& k) const
    {
      int ix = i + m_cbrtN*(j + m_cbrtN*k); // using ijk ordering
      return m_vel[d][ix];
    }

    /**
     * Const indexing into acc using the lattice indices
     */
    inline const double& acc(const int& d, const int& i, const int& j,
                             const int& k) const
    {
      int ix = i + m_cbrtN*(j + m_cbrtN*k); // using ijk ordering
      return m_acc[d][ix];
    }

    /**
     * Update the acceleration from "a=F/m", assuming constant mass
     * NOTE: pure virtual method from abstract base class ParticleList
     */
    virtual void calc_accel(const Vector3D& force)
    {
      assert(m_N == force.size());
      for (int d=0; d<force.dim(); ++d)
      {
        for (int i=0; i<m_N; ++i)
        {
          m_acc[d][i] = force[d][i]/m_vec_mass[i];
        }
      }
    }

  private:
    std::vector<double> m_vec_mass;
    int m_cbrtN; // size of the cubic lattice in each dim
};

