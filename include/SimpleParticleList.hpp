#pragma once

#include "ParticleList.hpp"

/**
 * A derived class for a simple constant mass particle list
 */
class SimpleParticleList : public ParticleList
{
  public:
    /**
     * Constructor: sets particle mass and calls ParticleList constructor
     */
    SimpleParticleList(const int& N, const double& m)
      : ParticleList(N)
    {
      assert(m > 0);
      m_mass=m;
    }

    /**
     * Mass of each particle
     */
    double get_mass() const
    {
      return m_mass;
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
          m_acc[d][i] = force[d][i]/m_mass;
        }
      }
    }

  private:
    double m_mass=-1;
};

