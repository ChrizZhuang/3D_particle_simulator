#pragma once

#include "Force.hpp"
#include "SimpleParticleList.hpp"
#include "LatticeParticleList.hpp"

/**
 * Class to represent gravity forces as a vector, g.
 * Force is calculated using particle mass
 */
class GravityForce : public Force
{
public:
  /**
   * Constructor: initialize the gravity force vector
   * @param g - 3-dimensional gravity vector
   */
  GravityForce(const int& N, const double (&g)[3]) : Force(N)
  {
    for (int d=0; d < 3; ++d)
      m_g[d] = g[d];
  }
  
  /**
   * Return the calculated force based on the ParticleList masses
   * @overload
   */
  virtual const Vector3D& calc_force(const ParticleList& pl)
  {
    assert(pl.get_N() == m_force.size());

    // Check if we have a simple particle list, just use it
    const SimpleParticleList* spl = 
      dynamic_cast<const SimpleParticleList*>(&pl);
    if (spl != nullptr) 
    {
      double mass = spl->get_mass();
      // Calculate the gravity force using a single particle mass
      for (int d=0; d < m_force.dim(); ++d)
        for (int i=0; i < m_force.size(); ++i)
          m_force[d][i] = mass*m_g[d];

      return m_force;
    }

    // Check if we have a lattice particle list, just use it
    const LatticeParticleList* lpl = 
      dynamic_cast<const LatticeParticleList*>(&pl);
    if (lpl != nullptr) 
    {
      const std::vector<double>& vecm = lpl->get_vec_mass();
      assert(vecm.size() == m_force.size());
      // Calculate the gravity force using a single particle mass
      for (int d=0; d < m_force.dim(); ++d)
        for (int i=0; i < m_force.size(); ++i)
          m_force[d][i] = vecm[i]*m_g[d];

      return m_force;
    }

    // Otherwise don't know how to calculate gravity force, return 0
    for (int d=0; d < m_force.dim(); ++d)
      for (int i=0; i < m_force.size(); ++i)
        m_force[d][i]=0.;

    return m_force;
  }

protected:
  double m_g[3]; // gravity vector, units "m/s2"
};

