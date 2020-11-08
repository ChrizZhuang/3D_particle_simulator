#pragma once

#include "Force.hpp"

/**
 * Class to represent friction force proportional to the velocity
 */
class FrictionForce : public Force
{
public:
  /**
   * Constructor: initialize the 1D spring force model
   * @param d - friction "drag" constant, units "kg/s"
   */
  FrictionForce(const int& N, const double& d) : Force(N)
  {
    assert(d >= 0);
    m_d = d;
  }
  
  /**
   * Return the calculated force based on the ParticleList velocity
   * @overload
   */
  virtual const Vector3D& calc_force(const ParticleList& pl)
  {
    assert(pl.get_N() == m_force.size());

    const Vector3D& vel = pl.get_vel();

    // Calculate the F = -d v
    for (int d=0; d < vel.dim(); ++d)
      for (int i=0; i < vel.size(); ++i)
        m_force[d][i] = -m_d*vel[d][i];

    return m_force;
  }

protected:
  double m_d; // friction "drag" constant, units "kg/s"
};

