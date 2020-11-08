#pragma once

#include "Force.hpp"

/**
 * Class to represent spring forces between particles on the x axis.
 * Forces are calculated only between neighbors, end point forces are zero.
 */
class Spring1DForce : public Force
{
public:
  /**
   * Constructor: initialize the 1D spring force model
   * @param c - spring constant, units "N/m"
   * @param a - spring equilibrium distance, units "m"
   */
  Spring1DForce(const int& N, const double& c, const double& a) : Force(N)
  {
    assert(c > 0);
    assert(a > 0);
    m_c = c;
    m_a = a;
  }
  
  /**
   * Return the calculated force based on the ParticleList
   * @overload
   */
  virtual const Vector3D& calc_force(const ParticleList& pl)
  {
    // TODO - implement the spring 1D force from HW5
    int N = pl.get_pos().size();
    return m_force;
  }

protected:
  double m_c; // spring constant, units "N/m"
  double m_a; // spring equilibrium distance, units "m"
};

