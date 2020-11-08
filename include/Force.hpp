#pragma once

#include <cassert>

#include "ParticleList.hpp"

/**
 * Abstract base class to represent forces on particles
 */
class Force
{
public:
  /**
   * Constructor: init the force, Luke
   */
  Force(const int& N) : m_force(N)
  {
  }
  
  /**
   * Return the calculated force based on the ParticleList
   * NOTE: pure virtual function, derived classes must implement
   */
  virtual const Vector3D& calc_force(const ParticleList&) = 0;

  /**
   * Update to be the sum of this and another force over all particles
   */
  virtual Force& operator+=(const Force& f)
  {
    // Make sure the Vector3D's are compatible
    const int dim=m_force.dim();
    assert(dim==f.m_force.dim());
    const int size=m_force.size();
    assert(size==f.m_force.size());

    // Add each of the components of all the vectors in m_force's
    for (int d=0; d < dim; ++d)
      for (int i=0; i < size; ++i)
        m_force[d][i] += f.m_force[d][i];

    return *this;
  }
 
protected:
  Vector3D m_force;
};

