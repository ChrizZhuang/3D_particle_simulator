#pragma once

#include <cassert>
#include <memory>

#include "Force.hpp"

/**
 * Derived class to represent a sum of forces on particles
 */
class ForceList : public Force
{
public:
  /**
   * Constructor: init the force, Luke
   */
  ForceList(const int& N) : Force(N)
  {
  }
  
  /**
   * Add a force shared pointer to our list of forces
   */
  void push_back(std::shared_ptr<Force> pf)
  {
    assert(pf != nullptr);
    m_forceList.push_back(pf);
  }

  /**
   * Return the sum of all forces given the ParticleList
   */
  virtual const Vector3D& calc_force(const ParticleList& pl)
  {
    // Set our force vectors to zero
    for (int d=0; d < m_force.dim(); ++d)
      for (int i=0; i < m_force.size(); ++i)
        m_force[d][i]=0.;
    
    // Loop through all the forces
    for (int f=0; f < m_forceList.size(); ++f)
    {
      // Update each force for this pl
      m_forceList[f]->calc_force(pl);
      // Accumulate its values in ours
      *this += *m_forceList[f];
    }

    return m_force;
  }
 
protected:
  std::vector<std::shared_ptr<Force> > m_forceList;
};

