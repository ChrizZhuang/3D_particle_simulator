#pragma once

#include <vector>

#include "Vector3D.hpp"

/**
 * Abstract base class for managing particle positions & velocities
 * NOTE: must implement the calc_accel() method in concrete derived class
 */
class ParticleList
{
  public:
    /**
     * Constructor for position, velocity, and acceleration, N particles
     */
    ParticleList(const int& N)
      : m_pos(N), m_vel(N), m_acc(N), m_N(N) { }

    virtual int get_N() const
    {
      return m_N;
    }

    /**
     * Return position vector
     */
    virtual Vector3D& get_pos()
    {
      return m_pos;
    }

    /**
     * Return velocity vector
     */
    virtual Vector3D& get_vel()
    {
      return m_vel;
    }

    /**
     * Return acceleration vector
     */
    virtual Vector3D& get_acc()
    {
      return m_acc;
    }

    /**
     * Const position vector
     */
    virtual const Vector3D& get_pos() const
    {
      return m_pos;
    }

    /**
     * Const velocity vector
     */
    virtual const Vector3D& get_vel() const
    {
      return m_vel;
    }

    /**
     * Const acceleration vector
     */
    virtual const Vector3D& get_acc() const
    {
      return m_acc;
    }
    
    /*
     * Update the acceleration from a given force
     * NOTE: pure virtual function, thus an abstract base class
     */
    virtual void calc_accel(const Vector3D& force) = 0;

  protected:
    Vector3D m_pos;
    Vector3D m_vel;
    Vector3D m_acc;
    int m_N=0;
};

