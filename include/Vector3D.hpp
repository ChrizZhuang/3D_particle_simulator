#pragma once

#include <vector>
#include <cassert>
#include <cmath>


/*
 * A class to wrap 3 vector<double> of a fixed size
 */
class Vector3D
{
public:
  /*
   * Constructor, initializes vectors to length N, filled with NaN's 
   */
  Vector3D(int N) 
  {
    assert(N >= 0);
    m_N=N;
    for (int d=0; d<s_dim; ++d)
      m_vec[d].resize(m_N, nan(""));
  }

  /*
   * Return the number of vectors, dim=3
   */
  int dim() const
  {
    return s_dim;
  }

  /*
   * Return the length of the vectors (for loops, for example)
   */
  int size() const
  {
    return m_N;
  }

  /*
   * Return the d-index vector (d=0,1,2 for x,y,z)
   * NOTE: do not change the vector length, an assert will fail
   */
  std::vector<double>& operator[](const int& dim)
  {
    assert((dim >= 0) && (dim < s_dim));
    // NOTE: this is a simple but expensive way to check vector length
    for (int d=0; d<s_dim; ++d)
      assert(m_vec[dim].size() == m_N);
    return m_vec[dim];
  }
    
  /*
   * Const version of operator[] above, cannot modify vectors
   */
  const std::vector<double>& operator[](const int& dim) const
  {
    assert((dim >= 0) && (dim < s_dim));
    // NOTE: this is a simple but expensive way to check vector length
    for (int d=0; d<s_dim; ++d)
      assert(m_vec[dim].size() == m_N);
    return m_vec[dim];
  }

 protected:
  const static int s_dim=3;
  int m_N=-1;
  std::vector<double> m_vec[s_dim]; // x, y, z components

};
