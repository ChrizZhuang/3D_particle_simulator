/**
 * Simple test harness that calls test functions
 */

#include <iostream>

// TODO - include any other tests you write here and in main()
void test_Vector3D();
void test_ParticleList();
void test_SimpleParticleList();
void test_LatticeParticleForce();
void test_ForceList();
void test_GravityForce();
void test_FrictionForce();
void test_func();

#ifdef TEST

int main()
{
  // Call all the tests - asserts will fail
  // test some objects
  test_Vector3D();
  test_ParticleList();
  test_SimpleParticleList();
  test_LatticeParticleForce();
  test_GravityForce();
  test_ForceList();
  test_FrictionForce();
  test_func();

  std::cout << "All tests passed!" << std::endl;
  return 0;
}

#endif // ifdef TEST
