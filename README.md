# HW7 / final project instructions

## Instructions
Goal: Implement velocity Verlet integrator for the three new forces
1. Revise / merge this template code with your hw6
2. Reuse your hw6 solution, refactoring any for your LatticeParticleList
   and force (gravity, friction, neighbor) calculations with the classes
   - **NOTE:** you will now calculate 3D solutions, and the boundary
     lattice particles should have no motion
3. Create tests for all your old and new code in `test`
4. Make sure your code, documentation, comments, etc. are good
5. Determine if the particles settle into "steady-state" for some values

## Grading
- 20% Getting reasonable results - prove it to us in your presentation!
- 20% Well-formatted, good naming, comments, documentation, instructions
- 20% Tests pass for both old and new files/functions/classes
- 20% Class design/encapsulation, avoiding cut-and-paste, good reuse
- 20% Final project presentation 
- (BONUS) 3D plots, extra refactoring of hw6, doxygen

### Sam's 3D visualization
- In `util/3D_visualization.ipynb`

### Build instructions
- In this directory, you can build your main executable with:
```
g++ --std=c++11 -Iinclude src/*.cpp -o ./run.exe; 
./run.exe -n 27 -nstep 10000 -f results/test1.txt -test shift -l 2 -time 2;
```
Please store your result.txt under the diectory results.

- To build and run all the tests:
```
cd test
g++ --std=c++11 -DTEST -I../include *.cpp -o run_tests.exe; 
./run_tests.exe;
```
- To visualize the results
```
cd util;
conda activate [your_python_environment];
jupyter notebook 3D_Visualization.ipynb
```
In the 3rd block of 3D_Visualization.ipynb, please specify the filename of your results in the results directory and the cbrtN which is the number of particle in each dimension(including the fixed ones).
Then run the whole code and you can see the results at the bottom of this file.

### Adding tests
* in `test/run_tests.cpp`, you can add any additional tests
* we have provided more new examples:
  `test_LatticeParticleForce.cpp, etc.`
* the build / run instructions above should then run your additional tests

