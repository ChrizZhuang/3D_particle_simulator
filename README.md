## 3D particle simulator

### Introduction 
This is a software package that can be used to simulate particles' movement within a cube.
Considered forces include:
1. Gravity force
1. Spring force
1. Drag force

### API
Tunable parameters include:
1. number of particles (should be an integer's cubic, int)
1. number of step (number of total step of Verlet integration, int)
1. filename (output filename, recommend to store inside the prebuilt *results* directory, str)
1. test (initial configuration you want to test, including 'equil', 'shift' and 'moving', str)
1. l (length of the cubic side, double)
1. time (total time you want to simulate, double)

*Illustration about 3 test configurations:
**'equil'** means that all the particles are placed uniformly, no change of position should happen.
**'shift'** means that all the particles are shift towards the positive x direction for half of the equil spacing.
**' moving'** means that all the particles have the initial velocity towards the positive x direction.

### Build instructions
- In this directory, you can build your main executable with:
```
g++ --std=c++11 -Iinclude src/*.cpp -o ./run.exe; 
./run.exe -n 27 -nstep 10000 -f results/test1.txt -test shift -l 2 -time 2;
```
Please store your result.txt under the diectory *results*.

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
In the 3rd block of 3D_Visualization.ipynb, please specify the filename of your results in the results directory and cbrtN which is the number of particle in each dimension(including the fixed ones).
Then run the whole code and you can see the results at the bottom of this file.


