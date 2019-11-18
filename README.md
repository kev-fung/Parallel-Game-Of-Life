# Conway's Game of Life Parallelised

A cellular automaton program, https://en.wikipedia.org/wiki/Conway%27s_Game_of_Life, which has been implemented with C++ and the Message Passing Interface (MPI) for computational optimisation. 

<p align="center">
  <img src="./misc/front_img.gif" alt="front_img" width="300">
</p>

Users may specify the domain size, time simulation, and number of computational cores for parallelisation.

## Requirements
* MPI
* C++
* Python

## How to run
1. Set up simulation parameters inside Main.cpp
2. Build files to generate the executable in /Debug
3. From terminal, run mpiexec in ./Debug/Assignment3.exe
4. After running the .exe the animation is generated.
5. To view the animation, run visualisation.ipynb

**NB: Game of life animation folder should be in the same location as the notebook**

## Software Analysis

A generated [example](https://github.com/kev-fung/Parallel-Game-Of-Life/blob/master/misc/output.html).

See code [document](https://github.com/kev-fung/Parallel-Game-Of-Life/blob/master/Code%20Analysis.pdf) for further information.
