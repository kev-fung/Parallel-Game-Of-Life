# Conway's Game of Life Parallelised

A cellular automaton program, https://en.wikipedia.org/wiki/Conway%27s_Game_of_Life, which has been implemented with C++ and the Message Passing Interface (MPI) for computational optimisation. Users may specify the domain size, time simulation, and number of computational cores for parallelisation.

<p align="center">
  <img src="./front_img.png" alt="front_img" width="300">
</p>

## How to run
1. Set up simulation parameters inside Main.cpp
2. Build files to generate the executable in /Debug
3. From terminal, run mpiexec in ./Debug/Assignment3.exe
4. If you haven't already you must install mpi!
5. After running the .exe the animation is generated.
6. To view the animation, run visualisation.ipynb

**NB: Game of life animation folder should be in the same location as the notebook!**

## Software Analysis
See code [document](https://github.com/kev-fung/Parallel-Game-Of-Life/blob/master/Code%20Analysis.pdf) for further information.
