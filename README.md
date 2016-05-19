# README #

### Parallel Programming Languages & Systems Exercise Sheet 2 ###

You will write a C/MPI program which implements a given parallel algorithm.

Since it is not feasible for everyone to have sole access to a real parallel computer, you will be
running your program on a normal DICE workstation, using Unix processes to emulate the
parallelism. This means that you will not get any feel for the real performance of the algorithm.
However, all the tricky conceptual issues of MPI programming will be very apparent to you
and dealing with these is the real purpose of the exercise.

Adaptive Quadrature with a Bag of Tasks
Adaptive Quadrature is a recursive algorithm that computes an approximation of the integral
of a function F(x), using static quadrature rules on adaptively refined sub-intervals of the
integration domain.

### How to run ###
[mymachine]mic: /usr/lib64/openmpi/bin/mpirun -c 5 aquadSolution
Area=7583461.801486