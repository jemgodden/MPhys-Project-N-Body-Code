# MPhys-Project-N-Body-Code
The N-Body code written to simulate interacting galaxies, for MPhys project.

Fourth generation of an N-Body code for modelling galaxy interactions. Each galaxy has a disk of test particles that do not generate a force on other bodies. 

Force on each body is found and used to give that body a new velocity and position, at each time step. Only calculates the acceleration of each particle once in an attempt to decrease runtime. Did not speed it up.
 
Uses velocity verlet leapfrog method for calculating velocities and positions during a time step.

This code runs a simulation and outputs particle positions at image time step and paths to text file to be plotted.

Interaction options in IntData.

Images plotted by separate function.
