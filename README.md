# Fock Dimers Limit Shapes

This is a repository accompanying the papers [1](https://arxiv.org/abs/2402.08798), [2](https://www.arxiv.org/abs/2407.19462). It consists of these essential parts:

1. Schottky uniformizations for definition of an M-curve and computation of corresponding differentials and theta functions. This uses the [jtem riemann project](https://www3.math.tu-berlin.de/geometrie/jtem/riemann/).
2. Creation of Z2 and Hex lattices with Fock weights corresponding to a choice of Harnack data.
3. Simulation and visualization of dimer configurations on a weighted graph. Simulation can run on a CUDA enabled GPU or multithreaded on CPU.
4. Calculation of limiting arctic curves and height functions based on Harnack data.

## Recommended installation procedure
1. Use VSCode as environment
2. Install the java plugin in VSCode
3. Install [Maven](https://maven.apache.org/). This is a dependency manager for Java that will help install necessary libraries.
4. Run a script that you want. 
    1. For dimer simulation scripts see the dimerSim folder E.g. [RunSimulations](dimerSim/RunSimulations.java)
    2. For fun with theta functions see e.g. the [Fay Identity verification](inUtil/FayIdentityVerification.java)