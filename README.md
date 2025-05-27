## Overview
The purpose of this repo is to test a parallel implementation of the computation of the algorithmic tangent operator as a step in the finite element method. The parallel computation will be performed in CUDA. 

## Finite Element Program

To set up the calculations for the tangent modulus, a finite element program was written for the simple problem of plane stress on a rectangular domain with rectangular elements. This finite element program provided a framework for the setting up of sample displacements and element geometries that could be used as an example input for the comparison of different methods for calculating the tangent operator. 

The finite element program was implemented with these parameters:
- Domain length=1 m (square)
- Number of elements in x direction=20
- Number of elements in y direction=20
- E=Young's modulus=200e3 MPa
- $E_T$=Hardening modulus MPa
- $\nu$=Poisson's ratio=0.3
- Thickness=0.01 m
- $\sigma_{y0}$=Initial yield stress=300 MPa
- $\beta$=Kinematic-isotropic interpolation factor=0 (kinematic only)

The program was implemented as stress-driven, which means that Newton-Raphson iterations were needed to compute the change in strain at each loading step. For testing, only a loading cycle (no unloading) was implemented, though the loading was split into several steps and the Newton-Raphson subroutine was run for each step. Approximately quadratic convergence was observed using Newton-Raphson, as expected. The elements were bidirectional quads (4-node quads). 

### Run
To run the program, run "script.mlx." Plots of the domain and graphs of some variables during the simulation will be displayed and saved. 

### Example Results

<img src="https://github.com/user-attachments/assets/1c8d1131-62d1-4f96-8331-967d5a24bd7e" alt="domain" width="600"/>
<img src="https://github.com/user-attachments/assets/77491f9e-ef72-436a-9f16-771e65c8a59e" alt="deformed" width="600"/>
<img src="https://github.com/user-attachments/assets/b0fb0669-c2a5-4f94-b274-9ac805b9ef2f" alt="stress" width="600"/>
