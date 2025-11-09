# Fortran ODE Solver 

This Fortran is a simple computational physics project that solves two classic ordinary differential equations (ODEs) using two different numerical methods. It compares the numerical approximations to the exact analytical solutions.

The program solves:
1. **Radioactive Decay** (1st-order ODE): $dN/dt = -\lambda N$
2. **Simple Harmonic Motion** (2nd-order ODE): $d^2y/dt^2 = -\omega^2 y$

Numerical methods:
1. **Euler Method** 
2. **Leapfrog Method** 

## Features 
* **Modular:** All numerical methods and physics functions are organised inside the `methods` module.
* **Double Precision:** Offers up to 15 decimal places of accuracy.
* **Data Output:** Generates easy to plot .dat files for both the true (analytical) and approximated (numerical) solutions.
* **Performance Logging:** Appends the time step `dt`, final error and computation time to log files, allowing for convergence and performance analysis. 

## How to Compile and Run
You will need a Fortran compiler, `gfortran` (GCC) is recommended.

1. Save the code as ode.f90 
2. Open your terminal and navigate to the directory containing the file. Run the following command: 

```bash
gfortran -O3 ode.f90 -o ode
```

3. After compiling, run the executable from your terminal:

```bash
./ode
```

4. The program will run all four simulations and print a summary of the results to your console. 

### Configuration
Simulation parameters such as time step, duration and intial conditions can be configured in the `program ode` block in `ode.f90`
