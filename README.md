# CCFULL_progame
CCFull computes fusion cross sections as a function of incident energy for heavy-ion nuclear reactions. The program implements four different theoretical approaches

# CCFull: Coupled-Channels Fusion Cross Section Calculator

A Fortran program for calculating fusion cross sections using various nuclear reaction models with coupled-channels formalism.

## Overview

CCFull computes fusion cross sections as a function of incident energy for heavy-ion nuclear reactions. The program implements four different theoretical approaches:

1. **Numerov (Direct) Calculation** 
2. **Discrete basis Method (DBM)**
3. **Modified Discrete basis Method (MDBM)** 
4. **R-matrix Theory** 

## Features

- Coupled-channels calculations for nuclear fusion reactions
- Multiple theoretical models for cross-validation
- Angular momentum distribution calculations

## Prerequisites

### Required Libraries
- **LAPACK** - Linear algebra package
- **BLAS** - Basic linear algebra subprograms
- **gfortran** - GNU Fortran compiler (version 7.0 or higher recommended)


## Program Structure
ccfull/
├── Makefile              # Build configuration
├── main.f90              # Main program
├── initialization.f90    # Initialization module
├── Interaction.f90       # Nuclear interaction potentials
├── CoulombWave.f90       # Coulomb wave functions
├── Numerov.f90          # Numerov integration method
├── DiscreteBasis.f90     # Discrete basis 
├── ModifiedDiscreteBasis.f90 # Modified Discrete basis
├── Rmatrix.f90          # R-matrix calculations
└── README.md            # This file

## Usage
compile: make


## Example Output
              E      sigma_num      sigma_dbm     sigma_mdbm  sigma_rmatrix      <L>(hbar)      <L>(hbar)      <L>(hbar)      <L>(hbar)      P (L = 0)      P (L = 0)      P (L = 0)      P (L = 0)
       55.00000    0.97448E-02    0.99270E-02    0.98093E-02    0.97437E-02        5.87031        5.87031        5.87031        5.87031    0.22629E-03    0.22967E-03    0.22705E-03    0.22625E-03
       56.00000        0.05489        0.05545        0.05497        0.05488        5.94333        5.94333        5.94333        5.94333    0.12592E-02    0.12643E-02    0.12534E-02    0.12590E-02
       57.00000        0.28583        0.28575        0.28441        0.28578        6.05134        6.05134        6.05134        6.05134    0.64144E-02    0.63493E-02    0.63474E-02    0.64135E-02
       58.00000        1.36500        1.35843        1.35333        1.36483        6.19272        6.19272        6.19272        6.19272    0.29721E-01    0.29472E-01    0.29433E-01    0.29718E-01
       59.00000        5.84375        5.80359        5.80532        5.84318        6.40451        6.40451        6.40451        6.40451    0.11981E+00    0.11981E+00    0.11957E+00    0.11980E+00
       60.00000       20.59856       20.61289       20.57591       20.59679        6.86092        6.86092        6.86092        6.86092    0.35491E+00    0.35802E+00    0.35676E+00    0.35488E+00
       61.00000       52.14434       52.45950       52.21295       52.13902        7.81887        7.81887        7.81887        7.81887    0.62839E+00    0.63535E+00    0.63023E+00    0.62830E+00
       62.00000       94.62476       95.05488       94.68164       94.61396        9.18913        9.18913        9.18913        9.18913    0.75138E+00    0.75517E+00    0.75072E+00    0.75126E+00
       63.00000      139.58987      140.09738      139.56203      139.57394       10.65032       10.65032       10.65032       10.65032    0.79530E+00    0.79767E+00    0.79439E+00    0.79518E+00
       64.00000      185.55959      186.08861      185.47784      185.53902       11.98384       11.98384       11.98384       11.98384    0.86342E+00    0.86742E+00    0.86342E+00    0.86330E+00
       65.00000      234.04527      234.57189      233.93054      234.02059       13.13045       13.13045       13.13045       13.13045    0.94486E+00    0.94765E+00    0.94444E+00    0.94474E+00
       66.00000      283.93526      284.44688      283.76725      283.90671       14.18620       14.18620       14.18620       14.18620    0.98364E+00    0.98670E+00    0.98280E+00    0.98351E+00
       67.00000      333.26115      333.72636      333.04373      333.22904       15.21129       15.21129       15.21129       15.21129    0.99507E+00    0.99845E+00    0.99440E+00    0.99493E+00
       68.00000      381.21016      381.62127      380.95965      381.17495       16.20563       16.20563       16.20563       16.20563    0.99820E+00    0.10017E+01    0.99777E+00    0.99806E+00
       69.00000      427.61803      427.95219      427.34838      427.58016       17.16333       17.16333       17.16333       17.16333    0.99921E+00    0.10030E+01    0.99890E+00    0.99906E+00
       70.00000      472.48080      472.74120      472.20209      472.44088       18.08211       18.08211       18.08211       18.08211    0.99962E+00    0.10032E+01    0.99936E+00    0.99947E+00
       71.00000      515.83672      515.98292      515.55608      515.79521       18.96273       18.96273       18.96273       18.96273    0.99981E+00    0.10035E+01    0.99956E+00    0.99966E+00
       72.00000      557.73621      557.73624      557.45847      557.69354       19.80734       19.80734       19.80734       19.80734    0.99990E+00    0.10035E+01    0.99964E+00    0.99975E+00
       73.00000      598.23608      598.09166      597.96377      598.19276       20.61860       20.61860       20.61860       20.61860    0.99995E+00    0.10035E+01    0.99968E+00    0.99979E+00
       74.00000      637.39727      637.06321      637.13062      637.35372       21.39926       21.39926       21.39926       21.39926    0.99997E+00    0.10035E+01    0.99969E+00    0.99981E+00
       75.00000      675.28418      674.72280      675.02140      675.24073       22.15197       22.15197       22.15197       22.15197    0.99999E+00    0.10035E+01    0.99970E+00    0.99982E+00



Citation
If you use this code in your research, please cite:

[Computer Physics Communications 123 (1999) 143–152]
[Phys. Rev. C 112, 024618 ]
