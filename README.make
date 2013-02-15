==================================================
           V P F F T + + 
==================================================

-----------
Description:  This is a C++ implementation of Ricardo Lebensohn's VPFFT, 
-----------   a crystal plasticity simulation code.  The original 
              implementation (fft3.for) is in Fortran 77.  This C++ VPFFT is
              developed to remove the requirement of material science specific
              knowledge from porting this code. 

-----------
Authors:      Frankie Li (li31@llnl.gov)
-----------

-----------
Requirements:  
               *Eigen  -- A templated linear algebra library designed for
                automatic vectorization.
               *FFTW   -- Fourier transform library.
               
-----------
Content:
-----------
For convinience, libXDM is included in this SVN.  Eigen, FFTW, and BOOST will
have to be downloaded and specified in the Makefile

Src/           Source code for VPFFT
               -- MaterialGrid.h/cpp  - implementation of VPFFT
               -- UnitTester.h/cpp    - example run case for single time step
               -- LinearAlgebra.h/cpp - Linear algebra functions - to
                                        be merged with the rest of the required
                                        functions from libXDM.
libXDM/        Source for libXDM


Makefile       make VPFFT++   - compiles optimized code
               make debug     - compiles debug code with -O0
-----------
TODO:
-----------
1.  Profile and compare with parent app.
2.  Put in mechanism and interface for dimensional-reduction/coarse-scaling
scheme.