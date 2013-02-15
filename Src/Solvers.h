//----------------------------------------------------------------
//
// Solvers
//
//
//----------------------------------------------------------------


#ifndef VPFFT_SOLVERS_H
#define VPFFT_SOLVERS_H

#include "LinearAlgebra.h"
#include "MaterialGrid.h"


namespace VPFFT
{
  namespace Solvers
  {
    //--------------------
    // classic solution
    //--------------------
    template <typename T> int sign(T val) {
      return (T(0) < val) - (val < T(0));
    }
   
    using namespace LinearAlgebra;
    using std::vector;


    //--------------------------------------------------------------------------------
    //  NewtonRaphson -
    //    template< class ResidualOperator< class VariableType >, class ResidualDerivativeOperator< class VariableType > >
    //    NewtonRaphson( ResidualOperator R_Op, const VariableType X )
    //
    //   ResidualOperator           <->  function object, which maps X :-> a vector of 1      * dim(X)
    //   ResidualDerivativeOperator <->  function object, which maps X :-> a matrix of dim(X) * dim(X)
    //
    //  - to be done later
    //--------------------------------------------------------------------------------
    
    //--------------------------------------------------------------------------------
    //  SolveConstitutiveEquations
    //    -  Solves the constitutive equation given the applied strain rate; i.e.,
    //       find the stress state that is compatible with the strain rate defined.
    //
    //
    //  Would probably be good if there's a generalized "Newton Raphson" or something that
    //  enforces the constrain needed...
    //--------------------------------------------------------------------------------
    EigenRep SolveConstitutiveEquations( const EigenRep           & InitialStressState,
                                         const vector<EigenRep>   & SchmidtTensors,
                                         const vector<Float>      & CRSS,
                                         const vector<Float>      & GammaDotBase,         // reference shear rate
                                         const vector<int>        & RateSensitivity,
                                         const EigenRep           & StrainRate,        
                                         Float                      EpsilonConvergence,
                                         Float                      MaxResolvedStress,
                                         int                        MaxNumIterations );



    //--------------------------------------------------------------------------------
    //  SolveConstrainedConstitutiveEquations
    //   Solve the constitutive equation with the inhomogeneity and all
    //
    //--------------------------------------------------------------------------------
    EigenRep SolveConstrainedConstitutiveEquations( const EigenRep           & InitialStressState,
                                                    const vector<EigenRep>   & SchmidtTensors,
                                                    const vector<Float>      & CRSS,
                                                    const vector<Float>      & GammaDotBase,         // reference shear rate
                                                    const vector<int>        & RateSensitivity,
                                                    const EigenRep           & MacroscopicStrainRate,
                                                    const EigenRep           & LagrangeMultiplier,
                                                    const EigenRep           & LocalDisplacementVariation,
                                                    const SMatrix5x5         & HomogeonousReference,
                                                    Float                      EpsilonConvergence,
                                                    Float                      MaxResolvedStress,
                                                    int                        MaxNumIterations,
                                                    Float                    *NR_Error_Out );



    
    //--------------------------------------------------------------------------------
    //
    //  SolveSachEquation
    //  Given the Schmidt tensor in the sample frame, return Sach's solution.
    //--------------------------------------------------------------------------------
    EigenRep SolveSachEquation( const vector<EigenRep> & SchmidtTensors,
                                const vector<Float>    & CRSS,
                                const EigenRep & StrainRate );

    
    //--------------------------------------------------------------------------------
    //
    //  SolveConstrainedConstitutiveEquations
    //
    //--------------------------------------------------------------------------------



    //--------------------------------------------------------------------------------
    //
    //  ApplyConstitutiveEquations
    //  -- Apply Constitutive Equation to a given stress state and return the
    //     strain rate as a result.
    //
    //
    //--------------------------------------------------------------------------------
    EigenRep ApplyConstitutiveEquations( const EigenRep           & StressState,
                                         const vector<EigenRep>   & SchmidtTensors,
                                         const vector<Float>      & CRSS,
                                         const vector<Float>      & GammaDotBase,         // reference shear rate
                                         const vector<int>      & RateSensitivity );


    //--------------------------------------------------------------------------------
    //
    //  CalculateShearRate
    //  -- Apply Constitutive Equation to a given stress state and return the
    //     shear rate as a result.  AntiSchmidtTensors are in the same frame as the schmidt tensor
    //
    //
    //--------------------------------------------------------------------------------
    SMatrix3x3 CalculateShearRate( const EigenRep           & StressState,
                                   const vector<EigenRep>   & SchmidtTensors,
                                   const vector<SMatrix3x3> & AntiSchmidtTensors,
                                   const vector<Float>      & CRSS,
                                   const vector<Float>      & GammaDotBase,         // reference shear rate
                                   const vector<int>        & RateSensitivity );
    
    //--------------------------------------------------------------------------------
    //
    //  CalculateHardening
    //     Current used for testing the other parts of the code.
    //
    //
    //  This is almost a direct translation from Ricardo's "HARDEN" routine from fft3.for
    //
    //--------------------------------------------------------------------------------
    vector<Float> CalculateHardening( Float & AccumulatedShear,
                                      const Float              & TimeStep,
                                      const EigenRep           & StressState,
                                      const vector<EigenRep>   & SchmidtTensors,
                                      const vector<Float>      & CRSS,
                                      const vector<Float>      & GammaDotBase,         // reference shear rate
                                      const vector<int>        & RateSensitivity,
                                      const vector< vector< Float > > & HardeningMatrix,
                                      const vector<Float> & Tau0,
                                      const vector<Float> & Tau1,
                                      const vector<Float> & Theta0,
                                      const vector<Float> & Theta1  );

  }
}


#endif
