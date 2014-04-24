#include "Solvers.h"
#include <cmath>

namespace VPFFT
{
  namespace Solvers
  {

    
    //--------------------------------------------------------------------------------
    //  SolveConstitutiveEquations
    //    -  Solves the constitutive equation given the applied strain rate; i.e.,
    //       find the stress state that is compatible with the strain rate defined.
    //
    //  Applying Newton Raphson to solve the the stress state given the strain state of
    //  the system.
    //
    //--------------------------------------------------------------------------------
    EigenRep SolveConstitutiveEquations( const EigenRep           & InitialStressState,   // initial guess, either from Sach's or previous iteration
                                         const vector<EigenRep>   & SchmidtTensors,
                                         const vector<Float>      & CRSS,
                                         const vector<Float>      & GammaDotBase,         // reference shear rate
                                         const vector<int>      & RateSensitivity,
                                         const EigenRep           & StrainRate,           // Current strain rate - constant from caller
                                         Float                      EpsilonConvergence,
                                         Float                      MaxResolvedStress,
                                         int                        MaxNumIterations )
    {
      const Float Coeff = 0.2; // global fudge factor
      
      EigenRep CurrentStressState = InitialStressState;
      int RemainingIterations     = MaxNumIterations;
      EigenRep NewStressState = InitialStressState;
      EigenRep SavedState(0, 0, 0, 0, 0);  // used to return to old state
      
      while( RemainingIterations > 0 )
      {
        
        bool AdmissibleStartPointFound = false;
        std::vector<Float> RSS( SchmidtTensors.size(), 0 );            //  This is really RSS / tau
        do    // refresh critical resolved shear stress.
              // check to see if it is outside of the yield
              // surface, and therefore inadmissible
        {
          AdmissibleStartPointFound = true;
          for( int i = 0; i < SchmidtTensors.size(); i ++ )
          {
            RSS[i] = InnerProduct( SchmidtTensors[i], CurrentStressState) / CRSS[i];

            if( std::fabs( RSS[i] ) < 1e-10 )
              RSS[i] = 0;
            if( std::fabs( RSS[i] ) > MaxResolvedStress )
              AdmissibleStartPointFound = false;
          }
          if( !AdmissibleStartPointFound )
            CurrentStressState = SavedState + ( CurrentStressState - SavedState ) * Coeff;
          
          RemainingIterations --;
          if( RemainingIterations < 0 )
            return NewStressState;
        } while ( !AdmissibleStartPointFound ); 
          
        std::vector<Float> GammaDot( SchmidtTensors.size(), 0 );            //  This is really RSS / tau
        for( int i = 0; i < SchmidtTensors.size(); i ++ )
        {
          if( RateSensitivity[i] -1 > 0 )
            GammaDot[i] = GammaDotBase[i] * std::pow( std::fabs( RSS[i] ), static_cast<int>( RateSensitivity[i] - 1 ) );
          else
            GammaDot[i] = GammaDotBase[i];          
        }

        // Construct residual vector R(Sigma), where Sigma is stress.  R(Sigma) is a 5d vector
        EigenRep Residual( 0, 0, 0, 0, 0 ); // current estimate
        SMatrix5x5 ResidualJacobian;
        ResidualJacobian.SetZero();

        for( int i = 0; i < SchmidtTensors.size(); i ++ )
        {
          Residual += SchmidtTensors[i] * GammaDot[i] * RSS[i];  // Residual = StrainRate - sum_{k} [ m^{k}_i * |r_ss/ crss|^n ]

          // Construct F', or the Jacobian
          ResidualJacobian -= OuterProduct( SchmidtTensors[i], SchmidtTensors[i] )
            * RateSensitivity[i] * GammaDot[i] / CRSS[i];

        }
        
        Residual =  Residual - StrainRate ;  // need the negative residual, instead of E - R
        SavedState            = CurrentStressState;
        EigenRep Delta_Stress = LU_Solver( ResidualJacobian, Residual );     // <----------- Need to hangle error from this
        NewStressState        = CurrentStressState + Delta_Stress;        
        Float RelativeError = static_cast<Float>(2) * Delta_Stress.Norm() / ( NewStressState + CurrentStressState ).Norm();
        

        CurrentStressState = NewStressState;
        if( RelativeError < EpsilonConvergence )
        {
          break;
        }
      } // end while
      return NewStressState;
    }


    //--------------------------------------------------------------------------------
    //  SolveConstrainedConstitutiveEquations
    //
    //  --  This is *REALLY* close to the normal SolveConstitutiveEquations.  One could
    //      replace both of these functions with a function that takes a functor...  but
    //      may not do this till later.
    //--------------------------------------------------------------------------------
    EigenRep SolveConstrainedConstitutiveEquations( const EigenRep           & InitialStressState,
                                                    const vector<EigenRep>   & SchmidtTensors,
                                                    const vector<Float>      & CRSS,
                                                    const vector<Float>      & GammaDotBase,         // reference shear rate
                                                    const vector<int>      & RateSensitivity,
                                                    const EigenRep           & MacroscopicStrainRate,
                                                    const EigenRep           & LagrangeMultiplier,      // also known as - lagrange multiplier, or lambda(x)
                                                    const EigenRep           & LocalDisplacementVariation,
                                                    const SMatrix5x5         & HomogeonousReference,
                                                    Float                      EpsilonConvergence,
                                                    Float                      MaxResolvedStress,
                                                    int                        MaxNumIterations,
                                                    Float *NR_Error_Out )
    {
      const Float Coeff = 0.2; // global fudge factor
      
      EigenRep CurrentStressState = InitialStressState;
      int RemainingIterations     = MaxNumIterations;
      EigenRep NewStressState( -1, -3, -7, -9, -11 );
      Float RelativeError = 1e5;
      EigenRep SavedState = InitialStressState; 


      //-------------------
      // DEBUG
//       std::cout << "MacroscopicStrainRate " << MacroscopicStrainRate << std::endl;
//       std::cout << "LagrangeMultiplier " << LagrangeMultiplier << std::endl;
//       std::cout << "Local d-dot " << LocalDisplacementVariation << std::endl;
//       std::cout << "Local L \n " << HomogeonousReference << std::endl;
      //-------------------
      while( RemainingIterations > 0 )
      {
        EigenRep LoopStartStressState = CurrentStressState;
        bool AdmissibleStartPointFound = false;
        std::vector<Float> RSS( SchmidtTensors.size(), 0 );            //  This is really RSS / tau
        
        do    // refresh critical resolved shear stress.
              // check to see if it is outside of the yield
              // surface, and therefore inadmissible
        {
          AdmissibleStartPointFound = true;
          Float MaxCRSS = 0;
          
          for( int i = 0; i < SchmidtTensors.size(); i ++ )
          {
            
            RSS[i] = InnerProduct( SchmidtTensors[i], CurrentStressState) / CRSS[i];
            MaxCRSS = std::max( std::fabs( RSS[i] ), MaxCRSS );
            if( std::fabs( RSS[i] ) > MaxResolvedStress )
              AdmissibleStartPointFound = false;
          }
          
          if( !AdmissibleStartPointFound )
            CurrentStressState =  SavedState + ( CurrentStressState - SavedState ) * Coeff;
          
          RemainingIterations --;
          if( RemainingIterations < 0 )
          {
            return CurrentStressState; //  InitialStressState;      // not failing gracefully at all
          }
     
        } while ( !AdmissibleStartPointFound );

        std::vector<Float> GammaDot( SchmidtTensors.size(), 0 );            //  This is really RSS / tau
        for( int i = 0; i < SchmidtTensors.size(); i ++ )
        {
          if( RateSensitivity[i] != 1 )
            GammaDot[i] = GammaDotBase[i] * std::pow( std::fabs( RSS[i] ), ( RateSensitivity[i] - 1 ) );
          else
            GammaDot[i] = GammaDotBase[i];          
        }

        // Construct residual vector R(Sigma), where Sigma is stress.  R(Sigma) is a 5d vector
        EigenRep ResolvedStressSum( 0, 0, 0, 0, 0 ); // current estimate
        SMatrix5x5 ResidualJacobian;
        ResidualJacobian.SetZero();

        for( int i = 0; i < SchmidtTensors.size(); i ++ )
        {
          ResolvedStressSum += SchmidtTensors[i] * GammaDot[i] * RSS[i];  // Residual = StrainRate - sum_{k} [ m^{k}_i * |r_ss/ crss|^n ]

          // Construct F', or the Jacobian
          ResidualJacobian += OuterProduct( SchmidtTensors[i], SchmidtTensors[i] ) * RateSensitivity[i] * GammaDot[i] / CRSS[i]; // may need a sign
        }

        // -------------------------------
        //  If we use Eigen for this section, the
        //  expression templating will significantly improve
        //  the efficiency of this section significantly via
        //  explicit vectorization.  
        // -------------------------------
        SMatrix5x5 Identity5x5;
        Identity5x5.SetIdentity();
        ResidualJacobian  =   -Identity5x5 - (HomogeonousReference * ResidualJacobian) ;   // was I -
        EigenRep Residual =   CurrentStressState  - LagrangeMultiplier
          +   HomogeonousReference
          * ( ResolvedStressSum - MacroscopicStrainRate - LocalDisplacementVariation );     // validated
        
        SavedState            = CurrentStressState;
        EigenRep Delta_Stress = LU_Solver( ResidualJacobian, Residual );     // <----------- Need to hangle error from this
        NewStressState        = CurrentStressState + Delta_Stress;           // was - Delta
        RelativeError         =  (NewStressState - SavedState).Norm() /( ( NewStressState + SavedState) * static_cast<Float>(0.5)).Norm();
     
        CurrentStressState = NewStressState;
        if( RelativeError < EpsilonConvergence )
          break;
      } // end while

      *NR_Error_Out = RelativeError;

      return NewStressState;
    }
    
    
    //--------------------------------------------------------------------------------
    //
    //  SolveSachEquation
    //  Given the Schmidt tensor in the sample frame, return Sach's solution.
    //--------------------------------------------------------------------------------
    EigenRep SolveSachEquation( const vector<EigenRep> & SchmidtTensors,
                                const vector<Float>    & CRSS,
                                const EigenRep & StrainRate )
    {
      Float MaxCRSS = 0;
      for ( int i = 0; i < SchmidtTensors.size(); i ++ )
      {
        Float RSS = LinearAlgebra::InnerProduct( StrainRate,
                                                 SchmidtTensors[i] )  / CRSS[i];
        MaxCRSS = std::max( std::fabs( RSS ), MaxCRSS );
      }
      EigenRep SachGuess( StrainRate / MaxCRSS );
      return SachGuess;
    }
    

    //--------------------------------------------------------------------------------
    //
    //  ApplyConstitutiveEquations
    //  -- Apply Constitutive Equation to a given stress state and return the
    //     strain rate as a result.
    //
    //     Current used for testing the other parts of the code.
    //
    //--------------------------------------------------------------------------------
    EigenRep ApplyConstitutiveEquations( const EigenRep           & StressState,
                                         const vector<EigenRep>   & SchmidtTensors,
                                         const vector<Float>      & CRSS,
                                         const vector<Float>      & GammaDotBase,         // reference shear rate
                                         const vector<int>      & RateSensitivity)
    {

      EigenRep StrainRate( 0, 0, 0, 0, 0);
      for( int i = 0; i < SchmidtTensors.size(); i ++ )
      {
        Float GammaDot;
        Float RSS = InnerProduct( SchmidtTensors[i], StressState) / CRSS[i];
        if( RateSensitivity[i] > 0 )
          GammaDot = GammaDotBase[i] * std::pow( std::fabs( RSS ), static_cast<int>( RateSensitivity[i] ) );
        else
          GammaDot = GammaDotBase[i];
        StrainRate += SchmidtTensors[i] * GammaDot * sign( RSS );
      }
      return StrainRate;
    }
    
    //--------------------------------------------------------------------------------
    //
    //  CalculateShearRate
    //     Current used for testing the other parts of the code.
    //
    //--------------------------------------------------------------------------------
    SMatrix3x3 CalculateShearRate( const EigenRep           & StressState,
                                   const vector<EigenRep>   & SchmidtTensors,
                                   const vector<SMatrix3x3> & AntiSchmidtTensors,
                                   const vector<Float>      & CRSS,
                                   const vector<Float>      & GammaDotBase,         // reference shear rate
                                   const vector<int>        & RateSensitivity )
    {
      SMatrix3x3 ShearRate;
      ShearRate.SetZero();
      for( int i = 0; i < SchmidtTensors.size(); i ++ )
      {
        Float GammaDot;
        Float RSS = InnerProduct( SchmidtTensors[i], StressState) / CRSS[i];
        if( RateSensitivity[i] > 0 )
          GammaDot = GammaDotBase[i] * std::pow( std::fabs( RSS ), static_cast<int>( RateSensitivity[i] ) );
        else
          GammaDot = GammaDotBase[i];
        ShearRate += AntiSchmidtTensors[i] * GammaDot * sign( RSS );
      }
      return ShearRate;
    }


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
                                      const vector<Float> & Theta1  )
    {

      vector<Float> UpdatedCRSS    = CRSS;
      Float LocalShearAccumulation = 0;
      std::vector<Float> GammaDot( SchmidtTensors.size(), 0 );            
      for( int i = 0; i < SchmidtTensors.size(); i ++ )
      {
        Float RSS = InnerProduct( SchmidtTensors[i], StressState) / CRSS[i];
        if( RateSensitivity[i] > 0 )
          GammaDot[i] = GammaDotBase[i] * std::pow( std::fabs( RSS ), static_cast<int>( RateSensitivity[i] ) );
        else
          GammaDot[i] = GammaDotBase[i];
        LocalShearAccumulation += std::fabs( GammaDot[i] ) * TimeStep;
      }



      for( int i = 0; i < SchmidtTensors.size(); i ++ )
      {
        Float DTau = 0;
        for( int j = 0; j < SchmidtTensors.size(); j ++ )
          DTau += HardeningMatrix[i][j] * std::fabs( GammaDot[j] ) * TimeStep;


        Float Tiny = 1e-4 * Tau0[i];
        Float VoceFactor = 0;
        if( std::fabs( Theta0[i] ) > Tiny  )
        {
          VoceFactor = Theta1[i] * LocalShearAccumulation;
          if( std::fabs( Tau1[i]) > Tiny )
          {
            Float fRatio = std::fabs( Theta0[i] / Tau1[i] );
            Float ExpIni = exp( -AccumulatedShear       * fRatio );
            Float ExpDel = exp( -LocalShearAccumulation * fRatio );

            VoceFactor = VoceFactor - (fRatio * Tau1[i] * Theta1[i])/ fRatio * ExpIni
              * ( ExpDel - static_cast<Float>(1) ) - Theta1[i]/fRatio * ExpIni
              * ( ExpDel * ( ( AccumulatedShear + LocalShearAccumulation)
                             * fRatio + static_cast<Float>(1) )
                  - ( AccumulatedShear * fRatio + static_cast<Float>(1) ) );
          }
        }
        UpdatedCRSS[i] += DTau * VoceFactor / LocalShearAccumulation;
      }

      AccumulatedShear += LocalShearAccumulation;
      return UpdatedCRSS;
    }


    
  } // namespace Solver

  
}  // namespace VPFFT
