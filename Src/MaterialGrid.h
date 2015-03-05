#ifndef VPFFT_MATERIAL_GRID_H
#define VPFFT_MATERIAL_GRID_H

#include "LinearAlgebra.h"
#include <vector>
#include "BurgerVectorBasis.h"
//#include "fftw3.h"
#include <string>
#include "Solvers.h"
#include <omp.h>
#include <fftw3-mpi.h>
#include <queue>
#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_num_threads() 0
    #define omp_get_thread_num() 0
#endif

//--------- For MPI debug
#ifndef VPFFT_MPI_TAG
    #define VPFFT_MPI_TAG 2
#endif

namespace VPFFT
{
  namespace DataStructures
  { 
    using std::vector;
    using VPFFT::LinearAlgebra::SMatrix3x3;
    using VPFFT::LinearAlgebra::EigenRep;
    using VPFFT::LinearAlgebra::SMatrix5x5;
    using VPFFT::LinearAlgebra::SVector3;
    //-----------------------------
    //  class MaterialGrid
    //
    //        Main data structure for VPFFT.
    // Note:        To parameterize VPFFT for different
    //              architectures, one may consider
    //              using templates for different types
    //              of data structures used for the MaterialGrid,
    //              which can then be coupled with different
    //              FFT and solver strategies.  
    //  *** Designed explicitly for FFTW *ONLY
    //
    //   Consider this to define the interface for which all
    //   displacement field classes must satistfy. (probably
    //   should create an abstract class to enforce this)
    //
    //
    //   NOTE:  Changing precision in FFTW requires changing
    //          the library to be linked at compiled time.
    //          One should be aware that changing the definition
    //          of "Float" here would require the change of
    //          library as well.  Note also that FFTW is *NOT*
    //          compatible with arbitary precision arithmetic types
    //          such as those found in GMP and LEDA.
    //-----------------------------
    class MaterialGrid
    {
    public:

      typedef fftw_complex FFTW_Complex;
      typedef FFTW_Complex * FFTW_ComplexPtr;
      typedef fftw_plan    FFTW_Plan;
      
      MaterialGrid()
        : DimX( 0 ),
          DimY( 0 ),
          DimZ( 0 ), bInitialized( false ) {}
    
      MaterialGrid( int DimX_, int DimY_, int DimZ_ )
        : DimX( DimX_ ),
          DimY( DimY_ ),
          DimZ( DimZ_ )
      {
        InitializeMemory( DimX_, DimY_, DimZ_ );
      }

      //----------------------------
      // ParallelInitialize
      //   Distribute the initial values as specified by Orientations and Stress
      //   to all of the N processors.
      //----------------------------
      void SendDataAndInitialize( SMatrix3x3 * Orientations,
				  EigenRep * Stress );

      //----------------------------
      //  All clients must recv data.
      //----------------------------
      void RecvDataAndInitialize( );
      
      
      //----------------------------
      //  SetVoxelLength
      //----------------------------
      void SetVoxelLength( Float _xLength, Float _yLength, Float _zLength )
      {
        XVoxelLength = _xLength;
        YVoxelLength = _yLength;
        ZVoxelLength = _zLength;
      }

      //----------------------------
      //  UpdateSchmidtTensors
      //
      //   TODO:  Change the structure of SMatrix3x3, then apply
      //          alias from eigen to take advanage of the SIMD
      //          optimization - probably a good number of factors
      //          faster.
      //----------------------------
      void UpdateSchmidtTensors( );
      
      
      //----------------------------
      //  ToIndex
      //     A consistent way to convert x, y, z index to 1D index
      //
      //     following fftw's convention (or C convention)
      //----------------------------
      inline int ToIndex( int x, int y, int z ) const
      {
        return x * DimYZ + y * DimZ + z;
      }

      //----------------------------
      // Orientation - accesor
      //----------------------------
      inline const SMatrix3x3 & Orientation( int x, int y, int z ) const
      { return OrientationField[ ToIndex(x,y,z)]; }

      //----------------------------
      // Orientation - mutator
      //----------------------------    
      inline SMatrix3x3 & Orientation( int x, int y, int z )
      { return OrientationField[ ToIndex(x,y,z)]; }
      
      //----------------------------
      // SchmidtTensor - accesor
      //----------------------------
      inline const EigenRep & SchmidtTensor( int x, int y, int z, int SlipSystem ) const
      { return LocalSchmidtTensor[ ToIndex(x, y, z) ][ SlipSystem ]; }

      //----------------------------
      // SchmidtTensor - mutator
      //----------------------------    
      inline EigenRep & SchmidtTensor( int x, int y, int z, int SlipSystem )
      { return LocalSchmidtTensor[ ToIndex(x, y, z) ][ SlipSystem ]; }

      //----------------------------
      // SchmidtTensorList - accesor
      //----------------------------
      inline const vector<EigenRep> & SchmidtTensorList( int x, int y, int z ) const
      { return LocalSchmidtTensor[ ToIndex(x, y, z) ]; }

      //----------------------------
      // SchmidtTensorList - mutator
      //----------------------------    
      inline vector<EigenRep> & SchmidtTensorList( int x, int y, int z )
      { return LocalSchmidtTensor[ ToIndex(x, y, z) ]; }

      
      //----------------------------
      // Stress - accesor
      //----------------------------
      inline  const EigenRep & Stress( int x, int y, int z ) const
      { return StressField[ ToIndex(x, y, z) ]; }

      //----------------------------
      // Stress - mutator
      //----------------------------    
      inline EigenRep & Stress( int x, int y, int z )
      { return StressField[ ToIndex(x, y, z) ]; }

      //----------------------------
      // LagrangeMultiplier - accesor
      //----------------------------
      inline  const EigenRep & LagrangeMultiplier( int x, int y, int z ) const
      { return LagrangeMultiplierField[ ToIndex(x, y, z) ]; }

      //----------------------------
      // LagrangeMultiplier - mutator
      //----------------------------    
      inline EigenRep & LagrangeMultiplier( int x, int y, int z )
      { return LagrangeMultiplierField[ ToIndex(x, y, z) ]; }


      //----------------------------
      // SetPolarizationField
      //----------------------------    
      void SetPolarizationField( int x, int y, int z, const EigenRep & s );
      void AddPolarizationField( int x, int y, int z, const EigenRep & s );
      

      //-----------------------------------
      //  InitializeData
      //      Input data into MaterialGrid 
      //
      //-----------------------------------
      void InitializeData( SMatrix3x3 * Orientations,
                           EigenRep  * Stress  );
      
      //----------------------------
      //  RunSingleStrainStep
      //
      //  Purpose:   A single strain step in the VPFFT code involves
      //             calculating iteratively the equilibrium stress
      //             state of the system at every spatial point.  This
      //             typically involves tens to hundreds of FFT and
      //             Newton-Raphson steps.
      //
      //  Parameters:  MacroscopicStrainRate - the imposed, macroscopic strain
      //               rate of the system.
      //
      //               NumMaxIterations   - maximum number of iterations to be
      //               ran.  One iteration requires 1 FFT and 1 inverse FFT step.
      //
      //               ConvergenceEpsilon -  deviation of stress and strain from
      //                                     the solution to be accepted as converged.
      //----------------------------
      EigenRep RunSingleStrainStep( const EigenRep & MacroscopicStrainRate,
                                    int NumMaxIterations, Float TimeStep,
                                    Float ConvergenceEpsilon );
      

      //----------------------------
      //  AverageError
      //----------------------------
      
      //----------------------------
      //  UpdateOrientation
      //----------------------------

      
      int NumX() const { return DimX; };
      int NumY() const { return DimY; };
      int NumZ() const { return DimZ; };


      //----------------------------
      //  PrintStrain
      //----------------------------
      void PrintStrain( const std::string & filename );
      void PrintStress( const std::string & filename );
      void PrintPolarization( const std::string & filename );
      void PrintLagrangeMultiplier( const std::string & filename );
      
    private:

      int DimX;
      int DimY;
      int DimZ;
      int DimYZ;

      ptrdiff_t DimXLocal;      // Currently using "slap" decomposition for FFT.  DimX gets modified. 
      ptrdiff_t DimXLocalStart; // the use of ptrdiff_t is to conform with fftw-mpi's spec
      
      Float XVoxelLength;
      Float YVoxelLength;
      Float ZVoxelLength;
      
      bool bInitialized;

      //-------------------------------------------------------
      //
      //  Material Propertes -- DEBUG --  - to be read from file
      // to be incorperated into the class... 
      //-------------------------------------------------------

      vector< vector<Float> > LocalCRSS;
      vector<Float>           GammaDotBase;    // reference rate for each slip systems
      vector<int>             RateSensitivity; // rate sensitivity for each system


      //-------------------------------------------------------
      //    Hardening parameters.
      //
      vector< vector<Float> > HardeningMatrix;
      vector<Float>  Tau0;
      vector<Float>  Tau1;
      vector<Float>  Theta0;
      vector<Float>  Theta1;
      //-------------------------------------------------------

      //-------------------------------------------------------
      // vector< SMatrix3x3 >          OrientationField;
      // vector< EigenRep >            StressField;
      
      SMatrix3x3 *                  OrientationField;
      EigenRep   *                  StressField;
      vector< EigenRep >            LagrangeMultiplierField;
      vector< vector<EigenRep> >    LocalSchmidtTensor;
      vector< Float >               AccumulatedShear;
      

      //-----------------------------------
      //  WARNING!!!!!
      //  Ordering of the FFT vectors should depend
      //  on the implementation of FFT algorithm. If
      //  the algorithm is insensitive to strides,
      //  then m should really be "ordered like"
      //  m[NumElements][5]
      //  to improve performance when accessing
      //  individual matrix.  The current arrangement
      //  should have *EXTREMELY BAD* (20% cache hit if
      //  we're lucky)  cache performance in all of
      //  the non-FFT parts.  However, I have
      //  yet to test how well the FFT parts perform 
      //  with stride, so I am not messing with this
      //  implementation right now.
      //
      //  The current
      //  implementation is optimized for transparancy
      //  and readability.
      //
      //  
      //
      //-----------------------------------

      FFTW_ComplexPtr _DisplacementVariation[5];   //  5 components, in-place FFT
      FFTW_ComplexPtr _PolarizationField    [5];   

      //-----------------------------------
      // create plans
      //  switch to the following plans if we are using real only option in FFTW 
      //         FFT_ForwardPlan  = fftw_plan_dft_r2c_3d();
      //         FFT_BackwardPlan = fftw_plan_dft_c2r_3d();
      //-----------------------------------
      FFTW_Plan       _DisplacementBackwardPlan    [5];
      FFTW_Plan       _PolarizationFieldForwardPlan[5];
      
      //-----------------------------------
      // Accessors
      //
      //  (x, y, z) -> spatial index
      //  (i, j)    -> indices of the matrix
      //-----------------------------------
      inline const FFTW_Complex & DisplacementVariation( int x, int y, int z,
                                                         int i ) const
      { return _DisplacementVariation[ i ][ ToIndex(x, y, z) ]; }

      //-----------------------------------
      //  DisplacementVariation
      //-----------------------------------
      inline FFTW_Complex &  DisplacementVariation( int x, int y, int z,
                                                    int i ) 
      { return _DisplacementVariation[ i ][ ToIndex( x, y, z ) ]; }



      //-----------------------------------
      //  DisplacementVariation
      //-----------------------------------
      inline EigenRep DisplacementVariation( int x, int y, int z ) const
      {
        EigenRep oRes;
        oRes(0) = DisplacementVariation( x, y, z, 0 )[0];
        oRes(1) = DisplacementVariation( x, y, z, 1 )[0];
        oRes(2) = DisplacementVariation( x, y, z, 2 )[0];
        oRes(3) = DisplacementVariation( x, y, z, 3 )[0];
        oRes(4) = DisplacementVariation( x, y, z, 4 )[0];
        return oRes;
      }

            //-----------------------------------
      //  DisplacementVariation
      //-----------------------------------
      inline EigenRep DisplacementVariation_Im( int x, int y, int z ) const
      {
        EigenRep oRes;
        oRes(0) = DisplacementVariation( x, y, z, 0 )[1];
        oRes(1) = DisplacementVariation( x, y, z, 1 )[1];
        oRes(2) = DisplacementVariation( x, y, z, 2 )[1];
        oRes(3) = DisplacementVariation( x, y, z, 3 )[1];
        oRes(4) = DisplacementVariation( x, y, z, 4 )[1];
        return oRes;
      }


      
      //-----------------------------------
      //  PolarizationField_Re
      //-----------------------------------
      inline EigenRep PolarizationField_Re( int x, int y, int z ) const
      {
        EigenRep oRes;
        oRes(0) = PolarizationField( x, y, z, 0 )[0];
        oRes(1) = PolarizationField( x, y, z, 1 )[0];
        oRes(2) = PolarizationField( x, y, z, 2 )[0];
        oRes(3) = PolarizationField( x, y, z, 3 )[0];
        oRes(4) = PolarizationField( x, y, z, 4 )[0];
        return oRes;
      }

      
      //-----------------------------------
      //  PolarizationField_Re
      //-----------------------------------
      inline EigenRep PolarizationField_Im( int x, int y, int z ) const
      {
        EigenRep oRes;
        oRes(0) = PolarizationField( x, y, z, 0 )[1];
        oRes(1) = PolarizationField( x, y, z, 1 )[1];
        oRes(2) = PolarizationField( x, y, z, 2 )[1];
        oRes(3) = PolarizationField( x, y, z, 3 )[1];
        oRes(4) = PolarizationField( x, y, z, 4 )[1];
        return oRes;
      }

      
      
      //-----------------------------------
      //   PolarizationField
      //-----------------------------------      
      inline  const FFTW_Complex & PolarizationField( int x, int y, int z,
                                                      int i ) const 
      { return _PolarizationField[ i ][ ToIndex( x, y, z )]; }

      //-----------------------------------
      //  PolarizationField
      //-----------------------------------
      inline  FFTW_Complex & PolarizationField( int x, int y, int z,
                                                   int i ) 
      { return _PolarizationField[ i ][ ToIndex( x, y, z ) ]; }


      //-----------------------------------
      // CRSS - gets the CRSS of a system.  Note
      // that a reference is returned directly.
      //-----------------------------------
      inline vector<Float> & CRSS( int x, int y, int z )
      {
        return LocalCRSS[ ToIndex( x, y, z ) ];
      }
      
      //-----------------------------------
      //  Initialize
      //      Uniform initialization
      //
      //-----------------------------------
      void InitializeMemory( int DimX_, int DimY_, int DimZ_ );


      //-----------------------------------
      //  InitializeComputation
      //    Run through the initialization steps required
      //    for the VPFFT algorithm.
      //
      //-----------------------------------
      void InitializeComputation( const EigenRep & MacroscopicStrainRate );  
      
      //-----------------------------------
      //  Cleanup
      //-----------------------------------
      void Cleanup();



      int DEBUG_INT;  // Use for DEBUG state saving only!!!!
      
      //------------------------- DEBUG PUBLIC: - to be made private later
    public:

      //--------------------------------------------------------------------------------
      //  BuildLinearRefMediumStiffness
      //
      //--------------------------------------------------------------------------------
      SMatrix5x5 BuildLinearRefMediumStiffness( );
      
      //-----------------------------------
      //  ConstructGreenOperator
      //
      //   WaveVector <-->  \vec{k}
      //-----------------------------------
      LinearAlgebra::Tensor4Rank ConstructGreensOperator( const SMatrix5x5 & LinearRefMediumStiffness,
                                                          const SVector3 & WaveVector ) const;
      
      //--------------------------------------------------------------------------------
      //  StrainRateEvolutionWithFFT
      //
      //  Purpose:  Update the new local strain rate deviation at each point by
      //            following the algorithm perscribed in Lebensohn, 2001.  This
      //            is also where FFT is applied, and therefore substatial amount
      //            of computation is contained within this function. 
      //--------------------------------------------------------------------------------
      void StrainRateEvolutionWithFFT( const SMatrix5x5 & LinearRefMediumStiffness );

      //--------------------------------------------------------------------------------
      //  OrientationEvolutionWithFFT
      //   - Exactly the same as StrainRateEvolutionWithFFT, but uses an anti-symmetrized
      //     Greens' operator.  Unfortunately, displacement variation will have to be
      //     used to store the orientation variation ( \delta \Omega )
      //--------------------------------------------------------------------------------
      void OrientationEvolutionWithFFT( const SMatrix5x5 & LinearRefMediumStiffness,
                                        const Float TimeStep );

      //--------------------------------------------------------------------------------
      //  CalculateShearRate
      //
      //  \dot{omega}^{p}_{ij} = \sum_k \alpha_{ij}^{k}( \vec{x} ) \dot{\gamma}^k{\vec{x}}
      //
      //  \alpha_{ij}^{k} -  Anti-symmetric schmidt tensor for the kth system.
      //
      //--------------------------------------------------------------------------------
      SMatrix3x3 ShearRate( int i, int j, int k);
      

      //----------------------------
      //  UpdateHardeningCRSS
      //
      //  Purpose:  Update the hardening models.  Can be used as a way 
      //            to plug into high resolution, multi-physics simulation.
      //----------------------------
      void UpdateHardeningCRSS( Float TimeStep );

      //----------------------------
      //  UpdateVoxelLength
      //----------------------------
      void UpdateVoxelLength( Float TimeStep, const SMatrix3x3 & MacroscopicStrainRate );
      
      //--------------------------------------------------------------------------------
      //  UpdateGlobalVoxelLength
      //
      //  UpdateGlobalGridSize - note that in the current VPFFT, we are assuming that the
      //  entire grid deforms uniformly due to imposed strain rate, which is clearly not
      //  true in general.
      //--------------------------------------------------------------------------------
      void UpdateGlobalVoxelLength( const EigenRep & MacroscopicStrainRate, Float dt );
      
      //--------------------------------------------------------------------------------
      //  UpdateLagrangeMultipler
      //
      //     LocalStrainRateVariation used here comes from the difference of
      //     the macroscopic strain and application of the constitutive equation
      //     on the new stress state.
      //--------------------------------------------------------------------------------
      void UpdateLagrangeMultiplier( int i, int j, int k,
                                     const SMatrix5x5 & LinearRefStiffness,
                                     const EigenRep & LocalStrainRate,
                                     const EigenRep & MacroscopicStrainRate,
                                     Float FudgeFactor );

      
      //--------------------------------------------------------------------------------
      //  ManageYeildSurface
      //
      //    Make sure that the stress state is no more than 10% outside of the yield surface.
      //    ( 110% from Richardo's code )
      //--------------------------------------------------------------------------------
      inline Float GetMaxRSS( const EigenRep & StressState,
                              const vector<EigenRep> & SchmidtTensors,
                              const vector<Float> & _LocalCRSS ) const
      {
        Float MaxRSS = 0;
        for ( int SlipSystem = 0; SlipSystem < SchmidtTensors.size(); SlipSystem ++ )
        { 
          Float RSS = LinearAlgebra::InnerProduct( SchmidtTensors[ SlipSystem ], StressState );
          MaxRSS = std::max( std::fabs( RSS / _LocalCRSS[SlipSystem] ), MaxRSS );
        }
        if( MaxRSS < 1e-10 )
        {
          std::cerr << "Exit at MaxRSS < 1e-10, GetMaxRSS " << std::endl;
          exit(0);
        }
        return MaxRSS;
      }
      
      //--------------------------------------------------------------------------------
      //  RunVPFFT
      //--------------------------------------------------------------------------------
      void RunVPFFT( const EigenRep & MacroscopicStrainRate,
                     int NumMaxIterations,
                     const int NumTimeSteps,
                     const Float TimeStep,
                     const Float ConvergenceEpsilon);        
      
    };
  }
} // namespace VPFFT

#endif
