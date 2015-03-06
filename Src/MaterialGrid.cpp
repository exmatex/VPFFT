#include "MaterialGrid.h"
#include <fstream>
#include <sstream>
#include <iomanip>
namespace VPFFT
{
  namespace DataStructures
  {
    using LinearAlgebra::SMatrix5x5;
    
    
    //--------------------------------------------------------------------------------
    //
    //--------------------------------------------------------------------------------
    void MaterialGrid::PrintStrain( const std::string & filename )
    {
      std::ofstream Outfile( filename.c_str() );
      for( int k = 0; k < DimZ; k ++ )
      {
        for( int j = 0; j < DimY; j ++ )
        {
          for( int i = 0; i < DimX; i ++ )
          {
            Outfile << i << " " << j << " " << k << " ";
            for( int n = 0; n < 5; n ++ )
            {
              Outfile << DisplacementVariation(i, j, k, n)[0] << " "
                      << DisplacementVariation(i, j, k, n)[1] << " ";
            }
            Outfile<< std::endl;
          }
        }
      }
      Outfile.close();
    }

    //--------------------------------------------------------------------------------
    //
    //--------------------------------------------------------------------------------
    void MaterialGrid::PrintStress( const std::string & filename )
    {
      std::ofstream Outfile( filename.c_str() );
      for( int k = 0; k < DimZ; k ++ )
      {
        for( int j = 0; j < DimY; j ++ )
        {
          for( int i = 0; i < DimX; i ++ )
          {
            Outfile << i << " " << j << " " << k << " "
                    << Stress(i, j, k) << std::endl;
          }
        }
      }
      Outfile.close();
    }

    //--------------------------------------------------------------------------------
    //
    //--------------------------------------------------------------------------------
    void MaterialGrid::PrintLagrangeMultiplier( const std::string & filename )
    {
      std::ofstream Outfile( filename.c_str() );
      for( int i = 0; i < DimX; i ++ )
      {
        for( int j = 0; j < DimY; j ++ )
        {
          for( int k = 0; k < DimZ; k ++ )
          {
            Outfile << i << " " << j << " " << k << " "
                    << LagrangeMultiplier(i, j, k) << std::endl;
          }
        }
      }
      Outfile.close();
    }

    //--------------------------------------------------------------------------------
    //
    //--------------------------------------------------------------------------------
    void MaterialGrid::PrintPolarization( const std::string & filename )
    {
      std::ofstream Outfile( filename.c_str() );
      for( int k = 0; k < DimZ; k ++ )
      {
        for( int j = 0; j < DimY; j ++ )
        {
          for( int i = 0; i < DimX; i ++ )
          {
            Outfile << i << " " << j << " " << k << " ";
            for( int n = 0; n < 5; n ++ )
            {
              Outfile << PolarizationField(i, j, k, n)[0] << " "
                      << PolarizationField(i, j, k, n)[1] << " ";
            }
            Outfile<< std::endl;
          }
        }
      }
      Outfile.close();
    }
    
    //--------------------------------------------------------------------------------
    //  UpdateSchmidtTensors
    //
    //
    //--------------------------------------------------------------------------------
    void  MaterialGrid::UpdateSchmidtTensors()
    {
      // get global basis
      const VPFFT
        ::FCC_CrystalTest
        ::FCC_SchmidtBasis & FCC_Schmidt = VPFFT::FCC_CrystalTest::FCC_SchmidtBasis::Get();
      
#pragma omp parallel for					
      for( int i = 0; i < LocalSchmidtTensor.size(); i ++ )
      {        
        SMatrix3x3 OrientationTranspose = OrientationField[i];
        OrientationTranspose.Transpose();
        for( int j = 0; j < LocalSchmidtTensor[i].size(); j ++ )
        {
          LocalSchmidtTensor[i][j] = EigenRep(OrientationField[i] * FCC_Schmidt( j ) * OrientationTranspose ); // THIS IS VALIDATED
        }
      }
    }

    //--------------------------------------------------------------------------------
    //  SetPolarizationField
    //
    //--------------------------------------------------------------------------------
    void MaterialGrid::SetPolarizationField( int x, int y, int z, const EigenRep & s )
    {
      for( int i = 0; i < 5; i ++ )
      {
        PolarizationField( x, y, z, i )[0] = s(i);
        PolarizationField( x, y, z, i )[1] = 0;
      }
    }
    
    //--------------------------------------------------------------------------------
    //  AddPolarizationField
    //
    //--------------------------------------------------------------------------------
    void MaterialGrid::AddPolarizationField( int x, int y, int z, const EigenRep & s )
    {
      for( int i = 0; i < 5; i ++ )
        PolarizationField( x, y, z, i )[0] += s(i);  // no imaginary part
    }

    //--------------------------------------------------------------------------------
    //  InitializeComputation
    //--------------------------------------------------------------------------------
    void MaterialGrid::InitializeComputation( const EigenRep & MacroscopicStrainRate )
    {
      //----------------------
      //  TO DEBUG - must be unified
      const Float NewtonIterationEpsilon = 5e-4;
      const int   MaxNewtonIterations    = 1000;

      int MinRateSensitivity    = *std::min_element( RateSensitivity.begin(), RateSensitivity.end() );
      Float ImposedMaxResolvedStress   = static_cast<Float>(2) * pow( InnerProduct( MacroscopicStrainRate, MacroscopicStrainRate ) / GammaDotBase[0],
                                                                      static_cast<Float>(1. / float(MinRateSensitivity) ) );

      ImposedMaxResolvedStress = std::max( static_cast<Float>(2.0), ImposedMaxResolvedStress );

#pragma omp parallel for                      
      for( int i = 0; i < DimX; i ++ )
      {
        for( int j = 0; j < DimY; j ++ )
        {
          for( int k = 0; k < DimZ; k ++ )
          {
            // only need this for the first time.
            EigenRep InitialGuess = Solvers::SolveSachEquation( SchmidtTensorList(i, j, k),
                                                                CRSS(i, j, k),
                                                                MacroscopicStrainRate );
            
            Float MaxResolvedStress = GetMaxRSS( InitialGuess,
                                                 SchmidtTensorList(i, j, k),
                                                 CRSS(i, j, k) );
            
            if( MaxResolvedStress > 1.1 )   // no more than 10% past YS
              InitialGuess /= MaxResolvedStress;
            
            if( MaxResolvedStress < 1e-10 )
            {
              std::cout << "Max Resolved Stress < 1e-10 " << std::endl;
              exit(0);
            }

            Stress(i, j, k) = Solvers::SolveConstitutiveEquations( InitialGuess,
                                                                   SchmidtTensorList(i, j, k),
                                                                   CRSS(i, j, k),
                                                                   GammaDotBase,
                                                                   RateSensitivity,
                                                                   MacroscopicStrainRate,
                                                                   NewtonIterationEpsilon,
                                                                   ImposedMaxResolvedStress,
                                                                   MaxNewtonIterations );

            LagrangeMultiplier(i, j, k) = Stress(i, j, k);
            SetPolarizationField( i, j, k, Stress(i, j, k) );   

          }
        }
      }
      //      std::cout << "Finisiehd Initialization [ press enter to continue ] " << std::endl;
      //      std::cin.get();
    }
    
    //--------------------------------------------------------------------------------
    //  BuildLinearRefMediumStiffness
    //
    //    Should consider some form of artificial normalization so that the sum
    //    will end up around order unity....  the error as a function from number
    //    of element isn't exactly that awesome...
    //--------------------------------------------------------------------------------
    SMatrix5x5 MaterialGrid::BuildLinearRefMediumStiffness(   )
    {
      SMatrix5x5 M;         // Linear Modouli
            
      M.SetZero();
      int NumDataPoints = static_cast<int>( LocalSchmidtTensor.size() );
#pragma omp parallel for					
      for( int i = 0; i < LocalSchmidtTensor.size(); i ++ )
      {
        SMatrix5x5 M_Local;
        M_Local.SetZero();
        for( int SlipSystem = 0; SlipSystem < LocalSchmidtTensor[i].size(); SlipSystem ++ )
        {
          Float GammaDot;
          Float RSS = InnerProduct( LocalSchmidtTensor[i][ SlipSystem ],
                                    StressField[ i ] ) / LocalCRSS[i][ SlipSystem ];
          if( RateSensitivity[ SlipSystem ] - 1 > 0 )
            GammaDot = GammaDotBase[ SlipSystem ]
              *  std::pow( std::fabs(RSS), static_cast<int>( RateSensitivity[ SlipSystem ] - 1 ) );
          else
            GammaDot = GammaDotBase[ SlipSystem ];         
          
          M_Local += ( OuterProduct( LocalSchmidtTensor[i][ SlipSystem ],
                                     LocalSchmidtTensor[i][ SlipSystem ] ) * GammaDot / LocalCRSS[i][ SlipSystem ] );
        }
        
        using namespace Eigen;
        typedef Matrix< Float, 5, 5 > Matrix5x5;
        Map< Matrix5x5 > M_Mapped( M_Local.m );
        M_Mapped = M_Mapped.inverse();
        M.MovingAverage( M_Local, i );
      }
      
      M.NormalizeValues();      
      return M;
    }

    //--------------------------------------------------------------------------------
    //  Cleanup
    //--------------------------------------------------------------------------------
    void MaterialGrid::Cleanup()
    {
      if( bInitialized )
      {
        for( int i = 0; i < 5; i ++ )
        {
          fftw_destroy_plan( _DisplacementBackwardPlan[i] );
          fftw_destroy_plan( _PolarizationFieldForwardPlan[i] );

          fftw_free( _DisplacementVariation[i] );          
          fftw_free( _PolarizationField[i] );
        }
      }
      fftw_cleanup_threads();
      // fftw_cleanup();
    }
    
    //--------------------------------------------------------------------------------
    //
    //  InitializeMemory
    //    Takes care of all memory allocation only
    //--------------------------------------------------------------------------------
    void MaterialGrid::InitializeMemory( int DimX_, int DimY_, int DimZ_ ) 
    {
      //-------------------------------------------------------
      //  FFTW OpenMP initialization
      //-------------------------------------------------------
      fftw_init_threads( );
      std::cout << "FFTW setup with " << NumThreads << " threads " << std::endl;
      fftw_plan_with_nthreads( NumThreads );
      //-------------------------------------------------------
      //
      //               D E B U G
      //
      // to be incorperated into the class... 
      //-------------------------------------------------------
      GammaDotBase    = std::vector<Float>(FCC_CrystalTest::NumSystems, 15.5);   // was 1 when debuggin
      RateSensitivity = std::vector<int>(FCC_CrystalTest::NumSystems, 10);
      //-------------------------------------------------------

      DimYZ = DimY * DimZ;
      
      int NumElements = DimX * DimY * DimZ;

      OrientationField.resize       ( NumElements );
      StressField.resize            ( NumElements );
      StrainRateField.resize        ( NumElements );
      LagrangeMultiplierField.resize( NumElements );
      LocalSchmidtTensor.resize     ( NumElements );
      LocalCRSS.resize              ( NumElements );
      AccumulatedShear.resize       ( NumElements );

      using FCC_CrystalTest::FCC_SchmidtBasis;
      //-------------------------------------------------------
      //   Hardening parameters, DEBUG values, even for hardening matrix
      Tau0.resize  ( FCC_SchmidtBasis::NumSlipSystems(), 11.0 );
      Tau1.resize  ( FCC_SchmidtBasis::NumSlipSystems(), 11.0 );
      Theta0.resize( FCC_SchmidtBasis::NumSlipSystems(), 430.0);
      Theta1.resize( FCC_SchmidtBasis::NumSlipSystems(), 110.0 );
      
      HardeningMatrix.resize( FCC_SchmidtBasis::NumSlipSystems() );
      for( int i = 0; i < FCC_SchmidtBasis::NumSlipSystems(); i ++ )
        HardeningMatrix[i].resize( FCC_SchmidtBasis::NumSlipSystems(), 1.0 );
      //-------------------------------------------------------
      
   
      for( int i = 0; i < NumElements; i ++ )
      {
        LocalSchmidtTensor[i].resize( FCC_SchmidtBasis::NumSlipSystems() );
        LocalCRSS[i] = vector<Float>( FCC_SchmidtBasis::NumSlipSystems(), 11 ); // hard coded
      }
      
      for( int i = 0; i < 5; i ++ )
      {
        _DisplacementVariation[i]
          = static_cast<FFTW_ComplexPtr>( fftw_malloc( sizeof(FFTW_Complex) * NumElements ) );
        
        
        _DisplacementBackwardPlan[i]  = fftw_plan_dft_3d( DimX, DimY, DimZ,
                                                          _DisplacementVariation[i],
                                                          _DisplacementVariation[i],
                                                          FFTW_BACKWARD,
                                                          FFTW_MEASURE | FFTW_DESTROY_INPUT );
        _PolarizationField[i]
          = static_cast<FFTW_ComplexPtr>( fftw_malloc( sizeof(FFTW_Complex) * NumElements ) );

        
        bzero( _PolarizationField[i], sizeof( FFTW_Complex ) * NumElements );


        _PolarizationFieldForwardPlan[i]  = fftw_plan_dft_3d( DimX, DimY, DimZ,
                                                              _PolarizationField[i],
                                                              _PolarizationField[i],
                                                              FFTW_FORWARD,
                                                              FFTW_MEASURE | FFTW_DESTROY_INPUT  );
      }
      
      bInitialized = true;
    }
 
    //--------------------------------------------------------------------------------
    //  ConstructGreensOperator
    //--------------------------------------------------------------------------------
    LinearAlgebra::Tensor4Rank MaterialGrid::ConstructGreensOperator( const SMatrix5x5 & LinearRefMediumStiffness,
                                                                      const SVector3 & WaveVector ) const
    {
      using Eigen::Matrix;
      using VPFFT::LinearAlgebra::Tensor4Rank;
      typedef Matrix< Float, 4, 4 > Matrix4x4;
      Matrix4x4 AcousticMatrix = Matrix4x4::Zero();

      Tensor4Rank L = LinearRefMediumStiffness.ToTensor4RankRep( );


      for( int k = 0; k < 3; k ++ )
        for( int l = 0; l < 3; l ++ )
          for( int i = 0; i < 3; i ++ )
            for( int j = 0; j < 3; j ++ )
              AcousticMatrix( i, k ) += WaveVector[l] * WaveVector[j] * L(i, j, k, l);
      
      for( int i = 0; i < 3; i ++ )
        AcousticMatrix(3, i) = AcousticMatrix(i, 3) = WaveVector[i];

      AcousticMatrix(3, 3) = 0;
      AcousticMatrix       = AcousticMatrix.inverse().eval();
      
      
      Tensor4Rank GreensOperator4;
      for( int i = 0; i < 3; i ++ )
        for( int j = 0; j < 3; j ++ )
          for( int k = 0; k < 3; k ++ )
            for( int l = 0; l < 3; l ++ )
              GreensOperator4( i, j, k, l ) = - WaveVector[j] * WaveVector[l] * AcousticMatrix(i, k);

      return GreensOperator4;
    }
    
    //--------------------------------------------------------------------------------
    //   StrainRateEvolutionWithFFT
    //
    //   WARNING!  Normalization must be applied for FFTW results
    //
    //   DISTANCE UNIT IS 1 right now - need scaling.
    //--------------------------------------------------------------------------------
    void MaterialGrid::StrainRateEvolutionWithFFT( const SMatrix5x5 & LinearRefMediumStiffness )
    {
      using LinearAlgebra::Tensor4Rank;

      SMatrix5x5 L_Inv = -LinearRefMediumStiffness;
      using namespace Eigen;
      typedef Matrix< Float, 5, 5 > Matrix5x5;
      Map< Matrix5x5 > M_Mapped( L_Inv.m );
      M_Mapped = M_Mapped.inverse();
      
      
      
      std::stringstream Suffix;
      Suffix << "." << DEBUG_INT;
      //      PrintLagrangeMultiplier( "LagrangeMultiplier.txt" + Suffix.str() );
      // PrintStress("Stress.txt"  + Suffix.str() );

      
      //      PrintPolarization( "PolarizationBefore.txt" + Suffix.str() );
      for( int i = 0; i < 5; i ++ )   
        fftw_execute( _PolarizationFieldForwardPlan[i] );
      // PrintPolarization( "PolarizationAfter.txt"  + Suffix.str() );

      int CutoffX = floor( DimX / 2 );
      int CutoffY = floor( DimY / 2 );
      int CutoffZ = floor( DimZ / 2 );

      bool IsEvenX = (DimX % 2 == 0 );
      bool IsEvenY = (DimY % 2 == 0 );
      bool IsEvenZ = (DimZ % 2 == 0 );
      
      const Float Two_PI = static_cast<Float>( 2 ) * PI;

      //---------------------------------------------
      // for all wave vectors -
      //
      //  WARNING!  This is specified for FFTW3 output
      //            only.  The definition of the output
      //            ordering from the DFT ran using FFTW
      //            can be found in section 4.8 of FFTW3.3.1
      //            document (pdf format, on fftw website.)
      //   
      //            This also means that the first "half" of
      //            each dimension is the normal wave vectors,
      //            and the second half is the "negatives"
      //
      //  NOTE:     Current implementation is designed for
      //            the normal DFT which does *NOT* assume
      //            the initial input to be the special case
      //            of being purely real instead of complex.
      //            Consequently, we are not taking advantage
      //            of the symmetric nature of the FFT output,
      //            which results in a factor of 2 increase
      //            in computation.  Note that an *adaptation*
      //            of this function is necessary to take
      //            advantage of the FFT special case.
      //
      //  WARNING:  FFTW computes *UNNORMALIZED* DFTs.
      //            Normalization step must e taken at some
      //            point.
      //
      //----------------------------------------------
      int DEBUG_NUM_CALLED = 0;
      for( int Nx = 0; Nx < DimX; Nx ++ )
      {
        for( int Ny = 0; Ny < DimY; Ny ++ )
        {
          for( int Nz = 0; Nz < DimZ; Nz ++ )
          {
            if( Nx != 0 || Ny != 0 || Nz != 0)
            {
              
              SVector3 WaveVector( Nx, Ny, Nz );
              
              if( Nx > CutoffX )
                WaveVector[0] -=  DimX;
              if( Ny > CutoffY )
                WaveVector[1] -=  DimY;
              if( Nz > CutoffZ )
                WaveVector[2] -=  DimZ;
              
              
              WaveVector[0] /= (DimX * XVoxelLength );
              WaveVector[1] /= (DimY * YVoxelLength );  // constant factor in front of the K vector doesn't matter - cancelled out in the
              WaveVector[2] /= (DimZ * ZVoxelLength );  // construction of green's operator  (that's why normalization doesn't change anything)
              WaveVector.Normalize();   
              
              Tensor4Rank GreensOperator4;
              SMatrix5x5  SymGreensOperator;
              
              
              if( (   IsEvenX && (Nx == CutoffX || Nx == CutoffX + 1) )
                  || (IsEvenY && (Ny == CutoffY || Ny == CutoffY + 1) )
                  || (IsEvenZ && (Nz == CutoffZ || Nz == CutoffZ + 1) ) )
                GreensOperator4   = L_Inv.ToTensor4RankRep();
              else
                GreensOperator4   = ConstructGreensOperator( LinearRefMediumStiffness, WaveVector );
            
              SymGreensOperator = SMatrix5x5( LinearAlgebra::Symmetrize( GreensOperator4 ) );  // this is the correct one, debugging
              
              // actually do the work
              for( int i = 0; i < 5; i ++ )
              {
                DisplacementVariation( Nx, Ny, Nz, i )[0] = 0;
                DisplacementVariation( Nx, Ny, Nz, i )[1] = 0;
                for( int j = 0; j < 5; j ++ )
                {
                  DisplacementVariation( Nx, Ny, Nz, i )[0] -= SymGreensOperator( i, j )
                                                             * PolarizationField( Nx, Ny, Nz, j )[0];
                  DisplacementVariation( Nx, Ny, Nz, i )[1] -= SymGreensOperator( i, j )
                                                             * PolarizationField( Nx, Ny, Nz, j )[1];
                }
              }
            }
            else   // Worry about the branching cost when we optimize - this is an example.
            {
              for( int i = 0; i < 5; i ++ )
              {
                DisplacementVariation( 0, 0, 0, i )[0] = 0;
                DisplacementVariation( 0, 0, 0, i )[1] = 0;
              }
            } // WaveNumber != (0, 0, 0)
          } // DimZ
        } // DimY
      } // DimX

      for(int i = 0; i < 5; i ++ )
        fftw_execute( _DisplacementBackwardPlan[i] );

      Float NormalizationFactor = static_cast<Float>( 1 ) / (DimX * DimY * DimZ);
      // All displacement will be real
      for( int i = 0; i < 5; i ++ )
        for( int j = 0; j < DimX * DimY * DimZ; j ++ )
          _DisplacementVariation[i][j][0] *= NormalizationFactor;
        
      
      //--------------------
      //   DEBUG
      //
      DEBUG_INT ++;
      //--------------------
      
    }  // end StrainRateEvolutionWithFFT

    //--------------------------------------------------------------------------------
    //  OrientationEvolutionWithFFT
    //   - Exactly the same as StrainRateEvolutionWithFFT, but uses an anti-symmetrized
    //     Greens' operator.  Unfortunately, displacement variation will have to be
    //     used to store the orientation variation ( \delta \Omega )
    //--------------------------------------------------------------------------------
    void MaterialGrid::OrientationEvolutionWithFFT( const SMatrix5x5 & LinearRefMediumStiffness,
                                                    const Float TimeStep )
    {
      using LinearAlgebra::Tensor4Rank;

      SMatrix5x5 L_Inv = -LinearRefMediumStiffness;
      using namespace Eigen;
      typedef Matrix< Float, 5, 5 > Matrix5x5;
      Map< Matrix5x5 > M_Mapped( L_Inv.m );
      M_Mapped = M_Mapped.inverse();
      for( int i = 0; i < 5; i ++ )
      {
        fftw_execute( _PolarizationFieldForwardPlan[i] );
      }
      
      int CutoffX = floor( DimX / 2 );
      int CutoffY = floor( DimY / 2 );
      int CutoffZ = floor( DimZ / 2 );

      bool IsEvenX = (DimX % 2 == 0 );
      bool IsEvenY = (DimY % 2 == 0 );
      bool IsEvenZ = (DimZ % 2 == 0 );
      
      const Float Two_PI = static_cast<Float>( 2 ) * PI;

      for( int Nx = 0; Nx < DimX; Nx ++ )
      {
        for( int Ny = 0; Ny < DimY; Ny ++ )
        {
          for( int Nz = 0; Nz < DimZ; Nz ++ )
          {
            if( Nx != 0 || Ny != 0 || Nz != 0)
            {
              SVector3 WaveVector( Nx, Ny, Nz );
              
              if( Nx > CutoffX )
                WaveVector[0] -=  DimX;
              if( Ny > CutoffY )
                WaveVector[1] -=  DimY;
              if( Nz > CutoffZ )
                WaveVector[2] -=  DimZ;


              WaveVector[0] /= (DimX * XVoxelLength );
              WaveVector[1] /= (DimY * YVoxelLength );  // constant factor in front of the K vector doesn't matter - cancelled out in the
              WaveVector[2] /= (DimZ * ZVoxelLength );  // construction of green's operator  (that's why normalization doesn't change anything)
              
              WaveVector.Normalize();   
              
              Tensor4Rank GreensOperator4;
              SMatrix5x5  SymGreensOperator;
                            
              if( (   IsEvenX && (Nx == CutoffX || Nx == CutoffX + 1) )
                  || (IsEvenY && (Ny == CutoffY || Ny == CutoffY + 1) )
                  || (IsEvenZ && (Nz == CutoffZ || Nz == CutoffZ + 1) ) )
                GreensOperator4   = L_Inv.ToTensor4RankRep();
              else
                GreensOperator4   = ConstructGreensOperator( LinearRefMediumStiffness, WaveVector );
            
              SymGreensOperator = SMatrix5x5(  GreensOperator4  );  // Not symmetrizing or anti-symmetrizing, waiting till the end
              
              // actually do the work
              for( int i = 0; i < 5; i ++ )
              {
                DisplacementVariation( Nx, Ny, Nz, i )[0] = 0;
                DisplacementVariation( Nx, Ny, Nz, i )[1] = 0;
                for( int j = 0; j < 5; j ++ )
                {
                  DisplacementVariation( Nx, Ny, Nz, i )[0] -= SymGreensOperator( i, j )
                                                             * PolarizationField( Nx, Ny, Nz, j )[0];
                  DisplacementVariation( Nx, Ny, Nz, i )[1] -= SymGreensOperator( i, j )
                                                             * PolarizationField( Nx, Ny, Nz, j )[1];
                }
              }
            }
            else   
            {
              for( int i = 0; i < 5; i ++ )
              {
                DisplacementVariation( 0, 0, 0, i )[0] = 0;
                DisplacementVariation( 0, 0, 0, i )[1] = 0;
              }
            } // WaveNumber != (0, 0, 0)
          } // DimZ
        } // DimY
      } // DimX

      for(int i = 0; i < 5; i ++ )
        fftw_execute( _DisplacementBackwardPlan[i] );

      Float NormalizationFactor = static_cast<Float>( 1 ) / (DimX * DimY * DimZ);
      // All displacement will be real
      for( int i = 0; i < 5; i ++ )
        for( int j = 0; j < DimX * DimY * DimZ; j ++ )
          _DisplacementVariation[i][j][0] *= NormalizationFactor;

      // The Displacement Variation is now the orientation 
      //  note that we may change the ordering so that we're running in single index instead of the 3-tuples.
      for( int Nx = 0; Nx < DimX; Nx ++ )
        for( int Ny = 0; Ny < DimY; Ny ++ )
          for( int Nz = 0; Nz < DimZ; Nz ++ )
          {
            Orientation( Nx, Ny, Nz ) += TimeStep * ( ( DisplacementVariation( Nx, Ny, Nz ).ToMatrixRep() ).AntiSymmetrize()
                                                      + ShearRate( Nx, Ny, Nz ) );
            
          }
      UpdateSchmidtTensors();
      
    }  // end EvoleOrientation

    
    //--------------------------------------------------------------------------------
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
    //--------------------------------------------------------------------------------
    EigenRep MaterialGrid::RunSingleStrainStep( const EigenRep & MacroscopicStrainRate,
                                                int NumMaxIterations, Float TimeStep,
                                                Float ConvergenceEpsilon )
    {
      InitializeComputation( MacroscopicStrainRate );
      
      //----------------------
      //  TO DEBUG - must be unified
      const Float NewtonIterationEpsilon = 5e-4;
      const Float InputMaxResolvedStress = 1e9; /// from Richardo's code
      const int   MaxNewtonIterations    = 100;
      //----------------------
      SMatrix5x5 LinearRefMediumStiffness = BuildLinearRefMediumStiffness();

      int   StepsTaken = 0;
      Float StrainRateErrorNorm = 0;
      Float StressErrorNorm     = 0;
      EigenRep TotalAveragedStress;
      
//-------------------- SETUP OpenMP
      //int NumThreads = 8;

      omp_set_num_threads( NumThreads );
      std::cout << "NumThreads " << NumThreads << std::endl;
      
      DEBUG_INT = 0;
      do
      {
        StrainRateErrorNorm = 0;
        StressErrorNorm     = 0;

        Float AverageVariationNorm = 0;
        Float AverageStressNorm    = 0;
        Float Max_NR_Error         = 0;
        //--------------------------------------
        //  for non-OpenMP
        //         EigenRep AveragedStress     ( 0, 0, 0, 0, 0 );
        //         EigenRep AveragedStrainRate ( 0, 0, 0, 0, 0 );
        //         EigenRep AveragedStrainRateVariation ( 0, 0, 0, 0, 0 );   // loop invariation = should be zero
        
 
        vector<EigenRep> AveragedStress             ( NumThreads );
        vector<EigenRep> AveragedStrainRate         ( NumThreads );
        vector<EigenRep> AveragedStrainRateVariation( NumThreads );   // loop invariation = should be zero
        for( int TId = 0; TId < NumThreads; TId ++ )
        {
          AveragedStress[TId] = EigenRep(0, 0, 0, 0, 0);
          AveragedStrainRate[TId] = EigenRep(0, 0, 0, 0, 0);
          AveragedStrainRateVariation[TId] = EigenRep(0, 0, 0, 0, 0);
        }
        
        StrainRateEvolutionWithFFT( LinearRefMediumStiffness );
        TotalAveragedStress = EigenRep(0, 0, 0, 0, 0);
        
#pragma omp parallel for                 \
  reduction( + : AverageVariationNorm )         \
  reduction( + : AverageStressNorm )            \
  reduction( + : Max_NR_Error )                 \
  reduction( + : StressErrorNorm )               \
  reduction( + : StrainRateErrorNorm )                 
#pragma omp nowait
        for( int i = 0; i < DimX; i ++ )
        {
//           if(i == 0)
//           {
//             std::cout << "Num Threads " << omp_get_num_threads() << std::endl;
//           }
          for( int j = 0; j < DimY; j ++ )
          {
            for( int k = 0; k < DimZ; k ++ )
            {
              EigenRep LocalStrainRateVariation = DisplacementVariation( i, j, k );                
              Float    NR_Error = 0;
              Float    MaxResolvedStress = GetMaxRSS( Stress(i, j, k),
                                                      SchmidtTensorList(i, j, k),
                                                      CRSS(i, j, k) );
              
              Stress(i, j, k) /= MaxResolvedStress;
              Stress(i, j, k)
                = Solvers::SolveConstrainedConstitutiveEquations( Stress(i, j, k),
                                                                  SchmidtTensorList(i, j, k),
                                                                  CRSS(i, j, k),
                                                                  GammaDotBase,
                                                                  RateSensitivity,
                                                                  MacroscopicStrainRate,
                                                                  LagrangeMultiplier(i, j, k),
                                                                  LocalStrainRateVariation,
                                                                  LinearRefMediumStiffness,
                                                                  NewtonIterationEpsilon,
                                                                  InputMaxResolvedStress,
                                                                  MaxNewtonIterations,
                                                                  &NR_Error );
               Max_NR_Error = std::max( NR_Error, Max_NR_Error );
               EigenRep LocalStrainRate 
                 = Solvers::ApplyConstitutiveEquations( Stress(i, j, k),
                                                        SchmidtTensorList(i, j, k),
                                                        CRSS(i, j, k), GammaDotBase, RateSensitivity );
              SetPolarizationField( i, j, k,
                                    Stress(i, j, k)
                                    -  LinearRefMediumStiffness * ( DisplacementVariation(i, j, k)  )   );
              UpdateLagrangeMultiplier( i, j, k, LinearRefMediumStiffness, LocalStrainRate, MacroscopicStrainRate, 1 );
              
              int ThreadID = omp_get_thread_num();
              
              // AveragedStress          += Stress(i, j, k);
//               AveragedStrainRate      += LocalStrainRate;

          
              AveragedStress[ThreadID]          += Stress(i, j, k);
              AveragedStrainRate[ThreadID]      += LocalStrainRate;

          
              EigenRep StrainRateError = LocalStrainRate - MacroscopicStrainRate - LocalStrainRateVariation ;  // \dot{epsilon}(x) - d(x); d(x) = \tilde(d) + E-dot
              EigenRep StressError     = Stress(i, j, k) - LagrangeMultiplier(i, j, k);
              StressErrorNorm      += std::sqrt( LinearAlgebra::InnerProduct( StressError,     StressError ) );
              StrainRateErrorNorm  += std::sqrt( LinearAlgebra::InnerProduct( StrainRateError, StrainRateError ) );
            }
          }
        }
        
        StressErrorNorm      /= static_cast<Float>( DimX * DimY * DimZ );
        StrainRateErrorNorm  /= static_cast<Float>( DimX * DimY * DimZ );
        
        EigenRep TotalAveragedStrainRate(0, 0, 0, 0, 0);

        for( int TId = 0; TId < NumThreads; TId ++ )
        {
          TotalAveragedStrainRate += AveragedStrainRate[TId];
          TotalAveragedStress += AveragedStress[TId];
        }
        
        TotalAveragedStress       /= static_cast<Float>( DimX * DimY * DimZ );  
        TotalAveragedStrainRate   /= static_cast<Float>( DimX * DimY * DimZ );
        
        StressErrorNorm     /= std::sqrt( LinearAlgebra::InnerProduct( TotalAveragedStress,     TotalAveragedStress) ) ;
        StrainRateErrorNorm /= std::sqrt( LinearAlgebra::InnerProduct( MacroscopicStrainRate, MacroscopicStrainRate ) );
        
        StepsTaken ++;
        // need convergence averages
        //        std::cout << " Max NR Error ============================================================== " << Max_NR_Error << std::endl;
   //      std::cout << " Averaged Stress           " << std::setw(10) << AveragedStress     << " || " << StressErrorNorm << std::endl;
//         std::cout << " Averaged Strain           " << std::setw(10) << AveragedStrainRate << " || " << StrainRateErrorNorm << std::endl;
      } while( StepsTaken < NumMaxIterations
               && ( StressErrorNorm > ConvergenceEpsilon
                    || StrainRateErrorNorm > ConvergenceEpsilon ) );

      OrientationEvolutionWithFFT( LinearRefMediumStiffness, TimeStep );  // move L out, and we can move this out
      UpdateVoxelLength  ( TimeStep, MacroscopicStrainRate.ToMatrixRep() );
      UpdateHardeningCRSS( TimeStep );
      
      
      return TotalAveragedStress;
    }
    
    //--------------------------------------------------------------------------------
    //  CalculateShearRate
    //
    //  \dot{omega}^{p}_{ij} = \sum_k \alpha_{ij}^{k}( \vec{x} ) \dot{\gamma}^k{\vec{x}}
    //
    //  \alpha_{ij}^{k} -  Anti-symmetric schmidt tensor for the kth system.
    //  \dot{\gamma} is the usual \gamma_0 \frac{(\m^{k} \sigma')^m}{\Tau}
    //
    //   Note that GammaDot can be precomputed here, then passed into CalculateShearRate 
    //   and CUpdateHardeningCRSS to reduce the number of computations.
    //
    //--------------------------------------------------------------------------------
    SMatrix3x3 MaterialGrid::ShearRate( int i, int j, int k )
    {
      // get global basis
      const VPFFT
        ::FCC_CrystalTest
        ::FCC_SchmidtBasis & FCC_Schmidt = VPFFT::FCC_CrystalTest::FCC_SchmidtBasis::Get();
      
      vector<SMatrix3x3> LocalAntiSchmidt( FCC_Schmidt.NumSlipSystems() );
      for( int n = 0; n < LocalAntiSchmidt.size(); n ++ )
      {        
        SMatrix3x3 OrientationTranspose = Orientation(i, j, k);
        OrientationTranspose.Transpose();
        LocalAntiSchmidt[n] = Orientation(i, j, k) * FCC_Schmidt.Anti( n ) * OrientationTranspose; // THIS IS VALIDATED
      }
      
      return Solvers::CalculateShearRate( Stress(i, j, k),
                                          SchmidtTensorList(i, j, k),
                                          LocalAntiSchmidt,
                                          CRSS(i, j, k),
                                          GammaDotBase, RateSensitivity );
      
    }

    //--------------------------------------------------------------------------------
    //  UpdateHardeningCRSS
    //
    //  Purpose:  Update the hardening models.  Can be used as a way 
    //            to plug into high resolution, multi-physics simulation.
    //--------------------------------------------------------------------------------
    void MaterialGrid::UpdateHardeningCRSS( Float TimeStep )
    {
      for( int i = 0; i < AccumulatedShear.size(); i ++ )
      {
/*
        if( i % 250 == 0 )
        {
          std::cout << "Before Hardening ------------------- " << std::endl;

          for( int n = 0; n < LocalCRSS[i].size(); n ++ )
            std::cout << LocalCRSS[i][n]   << " |>" << AccumulatedShear[i] << " <| ";
          std::cout << std::endl;
        }
*/
        LocalCRSS[i] = Solvers::CalculateHardening( AccumulatedShear[i], TimeStep,
                                                    StressField[i], 
                                                    LocalSchmidtTensor[i],
                                                    LocalCRSS[i],
                                                    GammaDotBase,
                                                    RateSensitivity,
                                                    HardeningMatrix,
                                                    Tau0, Tau1, Theta0, Theta1  );
/*        
        if( i % 250 == 0 )
        {
          std::cout << "After Hardening " << std::endl;

          for( int n = 0; n < LocalCRSS[i].size(); n ++ )
            std::cout << LocalCRSS[i][n]  << " |>" << AccumulatedShear[i] << " <| ";
          std::cout << "============================="  << std::endl;
        }
*/
      }
    }
    
    //--------------------------------------------------------------------------------
    //  UpdateLagrangeMultipler
    //
    //     LocalStrainRateVariation used here comes from the difference of
    //     the macroscopic strain and application of the constitutive equation
    //     on the new stress state.
    //
    //     This is clearly one of those functions that could be inlined and vectorized
    //     using SIMD. 
    //--------------------------------------------------------------------------------
    void MaterialGrid::UpdateLagrangeMultiplier( int i, int j, int k,
                                                 const SMatrix5x5 & LinearRefStiffness,
                                                 const EigenRep & LocalStrainRate,
                                                 const EigenRep & MacroscopicStrainRate,
                                                 Float FudgeFactor )
    {
      EigenRep T(0, 0, 0, 0, 0);
      //  only the real component of DisplacementVariation should contain values...
      T(0) = MacroscopicStrainRate(0) + DisplacementVariation( i, j, k, 0 )[0] - LocalStrainRate(0);
      T(1) = MacroscopicStrainRate(1) + DisplacementVariation( i, j, k, 1 )[0] - LocalStrainRate(1);
      T(2) = MacroscopicStrainRate(2) + DisplacementVariation( i, j, k, 2 )[0] - LocalStrainRate(2);
      T(3) = MacroscopicStrainRate(3) + DisplacementVariation( i, j, k, 3 )[0] - LocalStrainRate(3);
      T(4) = MacroscopicStrainRate(4) + DisplacementVariation( i, j, k, 4 )[0] - LocalStrainRate(4);
      
      LagrangeMultiplier( i, j, k ) += (LinearRefStiffness * T); 
    }

    
    //--------------------------------------------------------------------------------
    //  UpdateVoxelLength  -- only isotropic deformation is allowed
    //--------------------------------------------------------------------------------
    void MaterialGrid::UpdateVoxelLength( Float TimeStep, const SMatrix3x3 & MacroscopicStrainRate )
    {
      XVoxelLength *= ( 1.0 + MacroscopicStrainRate.m[0][0] * TimeStep);
      YVoxelLength *= ( 1.0 + MacroscopicStrainRate.m[1][1] * TimeStep);
      ZVoxelLength *= ( 1.0 + MacroscopicStrainRate.m[2][2] * TimeStep);
      
    }
    
    
    //--------------------------------------------------------------------------------
    //  RunVPFFT
    //--------------------------------------------------------------------------------
    void MaterialGrid::RunVPFFT( const EigenRep & MacroscopicStrainRate,
                                 int NumMaxIterations,
                                 const int NumTimeSteps,
                                 const Float TimeStep,
                                 const Float ConvergenceEpsilon)
    {
      std::ofstream os("Stress_Strain.txt");
      EigenRep MacroscopicStrain(0, 0, 0, 0, 0);
      for ( int i = 0; i < NumTimeSteps; i ++ )
      {
        EigenRep AveragedStress = RunSingleStrainStep( MacroscopicStrainRate, NumMaxIterations, TimeStep, ConvergenceEpsilon );
        MacroscopicStrain +=  MacroscopicStrainRate * TimeStep;       // update strain
        os << i << " " << MacroscopicStrain << " "
           << AveragedStress << " "
           << std::sqrt( LinearAlgebra::InnerProduct( MacroscopicStrain, MacroscopicStrain) ) << " "
           << std::sqrt( LinearAlgebra::InnerProduct( AveragedStress, AveragedStress) ) << " "
           << std::endl;
      }
    }
  } // DataStructures
  
}// namespace VPFFT

