#include "UnitTester.h"
#include <fstream>
#include "Config.h"

namespace VPFFT
{
  namespace UnitTester
  {

    //----------------------------------------------------------------------------------
    //  PrintGrid
    //    Write grid into a format that's readable by Ricardo's code, currently hard coded
    //    to follow exactly SetMaterialGrid
    //----------------------------------------------------------------------------------
    void PrintGrid( std::ofstream & os, const VPFFT::DataStructures::MaterialGrid & Grid )
    {      
      int XHalf = std::floor(Grid.NumX() / 2) + 1;
      int YHalf = std::floor(Grid.NumY() / 2) + 1;   // change due to indexing in Fortran for Ricardo's code

      VPFFT::LinearAlgebra::SVector3 LeftAngle( -40, 0, 0 );
      VPFFT::LinearAlgebra::SVector3 RightAngle(  0, 40 , 0 );
      for( int k = 1; k <= Grid.NumZ(); k ++ )
        for( int j = 1; j <= Grid.NumY(); j ++ )
          for( int i = 1; i <= Grid.NumX(); i ++ )
          {
            if( (i >= XHalf && j >= YHalf) || (i < XHalf && j < YHalf)  )
            {
              os << LeftAngle << " " << i << " " << j << " " << k
                 << " " << 1 << " " << 1 << std::endl;
            }
            else
            {
              os << RightAngle << " " << i << " " << j << " " << k
                 << " " << 2 << " " << 1 << std::endl;
            }
          }
    }
    
    //----------------------------------------------------------------------------------
    //
    //  SetMaterialGrid
    //   DEBUG
    //----------------------------------------------------------------------------------
    void SetMaterialGrid( VPFFT::DataStructures::MaterialGrid & Grid )
    {
      VPFFT::LinearAlgebra::SMatrix3x3 Left, Right;
      Left.BuildActiveEulerMatrix (  DEGREE_TO_RADIAN( -40 ), 0, 0 );
      Right.BuildActiveEulerMatrix( DEGREE_TO_RADIAN( 0), DEGREE_TO_RADIAN(40), 0 );
      VPFFT::LinearAlgebra::EigenRep StressState( 0, 0, 0, 0, 0 );


//       Left.SetIdentity();
//       Right.SetIdentity();
      
      int XHalf = std::floor(Grid.NumX() / 2);
      int YHalf = std::floor(Grid.NumY() / 2);
      
      for( int i = 0; i < Grid.NumX(); i ++ )
        for( int j = 0; j < Grid.NumY(); j ++ )
          for( int k = 0; k < Grid.NumZ(); k ++ )
          {
            if( (i >= XHalf && j >= YHalf) || (i < XHalf && j < YHalf)  )
            {
              Grid.Orientation(i, j, k) = Left;
              //Grid.Orientation(i, j, k) = Right;
              Grid.Stress     (i, j, k) = VPFFT::LinearAlgebra::EigenRep( 1, 1, 1, 1, 1 );       // this is only for debugging - values not used
            }
            else
            {
              //Grid.Orientation(i, j, k) = Right;
              Grid.Orientation(i, j, k) = Right;
              Grid.Stress     (i, j, k) = VPFFT::LinearAlgebra::EigenRep( -1, -1, -1, -1, -1 );   // this is only for debugging
              
            }
          }

      Grid.UpdateSchmidtTensors();
    }
    
    //----------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------
    void ReadRicardoData( const std::string & oFilename, int DimX, int DimY, int DimZ,
                          VPFFT::DataStructures::MaterialGrid & Grid )
    {
      std::ifstream is;
      is.open( oFilename.c_str() );
      
      int count = DimX * DimY * DimZ;
      
      for( int i = 0; i < count; i ++ )
      {
        int x, y, z, tmp;
        VPFFT::LinearAlgebra::SVector3 EulerAngles;
        
        is >> std::skipws
           >> EulerAngles.m_fX >> std::skipws
           >> EulerAngles.m_fY >> std::skipws
           >> EulerAngles.m_fZ >> std::skipws
           >> x  >> std::skipws
           >> y  >> std::skipws
           >> z  >> std::skipws
           >> tmp  >> std::skipws
           >> tmp >> std::skipws;

        EulerAngles = VPFFT::LinearAlgebra::DegreeToRadian( EulerAngles );
        Grid.Orientation( x-1, y-1, z-1 ).BuildActiveEulerMatrix( EulerAngles.m_fX, EulerAngles.m_fY, EulerAngles.m_fZ );
      }
      is.close();
      Grid.UpdateSchmidtTensors();
    }

    //----------------------------------------------------------------------------------
    //  WriteRandomMaterial
    //----------------------------------------------------------------------------------
    void WriteRandomMaterial ( VPFFT::LinearAlgebra::SMatrix3x3 * OrientList,
			       VPFFT::LinearAlgebra::EigenRep * Stress,
			       int NumX, int NumY, int NumZ,
			       const VPFFT::DataStructures::MaterialGrid & Grid,
			       const std::string & Filename )
    {

      std::ofstream os( Filename.c_str() );
      
      for( int i = 0; i < Grid.NumX(); i ++ )
        for( int j = 0; j < Grid.NumY(); j ++ )
          for( int k = 0; k < Grid.NumZ(); k ++ )
          {
            
	    VPFFT::LinearAlgebra::SMatrix3x3 Orient = OrientList[ Grid.ToIndex( i, j, k ) ];
	    
	    for( int m = 0; m < 3; m ++ )
	      for( int n = 0; n < 3; n ++ )
		os << Orient.m[m][n] << " ";
	    os << Stress[ Grid.ToIndex( i, j, k ) ] << std::endl; 
	  }
      os.close();
    }
    
    //----------------------------------------------------------------------------------
    //  ReadRandomMaterial
    //----------------------------------------------------------------------------------
    void ReadRandomMaterial ( VPFFT::LinearAlgebra::SMatrix3x3 * OrientList,
			      VPFFT::LinearAlgebra::EigenRep * Stress,
			      int NumX, int NumY, int NumZ,
			      const VPFFT::DataStructures::MaterialGrid & Grid,
			      const std::string & Filename )
    {
      
      std::ifstream is( Filename.c_str() );
      
      for( int i = 0; i < Grid.NumX(); i ++ )
        for( int j = 0; j < Grid.NumY(); j ++ )
          for( int k = 0; k < Grid.NumZ(); k ++ )
          {
	    VPFFT::LinearAlgebra::SMatrix3x3 & Orient = OrientList[ Grid.ToIndex(i, j, k ) ];
	   
	    for( int m = 0; m < 3; m ++ )
	      for( int n = 0; n < 3; n ++ )
		is >> std::skipws >> Orient.m[m][n];
	    
	    VPFFT::LinearAlgebra::EigenRep & S = Stress[ Grid.ToIndex( i, j, k ) ];
	    for( int m = 0; m < 5; m ++ )
	      is >> std::skipws >> S.m[ m ];
	  }
      is.close();
    }
    //----------------------------------------------------------------------------------
    //
    //  SetMaterialGrid
    //   DEBUG
    //----------------------------------------------------------------------------------
    void SetRandomMaterialGrid( int NumX, int NumY, int NumZ,
				VPFFT::DataStructures::MaterialGrid & Grid,
                                int NumSeedPoints, unsigned int RandSeed)
    {
      srand( RandSeed );

      // build seed points

      std::vector< VPFFT::LinearAlgebra::SMatrix3x3 > RandOrientationList;
      std::vector< VPFFT::LinearAlgebra::SVector3 >   CenterList;

      for( int i = 0; i < NumSeedPoints; i ++ )
      {
        VPFFT::LinearAlgebra::SMatrix3x3 m;

        m.BuildActiveEulerMatrix (  DEGREE_TO_RADIAN( rand() % 360 ),
                                    DEGREE_TO_RADIAN( rand() % 180 ),
                                    DEGREE_TO_RADIAN( rand() % 360 ) );  // this is not uniform sampling - just for debugging
	
        VPFFT::LinearAlgebra::SVector3 v;

        v.Set( rand() % Grid.NumX(),
               rand() % Grid.NumY(),
               rand() % Grid.NumZ() );

        RandOrientationList.push_back( m );
        CenterList.push_back( v );
      }
      
      VPFFT::LinearAlgebra::SMatrix3x3 Left, Right;
      VPFFT::LinearAlgebra::EigenRep StressState( 0, 0, 0, 0, 0 );
      
      for( int i = 0; i < Grid.NumX(); i ++ )
        for( int j = 0; j < Grid.NumY(); j ++ )
          for( int k = 0; k < Grid.NumZ(); k ++ )
          {
            int MinIndex = 0;
            float MinDist = std::numeric_limits<float>::max();
            
            VPFFT::LinearAlgebra::SVector3 CurLocation ( i, j, k );
            for( int n = 0; n < NumSeedPoints; n ++ )
            {
              float Dist = (CurLocation - CenterList[n]).GetLength();
              if( Dist < MinDist )
              {
                MinIndex = n;
                MinDist  = Dist;
              }
            }
	    Grid.Orientation( i, j, k ) = RandOrientationList[MinIndex];
	    Grid.Stress     ( i, j, k ) = VPFFT::LinearAlgebra::EigenRep( 0, 0, 0, 0, 0 );   // this is only for debugging
          }
      Grid.UpdateSchmidtTensors();
    }

    //----------------------------------------------------------------------------------
    //
    //   SimpleOpenMPTest()
    //      a hard coded test example
    //
    //----------------------------------------------------------------------------------
    void SimpleOpenMPTest()
    {
      using VPFFT::LinearAlgebra::EigenRep;
      using VPFFT::LinearAlgebra::SMatrix5x5;
      using VPFFT::LinearAlgebra::SMatrix3x3;
  
      VPFFT::FCC_CrystalTest::FCC_SchmidtBasis FCC_Schmidt = VPFFT::FCC_CrystalTest::FCC_SchmidtBasis::Get();

      
      std::cout << "Begin Test --------- " << std::endl;
      VPFFT::DataStructures::MaterialGrid SampleGrid( 32, 32, 32 );
      //   VPFFT::DataStructures::MaterialGrid SampleGrid( 128, 128, 128 );
      
      
      std::string Filename = "TXFFT";
      //ReadRicardoData( Filename, 64, 64, 64, SampleGrid );
      SetMaterialGrid( SampleGrid );
      VPFFT::LinearAlgebra::SMatrix3x3 StrainRateM;
      StrainRateM.SetZero();
      StrainRateM.m[2][2] = 1;
      StrainRateM.m[1][1] = -0.5;
      StrainRateM.m[0][0] = -0.5;

      SampleGrid.SetVoxelLength( 1, 1, 1);
      std::cout << "Input Strain Rate  " << std::endl;
      std::cout << StrainRateM << std::endl;
      std::cout << "|| ===============================" << std::endl;
      EigenRep MacroscopicStrain( StrainRateM );
      std::cout << "Macroscopic Strain " << MacroscopicStrain << std::endl;

      //   exit(0);
      //SampleGrid.RunSingleStrainStep( MacroscopicStrain, 50, 0.01, 1e-5);
      //     SampleGrid.RunVPFFT( MacroscopicStrain, 50, 4, 0.01, 1e-5);
      SampleGrid.RunVPFFT( MacroscopicStrain, 100, 100, 0.03, 1e-5);
      //SampleGrid.RunVPFFT( MacroscopicStrain, 50, 100, 2, 1e-5);
      
      std::ofstream outfile("TestFile.dx");
      PrintGrid( outfile, SampleGrid );
      outfile.close();
            
      std::cout << "|| ===============================" << std::endl;
      
    }


    //----------------------------------------------------------------------------------
    // OpenMPTest
    //    OpenMP test reads in a config file to run a VPFFT test
    //
    //----------------------------------------------------------------------------------
    void OpenMPTest( const std::string & ConfigFilename, int NumThreads )
    {
      Utilities::ConfigFile InputFile;
      
      InputFile.Read( ConfigFilename );
      
      using VPFFT::LinearAlgebra::EigenRep;
      using VPFFT::LinearAlgebra::SMatrix5x5;
      using VPFFT::LinearAlgebra::SMatrix3x3;
  
      VPFFT::FCC_CrystalTest::FCC_SchmidtBasis FCC_Schmidt = VPFFT::FCC_CrystalTest::FCC_SchmidtBasis::Get();

      
      std::cout << "Begin Test --------- " << std::endl;
      VPFFT::DataStructures::MaterialGrid SampleGrid( InputFile.NumX, 
						      InputFile.NumY, 
						      InputFile.NumZ,
						      NumThreads );
      const int NumDomains = 30;
      const int RandSeed = 37;
      SetRandomMaterialGrid( InputFile.NumX, 
			     InputFile.NumY, 
			     InputFile.NumZ,
			     SampleGrid,
			     NumDomains, RandSeed );
      std::cout << "Finished grid generation " << std::endl;
      
      VPFFT::LinearAlgebra::SMatrix3x3 StrainRateM = InputFile.InputStrainRate;


      SampleGrid.SetVoxelLength( 1, 1, 1);
      std::cout << "Input Strain Rate  " << std::endl;
      std::cout << StrainRateM << std::endl;
      std::cout << "|| ===============================" << std::endl;
      EigenRep MacroscopicStrain( StrainRateM );
      std::cout << "Macroscopic Strain " << MacroscopicStrain << std::endl;
      
      SampleGrid.RunVPFFT( MacroscopicStrain, 
			   InputFile.MaxNRIter, 
			   InputFile.NumTimeIter,
			   InputFile.TimeStep,
			   InputFile.NRTolerence );
      
      std::ofstream outfile("TestFile.dx");
      PrintGrid( outfile, SampleGrid );
      outfile.close();
            
      std::cout << "|| ===============================" << std::endl;
      
    }

    

  }
}
