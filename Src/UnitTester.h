#ifndef UNIT_TESTER_H
#define UNIT_TESTER_H


#include "BurgerVectorBasis.h"
#include "LinearAlgebra.h"
#include "Solvers.h"
#include <vector>

#include "MaterialGrid.h"



namespace VPFFT
{

  namespace UnitTester
  {

    void ReadRicardoData( const std::string & oFilename, int DimX, int DimY, int DimZ,
                          VPFFT::DataStructures::MaterialGrid & Grid );


    void MPI_VPFFT_Test( int argc, char* argv[] );   // kitchen sink 
    void BasicTestMain();   // kitchen sink 

    
    //--------------------------------------------------
    //  
    //--------------------------------------------------
    void SerialTestMain( int argc, char* argv[] );
    void OpenMPTestMain( int argc, char* argv[] );

    void GenerateExampleState( int argc, char* argv[] );
    
    
    
    void SetMaterialGrid( VPFFT::DataStructures::MaterialGrid & Grid);

    //-----------------------------------------
    //  SetMaterialGrid
    // 
    //-----------------------------------------
    void SetMaterialGrid( VPFFT::DataStructures::MaterialGrid & Grid,
			  VPFFT::LinearAlgebra::SMatrix3x3 * OrientList,
			  VPFFT::LinearAlgebra::EigenRep * Stress,
			  int NumX, int NumY, int NumZ );
    
    //-----------------------------------------
    //  GenerateRandomMicroStructure
    // 
    //  Generate a random microstructure, then output
    //  the result into a text file for future input.
    //-----------------------------------------
    void WriteRandomMaterial( VPFFT::LinearAlgebra::SMatrix3x3 * OrientList,
			      VPFFT::LinearAlgebra::EigenRep * Stress,
			      int NumX, int NumY, int NumZ,
			      const VPFFT::DataStructures::MaterialGrid & Grid,
			      const std::string & Filename );
    
    //-----------------------------------------
    //  ReadRandomMaterial
    //-----------------------------------------
    void ReadRandomMaterial( VPFFT::LinearAlgebra::SMatrix3x3 * OrientList,
			     VPFFT::LinearAlgebra::EigenRep * Stress,
			     int NumX, int NumY, int NumZ,
			     const VPFFT::DataStructures::MaterialGrid & Grid,
			     const std::string & Filename );
    


    void SetRandomMaterialGrid( VPFFT::LinearAlgebra::SMatrix3x3 * OrientList,
                                VPFFT::LinearAlgebra::EigenRep * Stress,
                                int NumX, int NumY, int NumZ,
                                const VPFFT::DataStructures::MaterialGrid & Grid,
                                int NumSeedPoints, unsigned int RandSeed);
  }
  
}

#endif
