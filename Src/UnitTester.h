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

    //-------------------------------------------
    //  Default, hard coded open MP test used
    //  to make sure that the compilation is correct
    //-------------------------------------------
    void SimpleOpenMPTest();   
    
    //-------------------------------------------
    //
    //-------------------------------------------
    void OpenMPTest( const std::string & ConfigFilename, int NumThreads );   
    
    void SetMaterialGrid( VPFFT::DataStructures::MaterialGrid & Grid);

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
    


    void SetRandomMaterialGrid( int NumX, int NumY, int NumZ,
                                const VPFFT::DataStructures::MaterialGrid & Grid,
                                int NumSeedPoints, unsigned int RandSeed);
    
  }
  
}

#endif
