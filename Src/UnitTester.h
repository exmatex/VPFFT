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


    
  }
  
}

#endif
