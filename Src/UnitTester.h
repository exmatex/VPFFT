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

    void SetMaterialGrid( VPFFT::DataStructures::MaterialGrid & Grid);

    void SetRandomMaterial( VPFFT::DataStructures::MaterialGrid & Grid,
                            int NumSeedPoints, unsigned int RandSeed );

    void SetRandomMaterialGrid( VPFFT::LinearAlgebra::SMatrix3x3 * OrientList,
                                VPFFT::LinearAlgebra::EigenRep * Stress,
                                int NumX, int NumY, int NumZ,
                                const VPFFT::DataStructures::MaterialGrid & Grid,
                                int NumSeedPoints, unsigned int RandSeed);
  }
  
}

#endif
