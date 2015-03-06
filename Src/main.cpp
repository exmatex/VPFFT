#include <iostream>

#include "BurgerVectorBasis.h"
#include "LinearAlgebra.h"
#include "Solvers.h"
#include <vector>

#include "MaterialGrid.h"
#include "UnitTester.h"


int main( int argc, char ** argv)
{
  
//  VPFFT::UnitTester::SimpleOpenMPTest();
  
  if( argc != 3 )
  {
    std::cout << "Usage: " << argv[0] << " <Input File> <Num OpenMP Threads>" << std::endl; 
    exit( 0 );
  }
  std::cout << " --- Number of threads " << atoi( argv[2] ) << std::endl;
  VPFFT::UnitTester::OpenMPTest( "Test.input", atoi( argv[2] ) );
  return 0;
}
