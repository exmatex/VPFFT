#ifndef VPFFT_CONFIG_H_
#define VPFFT_CONFIG_H

#include "Parser.h"
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include "LinearAlgebra.h"
#include "BurgerVectorBasis.h"

namespace Utilities
{
  using std::endl;
  using std::cout;
  using std::vector;
  using namespace VPFFT::LinearAlgebra;
  //----------------------------
  //  Configuration File
  //----------------------------
  class ConfigFile
  {
  public:
    enum RunType
      {
        eInitTwoGrain,
        eInitFile, // doesn't exist
      };

    

    int NumX;
    int NumY;
    int NumZ;

    
    SMatrix3x3 InputStrainRate;
    int        MaxNRIter;
    float      NRTolerence;
    float      TimeStep;
    int        NumTimeIter;
    RunType    InitMethod;


    std::string InputMicrostructureFilename;
    //-------------- Model parameters
    std::vector<Float> CRSS;
    std::vector<Float> GammaDotBase;
    std::vector<int> RateSensitivity;

    
    bool Read( const string & Filename );
  private:
    
    bool Parse( const string & StringBuffer );
    char* ReadFileToBuf( Size_Type & nBufferSize, std::string filename );

    //-----------------------
    //  return true if first token is space, i.e., this
    //  line is a continuation of the previous line
    //-----------------------
    bool HasLineContinuation( int CurrentLine,
                              int NumLineContinuation,
                              const vector< vector< string> > & vsTokens,
                              const char * KeywordStrings[],
                              int NumKeywords )
    {
      if ( CurrentLine + NumLineContinuation >= vsTokens.size() )
        return false;
      std::cout << "Check line continuation" << std::endl;

      for( int d = 1; d <= NumLineContinuation; d ++ )
      {
        std::cout << d << std::endl;
        for( int n = 0; n < NumKeywords; n ++ )
        {
          if( vsTokens[CurrentLine + d][0].find_first_of( KeywordStrings[ n ] ) == 0 )
            return false;
          
        }
      }
      return true;
    }
      
  };
}

#endif
