#include "Config.h"



namespace Utilities
{

  //-----------------------------------
  // Read
  //-----------------------------------
  bool ConfigFile::Read( const std::string & ConfigFilename )
  {
    Size_Type  nBufferSize = 0;
    char *pBuffer = ReadFileToBuf( nBufferSize, ConfigFilename);
    
    if( nBufferSize<= 0 || !pBuffer )
    {
      std::cerr << "ERROR: Unable to read config file" << endl; 
      return false;
    }
	
    if( !Parse( std::string( pBuffer, nBufferSize ) ) )
    {
      std::cerr << "ERROR: Unable to parse config file" << endl;
      return false;
    }
    
    if( pBuffer )
      delete[] pBuffer;
    
    return true;
  }
  
  //-----------------------------------
  //  Parser
  //-----------------------------------
  bool ConfigFile::Parse( const string & sBuf )
  {
 
    const char* KeywordStrings[] =
      {
        
        "#",  // comment

        "GridDimension",
        "InputStrainRate",
        "MaxNRIter",
        "TimeStep",
        "NumTimeIter",
        "NRTolerence",
        "TwoGrainInit",
        "FileInit",

        "CRSS",
        "GammaDotBase",
        "RateSensitivity"
      };
    
    enum EKeyword
    {
      // comment
      eComment,

      eGridDimension,
      eInputStrainRate,
      eMaxNRIter,
      eTimeStep,
      eNumTimeIter,
      eNRTolerence,
      eTwoGrainInit,
      eFileInit,
      
      eCRSS,
      eGammaDotBase,
      eRateSensitivity,
      
      //  error check on number of keywords
      eNumKeywords
    };

    vector< vector<string> > vsTokens;
    Tokenize( vsTokens, sBuf, " \t\n");
  
    
    bool *vInitializationCheck;
    bool *vRequirementCheck;
    vInitializationCheck = new bool[eNumKeywords];   // a check list to see which variable is not initialized
    vRequirementCheck    = new bool[eNumKeywords];   // a check list to see which variable is not initialized

    for( int i = 0; i < eNumKeywords; i ++ )
    {
      vInitializationCheck[i] = false;
      vRequirementCheck[i]    = true;
    }
    vRequirementCheck[eFileInit] = false;
    vRequirementCheck[eTwoGrainInit] = false;
    
    
    for(Size_Type i = 0; i < vsTokens.size(); i ++)
    {
      Size_Type iFirstToken;
      
      if ( vsTokens[i].size() == 0)
      {
        iFirstToken = eComment;
      }
      else if ( vsTokens[i][0].find_first_of( KeywordStrings[ eComment ] ) == 0) // if first token is comment
      {
        iFirstToken = eComment;
      }
      else
      {
        // Identify the keyword in the beginning of the line
        for( iFirstToken = 0; iFirstToken < eNumKeywords; iFirstToken ++ )
        {
          if(strcmp(KeywordStrings[iFirstToken], vsTokens[i][0].c_str()) == 0)
            break;
        }
         std::cout << i + 1 << " " << KeywordStrings[iFirstToken] << " "
                   << vsTokens[i][0].c_str() << " Token " << iFirstToken << endl; 
      }
      vInitializationCheck[ iFirstToken ] = true;

      switch( iFirstToken )   // look at first token of each line
      {
        
        case eComment:  // Comment
          break;
        case eGridDimension:
          RUNTIME_ASSERT( vsTokens[i].size() >= 4,
                          "GridDimension: Expect 3 parameters [NumX NumY NumZ]");
          
          NumX = atoi( vsTokens[i][1].c_str() );
          NumY = atoi( vsTokens[i][2].c_str() );
          NumZ = atoi( vsTokens[i][3].c_str() );
          break;
          
        case eInputStrainRate:
          {
            bool bSuccess =  HasLineContinuation( i, 3, vsTokens, KeywordStrings, eNumKeywords );
            RUNTIME_ASSERT( bSuccess, "InputStrainRate: Expect s00 s11 s22 \n <empty space> s10 s11 s12 \n <empty space> s20 s21 s22  " );

            for(int di = 0; di < 3; di ++ )
            {
              ++i;
              InputStrainRate.m[di][0] =  atof( vsTokens[i][0].c_str() );
              InputStrainRate.m[di][1] =  atof( vsTokens[i][1].c_str() );
              InputStrainRate.m[di][2] =  atof( vsTokens[i][2].c_str() );
            }
            std::cout << "Parsed input strain rate " << std::endl;
            std::cout << InputStrainRate << std::endl;
          }
          break;
        case eMaxNRIter:
          MaxNRIter = atoi( vsTokens[i][1].c_str() );
          break;

        case eTimeStep:
          TimeStep = atof( vsTokens[i][1].c_str() );
          break;

        case eNumTimeIter:
          NumTimeIter = atoi( vsTokens[i][1].c_str() );
          break;
          
        case eNRTolerence:
          NRTolerence = atof( vsTokens[i][1].c_str() );
          break;
          
        case eTwoGrainInit:
          InitMethod = eInitTwoGrain;
          break;
          
        case eFileInit:
          {
          RUNTIME_ASSERT( vsTokens[i].size() >= 2,
                          "FileInit:  Expect Filename for initialization" );
          
          InputMicrostructureFilename = vsTokens[i][1];
          std::ifstream fd;
          fd.open( InputMicrostructureFilename.c_str() );
          if( !fd  )
          {
            std::cout << "FileInit " << InputMicrostructureFilename << " Not found " << std::endl;
            exit(0);                      
          }
          fd.close();
          InitMethod = eInitFile;
          break;
          }
        case eCRSS:
          RUNTIME_ASSERT( vsTokens[i].size() >= VPFFT::FCC_CrystalTest::NumSystems,
                          "GridDimension: Expect the same number of parameters as number slip systems in the symmetry\n");

          CRSS.resize( VPFFT::FCC_CrystalTest::NumSystems );
          for( int n = 1; n <= VPFFT::FCC_CrystalTest::NumSystems; n ++ )
            CRSS[n] = atof( vsTokens[i][n].c_str() );
          
          break;
        case eGammaDotBase:
          RUNTIME_ASSERT( vsTokens[i].size() >= VPFFT::FCC_CrystalTest::NumSystems,
                          "GridDimension: Expect the same number of parameters as number slip systems in the symmetry\n");
          GammaDotBase.resize( VPFFT::FCC_CrystalTest::NumSystems );
          for( int n = 1; n <= VPFFT::FCC_CrystalTest::NumSystems; n ++ )
            GammaDotBase[n] = atof( vsTokens[i][n].c_str() );
          
          break;
        case eRateSensitivity:
          RUNTIME_ASSERT( vsTokens[i].size() >= VPFFT::FCC_CrystalTest::NumSystems,
                          "GridDimension: Expect the same number of parameters as number slip systems in the symmetry\n");
          RateSensitivity.resize( VPFFT::FCC_CrystalTest::NumSystems );
          for( int n = 1; n <= VPFFT::FCC_CrystalTest::NumSystems; n ++ )
            RateSensitivity[n] = atoi( vsTokens[i][n].c_str() );

          break;
          
        default:
          {
            std::cerr << "[ConfigFile] Error: syntax not recognized:  Line " << i
                      << " Keyword: " << vsTokens[i][0] << endl;
            return false;
          }
      }
    } // end loop over tokens
    
    bool Success = true;
    for( int i = 0; i < eNumKeywords; i ++ )
    {
      if( !vInitializationCheck[i] && vRequirementCheck[i] )
      {
        std::cerr << "[ConfigFile] Initialization Error: Missing parameter: "
                  << KeywordStrings[i]  << " not optional." << std::endl;
        Success = false;
      }
    }
    if( ! vInitializationCheck[eFileInit] &&
        ! vInitializationCheck[eTwoGrainInit]  )
    {
      std::cerr << "Initialization scheme (either FileInit or TwoGrain init) needs to be specified." << std::endl;
      Success = false;
    }
    
    if( !Success )
    {
      std::cerr << "[ConfigFile] Initialization Failed" << std::endl;
      exit(0);
    }
    delete [] vInitializationCheck;
    delete [] vRequirementCheck;

    std::cout << "GridDimension "
              << NumX << " "
              << NumY << " "
              << NumZ  << std::endl;
    std::cout << "InputStrainRate " << std::endl
              << InputStrainRate    << std::endl;
    std::cout << "MaxNRIter "  << MaxNRIter << std::endl;
    std::cout << "TimeStep "   << TimeStep << std::endl;
    std::cout << "InitMethod " << InitMethod << std::endl;
    
  }


  //--------------------------------------------
  //  Returns a Null terminated C-string
  //--------------------------------------------
  char* ConfigFile::ReadFileToBuf( Size_Type & nBufferSize, std::string filename)
  {
    std::ifstream infile;
    infile.open(filename.c_str());
    FILE *pFile;
    char *pBuffer;

    if ( (pFile = fopen(filename.c_str(), "r" )) == NULL ){
      cerr << "[ReadFileToBuf]Cannot Open File: " << filename.c_str()  << endl;
      return NULL;
    }

    // Get the size of the file
    fseek( pFile, 0L, SEEK_END ); // Position to end of file
    nBufferSize = ftell( pFile );           // Get file length
    rewind( pFile );                                // Back to start of file

    if( nBufferSize <= 0 )
      return NULL;
    
    // Read in the entire file and close the file handle
    pBuffer = new char[nBufferSize + 1];

    if ( fread( pBuffer, nBufferSize, 1, pFile ) <= 0 ){
      fclose( pFile );
      cerr << "[ReadFileToBuf]:File read error" << endl;
      return NULL;
    }

    fclose( pFile );

    pBuffer[ nBufferSize ] = '\0'; // NULL termination of string
    return pBuffer;
  }


  
}
