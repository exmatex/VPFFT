#ifndef BURGER_VECTOR_BASIS_H
#define BURGER_VECTOR_BASIS_H

#include <iostream>
#include <vector>
#include <cmath>
#include "LinearAlgebra.h"
namespace VPFFT
{

  namespace FCC_CrystalTest
  {
    static const int NumSystems = 12;

    // change to the same order as Ricardo's so that the debugging be more streamlined
    static const Float BurgerVector[NumSystems][3] =
      {
        {0,    1,    1},
        {1,    0,    1},
        {1,   -1,    0},
        {0,    1,   -1},
        {1,    0,    1},
        {1,    1,    0},
        {0,    1,    1},
        {1,    0,   -1},
        {1,    1,    0},
        {0,    1,   -1},
        {1,    0,   -1},
        {1,   -1,    0}
        
//         {-1,  1,  0},
//         { 1,  0, -1},
//         { 0, -1,  1},
//         {-1, -1,  0},
//         { 1,  0,  1},
//         { 0,  1, -1},
//         { 1,  1,  0},
//         {-1,  0,  1},
//         { 0, -1, -1},
//         { 1, -1,  0},
//         {-1,  0, -1},
//         { 0,  1,  1}
      };

    static const Float SlipPlaneNormal[NumSystems][3] =
      {

        {1,    1,   -1},
        {1,    1,   -1},
        {1,    1,   -1},
        {1,   -1,   -1},
        {1,   -1,   -1},
        {1,   -1,   -1},
        {1,   -1,    1},
        {1,   -1,    1},
        {1,   -1,    1},
        {1,    1,    1},
        {1,    1,    1},
        {1,    1,    1}


        
//         { 1,  1,  1},
//         { 1,  1,  1},
//         { 1,  1,  1},
//         { 1, -1, -1},
//         { 1, -1, -1},
//         { 1, -1, -1},
//         {-1,  1, -1},
//         {-1,  1, -1},
//         {-1,  1, -1},
//         {-1, -1,  1},
//         {-1, -1,  1},
//         {-1, -1,  1}        
      };


    using VPFFT::LinearAlgebra::SMatrix3x3;
    using std::vector;
    //---------------------
    //  FCC_SchmidtBasis
    //
    //
    //---------------------
    class FCC_SchmidtBasis
    {
    public:
      static int         NumSlipSystems(   )       { return NumSystems;  }
      const SMatrix3x3 & operator()( int n ) const { return SchmidtBasis[n]; }
      SMatrix3x3         operator()( int n )       { return SchmidtBasis[n]; }

      const SMatrix3x3 & Anti( int n ) const { return AntiSchmidtBasis[n]; }
      SMatrix3x3         Anti( int n )       { return AntiSchmidtBasis[n]; }
      
      static FCC_SchmidtBasis & Get() 
      {
        static FCC_SchmidtBasis  BasisTensors;
        return BasisTensors;
      }
    private:

      FCC_SchmidtBasis()
      {
        BuildBasisList();
        BuildAntiBasisList();
      }
      
      vector<SMatrix3x3> SchmidtBasis;
      vector<SMatrix3x3> AntiSchmidtBasis;
      
      //---------------------
      //  OuterProduct
      //---------------------
      SMatrix3x3 OuterProduct( const Float v1[3], const Float v2[3] )
      {
        SMatrix3x3 oRes;
        for( int i = 0; i < 3; i ++ )
          for( int j = 0; j < 3; j ++ )
            oRes.m[i][j]  = v1[i] * v2[j];
        return oRes;
      }

      //---------------------
      //  BuildBasisList
      //---------------------
      void BuildBasisList()
      {
        SchmidtBasis.clear();
        for( int i = 0; i < NumSystems; i ++ )
        {
          SMatrix3x3 m1 = OuterProduct( BurgerVector[i], SlipPlaneNormal[i] );
          SMatrix3x3 m2 = OuterProduct( SlipPlaneNormal[i], BurgerVector[i]  );
          SMatrix3x3 mSym = (m1 + m2) / (2.0 * std::sqrt( 6.0 ) );
          SchmidtBasis.push_back( mSym );
        }
      }

      //---------------------
      //  BuildAntiBasisList
      //---------------------
      void BuildAntiBasisList()
      {
        AntiSchmidtBasis.clear();
        for( int i = 0; i < NumSystems; i ++ )
        {
          SMatrix3x3 m1 = OuterProduct( BurgerVector[i], SlipPlaneNormal[i] );
          SMatrix3x3 m2 = OuterProduct( SlipPlaneNormal[i], BurgerVector[i]  );
          SMatrix3x3 mSym = (m1 - m2) / (2.0 * std::sqrt( 6.0 ) );
          AntiSchmidtBasis.push_back( mSym );
        }
      }
        
    };
    
    
  }
  
}// end VPFFT


#endif

