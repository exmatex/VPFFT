//---------------------------------------------------
//
//  LinearAlgebra.h - contains operators for simple 
//                    linear algebra operations in
//                    VPFFT.
//---------------------------------------------------

#ifndef _VPFFT_LINEAR_ALGEBRA_H
#define _VPFFT_LINEAR_ALGEBRA_H

#include "Types.h"
#include "Debug.h"
#include "Error.h"
#include <vector>
// #include "EigenBasis.h"

#include <Eigen/Dense>


#ifndef PI
#define PI ((FLOAT)  3.14159265358979323846)
#endif

#ifndef SQRT3
#define SQRT3 ((FLOAT)1.73205080756888f)
#endif

#define DEGREE_TO_RADIAN( degree ) ((degree) * (PI / (FLOAT)180.0f))
#define RADIAN_TO_DEGREE( radian ) ((radian) * ( (FLOAT) 180.0f / PI))





namespace VPFFT
{
  namespace LinearAlgebra
  {

    class SVector3;
  
    ////////////////////////////////////////////////////////////////////////////////////////////////

    class SMatrix3x3
    {
    public:
      SMatrix3x3();
      SMatrix3x3(FLOAT pMatrix[3][3]);
      void SetIdentity();
      void SetZero();

      //  BuildActiveSmalLRotation
      void BuildActiveSmallRotation(FLOAT fPhi, FLOAT fTheta, FLOAT fPsi);
    
      //   BuildActiveEulerMatrix - oEulerAngles are in radians 
      void BuildActiveEulerMatrix(FLOAT fPhi, FLOAT fTheta, FLOAT fPsi);

      void BuildPassiveEulerMatrix(FLOAT fPhi, FLOAT fTheta, FLOAT fPsi);

      // Builds matrix that rotates about the specified axis
      // Precondition:  oAxis must be a UNIT VECTOR
      void BuildRotationAboutAxis( const SVector3 &oAxis, FLOAT fAngle );

    
      //  Given a unit vector oProjDirection, build a projection matrix.
      void BuildProjectionMatrix(const SVector3 &oProjDirection );

      //----------------------------
      //  Symmetrize
      //----------------------------
      SMatrix3x3 Symmetrize()
      {
        SMatrix3x3 mRes;
        for( int i = 0; i < 3; i ++ )
          for( int j = 0; j < 3; j ++ )
            mRes.m[i][j] = 0.5 * ( m[i][j] + m[j][i] );
        return mRes;
      }

      //----------------------------
      //  AntiSymmetrize
      //----------------------------
      SMatrix3x3 AntiSymmetrize()
      {
        SMatrix3x3 mRes;
        for( int i = 0; i < 3; i ++ )
          for( int j = 0; j < 3; j ++ )
            mRes.m[i][j] = 0.5 * ( m[i][j] - m[j][i] );
        return mRes;
      }


      
      // return a vector of Euler angles in ZYZ convention in radians
      SVector3 GetEulerAngles() const;
    
      FLOAT Trace() const;
    
      void Transpose();

      SMatrix3x3 operator*( Float f ) const;
      SMatrix3x3 & operator*=( Float f );
      SMatrix3x3 operator/( Float f ) const;
      SMatrix3x3 & operator/=( Float f );
    
      SMatrix3x3 operator*( const SMatrix3x3& oRHS ) const;
      SMatrix3x3 operator+( const SMatrix3x3& oRHS ) const;
      SMatrix3x3 & operator+=( const SMatrix3x3& oRHS );
      SMatrix3x3 operator-( const SMatrix3x3& oRHS ) const;
      SVector3 operator*( const SVector3 &oRHS ) const;
    
      SMatrix3x3 Symmetrize(  ) const;
    
      FLOAT m[3][3];
    };

  
    std::ostream & operator<< ( std::ostream & os, const SMatrix3x3 &m ); 
  
    ////////////////////////////////////////////////////////////////////////////////////////////////

    class SVector3
    {
    public:
      SVector3();
      SVector3( FLOAT fX, FLOAT fY, FLOAT fZ );

      void Set( FLOAT fX, FLOAT fY, FLOAT fZ );

      FLOAT Get( UInt nAxis ) const;
      void Set( UInt nAxis, FLOAT fVal );
    
      FLOAT GetLength() const;
      
      // Normalize the vector
      void Normalize();

      // Assignment & comparisons
      SVector3 &operator=( const SVector3 &oRHS );
        
      // Add & subtract
      SVector3 operator+( const SVector3 &oRHS ) const;
      SVector3 &operator+=( const SVector3 &oRHS );
      SVector3 operator-( const SVector3& oRHS ) const;
      SVector3 operator-( ) const;
    
      SVector3 &operator-=( const SVector3 &oRHS );
    
      // Multiply & divide
      SVector3 operator*( FLOAT fRHS ) const;
      SVector3 &operator*=( FLOAT fRHS );
      SVector3 operator/( FLOAT fRHS ) const;
      SVector3 &operator/=( FLOAT fRHS );		

      // Accessor

      FLOAT operator[]( Int i) const;
      Float& operator[](Int i);
    
      FLOAT m_fX;
      FLOAT m_fY;
      FLOAT m_fZ;
   
    
    };
  
    SVector3 Cross( const SVector3 &oLHS, const SVector3 &oRHS );
    FLOAT Dot( const SVector3 &oLHS, const SVector3 &oRHS );
    SMatrix3x3 OuterProduct( const SVector3 & oLHS, const SVector3 & oRHS );
    SVector3 operator*( Float fLHS, const SVector3 &oRHS );
    std::ostream & operator<< ( std::ostream & os, const SVector3 &v ); 
  
    SMatrix3x3 operator*( Float fLHS, const SMatrix3x3 &oRHS );
    SVector3 DegreeToRadian( const SVector3 & oRHS );
    SVector3 RadianToDegree( const SVector3 & oRHS );





    //-----------------------------
    //  Singleton to represent eigen-basis
    //-----------------------------
    class EigenBasis
    {
    public:
    
      static EigenBasis & Get() 
      {
        static EigenBasis oBasis;
        return oBasis;
      }

      //-----------------------------
      //  Default constructor (also constructor functor)
      //-----------------------------
      EigenBasis() { BuildBasisList(); }

      // Accessor
      const std::vector<SMatrix3x3>  & GetBasisList() const { return BasisList;    }
      const SMatrix3x3 & operator()( int n ) const          { return BasisList[n]; }   // NO ASSERT
    
    private:

      std::vector<SMatrix3x3> BasisList;
    
      void BuildBasisList()
      {
        BasisList.resize(6);
            
        Float m[3][3];
        //      basis(:, :, 1) =  1/sqrt(2) .* [-1 0 0; 0  1 0; 0 0 0 ];
        m[0][0] = -1;   m[0][1] = 0;  m[0][2] = 0;
        m[1][0] =  0;   m[1][1] = 1;  m[1][2] = 0;
        m[2][0] =  0;   m[2][1] = 0;  m[2][2] = 0;
        BasisList[0] =  Float( 1 ) / sqrt( Float( 2 ) ) * SMatrix3x3( m ) ;
      
        //       basis(:, :, 2) =  1/sqrt(6) .* [-1 0 0; 0 -1 0; 0 0 2 ];
        m[0][0] = -1;   m[0][1] = 0;  m[0][2] = 0;
        m[1][0] =  0;   m[1][1] = -1; m[1][2] = 0;
        m[2][0] =  0;   m[2][1] = 0;  m[2][2] = 2;
        BasisList[1] = Float( 1 ) / sqrt( Float( 6 ) ) * SMatrix3x3( m );

        //       basis(:, :, 3) =  1/sqrt(2) .* [ 0 0 0; 0  0 1; 0 1 0 ];
        m[0][0] =  0;   m[0][1] = 0;  m[0][2] = 0;
        m[1][0] =  0;   m[1][1] = 0;  m[1][2] = 1;
        m[2][0] =  0;   m[2][1] = 1;  m[2][2] = 0;
        BasisList[2] = Float( 1 ) / sqrt( Float( 2 ) ) * SMatrix3x3( m );

        //       basis(:, :, 4) =  1/sqrt(2) .* [ 0 0 1; 0  0 0; 1 0 0 ];
        m[0][0] =  0;   m[0][1] = 0;  m[0][2] = 1;
        m[1][0] =  0;   m[1][1] = 0;  m[1][2] = 0;
        m[2][0] =  1;   m[2][1] = 0;  m[2][2] = 0;
        BasisList[3] = Float( 1 ) / sqrt( Float( 2 ) ) * SMatrix3x3( m );

        //       basis(:, :, 5) =  1/sqrt(2) .* [ 0 1 0; 1  0 0; 0 0 0 ];
        m[0][0] =  0;   m[0][1] = 1;  m[0][2] = 0;
        m[1][0] =  1;   m[1][1] = 0;  m[1][2] = 0;
        m[2][0] =  0;   m[2][1] = 0;  m[2][2] = 0;
        BasisList[4] = Float( 1 ) / sqrt( Float( 2 ) ) * SMatrix3x3( m );

        //       basis(:, :, 6) =  1/sqrt(3) .* [ 1 0 0; 0  1 0; 0 0 1 ];
        m[0][0] =  1;   m[0][1] = 0;  m[0][2] = 0;
        m[1][0] =  0;   m[1][1] = 1;  m[1][2] = 0;
        m[2][0] =  0;   m[2][1] = 0;  m[2][2] = 1;
        BasisList[5] = Float( 1 ) / sqrt( Float( 3 ) ) * SMatrix3x3( m );
      }
    };
    


    //----------------------------------------------------
#ifndef EIGEN_DIM
#define EIGEN_DIM 5 
#endif
    
    
    struct Tensor4Rank;  // Fourth rank tensor
    struct EigenRep;
    struct SMatrix5x5;
    
    
    //----------------------------
    //  4th Rank tensor
    //
    //  Could really improve the element ordering... maybe
    //  replace with a sparse vector
    //----------------------------
    struct Tensor4Rank
    {
      Float m[81];

      void SetZero() { memset( m, 0, sizeof(Float) * 81 ); }

      Float & operator()( int i, int j, int k, int l )
      { return m[ i * 27 + j * 9 + k * 3 + l ];  }
      
      Float  operator()( int i, int j, int k, int l ) const
      { return m[ i * 27 + j * 9 + k * 3 + l ];  }
    };

    //----------------------------
    //  InnerProduct
    //----------------------------
    Float       InnerProduct( const EigenRep    & LHS, const EigenRep    & RHS );
    Float       InnerProduct( const SMatrix3x3  & LHS, const SMatrix3x3  & RHS );
    Float       InnerProduct( const Tensor4Rank & LHS, const Tensor4Rank & RHS );
    
    SMatrix5x5  OuterProduct( const EigenRep   & LHS, const EigenRep   & RHS );
    Tensor4Rank OuterProduct( const SMatrix3x3 & LHS, const SMatrix3x3 & RHS );

    Tensor4Rank Symmetrize( const Tensor4Rank & T );
    Tensor4Rank AntiSymmetrize( const Tensor4Rank & T );
    
    //-----------------------------
    //  EigenRep - Eigen strain representation
    //-----------------------------
    struct EigenRep
    {
      Float m[EIGEN_DIM];

      EigenRep(){}

      //-----------------------------
      //  EigenRep
      //-----------------------------
      EigenRep( const SMatrix3x3 & Mat )
      {
        EigenBasis Basis = EigenBasis::Get();
        for ( int i = 0; i < EIGEN_DIM; i ++ )
          m[i] = InnerProduct( Mat, Basis( i ) );
      }

      //-----------------------------
      // Ctor
      //-----------------------------
      EigenRep( const std::vector<Float> &v )
      {
        RUNTIME_ASSERT( v.size() == EIGEN_DIM,
                        "Error: Constructor Does not take vector<Float>.size() != EIGEN_DIM");
        for( int i = 0; i < EIGEN_DIM; i ++ )
          m[i] = v[i];
      }

      //-----------------------------
      // Ctor
      //-----------------------------
      EigenRep( Float a, Float b, Float c, Float d, Float e )
      {
        m[0] = a;
        m[1] = b;
        m[2] = c;
        m[3] = d;
        m[4] = e;
      }

      inline Float operator()( int n ) const
      {
        DEBUG_ASSERT( n < EIGEN_DIM, "EigenRep::operator() ");
        return m[n];
      } 
      
      inline Float & operator()( int n )
      {
        DEBUG_ASSERT( n < EIGEN_DIM, "EigenRep::operator() ");
        return m[n];
      }

      //----------------------------
      //  ToMatrixRep
      //----------------------------
      SMatrix3x3 ToMatrixRep() const
      {
        SMatrix3x3 Result;
        EigenBasis Basis = EigenBasis::Get();
        Result.SetZero();
        for( int i = 0; i < EIGEN_DIM; i ++ )
          Result += m[i] * Basis( i ); 
        return Result;
      }

      //----------------------------
      //  Operator *=
      //----------------------------
      EigenRep & operator*= ( Float f )
      {
        m[0] *= f;
        m[1] *= f;
        m[2] *= f;
        m[3] *= f;
        m[4] *= f; 
        return *this;
      }

      //----------------------------
      //  Operator +=
      //----------------------------
      EigenRep & operator+= ( const EigenRep & RHS )
      {
        m[0] += RHS.m[0];
        m[1] += RHS.m[1];
        m[2] += RHS.m[2];
        m[3] += RHS.m[3];
        m[4] += RHS.m[4]; 
        return *this;
      }


      //----------------------------
      //  Operator -=
      //----------------------------
      EigenRep & operator-= ( const EigenRep & RHS )
      {
        m[0] -= RHS.m[0];
        m[1] -= RHS.m[1];
        m[2] -= RHS.m[2];
        m[3] -= RHS.m[3];
        m[4] -= RHS.m[4]; 
        return *this;
      }

      //----------------------------
      //  Operator +=
      //----------------------------
      EigenRep operator+ ( const EigenRep & RHS  ) const
      {
        EigenRep oRes;
        oRes.m[0] = m[0] + RHS.m[0];
        oRes.m[1] = m[1] + RHS.m[1];
        oRes.m[2] = m[2] + RHS.m[2];
        oRes.m[3] = m[3] + RHS.m[3];
        oRes.m[4] = m[4] + RHS.m[4]; 
        return oRes;
      }

      
      //----------------------------
      //  Operator +=
      //----------------------------
      EigenRep operator- ( const EigenRep & RHS  ) const
      {
        EigenRep oRes;
        oRes.m[0] = m[0] - RHS.m[0];
        oRes.m[1] = m[1] - RHS.m[1];
        oRes.m[2] = m[2] - RHS.m[2];
        oRes.m[3] = m[3] - RHS.m[3];
        oRes.m[4] = m[4] - RHS.m[4]; 
        return oRes;
      }

      //----------------------------
      //
      //----------------------------
      EigenRep operator- (   ) const
      {
        EigenRep oRes;
        oRes.m[0] = -m[0];
        oRes.m[1] = -m[1];
        oRes.m[2] = -m[2];
        oRes.m[3] = -m[3];
        oRes.m[4] = -m[4];
        return oRes;
      }
      //----------------------------
      //  Operator *
      //----------------------------
      EigenRep  operator* ( Float f ) const
      {
        EigenRep LHS;
        LHS.m[0] = m[0] * f;
        LHS.m[1] = m[1] * f;
        LHS.m[2] = m[2] * f;
        LHS.m[3] = m[3] * f;
        LHS.m[4] = m[4] * f;
        return LHS;
      }

      //----------------------------
      // operator/
      //----------------------------
      EigenRep  operator/ ( Float f ) const
      {
        EigenRep oRes = *this;
        oRes /= f;
        return oRes;
      }

      //----------------------------
      //
      //----------------------------
      EigenRep  & operator/= ( Float f ) 
      {
        m[0] /= f;
        m[1] /= f;
        m[2] /= f;
        m[3] /= f;
        m[4] /= f;
        return *this;
      }
      
      //----------------------------
      //  Symmetrize
      //----------------------------
      EigenRep Symmetrize()
      {
        SMatrix3x3 m = ToMatrixRep();
        SMatrix3x3 mRes;
        for( int i = 0; i < 3; i ++ )
          for( int j = 0; j < 3; j ++ )
            mRes.m[i][j] = 0.5 * ( m.m[i][j] + m.m[j][i] );
        
        return EigenRep( mRes );
      }
      
      //----------------------------
      //  AntiSymmetrize
      //----------------------------
      EigenRep AntiSymmetrize()
      {
        SMatrix3x3 m = ToMatrixRep();
        SMatrix3x3 mRes;
        for( int i = 0; i < 3; i ++ )
          for( int j = 0; j < 3; j ++ )
            mRes.m[i][j] = 0.5 * ( m.m[i][j] - m.m[j][i] );
        
        return EigenRep( mRes );
      }


      //----------------------------
      //  Norm
      //----------------------------
      Float Norm()
      {
        return std::sqrt(m[0] * m[0] +
                         m[1] * m[1] +
                         m[2] * m[2] +
                         m[3] * m[3] +
                         m[4] * m[4] ) ;
      }
    };
    

    //----------------------------
    //----------------------------
    std::ostream & operator<< ( std::ostream & os, const EigenRep &v );
    std::ostream & operator<< ( std::ostream & os, const SMatrix5x5 &o );
    std::ostream & operator<< ( std::ostream & os, const Tensor4Rank &T ); 

    //----------------------------
    //  SMatrix5x5 - directly from SMatrix3x3
    //
    //  5x5 comes from the fact that there are 5 eigen basis 
    //
    //----------------------------
    class SMatrix5x5
    {
    public:
      SMatrix5x5() {}

      SMatrix5x5( const Tensor4Rank & T )
      {
        EigenBasis Basis = EigenBasis::Get();
                
        SetZero();
        for( int Slip1 = 0; Slip1 < 5; Slip1 ++)
          for( int Slip2 = 0; Slip2 < 5; Slip2 ++)
            for( int i = 0; i < 3; i ++ )
              for( int j = 0; j < 3; j ++ )
                for( int k = 0; k < 3; k ++ )
                  for( int l = 0; l < 3; l ++ )
                    operator()( Slip1, Slip2 ) += T(i, j, k, l) 
                      * Basis(Slip1).m[i][j] * Basis(Slip2).m[k][l];
      }
      
      SMatrix5x5(Float pMatrix[5][5]);
      void SetIdentity();
      void SetZero();
      
      void  Transpose();
      
      SMatrix5x5 operator*( Float f ) const;
      SMatrix5x5 & operator*=( Float f );
      SMatrix5x5 operator/( Float f ) const;
      SMatrix5x5 & operator/=( Float f );
      
      SMatrix5x5 operator*( const SMatrix5x5& oRHS ) const;
      SMatrix5x5 operator+( const SMatrix5x5& oRHS ) const;
      SMatrix5x5 & operator+=( const SMatrix5x5& oRHS );
      SMatrix5x5 operator-( const SMatrix5x5& oRHS ) const;
      SMatrix5x5 operator-(  ) const;
      SMatrix5x5 & operator-=( const SMatrix5x5& oRHS );
      
      inline Float   operator() ( size_t i, size_t j ) const { return m[ j * 5 + i ];  }  // this ordering is necessary for Eigen to work
      inline Float & operator() ( size_t i, size_t j )       { return m[ j * 5 + i ];  }

      void MovingAverage( const SMatrix5x5 & RHS, int CurrentCount );
      Float m[25];

      //----------------------------
      //  ToTensor4RankRep()
      //----------------------------
      Tensor4Rank ToTensor4RankRep() const
      {
        EigenBasis Basis = EigenBasis::Get();

        Tensor4Rank oRes;
        oRes.SetZero();
        
        for( int i = 0; i < 3; i ++ )
          for( int j = 0; j < 3; j ++ )
            for( int k = 0; k < 3; k ++ )
              for( int l = 0; l < 3; l ++ )
                for( int Slip1 = 0; Slip1 < 5; Slip1 ++)
                  for( int Slip2 = 0; Slip2 < 5; Slip2 ++)
                    oRes(i, j, k, l) += operator()( Slip1, Slip2 )
                      * Basis(Slip1).m[i][j] * Basis(Slip2).m[k][l];
        return oRes;
      }

      //----------------------------
      //  NormalizeValues() --
      //   remove all over and underflows
      //----------------------------
      void NormalizeValues()
      {
        for ( int i = 0; i < 5; i ++ )
          for ( int j = 0; j < 5; j ++ )
          {
            if( std::fabs( operator()( i, j ) ) < 1e-10 )
              operator()(i, j) = 0;
          }
      }
    };

    EigenRep operator*( const SMatrix5x5 & LHS, const EigenRep & RHS);
    EigenRep operator*(const EigenRep & LHS,  const SMatrix5x5 & RHS);
    //----------------------------
    //  Linear Solvers
    //
    //----------------------------

    //----------------------------
    //  LU_Solver
    //
    //  Solving x in A x = b, given
    //  that A is a 5x5 matrix.  x and b are 5 x 1.
    //----------------------------
    EigenRep LU_Solver( const SMatrix5x5 & A, const EigenRep & b );
    
  }//  end LinearAlgebra namespace
  
} // end VPFFT namespace

#endif
