//---------------------------------------------------
//
//  LinearAlgebra.cpp
//
//
//---------------------------------------------------

#include "LinearAlgebra.h"

namespace VPFFT
{
  namespace LinearAlgebra
  {


    
  //----------------------------------------------------------------------------------------------
  // Public : SMatrix3x3
  //----------------------------------------------------------------------------------------------

  SMatrix3x3::SMatrix3x3()
  {}

  //----------------------------------------------------------------------------------------------
  // Public : SMatrix3x3
  //----------------------------------------------------------------------------------------------
  SMatrix3x3::SMatrix3x3(FLOAT pMatrix[3][3])
  {
    m[0][0] = pMatrix[0][0];
    m[1][0] = pMatrix[1][0];
    m[2][0] = pMatrix[2][0];

    m[0][1] = pMatrix[0][1];
    m[1][1] = pMatrix[1][1];
    m[2][1] = pMatrix[2][1];


    m[0][2] = pMatrix[0][2];
    m[1][2] = pMatrix[1][2];
    m[2][2] = pMatrix[2][2];

  }

  //----------------------------------------------------------------------------------------------
  //  Symmetrize
  //----------------------------------------------------------------------------------------------
  SMatrix3x3 SMatrix3x3::Symmetrize( ) const
  {
    SMatrix3x3 oRes;
    for( int i = 0; i < 3; i ++ )
      for( int j = 0; j < 3; j ++ )
        oRes.m[i][j] = static_cast<Float>(0.5) * ( m[i][j] + m[j][i] );
    return oRes;
  }

  //----------------------------------------------------------------------------------------------
  // Public : SetIdentity
  //----------------------------------------------------------------------------------------------
  void SMatrix3x3::SetIdentity()
  {
    memset( m, 0, sizeof(FLOAT) * 9 );
    m[0][0] = m[1][1] = m[2][2] = (FLOAT)1.0;
  }

  //----------------------------------------------------------------------------------------------
  // Public : SetZero
  //----------------------------------------------------------------------------------------------
  void SMatrix3x3::SetZero()
  {
    memset( m, 0, sizeof(FLOAT) * 9 );
  }

 //----------------------------------------------------------------------------------------------
  // Public : BuildRotationAboutAxis.
  //
  // Description:  An active transformation matrix that will perform the rotation in a positive sence.
  // (i.e., follows the right hand rule.)
  //
  // Precondition:  oAxis must be a UNIT VECTOR
  //----------------------------------------------------------------------------------------------
  void SMatrix3x3::BuildRotationAboutAxis( const SVector3 &oAxis, FLOAT fAngle )
  {
    // 1. Perform transformations which align rotation axis with one of coordinate axis (x, y, z).
    //	  This is essentially a transformation matrix with the given vector as one of the rows, and
    //    orthogonal vectors as the other rows.
    // 2. Perform rotation about the axis
    // 3. Do inverse of (1)

    // Alternatively, use the method described on page 79-80 of Eric Lengyel's "Mathematics for 3d game programming
    // and computer graphics". The core idea involves finding basis vectors and writing the rotation as a
    // linear combination of these vectors.  The equation used here is a transpose of the one derived in Lengyel's
    // book, since multiplication is done here with the vector on the LHS and the matrix on the RHS.

    // Need proper epsilon
    DEBUG_ASSERT( fabs( oAxis.GetLength() - Float(1.0) ) < 1E-6,
                  "[SMatrix3x3::BuildRotationAboutAxis]: ERROR, oAxis needs to be a UNIT VECTOR\n" );
    
    FLOAT fCosTheta = cosf( fAngle );
    FLOAT fSinTheta = sinf( fAngle );
    m[0][0] = fCosTheta + (1 - fCosTheta) * (oAxis.m_fX * oAxis.m_fX);
    m[1][0] = (1 - fCosTheta)*(oAxis.m_fX * oAxis.m_fY) + (fSinTheta * oAxis.m_fZ);
    m[2][0] = (1 - fCosTheta)*(oAxis.m_fX * oAxis.m_fZ) - (fSinTheta * oAxis.m_fY);

    m[0][1] = (1 - fCosTheta)*(oAxis.m_fX * oAxis.m_fY) - (fSinTheta * oAxis.m_fZ);
    m[1][1] = fCosTheta + (1 - fCosTheta) * (oAxis.m_fY * oAxis.m_fY);
    m[2][1] = (1 - fCosTheta)*(oAxis.m_fY * oAxis.m_fZ) + (fSinTheta * oAxis.m_fX);

    m[0][2] = (1 - fCosTheta)*(oAxis.m_fX * oAxis.m_fZ) + (fSinTheta * oAxis.m_fY);
    m[1][2] = (1 - fCosTheta)*(oAxis.m_fY * oAxis.m_fZ) - (fSinTheta * oAxis.m_fX);
    m[2][2] = fCosTheta + (1 - fCosTheta) * (oAxis.m_fZ * oAxis.m_fZ);
  }
  //----------------------------------------------------------------------------------------------
  // Public : BuildPassiveEulerMatrix
  // Range:  Phi = [0, 2Pi], Theta = [0, Pi], Psi = [0, 2Pi]
  //----------------------------------------------------------------------------------------------
  void SMatrix3x3::BuildPassiveEulerMatrix( FLOAT fPhi, FLOAT fTheta, FLOAT fPsi)
  {
    FLOAT fCosX = cosf( fPhi ), fCosY = cosf( fTheta ), fCosZ = cosf( fPsi );
    FLOAT fSinX = sinf( fPhi ), fSinY = sinf( fTheta ), fSinZ = sinf( fPsi );

    DEBUG_ASSERT(fTheta <= PI, "BuildEulerRotationMatrix: Theta out of range");
    
    m[0][0] = fCosZ * fCosX - fCosY * fSinX * fSinZ;
    m[1][0] = -fSinZ * fCosX - fCosY * fSinX * fCosZ;
    m[2][0] = fSinY * fSinX;

    m[0][1] = fCosZ * fSinX + fCosY * fCosX * fSinZ;
    m[1][1] = -fSinZ * fSinX + fCosY * fCosX * fCosZ;
    m[2][1] = -fSinY * fCosX;

    m[0][2] = fSinZ * fSinY;
    m[1][2] = fCosZ * fSinY;
    m[2][2] = fCosY;


  }

  //----------------------------------------------------------------------------------------------
  // BuildActiveSmallRotation
  //  -- this is made for infinitesimal approximation
  //----------------------------------------------------------------------------------------------
  void SMatrix3x3::BuildActiveSmallRotation(FLOAT fPhi, FLOAT fTheta, FLOAT fPsi)
  {
    m[0][0] =   Float( 1 );
    m[1][0] =   fPsi      ;
    m[2][0] = - fTheta    ;

    m[0][1] = - fPsi      ;
    m[1][1] =   Float( 1 );
    m[2][1] =   fPhi      ;

    m[0][2] =   fTheta    ;
    m[1][2] = - fPhi      ;
    m[2][2] =   Float( 1 );
  }
  
  //----------------------------------------------------------------------------------------------
  // Public : BuildActiveEulerMatrix
  // Range:  Phi = [0, 2Pi], Theta = [0, Pi], Psi = [0, 2Pi]
  //----------------------------------------------------------------------------------------------
  void SMatrix3x3::BuildActiveEulerMatrix(FLOAT fPhi, FLOAT fTheta, FLOAT fPsi)
  {

    FLOAT fCosPhi = cosf( fPhi ), fCosTheta = cosf( fTheta ), fCosPsi = cosf( fPsi );
    FLOAT fSinPhi = sinf( fPhi ), fSinTheta = sinf( fTheta ), fSinPsi = sinf( fPsi );

    DEBUG_ASSERT(fTheta <= PI, "BuildEulerRotationMatrix: Theta out of range");

    m[0][0] = fCosPhi * fCosPsi - fSinPhi * fCosTheta * fSinPsi;
    m[1][0] = fSinPhi * fCosPsi + fCosPhi * fCosTheta * fSinPsi;
    m[2][0] = fSinTheta * fSinPsi;

    m[0][1] = -fCosPhi * fSinPsi - fSinPhi * fCosTheta * fCosPsi;
    m[1][1] = -fSinPhi * fSinPsi + fCosPhi * fCosTheta * fCosPsi;
    m[2][1] = fSinTheta * fCosPsi;

    m[0][2] = fSinPhi * fSinTheta;
    m[1][2] = -fCosPhi * fSinTheta;
    m[2][2] = fCosTheta;

  }
  
  //----------------------------------------------------------------------------------------------
  // Public : GetEulerAngles
  // return a vector of Euler angles in ZYZ convention in radians
  //
  //  This is taken from Bob's code
  //----------------------------------------------------------------------------------------------
  SVector3 SMatrix3x3::GetEulerAngles() const
  {
    SVector3 oEulerAngles;

    Float fCosThresh = 0.999999; 
    

    if( m[2][2] > fCosThresh )	                    //	 is chi approx 0.?
    {
      oEulerAngles.m_fX = Float( 0 );				// set omega and chi to zero
      oEulerAngles.m_fY = Float( 0 );
      oEulerAngles.m_fZ = atan2( m[1][0], m[0][0] );
    }
    else if ( m[2][2] < - fCosThresh )              //  is chi approx pi?
    {
      oEulerAngles.m_fX = Float( 0 );
      oEulerAngles.m_fY = PI;
      oEulerAngles.m_fZ = atan2( m[0][1], m[0][0] );
    }
    else                                            //  chi is not zero or pi
    { 
      oEulerAngles.m_fX = atan2( m[0][2], -m[1][2] );
      oEulerAngles.m_fY = atan2( sqrt( m[2][0] * m[2][0] +m[2][1] * m[2][1] ), m[2][2] );
      oEulerAngles.m_fZ = atan2( m[2][0], m[2][1] );
    }

    // bring back to proper region
    for (Int i = 0; i < 3; i ++ ) //			  % atan2 returns in [-pi:pi], we want [0:2pi]
      if( oEulerAngles[i] < Float (0 ) )
        oEulerAngles[i]+= 2 * PI;
    
 
    
    return oEulerAngles;
  }
  
  //----------------------------------------------------------------------------------------------
  // Public : BuildProjectionMatrix  -- oProjDir is a unit vector 
  //----------------------------------------------------------------------------------------------  
  void SMatrix3x3::BuildProjectionMatrix(const SVector3 &oProjDir )
  {
    m[0][0] = oProjDir.m_fX * oProjDir.m_fX;
    m[1][1] = oProjDir.m_fY * oProjDir.m_fY;
    m[2][2] = oProjDir.m_fZ * oProjDir.m_fZ;

    m[0][1] = m[1][0] = oProjDir.m_fX * oProjDir.m_fY;
    m[0][2] = m[2][0] = oProjDir.m_fX * oProjDir.m_fZ;
    m[1][2] = m[2][1] = oProjDir.m_fY * oProjDir.m_fZ;
  }
  
  //----------------------------------------------------------------------------------------------
  // Public : Trace
  //----------------------------------------------------------------------------------------------  
  FLOAT SMatrix3x3::Trace() const
  {
    FLOAT fRet;
    fRet = m[0][0] + m[1][1] + m[2][2];
    return fRet;
  }

  //----------------------------------------------------------------------------------------------
  // Public : Transpose
  //----------------------------------------------------------------------------------------------
  void SMatrix3x3::Transpose()
  {
    // diagonals do not change, and we only have to swap
    // the upper triangle
    for ( int nRow = 0; nRow < 2; nRow++ )
    {
      for ( int nCol = nRow +1; nCol < 3; nCol++ )
      {
        // Swap
        FLOAT fTemp = m[nRow][nCol];
        m[nRow][nCol] = m[nCol][nRow];
        m[nCol][nRow] = fTemp;
      }
    }
  }

  //----------------------------------------------------------------------------------------------
  //  operator* float
  //----------------------------------------------------------------------------------------------
  SMatrix3x3 SMatrix3x3::operator*( Float f ) const
  {
    SMatrix3x3 oRes;

    oRes.m[0][0] =  m[0][0] * f;
    oRes.m[0][1] =  m[0][1] * f;
    oRes.m[0][2] =  m[0][2] * f;
    
    oRes.m[1][0] =  m[1][0] * f;
    oRes.m[1][1] =  m[1][1] * f;
    oRes.m[1][2] =  m[1][2] * f;
  
    oRes.m[2][0] =  m[2][0] * f;
    oRes.m[2][1] =  m[2][1] * f;
    oRes.m[2][2] =  m[2][2] * f;
    
    return oRes;
  }

  //----------------------------------------------------------------------------------------------
  //  operator *= float
  //----------------------------------------------------------------------------------------------
  SMatrix3x3 & SMatrix3x3::operator*= ( Float f )
  {
    m[0][0] *= f;
    m[0][1] *= f;
    m[0][2] *= f;
   
    m[1][0] *= f;
    m[1][1] *= f;
    m[1][2] *= f;
  
    m[2][0] *= f;
    m[2][1] *= f;
    m[2][2] *= f;
    
    return *this; 
  }

    //----------------------------------------------------------------------------------------------
  //  operator/ float
  //----------------------------------------------------------------------------------------------
  SMatrix3x3 SMatrix3x3::operator/( Float f ) const
  {
    SMatrix3x3 oRes;

    oRes.m[0][0] =  m[0][0] / f;
    oRes.m[0][1] =  m[0][1] / f;
    oRes.m[0][2] =  m[0][2] / f;
    
    oRes.m[1][0] =  m[1][0] / f;
    oRes.m[1][1] =  m[1][1] / f;
    oRes.m[1][2] =  m[1][2] / f;
 
    oRes.m[2][0] =  m[2][0] / f;
    oRes.m[2][1] =  m[2][1] / f;
    oRes.m[2][2] =  m[2][2] / f;
    
    return oRes;
  }

  //----------------------------------------------------------------------------------------------
  //  operator /= float
  //----------------------------------------------------------------------------------------------
  SMatrix3x3 & SMatrix3x3::operator/= ( Float f )
  {
    m[0][0] /= f;
    m[0][1] /= f;
    m[0][2] /= f;
   
    m[1][0] /= f;
    m[1][1] /= f;
    m[1][2] /= f;
  
    m[2][0] /= f;
    m[2][1] /= f;
    m[2][2] /= f;
    
    return *this; 
  }

  //----------------------------------------------------------------------------------------------
  // Public : operator+
  //----------------------------------------------------------------------------------------------
  SMatrix3x3 SMatrix3x3::operator+( const SMatrix3x3& oRHS ) const
  {
    SMatrix3x3 oRes;

    oRes.m[0][0] =  m[0][0] + oRHS.m[0][0];
    oRes.m[0][1] =  m[0][1] + oRHS.m[0][1];
    oRes.m[0][2] =  m[0][2] + oRHS.m[0][2];
   
    oRes.m[1][0] =  m[1][0] + oRHS.m[1][0];
    oRes.m[1][1] =  m[1][1] + oRHS.m[1][1];
    oRes.m[1][2] =  m[1][2] + oRHS.m[1][2];
  
    oRes.m[2][0] =  m[2][0] + oRHS.m[2][0];
    oRes.m[2][1] =  m[2][1] + oRHS.m[2][1];
    oRes.m[2][2] =  m[2][2] + oRHS.m[2][2];
    
    return oRes;
  }
  
  //----------------------------------------------------------------------------------------------
  // Public : operator+=
  //----------------------------------------------------------------------------------------------
  SMatrix3x3 & SMatrix3x3::operator+=( const SMatrix3x3& oRHS ) 
  {
    m[0][0] +=  oRHS.m[0][0];
    m[0][1] +=  oRHS.m[0][1];
    m[0][2] +=  oRHS.m[0][2];
   
    m[1][0] +=  oRHS.m[1][0];
    m[1][1] +=  oRHS.m[1][1];
    m[1][2] +=  oRHS.m[1][2];
  
    m[2][0] +=  oRHS.m[2][0];
    m[2][1] +=  oRHS.m[2][1];
    m[2][2] +=  oRHS.m[2][2];
    
    return *this;
  }
  
  //----------------------------------------------------------------------------------------------
  // Public : operator-
  //----------------------------------------------------------------------------------------------
  SMatrix3x3 SMatrix3x3::operator-( const SMatrix3x3& oRHS ) const
  {
    SMatrix3x3 oRes;

    oRes.m[0][0] =  m[0][0] - oRHS.m[0][0];
    oRes.m[0][1] =  m[0][1] - oRHS.m[0][1];
    oRes.m[0][2] =  m[0][2] - oRHS.m[0][2];
   
    oRes.m[1][0] =  m[1][0] - oRHS.m[1][0];
    oRes.m[1][1] =  m[1][1] - oRHS.m[1][1];
    oRes.m[1][2] =  m[1][2] - oRHS.m[1][2];
  
    oRes.m[2][0] =  m[2][0] - oRHS.m[2][0];
    oRes.m[2][1] =  m[2][1] - oRHS.m[2][1];
    oRes.m[2][2] =  m[2][2] - oRHS.m[2][2];
    
    return oRes;
  }

  //----------------------------------------------------------------------------------------------
  // Public : operator*
  //----------------------------------------------------------------------------------------------
  SMatrix3x3 SMatrix3x3::operator*( const SMatrix3x3& oRHS ) const
  {
    SMatrix3x3 oRes;

    // Matrix multiplication
    oRes.m[0][0] = m[0][0] * oRHS.m[0][0] + m[0][1] * oRHS.m[1][0] + m[0][2] * oRHS.m[2][0];
    oRes.m[0][1] = m[0][0] * oRHS.m[0][1] + m[0][1] * oRHS.m[1][1] + m[0][2] * oRHS.m[2][1];
    oRes.m[0][2] = m[0][0] * oRHS.m[0][2] + m[0][1] * oRHS.m[1][2] + m[0][2] * oRHS.m[2][2];
   

    oRes.m[1][0] = m[1][0] * oRHS.m[0][0] + m[1][1] * oRHS.m[1][0] + m[1][2] * oRHS.m[2][0];
    oRes.m[1][1] = m[1][0] * oRHS.m[0][1] + m[1][1] * oRHS.m[1][1] + m[1][2] * oRHS.m[2][1];
    oRes.m[1][2] = m[1][0] * oRHS.m[0][2] + m[1][1] * oRHS.m[1][2] + m[1][2] * oRHS.m[2][2];
  

    oRes.m[2][0] = m[2][0] * oRHS.m[0][0] + m[2][1] * oRHS.m[1][0] + m[2][2] * oRHS.m[2][0];
    oRes.m[2][1] = m[2][0] * oRHS.m[0][1] + m[2][1] * oRHS.m[1][1] + m[2][2] * oRHS.m[2][1];
    oRes.m[2][2] = m[2][0] * oRHS.m[0][2] + m[2][1] * oRHS.m[1][2] + m[2][2] * oRHS.m[2][2];
    
    return oRes;
  }

  //----------------------------------------------------------------------------------------------
  // Public : operator*
  //----------------------------------------------------------------------------------------------
  SVector3 SMatrix3x3::operator*( const SVector3 &oRHS ) const
  {
    SVector3 oRes;

    oRes.m_fX = m[0][0] * oRHS.m_fX + m[0][1] * oRHS.m_fY + m[0][2] * oRHS.m_fZ;
    oRes.m_fY = m[1][0] * oRHS.m_fX + m[1][1] * oRHS.m_fY + m[1][2] * oRHS.m_fZ;
    oRes.m_fZ = m[2][0] * oRHS.m_fX + m[2][1] * oRHS.m_fY + m[2][2] * oRHS.m_fZ;
    
    return oRes;
    
  }

  //----------------------------------------------------------------------------------------------
  // Public : operator* (left multiply)
  //----------------------------------------------------------------------------------------------
  SMatrix3x3 operator*( Float fLHS, const SMatrix3x3 &oRHS )
  {
    return oRHS * fLHS;
  }
  
  //----------------------------------------------------------------------------------------------
  // operator<<
  //----------------------------------------------------------------------------------------------
  std::ostream & operator<< ( std::ostream & os, const SMatrix3x3 &m )
  {
    for ( int i = 0; i < 3; i ++ )
    {
      for ( int j = 0; j < 3; j ++ )
      {
        if( fabs( m.m[i][j] ) > 0.00001 )
          os << m.m[i][j] << " ";
        else
          os << 0 << " ";
      }
      os << std::endl;
    }
    return os;
  }
  
  //----------------------------------------------------------------------------------------------
  // operator<<
  //----------------------------------------------------------------------------------------------
  std::ostream & operator<< ( std::ostream & os, const SVector3 &v )
  {
    os << v.m_fX << " " << v.m_fY << " " << v.m_fZ;
    return os;
  }
  
  //----------------------------------------------------------------------------------------------
  // Public : SVector3
  //----------------------------------------------------------------------------------------------
  SVector3::SVector3()
  {}

  //----------------------------------------------------------------------------------------------
  // Public : SVector3
  //----------------------------------------------------------------------------------------------
  SVector3::SVector3( FLOAT fX, FLOAT fY, FLOAT fZ )
  {
    Set( fX, fY, fZ );
  }

  //----------------------------------------------------------------------------------------------
  // Public : Set
  //----------------------------------------------------------------------------------------------
  void SVector3::Set( FLOAT fX, FLOAT fY, FLOAT fZ )
  {
    this->m_fX = fX;
    this->m_fY = fY;
    this->m_fZ = fZ;
  }
  
  //----------------------------------------------------------------------------------------------
  // Public : Get
  //----------------------------------------------------------------------------------------------
  FLOAT SVector3::Get( UInt nAxis ) const
  {
    if ( nAxis == 0 )
      return this->m_fX;
    else if ( nAxis == 1 )
      return this->m_fY;
    else
      return this->m_fZ;
  }

  //----------------------------------------------------------------------------------------------
  // Public : Set
  //----------------------------------------------------------------------------------------------
  void SVector3::Set( UInt nAxis, FLOAT fVal )
  {
    if ( nAxis == 0 )
      this->m_fX = fVal;
    else if ( nAxis == 1 )
      this->m_fY = fVal;
    else
      this->m_fZ = fVal;
  }
  
  //----------------------------------------------------------------------------------------------
  // Public : GetLength
  //----------------------------------------------------------------------------------------------
  FLOAT SVector3::GetLength() const
  {
    return sqrtf(m_fX * m_fX + m_fY * m_fY + m_fZ * m_fZ);
  }
  
  //----------------------------------------------------------------------------------------------
  // Public : Normalize
  //----------------------------------------------------------------------------------------------
  void SVector3::Normalize()
  {
    FLOAT fLength = GetLength();
    m_fX /= fLength;
    m_fY /= fLength;
    m_fZ /= fLength;
  }

  //----------------------------------------------------------------------------------------------
  // Public : operator=
  //----------------------------------------------------------------------------------------------
  SVector3 &SVector3::operator=( const SVector3 &oRHS )
  {
    m_fX = oRHS.m_fX;
    m_fY = oRHS.m_fY;
    m_fZ = oRHS.m_fZ;
    return *this;
  }

  //----------------------------------------------------------------------------------------------
  //  Public : operator+
  //----------------------------------------------------------------------------------------------
  SVector3 SVector3::operator+( const SVector3& oRHS ) const
  {
    SVector3 oRetObj;
    oRetObj.m_fX = m_fX + oRHS.m_fX;
    oRetObj.m_fY = m_fY + oRHS.m_fY;
    oRetObj.m_fZ = m_fZ + oRHS.m_fZ;
    return oRetObj;
  }

  //----------------------------------------------------------------------------------------------
  //  Public : operator+= 
  //----------------------------------------------------------------------------------------------
  SVector3 &SVector3::operator+=( const SVector3 &oRHS )
  {
    m_fX += oRHS.m_fX;
    m_fY += oRHS.m_fY;
    m_fZ += oRHS.m_fZ;
    return *this;
  }
  
  //----------------------------------------------------------------------------------------------
  //  operator -
  //----------------------------------------------------------------------------------------------
  SVector3 SVector3::operator-( ) const
  {
    SVector3 oRetObj;
    oRetObj.m_fX = - m_fX;
    oRetObj.m_fY = - m_fY;
    oRetObj.m_fZ = - m_fZ;
    return oRetObj;
  }
  //----------------------------------------------------------------------------------------------
  // Public : operator-
  //----------------------------------------------------------------------------------------------
  SVector3 SVector3::operator-( const SVector3& oRHS ) const
  {
    SVector3 oRetObj;
    oRetObj.m_fX = m_fX - oRHS.m_fX;
    oRetObj.m_fY = m_fY - oRHS.m_fY;
    oRetObj.m_fZ = m_fZ - oRHS.m_fZ;
    return oRetObj;
  }

  //----------------------------------------------------------------------------------------------
  //  Public : operator-= 
  //----------------------------------------------------------------------------------------------
  SVector3 &SVector3::operator-=( const SVector3 &oRHS )
  {
    m_fX -= oRHS.m_fX;
    m_fY -= oRHS.m_fY;
    m_fZ -= oRHS.m_fZ;
    return *this;
  }

  //----------------------------------------------------------------------------------------------
  // Public : operator*
  //----------------------------------------------------------------------------------------------
  SVector3 SVector3::operator*( FLOAT fRHS ) const
  {
    SVector3 oRetObj;
    oRetObj.m_fX = m_fX * fRHS;
    oRetObj.m_fY = m_fY * fRHS;
    oRetObj.m_fZ = m_fZ * fRHS;
    return oRetObj;
  }
  
  //----------------------------------------------------------------------------------------------
  // Public : operator*  (left multiply)
  //----------------------------------------------------------------------------------------------
  SVector3 operator*( FLOAT fLHS, const SVector3 &oRHS )
  {
    SVector3 oRetObj;
    oRetObj.m_fX = fLHS * oRHS.m_fX;
    oRetObj.m_fY = fLHS * oRHS.m_fY;
    oRetObj.m_fZ = fLHS * oRHS.m_fZ;
    return oRetObj;
  }


  //----------------------------------------------------------------------------------------------
  // Public : operator*=
  //----------------------------------------------------------------------------------------------
  SVector3 &SVector3::operator*=( FLOAT fRHS )
  {
    m_fX *= fRHS;
    m_fY*= fRHS;
    m_fZ *= fRHS;
    return *this;
  }


  //----------------------------------------------------------------------------------------------
  // Public : operator/
  //----------------------------------------------------------------------------------------------
  SVector3 SVector3::operator/( FLOAT fRHS ) const
  {
    SVector3 oRetObj;
    oRetObj.m_fX = m_fX / fRHS;
    oRetObj.m_fY = m_fY / fRHS;
    oRetObj.m_fZ = m_fZ / fRHS;
    return oRetObj;
  }

  //----------------------------------------------------------------------------------------------
  // Public : operator/=
  //----------------------------------------------------------------------------------------------
  SVector3 &SVector3::operator/=( FLOAT fRHS )
  {
    m_fX /= fRHS;
    m_fY /= fRHS;
    m_fZ /= fRHS;
    return *this;
  }

  //----------------------------------------------------------------------------------------------
  // Public : operator[]
  //----------------------------------------------------------------------------------------------
  FLOAT SVector3::operator[]( Int i) const
  {

    DEBUG_ASSERT( i >= 0 && i < 3, "ERROR:  index out of range for SVector3\n");
    switch(i)
    {
      case 0:
        return m_fX;
      case 1:
        return m_fY;
      case 2:
        return m_fZ;
    };
    // error branch  (throw exception in the future)
    return std::numeric_limits<Float>::signaling_NaN();
  }
  
  //----------------------------------------------------------------------------------------------
  // Public : operator[]
  //----------------------------------------------------------------------------------------------
  FLOAT& SVector3::operator[](Int i)
  {
    DEBUG_ASSERT( i >= 0 && i < 3, "ERROR:  index out of range for SVector3\n");
    switch(i)
    {
      case 0:
        return m_fX;
      case 1:
        return m_fY;
      case 2:
        return m_fZ;
    };
    // error branch  (throw exception in the future)
    return m_fX;
  }

  //----------------------------------------------------------------------------------------------
  // Cross
  //----------------------------------------------------------------------------------------------
  SVector3 Cross( const SVector3 &oLHS, const SVector3 &oRHS )
  {
    SVector3 oRes;

    oRes.m_fX = oLHS.m_fY * oRHS.m_fZ - oLHS.m_fZ * oRHS.m_fY;
    oRes.m_fY = oLHS.m_fZ * oRHS.m_fX - oLHS.m_fX * oRHS.m_fZ;
    oRes.m_fZ = oLHS.m_fX * oRHS.m_fY - oLHS.m_fY * oRHS.m_fX;

    return oRes;
  }

  //----------------------------------------------------------------------------------------------
  // Dot
  //----------------------------------------------------------------------------------------------
  FLOAT Dot( const SVector3 &oLHS, const SVector3 &oRHS )
  {
    return oLHS.m_fX * oRHS.m_fX + oLHS.m_fY * oRHS.m_fY + oLHS.m_fZ * oRHS.m_fZ;
  }

  //----------------------------------------------------------------------------------------------
  //  OuterProduct
  //----------------------------------------------------------------------------------------------
  SMatrix3x3 OuterProduct( const SVector3 & oLHS, const SVector3 & oRHS )
  {
    SMatrix3x3 oRes;
    oRes.m[0][0] = oLHS.m_fX * oRHS.m_fX;
    oRes.m[0][1] = oLHS.m_fX * oRHS.m_fY;
    oRes.m[0][2] = oLHS.m_fX * oRHS.m_fZ;

    oRes.m[1][0] = oLHS.m_fY * oRHS.m_fX;
    oRes.m[1][1] = oLHS.m_fY * oRHS.m_fY;
    oRes.m[1][2] = oLHS.m_fY * oRHS.m_fZ;
    
    oRes.m[2][0] = oLHS.m_fZ * oRHS.m_fX;
    oRes.m[2][1] = oLHS.m_fZ * oRHS.m_fY;
    oRes.m[2][2] = oLHS.m_fZ * oRHS.m_fZ;

    return oRes;
  }
  
  //----------------------------------------------------------------------------------------------
  // DegreeToRadian
  //----------------------------------------------------------------------------------------------
  SVector3 DegreeToRadian( const SVector3 & oRHS )
  {
    SVector3 oRes;
    oRes.m_fX =  DEGREE_TO_RADIAN( oRHS.m_fX );
    oRes.m_fY =  DEGREE_TO_RADIAN( oRHS.m_fY );
    oRes.m_fZ =  DEGREE_TO_RADIAN( oRHS.m_fZ );
    return oRes;
  }
  
  //----------------------------------------------------------------------------------------------
  // RadianToDegree
  //----------------------------------------------------------------------------------------------
  SVector3 RadianToDegree( const SVector3 & oRHS )
  {
    SVector3 oRes;
    oRes.m_fX =  RADIAN_TO_DEGREE( oRHS.m_fX );
    oRes.m_fY =  RADIAN_TO_DEGREE( oRHS.m_fY );
    oRes.m_fZ =  RADIAN_TO_DEGREE( oRHS.m_fZ );
    return oRes;
  }
  

    
    //----------------------------
    //  InnerProduct
    //----------------------------
    Float InnerProduct( const EigenRep  & LHS, const EigenRep  & RHS )
    {
      Float Result = 0;
      for( int i = 0; i < EIGEN_DIM; i ++ )
        Result += LHS(i) * RHS(i);
      
      return Result;
    }
    
    //----------------------------
    //  InnerProduct
    //----------------------------
    Float InnerProduct( const SMatrix3x3 & LHS, const SMatrix3x3 & RHS )
    {
      Float Result = 0;
      for( int i = 0; i < 3; i ++ )
        for( int j = 0; j < 3; j ++ )
          Result += LHS.m[i][j] * RHS.m[i][j];

      return Result;
    }

    //----------------------------
    //  InnerProduct
    //----------------------------
    Float  InnerProduct( const Tensor4Rank & LHS, const Tensor4Rank & RHS )
    {
      Float Result = 0;
      for( int i = 0; i < 3; i ++ )
        for( int j = 0; j < 3; j ++ )
          for( int k = 0; k < 3; k ++ )
            for( int l = 0; l < 3; l ++ )
              Result += LHS(i, j, k, l) * RHS(i, j, k, l);
      return Result;
    }
    
    //----------------------------
    // OuterProduct
    //----------------------------
    SMatrix5x5 OuterProduct( const EigenRep & LHS,   const EigenRep & RHS )
    {
      SMatrix5x5 oRes;

      for( int i = 0; i < EIGEN_DIM; i ++ )
        for( int j = 0; j < EIGEN_DIM; j ++ )
          oRes(i, j) = LHS(i) * RHS(j);
      
      return oRes;
    }

    //----------------------------
    //  OuterProduct
    //----------------------------
    Tensor4Rank OuterProduct( const SMatrix3x3 & LHS, const SMatrix3x3 & RHS )
    {
      Tensor4Rank oRes;
      for( int i = 0; i < 3; i ++ )
        for( int j = 0; j < 3; j ++ )
          for( int k = 0; k < 3; k ++ )
            for( int l = 0; l < 3; l ++ )
              oRes(i, j, k, l) = LHS.m[i][j] * RHS.m[k][l];
      
      return oRes;
    }
    
    //----------------------------
    //  Symmetrize
    //----------------------------
    Tensor4Rank Symmetrize( const Tensor4Rank & T )
    {
      Tensor4Rank oRes;
      oRes.SetZero();
      for( int i = 0; i < 3; i ++ )
        for( int j = 0; j < 3; j ++ )
          for( int k = 0; k < 3; k ++ )
            for( int l = 0; l < 3; l ++ )
              oRes(i, j, k, l) = ( T(i, j, k, l) + T(j, i, k, l)
                               +   T(i, j, l, k) + T(j, i, l, k) ) / static_cast<Float>( 4 ); 
      return oRes;
    }


    //----------------------------
    //  AntiSymmetrize
    //----------------------------
    Tensor4Rank AntiSymmetrize( const Tensor4Rank & T )
    {
      RUNTIME_ASSERT( 0, "AntiSymmetrizer not yet implemented" );
      Tensor4Rank oRes;
      oRes.SetZero();
      for( int i = 0; i < 3; i ++ )
        for( int j = 0; j < 3; j ++ )
          for( int k = 0; k < 3; k ++ )
            for( int l = 0; l < 3; l ++ )
              oRes(i, j, k, l) = ( T(i, j, k, l) + T(j, i, k, l)
                               +   T(i, j, l, k) + T(j, i, l, k) ) / static_cast<Float>( 4 ); 
      return oRes;
    }

    //----------------------------
    //  operator << 
    //----------------------------
    std::ostream & operator<< ( std::ostream & os, const EigenRep &v )
    {
      os << v.m[0] << " "
         << v.m[1] << " "
         << v.m[2] << " "
         << v.m[3] << " "
         << v.m[4] << " ";
      return os;
    }
    std::ostream & operator<< ( std::ostream & os, const Tensor4Rank &T )
    {
      for( int i = 0; i < 3; i ++ )
        for( int j = 0; j < 3; j ++ )
        {
          for( int k = 0; k < 3; k ++ )
          {
            os << "| ";
            for( int l = 0; l < 3; l ++ )
            {
              os << T(i, j, k, l) << " ";
            }
          }
          os << std::endl;
        }
      
      return os;
    }

    //----------------------------------------------------------------------------------------------
    //  Matrix Class
    //
    //----------------------------------------------------------------------------------------------
    
    
    
    //----------------------------------------------------------------------------------------------
    // Public : SMatrix5x5
    //----------------------------------------------------------------------------------------------
    SMatrix5x5::SMatrix5x5( Float  pMatrix[5][5])
    {
      memcpy( m, pMatrix, sizeof( Float) * 25 );
    }

    //----------------------------------------------------------------------------------------------
    // Public : SetZero
    //----------------------------------------------------------------------------------------------
    void SMatrix5x5::SetZero()
    {
      memset( m, 0, sizeof(Float) * 25 );
    }

    //----------------------------------------------------------------------------------------------
    // Public:  SetIdentity
    //----------------------------------------------------------------------------------------------
    void SMatrix5x5::SetIdentity()
    {
      memset( m, 0, sizeof(Float) * 25 );
      operator()( 0, 0 ) = 1;
      operator()( 1, 1 ) = 1;
      operator()( 2, 2 ) = 1;
      operator()( 3, 3 ) = 1;
      operator()( 4, 4 ) = 1;
    }

    //----------------------------------------------------------------------------------------------
    // Public : Transpose
    //----------------------------------------------------------------------------------------------
    void SMatrix5x5::Transpose()
    {
      // diagonals do not change, and we only have to swap
      // the upper triangle
      for ( int nRow = 0; nRow < 2; nRow++ )
        for ( int nCol = nRow +1; nCol < 5; nCol++ )
          std::swap( this->operator()(nRow, nCol), this->operator()(nCol, nRow) );
    }

    //----------------------------------------------------------------------------------------------
    // Public : Transpose
    //----------------------------------------------------------------------------------------------
    void SMatrix5x5::MovingAverage( const SMatrix5x5 & RHS, int CurrentCount )
    {
      for( int i = 0; i < 5; i ++ )
        for( int j = 0; j < 5; j ++ )
        {
          operator()( i, j ) += ( RHS(i, j) - operator()(i, j) ) / static_cast<Float>( CurrentCount + 1 );
        }
    }
    
    //----------------------------------------------------------------------------------------------
    //  operator* float
    //----------------------------------------------------------------------------------------------
    SMatrix5x5 SMatrix5x5::operator*( Float f ) const
    {
      SMatrix5x5 oRes;

      for ( int i = 0; i < 5; i ++ )
        for ( int j = 0; j < 5; j ++ )
          oRes( i, j ) = operator()( i, j ) * f;      
      return oRes;
    }

    //----------------------------------------------------------------------------------------------
    //  operator *= float
    //----------------------------------------------------------------------------------------------
    SMatrix5x5 & SMatrix5x5::operator*= ( Float f )
    {
      for ( int i = 0; i < 5; i ++ )
        for ( int j = 0; j < 5; j ++ )
          operator()(i, j) *= f;
      return *this; 
    }

    //----------------------------------------------------------------------------------------------
    //  operator/ float
    //----------------------------------------------------------------------------------------------
    SMatrix5x5 SMatrix5x5::operator/( Float f ) const
    {
//       if( fabs( f ) < 1e-4 )
//       {
//         std::cout << "ERROR " << std::endl;
//       }
      SMatrix5x5 oRes;
      for ( int i = 0; i < 5; i ++ )
        for ( int j = 0; j < 5; j ++ )
          oRes( i, j ) = operator()( i, j ) / f;
      return oRes;
    }

    //----------------------------------------------------------------------------------------------
    //  operator /= float
    //----------------------------------------------------------------------------------------------
    SMatrix5x5 & SMatrix5x5::operator/= ( Float f )
    {
      for ( int i = 0; i < 5; i ++ )
        for ( int j = 0; j < 5; j ++ )
          operator()( i, j ) /= f;
      
      return *this; 
    }

    //----------------------------------------------------------------------------------------------
    // Public : operator+
    //----------------------------------------------------------------------------------------------
    SMatrix5x5 SMatrix5x5::operator+( const SMatrix5x5& oRHS ) const
    {
      SMatrix5x5 oRes;
      for ( int i = 0; i < 5; i ++ )
        for ( int j = 0; j < 5; j ++ )
          oRes(i, j) = operator()( i, j ) + oRHS(i, j ); 
      
      return oRes;
    }
  
    //----------------------------------------------------------------------------------------------
    // Public : operator+=
    //----------------------------------------------------------------------------------------------
    SMatrix5x5 & SMatrix5x5::operator+=( const SMatrix5x5& oRHS ) 
    {
      for ( int i = 0; i < 5; i ++ )
        for ( int j = 0; j < 5; j ++ )
          operator()( i, j ) += oRHS(i, j ); 
      return *this;
    }
  
    //----------------------------------------------------------------------------------------------
    // Public : operator-
    //----------------------------------------------------------------------------------------------
    SMatrix5x5 SMatrix5x5::operator-( const SMatrix5x5& oRHS ) const
    {
      SMatrix5x5 oRes;
      for ( int i = 0; i < 5; i ++ )
        for ( int j = 0; j < 5; j ++ )
          oRes(i, j) = operator()( i, j ) - oRHS(i, j ); 
      return oRes;
    }


    //----------------------------------------------------------------------------------------------
    // Public : operator-
    //----------------------------------------------------------------------------------------------
    SMatrix5x5 SMatrix5x5::operator-( ) const
    {
      SMatrix5x5 oRes;
      for ( int i = 0; i < 5; i ++ )
        for ( int j = 0; j < 5; j ++ )
          oRes(i, j) = -operator()( i, j );
      return oRes;
    }

    
    //----------------------------------------------------------------------------------------------
    // Public : operator-=
    //----------------------------------------------------------------------------------------------
    SMatrix5x5 & SMatrix5x5::operator-=( const SMatrix5x5& oRHS ) 
    {
      for ( int i = 0; i < 5; i ++ )
        for ( int j = 0; j < 5; j ++ )
          operator()( i, j ) -= oRHS(i, j ); 
      return *this;
    }

    //----------------------------------------------------------------------------------------------
    // Public : operator*  -- clearly not made to be fast --
    //----------------------------------------------------------------------------------------------
    SMatrix5x5 SMatrix5x5::operator*( const SMatrix5x5& oRHS ) const
    {
      SMatrix5x5 oRes;
      oRes.SetZero();
      for( int i = 0; i < 5; i ++ )
        for( int j = 0; j < 5; j ++ )
          for( int k = 0; k < 5; k ++ )
            oRes(i, j) += operator()( i, k ) * oRHS(k, j ); 
      return oRes;
    }

    //----------------------------------------------------------------------------------------------
    // Public : operator* (left multiply)
    //----------------------------------------------------------------------------------------------
    EigenRep operator*( const SMatrix5x5 & LHS, const EigenRep & RHS)
    {
      EigenRep oRes( 0, 0, 0, 0, 0 );
      for( int i = 0; i < 5; i ++ )
        for( int j = 0; j < 5; j ++ )
          oRes(i) += LHS( i, j ) * RHS(j); 
      return oRes;
    }


    //----------------------------------------------------------------------------------------------
    // Public : operator* (left multiply)
    //----------------------------------------------------------------------------------------------
    EigenRep operator*(const EigenRep & LHS,  const SMatrix5x5 & RHS)
    {
      EigenRep oRes( 0, 0, 0, 0, 0 );
      for( int i = 0; i < 5; i ++ )
        for( int j = 0; j < 5; j ++ )
          oRes(j) += LHS( i ) * RHS(i, j); 
      return oRes;
    }
   
    //----------------------------------------------------------------------------------------------
    // Public : operator* (left multiply)
    //----------------------------------------------------------------------------------------------
    SMatrix5x5 operator*( Float fLHS, const SMatrix5x5 &oRHS )
    {
      return oRHS * fLHS;
    }
    
    //----------------------------------------------------------------------------------------------
    // operator<<
    //----------------------------------------------------------------------------------------------
    std::ostream & operator<< ( std::ostream & os, const SMatrix5x5 &m )
    {
      for ( int i = 0; i < 5; i ++ )
      {
        for ( int j = 0; j < 5; j ++ )
          os << m(i, j) << " ";
        os << std::endl;
      }
      return os;
    }

    
    //----------------------------------------------------------------------------------------------
    //  LU_Solver
    //
    //  Solving x in A x = b, given
    //  that A is a 5x5 matrix.  x and b are 5 x 1.
    //----------------------------------------------------------------------------------------------
    EigenRep LU_Solver( const SMatrix5x5 & A, const EigenRep & b )
    {
      using namespace Eigen;
      typedef Matrix< Float, 5, 5 > Matrix5x5;
      typedef Matrix< Float, 5, 1 > Vector5;
      
      Map< const Matrix5x5> A_Mapped( A.m );
      Map< const Vector5>   b_Mapped( b.m );

      EigenRep Solution;
      Map<Vector5>   x_Mapped( Solution.m );
      x_Mapped = A_Mapped.fullPivLu().solve( b_Mapped );
      Float relative_error = (A_Mapped * x_Mapped - b_Mapped).norm() / b_Mapped.norm(); // norm() is L2 norm
      
      //     std::cout << "LU Relative Error " << relative_error << std::endl;
      return Solution;
    }
    
  }  // namespace LinearAlgebra
}
