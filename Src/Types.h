#ifndef TYPES_H_
#define TYPES_H_

#include <string>
#include <cstdlib>
#include <cstring>
#include <mpi.h>
using std::string;

typedef int Int;
typedef unsigned int UInt;
typedef string::size_type Size_Type;

//-----------------------------------------------
typedef double Float; //!!!  <---------  CHECK HARD CODING OF MPI_REDUCE_SCATTER! (SEE BELOW AS WELL)

#ifndef VPFFT_MPI_FLOAT_TYPE
#define VPFFT_MPI_FLOAT_TYPE MPI_DOUBLE
#endif
//-----------------------------------------------

typedef double Double;
//typedef short int IntensityT;
typedef Float IntensityT;
typedef Float FLOAT;  // for 3dMath.h

typedef long int LInt;
typedef unsigned char U8;
typedef unsigned short U16;
typedef unsigned int   U32;
typedef float F32;
typedef bool Bool;
#ifndef NULL
#define NULL 0
#endif

#ifndef MAX_FLOAT
#define MAX_FLOAT 1.e38
#endif

#ifndef MIN_FLOAT
#define MIN_FLOAT -1.e37
#endif

#ifndef MAX_INT
#define MAX_INT (int)(2147483647)
#endif

#ifndef MIN_INT
#define MIN_INT (-2147483648)
#endif

#ifndef EPSILON
#define EPSILON  1.e-10
#endif










#endif
