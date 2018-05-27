
#ifndef REALMACRO_H_
#define REALMACRO_H_

#include <CL/opencl.h>
//#include <CL/cl_platform.h>
// floating point type macros : prefer 32bit for most hardware.
#define REAL float	//double
#define REAL2 cl_float2	//double2
#define REAL3 cl_float3	//double3
#define REAL4 cl_float4	//double4

#endif
