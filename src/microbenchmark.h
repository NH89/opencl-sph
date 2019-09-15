#include <sys/time.h>
#include <stdio.h>

#include "macros.h"

// Usage (see the "microbenchmark" branch of the Git repository.
// Applied to kernel launches in src/opencl/particle_system_host.c
//
// #define OPENCL_SPH_MICROBENCHMARK
// #include "../microbenchmark.h" 
// ...
// void My_function (args...){
//  timeInit(My_function) ;
//  timeStart(My_function) ;
//  ... 
//  function code   
//  ...
//  timeStop(My_function)   OR   timeStopS(My_function,  some_string ) 
// }

#ifdef OPENCL_SPH_MICROBENCHMARK
#define timeInit(marker)    struct timeval t ## marker ## 1, t ## marker ## 2;
        
#define timeStart(marker)   gettimeofday(&t ## marker ## 1,0); 
        
#define timeStop(marker)    gettimeofday(&t ## marker ## 2,0); \
		{ \
			printf("benchmark-%s\t%f\n", STRINGIFY(marker) \
				, ((float)(t ## marker ## 2 . tv_sec - t ## marker ## 1 . tv_sec) \
				+ ((float)(t ## marker ## 2 . tv_usec - t ## marker ## 1 . tv_usec))*1e-6) \
				); \
		}
		
#define timeStopS(marker,arg) gettimeofday(&t ## marker ## 2,0); \
		{ \
			printf("benchmark-%s\t%f\n", arg \
				, ((float)(t ## marker ## 2 . tv_sec - t ## marker ## 1 . tv_sec) \
				+ ((float)(t ## marker ## 2 . tv_usec - t ## marker ## 1 . tv_usec))*1e-6) \
				); \
		}
#else
#define timeInit(marker)
#define timeStart(marker)
#define timeStop(marker)
#define timeStopS(marker,str)
#endif
