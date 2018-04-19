* Data structure uses integer offsets, may cause an overflow with large memory allocations
	* Suggest using **cl_long** or **cl_ulong** for offsets instead of unsigned int 
* JSON configuration files for runtime parameters
* Get cumulative times from OpenCL function calls **call_for_all_particles_device_opencl**
* Fix indentation issues in testopencl.c
* Get target platform and device id from command line or configuration file
