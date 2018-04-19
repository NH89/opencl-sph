* Data structure uses integer offsets, may cause an overflow with large memory allocations
	* Suggest using **cl_long** or **cl_ulong** for offsets instead of unsigned int 
* JSON configuration files for runtime parameters
* Get cumulative times from OpenCL function calls **call_for_all_particles_device_opencl**
* Fix indentation issues in testopencl.c
* Get target platform and device id from command line or configuration file
* Use HDF5 for input and output


##Directory Structure

I (Jeffrey Kelling) suggest implementing something like the following directory structure later:

/include/opencl-sph/
        for all the headers and code which are now compiled into the static library
        set ${CMAKE_SOURCE_SIR}/include as include path
        use relative paths for includes in files in there
        the headers in here can then also be part of the install target defined in
        cmake
/src/
        for the separate (main) programs
        if a main program shall contain code which is not shared (thus not in
        /inlude/....) put the other source files and headers belonging to it in a
        subdirectory (e.g.: /src/test/ for /src/test.cpp
        #include code from /include/opencl-sph/ using angle brackets:
        #include <spencl-sph/example.h>
/share/kernels/
        for the opencl kernel code
