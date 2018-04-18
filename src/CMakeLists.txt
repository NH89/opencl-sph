set(opencl_sph_SRCS
		build_psdata.c
		config.c
		note.c
		particle_system.c
		stringly.c
		opencl/clerror.c
		opencl/particle_system_host.c
		opencl/platforminfo.c
		)

add_library(opencl_sph_lib STATIC ${opencl_sph_SRCS})

add_executable(opencl_sph  testopencl.c)
target_link_libraries(opencl_sph opencl_sph_lib ${OpenCL_LIBRARY} m) # m is a math library

add_executable(testStretchSolids  testStretchSolid.c)
target_link_libraries(testStretchSolids opencl_sph_lib ${OpenCL_LIBRARY} m) # m is a math library

add_executable(testFluid  testFluid.c)
target_link_libraries(testFluid opencl_sph_lib ${OpenCL_LIBRARY} m) # m is a math library

install(TARGETS opencl_sph testStretchSolids testFluid DESTINATION bin/)
