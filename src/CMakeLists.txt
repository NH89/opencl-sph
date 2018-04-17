file(GLOB opencl_SRCS opencl/*.c)

add_library(opencl_sph_lib STATIC ${opencl_SRCS})

set(opencl_sph_SRCS
	build_psdata.c
	config.c
	note.c
	particle_system.c
	stringly.c
	testopencl.c
)

add_executable(opencl_sph  ${opencl_sph_SRCS})
target_link_libraries(opencl_sph opencl_sph_lib ${OpenCL_LIBRARY} m) # m is a math library
install(TARGETS opencl_sph DESTINATION bin/)
