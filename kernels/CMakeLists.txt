file(GLOB kernels_FILES *.cl)

install(FILES ${kernels_FILES} DESTINATION share/opencl_sph/kernels/)
