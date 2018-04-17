#ifndef PLATFORM_INFO_H_
#define PLATFORM_INFO_H_

#define CL_USE_DEPRECATED_OPENCL_1_2_APIS   // needed to use clCreateCommandQueue on OpenCL-2.0
#include <CL/opencl.h>

typedef struct {
    cl_device_id   id;
    char *         name;
    cl_device_type type;
    char *         opencl_version;
    char *         driver_version;
    char *         extensions;
    cl_bool        compiler_available;
    cl_bool        linker_available;
    char *         profile;
    cl_uint        max_compute_units;
    cl_uint        max_workgroup_size;
    cl_uint        max_work_item_dimensions;
    size_t *       max_work_item_sizes;
    cl_uint        pref_vector_width;
    cl_bool        double_supported;
} Device;

typedef struct {
    cl_platform_id  id;
    char *          name;
    Device const *  devices;
    unsigned int    num_devices;
} Platform;

void get_opencl_platform_info(Platform const ** const pptr, unsigned int * const npptr);
void free_opencl_platform_info();

#endif
