#ifndef PARTICLE_SYSTEM_HOST_H_
#define PARTICLE_SYSTEM_HOST_H_

#define CL_USE_DEPRECATED_OPENCL_1_2_APIS   // needed to use clCreateCommandQueue on OpenCL-2.0
#include <CL/opencl.h>
#include "../particle_system.h"

typedef struct {
    psdata host_psdata;

    unsigned int num_fields;
    cl_mem names;               //Memory objects are described by a cl_mem object. Kernels take memory objects as input, and output to one or more memory objects.
    cl_mem names_offsets;
    cl_mem dimensions;
    cl_mem num_dimensions;
    cl_mem dimensions_offsets;
    cl_mem entry_sizes;

    cl_mem data;
    cl_mem data_sizes;
    cl_mem data_offsets;

    /* Device only buffers
     *
     * The gridcount stuff should probably go here too but I'll leave
     * it for now because lazy
     */

    cl_mem block_totals;
    cl_mem backup_prefix_sum;

    /* Program stuff, it's per psdata as the program changes based on data array layout */

    size_t num_kernels;
    const char * const * kernel_names;
    const cl_kernel * kernels;

    cl_program ps_prog;

    /* System dependent variables, also computed on buffer creation
     *
     * These are metadata for the various operations the device will perform
     */
    unsigned int num_grid_cells;
    size_t po2_workgroup_size;
    unsigned int num_blocks;
} psdata_opencl;

#ifdef MATLAB_MEX_FILE
psdata_opencl * get_stored_psdata_opencl();
void free_stored_psdata_opencl();
void opencl_use_buflist(psdata_opencl pso);
#endif

void init_opencl();
psdata_opencl create_psdata_opencl(psdata * data, const char * file_list);
void assign_pso_kernel_args(psdata_opencl pso);
void set_kernel_args_to_pso(psdata_opencl pso, cl_kernel kernel);
void free_psdata_opencl(psdata_opencl * pso);
void terminate_opencl();

void sync_psdata_device_to_host(psdata data, psdata_opencl pso);
void sync_psdata_fields_device_to_host(psdata data, psdata_opencl pso, size_t num_fields, const char * const * const field_names);
void sync_psdata_host_to_device(psdata data, psdata_opencl pso, int full);
void sync_psdata_fields_host_to_device(psdata data, psdata_opencl pso, size_t num_fields, const char * const * const field_names);

void populate_position_cuboid_device_opencl(psdata_opencl pso,
                                            REAL x1, REAL y1, REAL z1,
                                            REAL x2, REAL y2, REAL z2,
                                            unsigned int xsize,
                                            unsigned int ysize,
                                            unsigned int zsize);
void rotate_particles_device_opencl(psdata_opencl, REAL angle_x, REAL angle_y, REAL angle_z);
void call_for_all_particles_device_opencl(psdata_opencl, const char * kernel_name);

void compute_particle_bins_device_opencl(psdata_opencl);

void compute_density_device_opencl(psdata_opencl);
void compute_forces_device_opencl(psdata_opencl);
void step_forward_device_opencl(psdata_opencl);

#endif
