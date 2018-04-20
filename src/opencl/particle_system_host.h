#ifndef PARTICLE_SYSTEM_HOST_H_
#define PARTICLE_SYSTEM_HOST_H_

#define CL_USE_DEPRECATED_OPENCL_1_2_APIS   // needed to use clCreateCommandQueue on OpenCL-2.0
#include <CL/opencl.h>
#include "../particle_system.h"

typedef struct {
    unsigned int num_fields;
    unsigned int device_id;
    unsigned int device_ind;
    
    /* necessary stuff for allocate data location */ 
    cl_mem names;
    cl_mem names_offsets;
    cl_mem dimensions;
    cl_mem num_dimensions;
    cl_mem dimensions_offsets;
    cl_mem entry_sizes;
    cl_uint unit_length; // Length of just one item 
    
    /* This is every kind of data all together   */ 
    cl_mem data;
    cl_uint pnum;  // Current num of particles being tracked     
    cl_uint num_of_offset; // Ring buffer front start from 0
    cl_uint initial_offset; // Ring buffer front initialized at middle of buffer 

    /* Memory reserve for communication */
    cl_mem left_rec;
    cl_mem left_send;
    cl_mem right_rec;
    cl_mem right_send;
    
    /* Device only buffers
     *
     */
    cl_mem gridcount;
    cl_mem gridcell;
    cl_mem cellparticles;
    cl_mem celloffset;
    cl_mem gridres;
    cl_mem block_totals;
    cl_mem backup_prefix_sum;

     /* System dependent variables, also computed on buffer creation
     *
     * These are metadata for the various operations the device will perform
     */
    unsigned int num_grid_cells;
    size_t po2_workgroup_size;
    unsigned int num_blocks;
    
    /* Device specific kernels   */ 
    const cl_kernel * kernels;
    size_t num_kernels;
    const char * const * kernel_names;

  
    
    
} device_mem;


typedef struct {
    psdata host_psdata;

    unsigned int num_devices;
    unsigned int device_id_list[];
    
    // The memory structure of the devices
    device_mem device_mems[]; 
    

    /* One Program for complete context */ 
    cl_program ps_prog;
    const char * const * kernel_names;
    size_t num_kernels;
} psdata_opencl;

#ifdef MATLAB_MEX_FILE
psdata_opencl * get_stored_psdata_opencl();
void free_stored_psdata_opencl();
void opencl_use_buflist(psdata_opencl pso);
#endif

void init_opencl();
psdata_opencl create_psdata_opencl(psdata * data, const char * file_list,const uint & num_dev);
device_mem create_device_buffer(psdata *data, const uint device_id);
void assign_device_kernel_args(device_mem device);
void set_kernel_args_to_device(device_mem device, cl_kernel kernel);
void free_psdata_opencl(psdata_opencl * pso);
void free_psdata_device(device_mem * device);
void terminate_opencl();

// The only sync func used in c TODO 
void sync_psdata_device_to_host(psdata data, psdata_opencl pso);

// mex only 
void sync_psdata_fields_device_to_host(psdata data, psdata_opencl pso, size_t num_fields, const char * const * const field_names);
void sync_psdata_host_to_device(psdata data, psdata_opencl pso, int full);
void sync_psdata_fields_host_to_device(psdata data, psdata_opencl pso, size_t num_fields, const char * const * const field_names);

// Init TODO
void populate_position_cuboid_device_opencl(psdata_opencl pso,
                                            REAL x1, REAL y1, REAL z1,
                                            REAL x2, REAL y2, REAL z2,
                                            unsigned int xsize,
                                            unsigned int ysize,
                                            unsigned int zsize);


// Rot TODO
void rotate_particles_device_opencl(psdata_opencl, REAL angle_x, REAL angle_y, REAL angle_z);

// Call TODO 
void call_for_all_particles_device_opencl(psdata_opencl, const char * kernel_name);

// Include Comm TODO 
void compute_particle_bins_device_opencl(psdata_opencl);

//TODO
void compute_density_device_opencl(psdata_opencl);

//TODO 
void compute_forces_device_opencl(psdata_opencl);

//TODO 
void step_forward_device_opencl(psdata_opencl);

#endif
