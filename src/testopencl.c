#include <stdio.h>
#include <string.h>
#include <stdbool.h>

#include "build_psdata.h"
#include "config.h"
#include "note.h"
#include "particle_system.h"
#include "opencl/particle_system_host.h"
#include "3rdparty/whereami.h"
#include <CL/opencl.h>  ///added for kdevelop  //cl_ext.h
#define PI 3.1415926535

void printFieldOffsets(psdata data) {
    for (size_t i = 0; i < data.num_fields; ++i) {
        note(2, "%s: %u\n", data.names + data.names_offsets[i], data.data_offsets[i]);
    }
}

int main(int argc, char *argv[])
{
    bool verbose = false;
    if (argc >= 2){
      if (strcmp(argv[1], "-v")) verbose = true;
    } 
    verbose = true; /*  force debug output */                                                         if (verbose) printf("\nchk0 ");
    set_log_level(1);
    psdata data;                                                                                    if (verbose) printf("\nchk1 ");    
    load_config(CMAKE_SOURCE_DIR "/conf/solid.conf");                                               if (verbose) printf("\nchk2 ");  
    build_psdata_from_string(&data, get_config_section("psdata_specification"));                    if (verbose) printf("\nchk3 ");  
    REAL * position;
    REAL * originalpos;
    REAL * density0;
    REAL * rotation;
    unsigned int numSteps = 200; //2000
    PS_GET_FIELD(data, "position", REAL, &position);
    PS_GET_FIELD(data, "originalpos", REAL, &originalpos);
    PS_GET_FIELD(data, "density0", REAL, &density0);
    PS_GET_FIELD(data, "rotation", REAL, &rotation);                                                if (verbose) printf("\nchk4 ");

    init_opencl();                                                                                  if (verbose) printf("\nchk5 ");  
    psdata_opencl pso = create_psdata_opencl(&data, get_config_section("opencl_kernel_files"));     if (verbose) printf("chk5.1 ");     // NB psdata data , psdata_opencl pso (struct holding cl_mem buffers)    
    populate_position_cuboid_device_opencl(pso, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 6, 6, 6);          if (verbose) printf("chk5.2 ");     // first kernel call  "populate_position_cuboid"
                                                                                                    //work group size=0 (auto) but 512 for others
    call_for_all_particles_device_opencl(pso, "init_original_position");                            if (verbose) printf("chk5.3 ");     
    rotate_particles_device_opencl(pso, PI/4, 0, PI/6);                                             if (verbose) printf("\nchk6 ");  
    printf("\nNumber of time steps = %i\n",numSteps);
    
    for(int i = 0; i<numSteps;i++){
        compute_particle_bins_device_opencl(pso);                                                   if (verbose) printf("\nchk7 ");  
        cl_mem temp = pso.data;
        pso.data = pso.tempdata ;
        pso.tempdata = temp;                                                                        // pointer swap:  &pso.data <=> &pso.tempdata
        
        call_for_all_particles_device_opencl(pso, "compute_original_density");
        call_for_all_particles_device_opencl(pso, "compute_density");                               if (verbose) printf("\nchk8 "); 
        call_for_all_particles_device_opencl(pso, "compute_rotations_and_strains");                 if (verbose) printf("\nchk9 "); 
        call_for_all_particles_device_opencl(pso, "compute_stresses");
        call_for_all_particles_device_opencl(pso, "compute_forces_solids");                         if (verbose) printf("\nchk10 "); 
        call_for_all_particles_device_opencl(pso, "step_forward");
        
        sync_psdata_device_to_host(data, pso);
	    //write_psdata(data, i, "solid");
        if (i%5 == 0) write_psdata_ply(data, i, "solid");                                           // save every 5th frame.
    }
    free_psdata_opencl(&pso);                                                                       if (verbose) printf("\nchk13 "); 
    terminate_opencl();                                                                             if (verbose) printf("\nchk14 "); 
    unload_config();
    //display_psdata(data, NULL);
    free_psdata(&data);
    return 0;
}

///  Note sequence for fluid simulation is (NB req fluid.conf &=> fluid.cl:
//     compute_particle_bins_device_opencl(*pso);
//     call_for_all_particles_device_opencl(*pso, "compute_density");
//     call_for_all_particles_device_opencl(*pso, "compute_forces_fluids");
//     call_for_all_particles_device_opencl(*pso, "step_forward");

