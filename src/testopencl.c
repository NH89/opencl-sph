#include <stdio.h>
#include <string.h>
#include <stdbool.h>

#include "build_psdata.h"
#include "config.h"
#include "note.h"
#include "particle_system.h"
#include "opencl/particle_system_host.h"
#include "3rdparty/whereami.h"

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
    
    if (verbose) printf("chk0 ");
    set_log_level(1);
    psdata data;
    if (verbose) printf("chk1 ");    
    load_config(CMAKE_SOURCE_DIR "/conf/solid.conf");
    if (verbose) printf("chk2 ");  
    build_psdata_from_string(&data, get_config_section("psdata_specification"));
    if (verbose) printf("chk3 ");  
    double * position;
    double * originalpos;
    double * density0;

    double * rotation;
    uint numSteps = 2000;

    PS_GET_FIELD(data, "position", double, &position);
    PS_GET_FIELD(data, "originalpos", double, &originalpos);
    PS_GET_FIELD(data, "density0", double, &density0);

    PS_GET_FIELD(data, "rotation", double, &rotation);
    if (verbose) printf("chk4 ");  
    init_opencl();
    if (verbose) printf("chk5 ");  
        psdata_opencl pso = create_psdata_opencl(&data, get_config_section("opencl_kernel_files"));
    if (verbose) printf("chk5.1 ");         
            populate_position_cuboid_device_opencl(pso, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 6, 6, 6);
    if (verbose) printf("chk5.2 "); // for some reason work group size is zero after this, but 512 for others later.
            call_for_all_particles_device_opencl(pso, "init_original_position");
    if (verbose) printf("chk5.3 "); 
            rotate_particles_device_opencl(pso, PI/4, 0, PI/6);
    if (verbose) printf(" chk6 ");  
for(int i = 0; i<numSteps;i++)
   {
            compute_particle_bins_device_opencl(pso);
    if (verbose) printf(" chk7 ");  
            call_for_all_particles_device_opencl(pso, "compute_original_density");
            call_for_all_particles_device_opencl(pso, "compute_density");
    if (verbose) printf(" chk8 "); 
            call_for_all_particles_device_opencl(pso, "compute_rotations_and_strains");
    if (verbose) printf(" chk9 "); 
            call_for_all_particles_device_opencl(pso, "compute_stresses");
            call_for_all_particles_device_opencl(pso, "compute_forces_solids");
    if (verbose) printf(" chk10 "); 
            call_for_all_particles_device_opencl(pso, "step_forward");
        
            sync_psdata_device_to_host(data, pso);
	    write_psdata(data, i, "solid");
}
        free_psdata_opencl(&pso);
    if (verbose) printf(" chk13 "); 
    terminate_opencl();
    if (verbose) printf(" chk14 "); 
    unload_config();
    //display_psdata(data, NULL);

    free_psdata(&data);

    return 0;
}
