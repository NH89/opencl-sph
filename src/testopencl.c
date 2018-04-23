#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include<time.h>

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
    

    set_log_level(1);
    psdata data;

    load_config(CMAKE_SOURCE_DIR "/conf/solid.conf");

    build_psdata_from_string(&data, get_config_section("psdata_specification"));


    uint numSteps = 3000;

    init_opencl();

        psdata_opencl pso = create_psdata_opencl(&data, get_config_section("opencl_kernel_files"));
            populate_position_cuboid_device_opencl(pso, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 10, 10, 10);
            call_for_all_particles_device_opencl(pso, "init_original_position");
            rotate_particles_device_opencl(pso, PI/4, 0, PI/6);

    clock_t t1=clock();

    for(int i = 0; i<numSteps;i++)
    {
            compute_particle_bins_device_opencl(pso);

            call_for_all_particles_device_opencl(pso, "compute_original_density");
            call_for_all_particles_device_opencl(pso, "compute_density");
            call_for_all_particles_device_opencl(pso, "compute_rotations_and_strains");
            call_for_all_particles_device_opencl(pso, "compute_stresses");
            call_for_all_particles_device_opencl(pso, "compute_forces_solids");
            call_for_all_particles_device_opencl(pso, "step_forward");
        
            sync_psdata_device_to_host(data, pso);
            write_psdata(data, i, "solid");
    }
    clock_t t2=clock();

    printf("The time taken is.. %g ", (t2-t1));

    free_psdata_opencl(&pso);
    terminate_opencl();
    unload_config();
    free_psdata(&data);

    return 0;
}