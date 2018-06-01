#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <boost/program_options.hpp>                                                        // NB writen for Boost 1.67.0 , 1.58 will not work.
#include <string>
#include <vector>
#include <algorithm>
#include <iterator>
#include <iostream>

#include "build_psdata.h"
#include "config.h"
#include "note.h"
#include "particle_system.h"
#include "opencl/particle_system_host.h"

int main(int argc, char *argv[])
{
    namespace po = boost::program_options;                                                  // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("input-file", boost::program_options::value< std::vector<std::string> >(), "input file")
    ;
    int opt;
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
//        ("optimization", po::value<int>(&opt)->default_value(10), "optimization level")   // for device program compilation // move to json file ?
//        ("include-path,I", po::value< std::vector<std::string> >(), "include path")       // for device program compilation
        ("input-file", po::value< std::vector<std::string> >(), "input file")               // json config file, may reference an Hdf5 file of particles
    ;                                                                                       // json file to hold an instance of (the) host 'simulation' class.
    set_log_level(1);                                                                       // controls the behaviour of note(...)
    psdata data;
    load_config(CMAKE_SOURCE_DIR "/conf/solid.conf");                                       // replace with importing json file
    build_psdata_from_string(&data, get_config_section("psdata_specification"));            // import particle system data
    REAL * position;
    REAL * originalpos;
    REAL * density0;
    REAL * rotation;                                                                        // pointer to REAL data type - by default REAL=foat
    uint numSteps = 2000;
    PS_GET_FIELD(data, "position", REAL, &position);                                        // find start of each sub-array within 'data' array
    PS_GET_FIELD(data, "originalpos", REAL, &originalpos);
    PS_GET_FIELD(data, "density0", REAL, &density0);
    PS_GET_FIELD(data, "rotation", REAL, &rotation);                                        // end of data import
    init_opencl();
    psdata_opencl pso = create_psdata_opencl(&data, get_config_section("opencl_kernel_files"));       // creates cl_buffer objects
    populate_position_cuboid_device_opencl(pso, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 6, 6, 6);  //# generate particles, uses kernel "populate_position_cuboid"
                                                                                            //# replace this with external Hdf5 particle model file, and clEnqueueWriteBuffer(...)
    call_for_all_particles_device_opencl(pso, "init_original_position");                    //#
    rotate_particles_device_opencl(pso, M_PI_4, 0, M_PI/6);                                 //# end of model initialization
    
    for(int i = 0; i<numSteps;i++){                                                         // simulation loop
        compute_particle_bins_device_opencl(pso);
        call_for_all_particles_device_opencl(pso, "compute_original_density");
        call_for_all_particles_device_opencl(pso, "compute_density");
        call_for_all_particles_device_opencl(pso, "compute_rotations_and_strains");
        call_for_all_particles_device_opencl(pso, "compute_stresses");
        call_for_all_particles_device_opencl(pso, "compute_forces_solids");
        call_for_all_particles_device_opencl(pso, "step_forward");
        sync_psdata_device_to_host(data, pso);                                              // collect data from device
	    write_psdata(data, i, "solid");                                                     // write data to file
    }
    free_psdata_opencl(&pso);                                                               // free memory & tidy up
    terminate_opencl();
    unload_config();
    free_psdata(&data);
    return 0;
}
