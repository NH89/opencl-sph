#include "mex.h"

#include "../build_psdata.h"
#include "../config.h"
#include "../particle_system.h"
#include "../opencl/particle_system_host.h"

void onExit() {
    free_stored_psdata_opencl();
    terminate_opencl();

    free_stored_psdata();
}

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]) {
    
    if (nlhs != 0) mexErrMsgIdAndTxt("Init:ReturnError", "This function does not return anything.");

    psdata * data = get_stored_psdata();
    
    
    if (data != NULL) mexErrMsgIdAndTxt("Init:InitError", "Module already initialised");

    if (nrhs == 0 || !mxIsChar(prhs[0])) mexErrMsgIdAndTxt("Init:ArgError", "No config path given");

    size_t config_path_len = mxGetN(prhs[0]) + 1;
    char * config_path = mxMalloc(config_path_len);

    mxGetString(prhs[0], config_path, config_path_len);

    load_config(config_path);

    mxFree(config_path);

    
    create_stored_psdata_from_string(get_config_section("psdata_specification"));

    data = get_stored_psdata();

    init_opencl();

    psdata_opencl pso = create_psdata_opencl(data, get_config_section("opencl_kernel_files"));

    opencl_use_buflist(pso);

    unload_config();

    mexAtExit(onExit);
     
}
