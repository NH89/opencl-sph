#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS   // needed to use clCreateCommandQueue on OpenCL-2.0
#include <CL/opencl.h>

#include "platforminfo.h"
#include "../note.h"

static Platform const * _platforms = NULL; /* Const as fields should not change value after init */
static unsigned int     _nplatforms = 0;

void   get_opencl_platform_info(Platform const ** const pptr, unsigned int * npptr);
void   free_opencl_platform_info();

char * getDeviceInfoString(cl_device_id, cl_device_info);
void   printDeviceType(cl_device_type);
void   printDeviceBool(cl_bool);

void get_opencl_platform_info(Platform const ** const pptr, unsigned int * const npptr) {
    if (pptr == NULL) return;

    if (_platforms != NULL) {
        *pptr = _platforms;
        *npptr = _nplatforms;
    }

    cl_platform_id * platform_ids;
    cl_uint          num_platforms;

    clGetPlatformIDs(0, NULL, &num_platforms);
    platform_ids = malloc(num_platforms*sizeof(cl_platform_id));
    clGetPlatformIDs(num_platforms, platform_ids, NULL);

    _nplatforms = (unsigned int) num_platforms;
    if (npptr != NULL) *npptr = _nplatforms;
    
    note(1, "%d platform%s present.\n\n", num_platforms, num_platforms == 1 ? "" : "s");

    _platforms = malloc(num_platforms*sizeof(Platform));
    *pptr = _platforms;                            /* Point to our static global */
    Platform * platforms = (Platform*) _platforms; /* Open for editing */

    int p;

    /* Fill platforms */
    for (p = 0; p < num_platforms; ++p) {
        size_t         platform_name_size;
        char *         platform_name;
        cl_device_id * device_ids;
        cl_uint        num_devices;

        platforms[p].id = platform_ids[p];

        clGetPlatformInfo(platforms[p].id, CL_PLATFORM_NAME, 0, NULL, &platform_name_size);
        platform_name = (char*) malloc(platform_name_size);
        clGetPlatformInfo(platforms[p].id, CL_PLATFORM_NAME, platform_name_size, platform_name, NULL);
        note(1, "%s\n", platform_name);

        platforms[p].name = platform_name;

        clGetDeviceIDs(platforms[p].id, 
                       CL_DEVICE_TYPE_ALL,
                       0,
                       NULL,
                       &num_devices);
        device_ids = (cl_device_id*) malloc(num_devices*sizeof(cl_device_id));
        clGetDeviceIDs(platforms[p].id,
                       CL_DEVICE_TYPE_ALL,
                       num_devices,
                       device_ids,
                       NULL);

        note(1, "%d devices on this platform:\n", num_devices);

        platforms[p].devices = malloc(num_devices*sizeof(Device));
        platforms[p].num_devices = num_devices;

        int d;

        /* Fill devices for this platform */
        for (d = 0; d < num_devices; ++d) {
            Device * device = (Device*) &platforms[p].devices[d];

            device->id = device_ids[d];

            device->name = getDeviceInfoString(device->id, CL_DEVICE_NAME);
            note(1, "Name: %s\n", device->name);
            clGetDeviceInfo(device->id, CL_DEVICE_TYPE, sizeof(cl_device_type), &device->type, NULL);
            note(1, " > Type: "); printDeviceType(device->type); note(1, "\n");
            device->opencl_version = getDeviceInfoString(device->id, CL_DEVICE_VERSION);
            note(1, " > Version: %s\n", device->opencl_version);
            device->driver_version = getDeviceInfoString(device->id, CL_DRIVER_VERSION);
            note(1, " > Driver version: %s\n", device->driver_version);

            device->extensions = getDeviceInfoString(device->id, CL_DEVICE_EXTENSIONS);
            note(1, " > Extensions: %s\n", device->extensions);
            if (strstr(device->extensions, "cl_khr_fp64") != NULL) {
                note(1, "   ├ cl_khr_fp64 found -> double precision enabled\n");
                device->double_supported = 1;
            } else {
                note(1, "   ├ cl_khr_fp64 not found -> double precision disabled\n");
                device->double_supported = 0;
            }
            note(1, "\n");

            clGetDeviceInfo
                (device->id, CL_DEVICE_COMPILER_AVAILABLE, sizeof(cl_bool),
                 &device->compiler_available, NULL);
            note(1, " > Compiler available: "); printDeviceBool(device->compiler_available);
            note(1, "\n");
            clGetDeviceInfo
                (device->id, CL_DEVICE_LINKER_AVAILABLE, sizeof(cl_bool),
                 &device->linker_available, NULL);
            note(1, " > Linker available: "); printDeviceBool(device->linker_available);
            note(1, "\n");
            device->profile = getDeviceInfoString(device->id, CL_DEVICE_PROFILE);
            note(1, " > Profile: %s\n", device->profile);
            note(1, "\n");

            clGetDeviceInfo
                (device->id, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint),
                 &device->max_compute_units, NULL);
            note(1, " > Number of compute units: %u\n", device->max_compute_units);
            clGetDeviceInfo
                (device->id, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t),
                 &device->max_workgroup_size, NULL);
            note(1, " > Max work-group size: %zu\n", device->max_workgroup_size);

            clGetDeviceInfo
                (device->id, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(cl_uint),
                 &device->max_work_item_dimensions, NULL);
            note(1, " > Max work item dimensions: %u\n", device->max_work_item_dimensions);
            device->max_work_item_sizes = malloc(device->max_work_item_dimensions*sizeof(size_t));
            clGetDeviceInfo
                (device->id, CL_DEVICE_MAX_WORK_ITEM_SIZES,
                 device->max_work_item_dimensions*sizeof(size_t),
                 device->max_work_item_sizes, NULL);
            note(1, " > Max work item sizes:");
            int i;
            for (i = 0; i < device->max_work_item_dimensions-1; ++i) {
                note(1, " %zu,", device->max_work_item_sizes[i]);
            }
            note(1, " %zu\n", device->max_work_item_sizes[i]);

            if (device->double_supported)
                clGetDeviceInfo
                    (device->id, CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE, sizeof(cl_uint),
                     &device->pref_vector_width, NULL);
            else
                clGetDeviceInfo
                    (device->id, CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT, sizeof(cl_uint),
                     &device->pref_vector_width, NULL);
            note(1, " > Preferred vector width (floating point): %u\n", device->pref_vector_width);

            note(1, "\n");
        }

        free(device_ids);
    }

    free(platform_ids);
}

void free_opencl_platform_info() {
    if (_platforms == NULL) return;

    int p;

    for (p = 0; p < _nplatforms; ++p) {
        free(_platforms[p].name);

        int d;

        for (d = 0; d < _platforms[p].num_devices; ++d) {
            Device * device = (Device*) &_platforms[p].devices[d];

            free(device->name);
            free(device->opencl_version);
            free(device->driver_version);
            free(device->extensions);
            free(device->profile);
            free(device->max_work_item_sizes);
        }

        free((void*)_platforms[p].devices);
    }

    free((void*)_platforms);
}

/* Heap allocates the string */
char * getDeviceInfoString(cl_device_id device, cl_device_info param_name) {
    size_t device_info_size;
    char * device_info;

    clGetDeviceInfo(device, param_name, 0, NULL, &device_info_size);
    device_info = (char*) malloc(device_info_size);
    clGetDeviceInfo(device, param_name, device_info_size, device_info, NULL);

    return device_info;
}

void printDeviceBool(cl_bool b) {
    if (b == CL_TRUE) note(1, "yes");
    else note(1, "no");
}

void printDeviceType(cl_device_type type) {
    switch (type) {
        case CL_DEVICE_TYPE_CPU:
            note(1, "CPU"); break;
        case CL_DEVICE_TYPE_GPU:
            note(1, "GPU"); break;
        case CL_DEVICE_TYPE_ACCELERATOR:
            note(1, "Accelerator"); break;
        case CL_DEVICE_TYPE_DEFAULT:
            note(1, "Default"); break;
        default:
            break;
    }
}
