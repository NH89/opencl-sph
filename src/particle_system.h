#ifndef PARTICLE_SYSTEM_H_
#define PARTICLE_SYSTEM_H_

#include "macros.h"
#include "types.h"

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif

/* Unnecessary but convenient voodoo, type should match pointer type */
#define PS_P_PTR(d, position, type)  ((type*)  ((char*) (d).data + (d).data_offsets[position] )) //  used line 17 below

#define PS_GET_FIELD(data_, name, type, pointer_to_pointer) {                \
                    int __fldpos = get_field_psdata(data, name);            \
                    ASSERT(__fldpos != -1);                                 \
                    *pointer_to_pointer = PS_P_PTR(data, __fldpos, type);}                       // used in  'opencl/particle_system_host.c'  and  'testsolid.cpp'
                    //PS_P_PTR(d, position, type)
                    //*pointer_to_pointer = ((type*)  ( (char*) (data_).data + (data_).data_offsets[__fldpos] ) )

typedef struct {
    unsigned int num_fields;
    const char * names;
    unsigned int * names_offsets;
    unsigned int * dimensions;                  /* e.g. { 30, 3 } */
    unsigned int * num_dimensions;              /* 2 for example above */
    unsigned int * dimensions_offsets;
    unsigned int * entry_sizes;

    void * data;                                // Numerical data in contiguous array, facilitate uploading to GPU.
    unsigned int * data_sizes;
    unsigned int * data_offsets;

    /* Fields not needed in computation - won't be sent to GPU */
    unsigned int num_host_fields;
    char ** host_names;
    void ** host_data;
    unsigned int * host_data_size;
} psdata;

void display_psdata(psdata, const char * const * mask);
void write_psdata(psdata, int number, const char * Case);

void init_psdata_fluid( psdata * data, int pnum, REAL mass, REAL timestep, REAL smoothingradius,
       REAL xbound1, REAL ybound1, REAL zbound1, REAL xbound2, REAL ybound2, REAL zbound2 );

int get_field_psdata( psdata data, const char * name );
void set_field_psdata( psdata * data, const char * name, void * field, unsigned int size, unsigned int offset );

unsigned int psdata_names_size( psdata data );
unsigned int psdata_dimensions_size( psdata data );
unsigned int psdata_data_size( psdata data );

int create_host_field_psdata( psdata * data, const char * name, void * field, unsigned int size );
int get_host_field_psdata( psdata * data, const char * name );

void free_psdata( psdata * data );

#endif
