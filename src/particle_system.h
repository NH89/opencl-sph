#ifndef PARTICLE_SYSTEM_H_
#define PARTICLE_SYSTEM_H_

#include <assert.h>//#include "macros.h"
#include "types.h"

/* Unnecessary but convenient voodoo, type should match pointer type */
#define PS_P_PTR(d, position, type) ((type*) ((char*) (d).data + (d).data_offsets[position]))
#define PS_GET_FIELD(data, name, type, pointer_to_pointer) {\
    int __fldpos = get_field_psdata(data, name); \
    assert(__fldpos != -1); \
    *pointer_to_pointer = PS_P_PTR(data, __fldpos, type);}

typedef struct {
    /**
     * The numerical data is all stored in a contiguous array at data, to
     * facilitate uploading to the GPU.
     */
    unsigned int num_fields;
    const char * names;
    unsigned int * names_offsets;
    unsigned int * dimensions; /* e.g. { 30, 3 } */
    unsigned int * num_dimensions; /* 2 for example above */
    unsigned int * dimensions_offsets;
    unsigned int * entry_sizes;

    void * data;
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
void write_psdata_ply(psdata, int number, const char * Case);


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
