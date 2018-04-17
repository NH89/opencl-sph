#ifndef PARTICLE_SYSTEM_H_
#define PARTICLE_SYSTEM_H_

#include "macros.h"



#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif

/* Unnecessary but convenient voodoo, type should match pointer type */
#define PS_P_PTR(d, position, type) ((type*) ((char*) (d).data + (d).data_offsets[position]))
#define PS_GET_FIELD(data, name, type, pointer_to_pointer) {\
    int __fldpos = get_field_psdata(data, name); \
    ASSERT(__fldpos != -1); \
    *pointer_to_pointer = PS_P_PTR(data, __fldpos, type);}

typedef struct {
    /**
     * The numerical data is all stored in a contiguous array at data, to
     * facilitate uploading to the GPU.
     *
     * When data needs to be used by MATLAB, sync_to_mex should be used to
     * copy the flat stored data to the MEX managed arrays in host_data.
     * These are kept host-side throughout as I don't know anything about
     * their internal structure.
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

#ifdef MATLAB_MEX_FILE
psdata * create_stored_psdata_from_string(const char * string);
psdata * get_stored_psdata();
void free_stored_psdata();
#endif

void display_psdata(psdata, const char * const * mask);
void write_psdata(psdata, int number, const char * Case);

void init_psdata_fluid( psdata * data, int pnum, double mass, double timestep, double smoothingradius,
       double xbound1, double ybound1, double zbound1, double xbound2, double ybound2, double zbound2 );

int get_field_psdata( psdata data, const char * name );
void set_field_psdata( psdata * data, const char * name, void * field, unsigned int size, unsigned int offset );

unsigned int psdata_names_size( psdata data );
unsigned int psdata_dimensions_size( psdata data );
unsigned int psdata_data_size( psdata data );

int create_host_field_psdata( psdata * data, const char * name, void * field, unsigned int size );
int get_host_field_psdata( psdata * data, const char * name );

#ifdef MATLAB_MEX_FILE
void sync_to_mex( psdata * data );
int is_mex_field( const char * name );
void get_mex_field_name(const char * name, /* Output */ size_t * pReturnLength, char * mex_name);
#endif

void free_psdata( psdata * data );

#endif
