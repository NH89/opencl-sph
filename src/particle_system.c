#include <math.h>
#include <string.h>
#include <stdio.h>

#ifndef MATLAB_MEX_FILE
    #include <stdlib.h>
#else
    #include "mex.h"
#endif


#include "note.h"
#include "particle_system.h"
#include "build_psdata.h"
#include "opencl/particle_system_host.h"


static psdata * ps_instance = NULL;

#ifdef MATLAB_MEX_FILE
psdata * create_stored_psdata_from_string(const char * string) {
    ps_instance = malloc(sizeof(psdata));

    ps_instance->num_fields = 0;

    build_psdata_from_string(ps_instance, string);

    return ps_instance;
}

psdata * get_stored_psdata() {
    return ps_instance;
}

void free_stored_psdata() {
    free_psdata(ps_instance);

    free(ps_instance);
}
#endif

void display_entry(psdata data, size_t offset, size_t size) {
    if (size == 4) {
        note(2, "%f, ", *((REAL*)((char*)data.data + offset)));
        note(2, "%u, ", *((unsigned int*)((char*)data.data + offset)));
    }
    /*else if (size == 4) {
        note(2, "%u, ", *((unsigned int*)((char*)data.data + offset)));
    }
    */
}

void display_psdata(psdata data, const char * const * mask) {
    if (mask == NULL) {
        for (size_t field = 0; field < data.num_fields; ++field) {
            note(2, "%s:\n\n", data.names + data.names_offsets[field]);
            if (data.num_dimensions[field] == 1) {
                unsigned int d0 = (data.dimensions + data.dimensions_offsets[field])[0];
                for (size_t i = 0; i < d0; ++i) {
                    display_entry(data, data.data_offsets[field] + i*data.entry_sizes[field], data.entry_sizes[field]);
                }
            } else if (data.num_dimensions[field] == 2) {
                unsigned int d0 = (data.dimensions + data.dimensions_offsets[field])[0];
                unsigned int d1 = (data.dimensions + data.dimensions_offsets[field])[1];
                for (unsigned int i = 0; i < d1; ++i) {
                    for (unsigned int j = 0; j < d0; ++j) {
                        display_entry(data, data.data_offsets[field] + (i*d0+j)*data.entry_sizes[field], data.entry_sizes[field]);
                    }
                    note(2, "\n");
                }
            } else if (data.num_dimensions[field] == 3) {
                unsigned int d0 = (data.dimensions + data.dimensions_offsets[field])[0];
                unsigned int d1 = (data.dimensions + data.dimensions_offsets[field])[1];
                unsigned int d2 = (data.dimensions + data.dimensions_offsets[field])[2];
                for (unsigned int i = 0; i < d2; ++i) {
                    for (unsigned int j = 0; j < d1; ++j) {
                        for (unsigned int k = 0; k < d0; ++k) {
                            display_entry(data, data.data_offsets[field] + (i*d1*d0+j*d0+k)*data.entry_sizes[field], data.entry_sizes[field]);
                        }
                        note(2, "\n");
                    }
                    note(2, "\n");
                }
            }
            note(2, "\n\n");
        }
    }
}

void set_constraints(psdata data)
{

    //Constraints are set the similar way which set in the go_test.m file
    for (size_t field = 0; field < data.num_fields; ++field) {

        char *name = data.names + data.names_offsets[field];
        double* position;

        if (strcmp(name, "plane_constraints_particles") == 0) {
            unsigned int d0 = (data.dimensions + data.dimensions_offsets[field])[0];
            for (size_t i = 0; i < d0; ++i)
            {
                unsigned int value = *(unsigned int*)((char*)data.data + data.data_offsets[field] + i*data.entry_sizes[field]);

                PS_GET_FIELD(data, "position", double, &position);

                int position_field = get_field_psdata(data, "position");

                REAL positiony  = *(REAL *)((char*)data.data + data.data_offsets[position_field] + (3*i + 1)*data.entry_sizes[position_field]);

                if(positiony == -1)
                {
                    *(unsigned int *) ((char *) data.data + data.data_offsets[field] + i * data.entry_sizes[field]) = 1;
                }
                else if (positiony == 1)
                {
                    *(unsigned int *) ((char *) data.data + data.data_offsets[field] + i * data.entry_sizes[field]) = 2;
                }
            }
        }
    }
}

void write_psdata(psdata data, int number, const char* Case)
{

    for (size_t field = 0; field < data.num_fields; ++field)
    {

        char * name = data.names + data.names_offsets[field];


        char snum[5];
        sprintf(snum, "%d", number);


        if (strcmp(name, "position") == 0)
        {

            char filename[150];
            strcpy(filename,CMAKE_SOURCE_DIR);
            strcat(filename, "/output/");
            strcat(filename, Case);
            strcat(filename, "/position_");
            strcat(filename, snum);
            strcat(filename, ".csv");

            FILE *f = fopen(filename, "w");

            unsigned int d0 = (data.dimensions + data.dimensions_offsets[field])[0];
            unsigned int d1 = (data.dimensions + data.dimensions_offsets[field])[1];

            for (unsigned int i = 0; i < d1; ++i) {
                for (unsigned int j = 0; j < d0; ++j) {
                    fprintf(f, "%0.3lf,", *((REAL*)((char*)data.data + data.data_offsets[field] + (i*d0+j)*data.entry_sizes[field])));
                }
                fprintf(f,"\n");

            }
            fclose(f);
        }
    }

}

void init_psdata_fluid( psdata * data, int pnum, REAL mass, REAL timestep, REAL smoothingradius,
                        REAL xbound1, REAL ybound1, REAL zbound1,
                        REAL xbound2, REAL ybound2, REAL zbound2 )
{
    { /* Names */
    const char * names_ref[] = {
        "pnum",
        "n",
        "mass",
        "timestep",
        "smoothingradius",

        "position",
        "posnext",
        "velocity",
        "veleval",
        "velnext",
        "acceleration",
        "force",
        "density",
        "volume",

        "gridbounds",
        "gridres",
        "gridcell",
        "gridcount",
        "celloffset",
        "cellparticles"
    };

    unsigned int i;
    unsigned int sum = 0;

    data->num_fields = sizeof(names_ref)/sizeof(char*);
    data->names_offsets = malloc(data->num_fields*sizeof(unsigned int));

    for (i = 0; i < data->num_fields; ++i) {
        data->names_offsets[i] = sum;
        sum += strlen(names_ref[i]) + 1;
    }

    data->names = malloc(sum*sizeof(char));

    char * name_ptr;
    for (i = 0; i < data->num_fields; ++i) {
        name_ptr = (char*) data->names + data->names_offsets[i];
        strcpy(name_ptr, names_ref[i]);
    }
    }

    unsigned int num_gridcells_x = (unsigned int) (fabs(xbound2 - xbound1) / smoothingradius);
    unsigned int num_gridcells_y = (unsigned int) (fabs(ybound2 - ybound1) / smoothingradius);
    unsigned int num_gridcells_z = (unsigned int) (fabs(zbound2 - zbound1) / smoothingradius);

    { /* Dimensions */
    unsigned int dimensions[] = {
        1,
        1,
        1,
        1,
        1,

        pnum, 3,
        pnum, 3,
        pnum, 3,
        pnum, 3,
        pnum, 3,
        pnum, 3,
        pnum, 3,
        pnum,
        pnum,

        2, 3,
        3,
        pnum,
        num_gridcells_x, num_gridcells_y, num_gridcells_z,
        num_gridcells_x, num_gridcells_y, num_gridcells_z,
        pnum
    };

    unsigned int num_dimensions[] = { 1, 1, 1, 1, 1,  2, 2, 2, 2, 2, 2, 2, 1, 1,  2, 1, 1, 3, 3, 1 };

    unsigned int entry_sizes[] = {
        sizeof(int),
        sizeof(int),
        sizeof(REAL),
        sizeof(REAL),
        sizeof(REAL),

        sizeof(REAL),
        sizeof(REAL),
        sizeof(REAL),
        sizeof(REAL),
        sizeof(REAL),
        sizeof(REAL),
        sizeof(REAL),
        sizeof(REAL),
        sizeof(REAL),

        sizeof(REAL),
        sizeof(unsigned int),
        sizeof(unsigned int),
        sizeof(unsigned int),
        sizeof(unsigned int),
        sizeof(unsigned int)
    };

    data->num_dimensions     = malloc(data->num_fields*sizeof(unsigned int));
    data->dimensions_offsets = malloc(data->num_fields*sizeof(unsigned int));
    data->entry_sizes        = malloc(data->num_fields*sizeof(unsigned int));

    memcpy(data->num_dimensions, num_dimensions, data->num_fields*sizeof(unsigned int));
    memcpy(data->entry_sizes,    entry_sizes,    data->num_fields*sizeof(unsigned int));

    unsigned int i;
    unsigned int sum = 0;
    for (i = 0; i < data->num_fields; ++i) {
        data->dimensions_offsets[i] = sum;
        sum += num_dimensions[i];
    }

    data->dimensions = malloc(sum*sizeof(unsigned int));
    
    memcpy(data->dimensions, dimensions, sum*sizeof(unsigned int));
    }

    { /* Data */
    unsigned int data_size = 0, field_data_size;

    data->data_sizes = malloc(data->num_fields*sizeof(unsigned int));
    data->data_offsets = malloc(data->num_fields*sizeof(unsigned int));

    unsigned int i, j;
    for (i = 0; i < data->num_fields; ++i) {
        field_data_size = data->entry_sizes[i];

        for (j = 0; j < data->num_dimensions[i]; ++j) {
            field_data_size *= data->dimensions[data->dimensions_offsets[i]+j];
        }

        data->data_sizes[i] = field_data_size;
        data->data_offsets[i] = data_size;

        data_size += field_data_size;
    }

    data->data = calloc(data_size, 1);
    }

    { /* Assign some values */
    set_field_psdata(data, "pnum", &pnum, sizeof(int), 0);
    set_field_psdata(data, "mass", &mass, sizeof(REAL), 0);
    set_field_psdata(data, "timestep", &timestep, sizeof(REAL), 0);
    set_field_psdata(data, "smoothingradius", &smoothingradius, sizeof(REAL), 0);

    REAL gridbounds[] = { xbound1, xbound2, ybound1, ybound2, zbound1, zbound2 };
    set_field_psdata(data, "gridbounds", gridbounds, sizeof(gridbounds), 0);

    unsigned int gridres[] = { num_gridcells_x, num_gridcells_y, num_gridcells_z };
    set_field_psdata(data, "gridres", gridres, sizeof(gridres), 0);

    data->num_host_fields = 0;
    }

    { /* MATLAB bollocks */
#ifdef MATLAB_MEX_FILE
    mxArray * position_mex     = mxCreateNumericArray( data->num_dimensions[5],
                                                       data->dimensions+data->dimensions_offsets[5],
                                                       mxDOUBLE_CLASS, mxREAL );
    mxArray * posnext_mex      = mxCreateNumericArray( data->num_dimensions[6],
                                                       data->dimensions+data->dimensions_offsets[6],
                                                       mxDOUBLE_CLASS, mxREAL );
    mxArray * velocity_mex     = mxCreateNumericArray( data->num_dimensions[7],
                                                       data->dimensions+data->dimensions_offsets[7],
                                                       mxDOUBLE_CLASS, mxREAL );
    mxArray * veleval_mex      = mxCreateNumericArray( data->num_dimensions[8],
                                                       data->dimensions+data->dimensions_offsets[8],
                                                       mxDOUBLE_CLASS, mxREAL );
    mxArray * velnext_mex      = mxCreateNumericArray( data->num_dimensions[9],
                                                       data->dimensions+data->dimensions_offsets[9],
                                                       mxDOUBLE_CLASS, mxREAL );
    mxArray * acceleration_mex = mxCreateNumericArray( data->num_dimensions[10],
                                                       data->dimensions+data->dimensions_offsets[10],
                                                       mxDOUBLE_CLASS, mxREAL );
    mxArray * force_mex        = mxCreateNumericArray( data->num_dimensions[11],
                                                       data->dimensions+data->dimensions_offsets[11],
                                                       mxDOUBLE_CLASS, mxREAL );
    mxArray * density_mex      = mxCreateNumericArray( data->num_dimensions[12],
                                                       data->dimensions+data->dimensions_offsets[12],
                                                       mxDOUBLE_CLASS, mxREAL );
    mxArray * volume_mex       = mxCreateNumericArray( data->num_dimensions[13],
                                                       data->dimensions+data->dimensions_offsets[13],
                                                       mxDOUBLE_CLASS, mxREAL );

    mexMakeArrayPersistent(position_mex);
    mexMakeArrayPersistent(posnext_mex);
    mexMakeArrayPersistent(velocity_mex);
    mexMakeArrayPersistent(veleval_mex);
    mexMakeArrayPersistent(velnext_mex);
    mexMakeArrayPersistent(acceleration_mex);
    mexMakeArrayPersistent(force_mex);
    mexMakeArrayPersistent(density_mex);
    mexMakeArrayPersistent(volume_mex);

    create_host_field_psdata(data, "position_mex",     position_mex,     sizeof(mxArray*));
    create_host_field_psdata(data, "posnext_mex",      posnext_mex,      sizeof(mxArray*));
    create_host_field_psdata(data, "velocity_mex",     velocity_mex,     sizeof(mxArray*));
    create_host_field_psdata(data, "veleval_mex",      veleval_mex,      sizeof(mxArray*));
    create_host_field_psdata(data, "velnext_mex",      velnext_mex,      sizeof(mxArray*));
    create_host_field_psdata(data, "acceleration_mex", acceleration_mex, sizeof(mxArray*));
    create_host_field_psdata(data, "force_mex",        force_mex,        sizeof(mxArray*));
    create_host_field_psdata(data, "density_mex",      density_mex,      sizeof(mxArray*));
    create_host_field_psdata(data, "volume_mex",       volume_mex,       sizeof(mxArray*));

    sync_to_mex(data);
#endif
    }
}

int get_field_psdata(psdata data, const char * name) {
    unsigned int i;
    for (i = 0; i < data.num_fields; ++i) {
        if (strcmp(data.names + data.names_offsets[i], name) == 0) return i;
    }

    return -1;
}

void set_field_psdata(psdata * data, const char * name, void * field, unsigned int size, unsigned int offset) {
    int i = get_field_psdata(*data, name);

    if (i != -1) memcpy(data->data + data->data_offsets[i] + offset, field, size);
}

unsigned int psdata_names_size(psdata data) {
    char * name = (char*) data.names + data.names_offsets[data.num_fields-1];
    return (data.names_offsets[data.num_fields-1]
          + strlen(name) + 1) * sizeof(char);
}

unsigned int psdata_dimensions_size(psdata data) {
    unsigned int nf = data.num_fields;
    return (data.dimensions_offsets[nf-1] + data.num_dimensions[nf-1]) * sizeof(unsigned int);
}

unsigned int psdata_data_size(psdata data) {
    return data.data_offsets[data.num_fields-1] + data.data_sizes[data.num_fields-1];
}

/* Data needs to be heap allocated but it should be anyway, right? RIGHT? */
/* Don't worry I'll allocate that string for you.. */
int create_host_field_psdata(psdata * data, const char * name, void * field, unsigned int size) {
    unsigned int new_num_host_fields = data->num_host_fields + 1;

    /* Allocate */
    char ** new_host_names = malloc(new_num_host_fields*sizeof(char*));
    void ** new_host_data = malloc(new_num_host_fields*sizeof(void*));
    unsigned int * new_host_data_size = malloc(new_num_host_fields*sizeof(unsigned int));

    /* Copy */
    unsigned int i;
    for (i = 0; i < data->num_host_fields; ++i) {
        new_host_names[i] = data->host_names[i];
        new_host_data[i] = data->host_data[i];
        new_host_data_size[i] = data->host_data_size[i];
    }

    /* Free old stuff */
    if (data->num_host_fields > 0) {
        free(data->host_names);
        free(data->host_data);
        free(data->host_data_size);
    }

    char * name_alloc = malloc((strlen(name)+1)*sizeof(char));
    strcpy(name_alloc, name);

    /* Add new field */
    new_host_names[i] = name_alloc;
    new_host_data[i] = field;
    new_host_data_size[i] = size;

    data->num_host_fields = new_num_host_fields;

    /* Attach to struct */
    data->host_names = new_host_names;
    data->host_data = new_host_data;
    data->host_data_size = new_host_data_size;

    return i;
}

int get_host_field_psdata(psdata * data, const char * name) {
    unsigned int i;
    for (i = 0; i < data->num_host_fields; ++i) {
        if (strcmp(data->host_names[i], name) == 0) return i;
    }

    return -1;
}

#ifdef MATLAB_MEX_FILE
void sync_to_mex(psdata * data) {
    void * field;
    mxArray * mex_field;

    int i, f;
    for (i = 0; i < data->num_host_fields; ++i) {
        if (!is_mex_field(data->host_names[i])) continue;

        mex_field = (mxArray*) data->host_data[i];

        int name_length = strlen(data->host_names[i]) - 4;
        char * field_name = malloc((name_length + 1)*sizeof(char));
        strncpy(field_name, data->host_names[i], name_length);
        field_name[name_length] = '\0';
        f = get_field_psdata(*data, field_name);

        free(field_name);

        if (f == -1) continue;

        void * mex_field_ptr = mxGetData(mex_field);
        
        memcpy(mex_field_ptr, ((char*) data->data) + data->data_offsets[f], data->data_sizes[f]);
        /* TODO: This is shit, fix it
        int j;
        for (j = 0; j < data->data_sizes[f]/sizeof(REAL); ++j) {
            REAL val = ((REAL*) (data->data + data->data_offsets[f]))[j];
            mex_field_ptr[j] = val;
        }
        */
    }
}

int is_mex_field(const char * name) {
    const char * mex_cmp = "_mex";
    size_t mex_len = strlen(mex_cmp);

    size_t name_len = strlen(name);
    if (name_len < mex_len) return 0;

    const char * name_end = name + name_len - mex_len;
    if (strcmp(name_end, mex_cmp) == 0) return 1;
    
    return 0;
}

void get_mex_field_name(const char * name, size_t * pReturnLength, char * mex_name) {
    const char * mext = "_mex";

    if (pReturnLength != NULL) {
        *pReturnLength = strlen(name) + strlen(mext);
        return;
    } else {
        sprintf(mex_name, "%s%s", name, mext);
    }
}
#endif

void free_psdata( psdata * data ) {
    free((char*) data->names);
    free(data->names_offsets);
    
    free(data->dimensions);
    free(data->num_dimensions);
    free(data->dimensions_offsets);
    free(data->entry_sizes);

#ifdef MATLAB_MEX_FILE
    mxFree(data->data);
#else
    free(data->data);
#endif
    free(data->data_sizes);
    free(data->data_offsets);

    unsigned int i;
    for (i = 0; i < data->num_host_fields; ++i) {
#ifdef MATLAB_MEX_FILE
        if (is_mex_field(data->host_names[i])) mxDestroyArray(data->host_data[i]);
        else free(data->host_data[i]);
#else
        free(data->host_data[i]);
#endif

        free(data->host_names[i]);
    }
    note(2, "freed host data\n");

    if (data->num_host_fields > 0) {
        free(data->host_names);
        free(data->host_data);
        free(data->host_data_size);
    }
}
