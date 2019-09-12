#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "note.h"
#include "particle_system.h"
#include "build_psdata.h"
#include "opencl/particle_system_host.h"

static psdata * ps_instance = NULL;

void display_entry(psdata data, size_t offset, size_t size) {
    if (size == 8) {
        note(2, "%g, ", *((REAL*)((char*)data.data + offset)));
    } else if (size == 4) {
        note(2, "%u, ", *((unsigned int*)((char*)data.data + offset)));
    }
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

void write_psdata(psdata data, int number, const char* Case){
    for (size_t field = 0; field < data.num_fields; ++field){
        char * name = data.names + data.names_offsets[field];
        char snum[5];
        sprintf(snum, "%d", number);
        if (strcmp(name, "position") == 0){
            char filename[150];
            //strcpy(filename, "/media/aslab/data/hackthon_data/");
            //strcat(filename, Case);
            strcpy(filename, "position_");
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
void write_psdata_ply(psdata data, int number, const char* Case){
    for (size_t field = 0; field < data.num_fields; ++field){
        char * name = data.names + data.names_offsets[field];
        char snum[5];
        sprintf(snum, "%d", number);
        if (strcmp(name, "position") == 0){
            char filename[150];
            //strcpy(filename, "/media/aslab/data/hackthon_data/");
            //strcat(filename, Case);
            strcpy(filename, "position_");
            strcat(filename, snum);
            strcat(filename, ".ply");
            FILE *f = fopen(filename, "w");
            unsigned int d0 = (data.dimensions + data.dimensions_offsets[field])[0];
            unsigned int d1 = (data.dimensions + data.dimensions_offsets[field])[1];
            fprintf(f, "ply \n format ascii 1.0\n comment particle cloud from Fluids_v4\n element vertex %i\n", d1 );
            fprintf(f, "property float x\nproperty float y\nproperty float z\n");
            fprintf(f, "end_header\n");
            for (unsigned int i = 0; i < d1; ++i) {
                for (unsigned int j = 0; j < d0; ++j) {
                    fprintf(f, "%0.3lf ", *((REAL*)((char*)data.data + data.data_offsets[field] + (i*d0+j)*data.entry_sizes[field])));
                }
                fprintf(f,"\n");
            }
            fclose(f);
            fflush (f);
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

void free_psdata( psdata * data ) {
    free((char*) data->names);
    free(data->names_offsets);
    free(data->dimensions);
    free(data->num_dimensions);
    free(data->dimensions_offsets);
    free(data->entry_sizes);
    free(data->data);
    free(data->data_sizes);
    free(data->data_offsets);

    unsigned int i;
    for (i = 0; i < data->num_host_fields; ++i) {
        free(data->host_data[i]);
        free(data->host_names[i]);
    }
    note(2, "freed host data\n");

    if (data->num_host_fields > 0) {
        free(data->host_names);
        free(data->host_data);
        free(data->host_data_size);
    }
}
