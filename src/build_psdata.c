#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

#include "macros.h"
#include "build_psdata.h"
#include "note.h"
#include "particle_system.h"
#include "stringly.h"
#ifdef MATLAB_MEX_FILE
    #include "mex.h"
#endif

typedef struct psdata_field_spec psdata_field_spec;

struct psdata_field_spec
{
    char * name;
    unsigned int num_dimensions;
    unsigned int * dimensions;
    char * type; // Only REAL and unsigned int for now
    void * data; /* Optional */

    psdata_field_spec * next;
};

static psdata_field_spec * get_field_spec_by_name(psdata_field_spec * list, const char * name) {
    psdata_field_spec * field_cursor = list;

    while (field_cursor != NULL) {
        if (strcmp(field_cursor->name, name) == 0) {
            return field_cursor;
        }

        field_cursor = field_cursor->next;
    }

    return NULL;
}
static unsigned int get_value(psdata_field_spec * list, const char * name) { // Only for dimensions really
    psdata_field_spec * field_spec = get_field_spec_by_name(list, name);

    if (field_spec != NULL) {
        return *((unsigned int *)field_spec->data);
    }

    return UINT_MAX;
}

/* Depends on smoothingradius, gridbounds & pnum */
static void add_grid_arrays(psdata_field_spec * list) {
    psdata_field_spec * last = list;

    while (last != NULL && last->next != NULL) last = last->next;

    psdata_field_spec * smoothingradius_fs = get_field_spec_by_name(list, "smoothingradius");
    if (smoothingradius_fs->data == NULL) {
        note(1, "No value for smoothingradius\n");
        return;
    }
    REAL smoothingradius = *((REAL*)smoothingradius_fs->data);

    psdata_field_spec * gridbounds_fs = get_field_spec_by_name(list, "gridbounds");
    if (gridbounds_fs->data == NULL) {
        note(1, "No value for gridbounds\n");
        return;
    }
    REAL * gridbounds = gridbounds_fs->data;

    unsigned int pnum = get_value(list, "pnum");
    if (pnum == UINT_MAX) {
        note(1, "No value for pnum\n");
        return;
    }

    unsigned int num_gridcells_x = (unsigned int) (fabs(gridbounds[1] - gridbounds[0])/smoothingradius);
    unsigned int num_gridcells_y = (unsigned int) (fabs(gridbounds[3] - gridbounds[2])/smoothingradius);
    unsigned int num_gridcells_z = (unsigned int) (fabs(gridbounds[5] - gridbounds[4])/smoothingradius);

    char * name;
    unsigned int dimensions[5];
    char * type;

    // Gridres
    psdata_field_spec * gridres = malloc(sizeof(psdata_field_spec));
    last->next = gridres;

    name = "gridres";
    gridres->num_dimensions = 1;
    dimensions[0] = 3;
    type = "unsigned int";

    gridres->name = malloc(strlen(name)+1); strcpy(gridres->name, name);
    gridres->dimensions = malloc(gridres->num_dimensions * sizeof(unsigned int));
    memcpy(gridres->dimensions, dimensions, gridres->num_dimensions * sizeof(unsigned int));
    gridres->type = malloc(strlen(type)+1); strcpy(gridres->type, type);
    gridres->data = malloc(dimensions[0] * sizeof(unsigned int));
    ((unsigned int *) gridres->data)[0] = num_gridcells_x;
    ((unsigned int *) gridres->data)[1] = num_gridcells_y;
    ((unsigned int *) gridres->data)[2] = num_gridcells_z;

    // Gridcell
    psdata_field_spec * gridcell = malloc(sizeof(psdata_field_spec));
    gridres->next = gridcell;

    name = "gridcell";
    gridcell->num_dimensions = 1;
    dimensions[0] = pnum;
    type = "unsigned int";

    gridcell->name = malloc(strlen(name)+1); strcpy(gridcell->name, name);
    gridcell->dimensions = malloc(gridcell->num_dimensions * sizeof(unsigned int));
    memcpy(gridcell->dimensions, dimensions, gridcell->num_dimensions * sizeof(unsigned int));
    gridcell->type = malloc(strlen(type)+1); strcpy(gridcell->type, type);
    gridcell->data = NULL;

    // Gridcount
    psdata_field_spec * gridcount = malloc(sizeof(psdata_field_spec));
    gridcell->next = gridcount;

    name = "gridcount";
    gridcount->num_dimensions = 3;
    dimensions[0] = num_gridcells_x; dimensions[1] = num_gridcells_y; dimensions[2] = num_gridcells_z;
    type = "unsigned int";

    gridcount->name = malloc(strlen(name)+1); strcpy(gridcount->name, name);
    gridcount->dimensions = malloc(gridcount->num_dimensions * sizeof(unsigned int));
    memcpy(gridcount->dimensions, dimensions, gridcount->num_dimensions * sizeof(unsigned int));
    gridcount->type = malloc(strlen(type)+1); strcpy(gridcount->type, type);
    gridcount->data = NULL;

    // Celloffset
    psdata_field_spec * celloffset = malloc(sizeof(psdata_field_spec));
    gridcount->next = celloffset;

    name = "celloffset";
    celloffset->num_dimensions = 3;
    type = "unsigned int";

    celloffset->name = malloc(strlen(name)+1); strcpy(celloffset->name, name);
    celloffset->dimensions = malloc(celloffset->num_dimensions * sizeof(unsigned int));
    memcpy(celloffset->dimensions, dimensions, celloffset->num_dimensions * sizeof(unsigned int));
    celloffset->type = malloc(strlen(type)+1); strcpy(celloffset->type, type);
    celloffset->data = NULL;

    // Cellparticle
    psdata_field_spec * cellparticles = malloc(sizeof(psdata_field_spec));
    celloffset->next = cellparticles;

    name = "cellparticles";
    cellparticles->num_dimensions = 1;
    dimensions[0] = pnum;
    type = "unsigned int";

    cellparticles->name = malloc(strlen(name)+1); strcpy(cellparticles->name, name);
    cellparticles->dimensions = malloc(cellparticles->num_dimensions * sizeof(unsigned int));
    memcpy(cellparticles->dimensions, dimensions, cellparticles->num_dimensions * sizeof(unsigned int));
    cellparticles->type = malloc(strlen(type)+1); strcpy(cellparticles->type, type);
    cellparticles->data = NULL;

    cellparticles->next = NULL;
}
static psdata_field_spec * create_psdata_field_spec(char * line, psdata_field_spec * previous, psdata_field_spec * list) {
#define DIMENSIONS_PAD_LENGTH 5
#define TYPE_PAD_LENGTH 64
#define DATA_PAD_LENGTH 256

    psdata_field_spec * field_spec = malloc(sizeof(psdata_field_spec));

    char * endline = line + strlen(line);
    char * endsegment = line - 1;

    unsigned int segmentno = 0;

    unsigned int dimensions_pad[DIMENSIONS_PAD_LENGTH];
    field_spec->num_dimensions = 0;

    char type_pad[TYPE_PAD_LENGTH];
    type_pad[0] = '\0';

    char data_pad[DATA_PAD_LENGTH];
    char * data_ptr = data_pad;
    field_spec->data = NULL;

    while (endsegment < endline) {
        char * segment = strtok(endsegment + 1, "|");
        if (segment == NULL) break;

        endsegment = segment + strlen(segment);

        unsigned int wordno = 0;

        char * word = strtok(segment, ", ");
        while (word != NULL) {
            if (segmentno == 0 && wordno == 0)
            {
                field_spec->name = malloc(sizeof(char)*(strlen(word)+1));
                strcpy(field_spec->name, word);
            }
            else if (segmentno == 1)
            {
                if (field_spec->num_dimensions == DIMENSIONS_PAD_LENGTH) {
                    note(1, "%s: Over maximum number of dimensions\n", field_spec->name);
                    break;
                }

                unsigned int dim = strtol(word, NULL, 0);

                if (dim == 0 && list != NULL) {
                    dim = get_value(list, word);
                    if (dim == UINT_MAX) dim = 0;
                }

                dimensions_pad[field_spec->num_dimensions] = dim;
                ++field_spec->num_dimensions;
            }
            else if (segmentno == 2)
            {
                size_t prevlength = strlen(type_pad);
                size_t wordlength = strlen(word);

                if (prevlength + wordlength + 1 >= TYPE_PAD_LENGTH) {
                    note(1, "%s: Over maximum type string length\n", field_spec->name);
                    break;
                }

                if (prevlength > 0) {
                    type_pad[prevlength] = ' ';
                    ++prevlength;
                }
                
                strcpy(type_pad + prevlength, word);
            }
            else if (segmentno == 3)
            {
                if (strcmp(type_pad, "REAL") == 0) {
                    REAL val = strtod(word, NULL);
                    if (data_ptr + sizeof(REAL) - data_pad > DATA_PAD_LENGTH) {
                        note(1, "%s: Over maximum default value length\n", field_spec->name);
                        break;
                    }

                    memcpy(data_ptr, &val, sizeof(REAL));
                    data_ptr += sizeof(REAL);
                } else if (strcmp(type_pad, "unsigned int") == 0) {
                    unsigned int val = strtol(word, NULL, 0);
                    if (data_ptr + sizeof(unsigned int) - data_pad > DATA_PAD_LENGTH) {
                        note(1, "%s: Over maximum default value length\n", field_spec->name);
                        break;
                    }

                    memcpy(data_ptr, &val, sizeof(unsigned int));
                    data_ptr += sizeof(unsigned int);
                } else if (strcmp(type_pad, "int") == 0) {
                    int val = strtol(word, NULL, 0);
                    if (data_ptr + sizeof(int) - data_pad > DATA_PAD_LENGTH) {
                        note(1, "%s: Over maximum default value length\n", field_spec->name);
                        break;
                    }

                    memcpy(data_ptr, &val, sizeof(int));
                    data_ptr += sizeof(int);
                }
            }

            word = strtok(NULL, ", ");

            ++wordno;
        }

        ++segmentno;
    }

    field_spec->dimensions = malloc(field_spec->num_dimensions * sizeof(unsigned int));
    memcpy(field_spec->dimensions, dimensions_pad, field_spec->num_dimensions * sizeof(unsigned int));

    field_spec->type = malloc((strlen(type_pad)+1) * sizeof(char));
    strcpy(field_spec->type, type_pad);

    if (data_ptr - data_pad > 0) {
        field_spec->data = malloc(data_ptr - data_pad);
        memcpy(field_spec->data, data_pad, data_ptr - data_pad);
    } else {
        field_spec->data = NULL;
    }

    field_spec->next = NULL;

    if (previous != NULL) previous->next = field_spec;

    return field_spec;

#undef DIMENSIONS_PAD_LENGTH
#undef TYPE_PAD_LENGTH
#undef DATA_PAD_LENGTH
}
static void populate_psdata(psdata * data, psdata_field_spec * list) {
    psdata_field_spec * field_cursor = list;

    unsigned int num_fields = 0;

    while (field_cursor != NULL) {
        field_cursor = field_cursor->next;
        ++num_fields;
    }

    data->num_fields = num_fields;

    data->names_offsets = malloc(num_fields * sizeof(unsigned int));

    data->num_dimensions = malloc(num_fields * sizeof(unsigned int));
    data->dimensions_offsets = malloc(num_fields * sizeof(unsigned int));
    data->entry_sizes = malloc(num_fields * sizeof(unsigned int));

    data->data_sizes = malloc(num_fields * sizeof(unsigned int));
    data->data_offsets = malloc(num_fields * sizeof(unsigned int));

    // Get sizes and offsets

    size_t names_length = 0, dimensions_length = 0, data_size = 0;

    field_cursor = list;

    int f = 0;
    while (field_cursor != NULL) {
        data->names_offsets[f] = names_length;
        names_length += strlen(field_cursor->name)+1;

        data->num_dimensions[f] = field_cursor->num_dimensions;
        data->dimensions_offsets[f] = dimensions_length;
        data->entry_sizes[f] = sizeof_string_type(field_cursor->type);
        dimensions_length += field_cursor->num_dimensions;

        unsigned int this_data_size = data->entry_sizes[f];

        unsigned int d;
        for (d = 0; d < field_cursor->num_dimensions; ++d) this_data_size *= field_cursor->dimensions[d];

        data->data_sizes[f] = this_data_size;
        data->data_offsets[f] = data_size;
        data_size += this_data_size;

        field_cursor = field_cursor->next;
        ++f;
    }

    data->names = malloc(names_length * sizeof(char));
    data->dimensions = malloc(dimensions_length * sizeof(unsigned int));
#ifndef MATLAB_MEX_FILE
    data->data = calloc(data_size, 1);
#else
    data->data = mxCalloc(data_size, 1);

    mexMakeMemoryPersistent(data->data);
#endif

    // Copy in data

    field_cursor = list;

    f = 0;
    while (field_cursor != NULL) {
        strcpy(((char*)data->names) + data->names_offsets[f], field_cursor->name);

        memcpy(data->dimensions + data->dimensions_offsets[f], field_cursor->dimensions, field_cursor->num_dimensions * sizeof(unsigned int));

        if (field_cursor->data != NULL) {
            memcpy(((char*)data->data) + data->data_offsets[f], field_cursor->data, data->data_sizes[f]);
        }

        field_cursor = field_cursor->next;
        ++f;
    }

    data->num_host_fields = 0;

#ifdef MATLAB_MEX_FILE
    const char * mex_suffix = "_mex";
    size_t mex_suffix_length = strlen(mex_suffix);
    
    unsigned int * zerodims = { 0 };

    field_cursor = list;

    f = 0;
    while (field_cursor != NULL) {
        int use_field = 1;

        char * host_name = malloc(strlen(field_cursor->name) + mex_suffix_length + 1);

        sprintf(host_name, "%s%s", field_cursor->name, mex_suffix);

        mxArray * data_mex;

        if (strcmp(field_cursor->type, "REAL") == 0) {

            data_mex = mxCreateNumericArray(field_cursor->num_dimensions,
                                            field_cursor->dimensions, mxDOUBLE_CLASS, mxREAL);

        } else if (strcmp(field_cursor->type, "int") == 0) {

            data_mex = mxCreateNumericArray(field_cursor->num_dimensions,
                                            field_cursor->dimensions, mxINT32_CLASS, mxREAL);

        } else if (strcmp(field_cursor->type, "unsigned int") == 0) {

            data_mex = mxCreateNumericArray(field_cursor->num_dimensions,
                                            field_cursor->dimensions, mxUINT32_CLASS, mxREAL);

        } else {
            use_field = 0;
        }

        if (use_field) {
            mexMakeArrayPersistent(data_mex);

            create_host_field_psdata(data, host_name, data_mex, sizeof(mxArray*));
        }

        free(host_name);

        field_cursor = field_cursor->next;
        ++f;
    }

    sync_to_mex(data);
#endif
}
static void free_psdata_field_spec_list(psdata_field_spec * list) {
    free(list->name);
    free(list->dimensions);
    free(list->type);
    if (list->data != NULL) free(list->data);

    if (list->next != NULL) free_psdata_field_spec_list(list->next);

    free(list);
}
static void print_field_spec_list_names(psdata_field_spec * list) {
    size_t i = 0;

    while (list != NULL) {
        note(2, "%u: %s\n", i, list->name);
        list = list->next;
        ++i;
    }
}
static size_t length_field_spec_list(psdata_field_spec * list) {
    size_t length = 1;

    while (list->next != NULL) {
        ++length;

        list = list->next;
    }

    return length;
}
// Guarantees to link up all the entries in the given list 
// so they can be freed with the new root node. 
// Sorting before creating the single large array prevents alignment issues
static psdata_field_spec * sort_field_spec_list_by_type_size_descending(psdata_field_spec * list) {
    if (list == NULL) return NULL;

    size_t list_length = length_field_spec_list(list);
    psdata_field_spec ** entry_array = malloc(list_length * sizeof(psdata_field_spec*));
    char * used = calloc(list_length, sizeof(char));

    const char * type_order[] = {
        "REAL", "unsigned int"
    };

    size_t num_types = sizeof type_order / sizeof(char*);

    size_t type_number = 0;

    psdata_field_spec * list_ptr = list;
    size_t field_pos = 0;

    for (size_t i = 0; i < list_length; ++i) {
        // Scroll through non-matching types
        while ((type_number < num_types && strcmp(list_ptr->type, type_order[type_number]) != 0) ||
               // or if we ran out of types, through anything we already used - this'll cause 
               // problems later but I guess it's not our concern here
               used[field_pos]) {
            if (list_ptr->next == NULL) {
                ++type_number;
                list_ptr = list;
                field_pos = 0;
            } else {
                list_ptr = list_ptr->next;
                ++field_pos;
            }
        }

        entry_array[i] = list_ptr;
        used[field_pos] = 1;

        if (list_ptr->next) {
            list_ptr = list_ptr->next;
            ++field_pos;
        } else if (type_number >= num_types) {
            ASSERT(i == list_length - 1);
        }
    }

    for (size_t i = 0; i < list_length-1; ++i) {
        entry_array[i]->next = entry_array[i+1];
    }

    entry_array[list_length-1]->next = NULL;

    psdata_field_spec * reordered = entry_array[0];

    free(entry_array);
    free(used);

    return reordered;
}



void build_psdata_from_string(psdata * data, const char * string) {
    char * string_copy = malloc((strlen(string)+1) * sizeof(char));
    strcpy(string_copy, string);

    char * endstring = string_copy + strlen(string_copy);
    char * endline = string_copy - 1;

    psdata_field_spec * first = NULL;
    psdata_field_spec * current = NULL;

    while (endline < endstring) {
        char * line = strtok(endline + 1, "\n");
        if (line == NULL) break;

        endline = line + strlen(line);

        if (endline > line && line[0] == '#') continue;

        current = create_psdata_field_spec(line, current, first);

        if (first == NULL) first = current;
    }

    free(string_copy);

    add_grid_arrays(first);

    psdata_field_spec * reordered = sort_field_spec_list_by_type_size_descending(first);

    populate_psdata(data, reordered);

    free_psdata_field_spec_list(reordered);
}


void build_psdata(psdata * data, const char * path) {
    FILE * conf = fopen(path, "rb");

    if (conf == NULL) {
        note(1, "Could not read file %s\n", path);
        ASSERT(0);
    }

    fseek(conf, 0, SEEK_END);
    size_t conf_end = ftell(conf);

    char * contents = malloc(conf_end);

    fseek(conf, 0, SEEK_SET);
    fread(contents, 1, conf_end, conf);
    fclose(conf);

    contents[conf_end-1] = '\0';

    build_psdata_from_string(data, contents);

    free(contents);
}
