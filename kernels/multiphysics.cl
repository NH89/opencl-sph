#ifndef OPENCL_SPH_REAL_TYPE
#define OPENCL_SPH_REAL_TYPE float
#endif

typedef OPENCL_SPH_REAL_TYPE REAL;
#if OPENCL_SPH_REAL_TYPE == float
typedef float2 REAL2;
typedef float3 REAL3;
typedef float4 REAL4;
#elif OPENCL_SPH_REAL_TYPE == double
typedef double2 REAL2;
typedef double3 REAL3;
typedef double4 REAL4;
#else
#error "OPENCL_SPH_REAL_TYPE must be either float or double."
#endif

kernel void compute_forces_multiphysics (PSO_ARGS) {
    USE_GRID_PROPS

    USE_FIELD_FIRST_VALUE(n, uint) USE_FIELD_FIRST_VALUE(mass, REAL)
    USE_FIELD_FIRST_VALUE(restdens, REAL) USE_FIELD_FIRST_VALUE(stiffness, REAL)

    USE_FIELD(originalpos, REAL) USE_FIELD_FIRST_VALUE(smoothingradius, REAL)
    USE_FIELD(stress, REAL) USE_FIELD_FIRST_VALUE(mass, REAL) USE_FIELD(density0, REAL)
    USE_FIELD(rotation, REAL) USE_FIELD(force, REAL) USE_FIELD(velocity, REAL)
    USE_FIELD_FIRST_VALUE(viscosity, REAL) USE_FIELD(density, REAL) USE_FIELD(position, REAL)

    USE_FIELD(particletype, uint)

    unsigned int i = get_global_id(0);

    if (i >= n) return;

    REAL3 ipos = vload3(i, position);
    REAL3 ivel = vload3(i, velocity);

    INIT_FLUIDS_FORCE_COMPUTATION
    INIT_SOLIDS_FORCE_COMPUTATION

    REAL3 f_ext = (REAL3)(0, 0, 0);

    /*FOR_PARTICLES_IN_RANGE(i, j,
        if (j == i) continue;

        if (particletype[i] == 1 && particletype[j] == 1) {
            SOLIDS_FORCE_COMPUTATION
        } else if (particletype[i] == 2 && particletype[j] == 2) {
            FLUIDS_FORCE_COMPUTATION
        } else {
            // f_ext += Something something pressure gradient
        }
    )*/    
    {
    int gx, gy, gz;
    for (gz = -1; gz <= 1; ++gz)
        for (gy = -1; gy <= 1; ++gy){
            for (gx = -1; gx <= 1; ++gx){                                               // loop through 27 bins for all adjacent particles
                int cell = get_cell_at_offset(gridres, gridcell[i], gx, gy, gz);
                if (cell == -1) continue;                                               // if bin is not valid
                uint offset = celloffset[cell];
                for (uint jp = offset; jp < offset + gridcount[cell]; ++jp) {
                    uint j = cellparticles[jp];
                    if (j == i) continue;                                               // if same particle
                    if (particletype[i] == 1 && particletype[j] == 1) {
                        SOLIDS_FORCE_COMPUTATION                                        // solid-sold interaction 
                    } else if (particletype[i] == 2 && particletype[j] == 2) {
                        FLUIDS_FORCE_COMPUTATION                                        // fluid-fluid
                    } else {
                        // f_ext += Something something pressure gradient               // solid-fluid
                    }
                }  
            }
        }
    }

    FINALISE_SOLIDS_FORCE_COMPUTATION
    FINALISE_FLUIDS_FORCE_COMPUTATION

    REAL3 f = f_i + f_e + f_v + f_ext;

    APPLY_CUBE_BOUNDS(ipos, f, -2.0, 2.0)

    vstore3(f, i, force);
}
