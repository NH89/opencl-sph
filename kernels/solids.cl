// solids.cl ///////////////////////////////////////////////////////////////////////////////////////

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

////////// Physics stuff //////////

inline void contributeApq (REAL3 p, REAL3 q, REAL weight, REAL apq[9]) {
    REAL3 wp = weight*p;

    apq[I_00] += wp.x*q.x; apq[I_01] += wp.x*q.y; apq[I_02] += wp.x*q.z;
    apq[I_10] += wp.y*q.x; apq[I_11] += wp.y*q.y; apq[I_12] += wp.y*q.z;
    apq[I_20] += wp.z*q.x; apq[I_21] += wp.z*q.y; apq[I_22] += wp.z*q.z;
}
inline void computeCauchyStrain (REAL deformation[9], global REAL strain[6]) {
    strain[S_00] = deformation[I_00];
    strain[S_01] = 0.5*(deformation[I_01] + deformation[I_10]);
    strain[S_02] = 0.5*(deformation[I_02] + deformation[I_20]);
    strain[S_11] = deformation[I_11];
    strain[S_12] = 0.5*(deformation[I_12] + deformation[I_21]);
    strain[S_22] = deformation[I_22];
}
inline void computeStress (global REAL strain[6], REAL bulk_modulus, REAL shear_modulus, global REAL stress[6]) {
    REAL lame_lambda = bulk_modulus - (2/3) * shear_modulus;
    REAL lame_mu = shear_modulus;

    REAL c[36] = {
        lame_lambda + 2 * lame_mu, lame_lambda              , lame_lambda              , 0          , 0          , 0          ,
        lame_lambda              , lame_lambda + 2 * lame_mu, lame_lambda              , 0          , 0          , 0          ,
        lame_lambda              , lame_lambda              , lame_lambda + 2 * lame_mu, 0          , 0          , 0          ,
        0                        , 0                        , 0                        , 2 * lame_mu, 0          , 0          ,
        0                        , 0                        , 0                        , 0          , 2 * lame_mu, 0          ,
        0                        , 0                        , 0                        , 0          , 0          , 2 * lame_mu,
    };

    for (uint i = 0; i < 6; ++i) {
        uint ir = i*6;
        stress[i] = c[ir]*strain[0] + c[ir+1]*strain[1] + c[ir+2]*strain[2] + c[ir+3]*strain[3] + c[ir+4]*strain[4] + c[ir+5]*strain[5];
    }
}

////////// Particle creation and transformation //////////

kernel void init_original_position (PSO_ARGS) {
    USE_FIELD(position, REAL) 
    USE_FIELD(originalpos, REAL)
    
    USE_FIELD_FIRST_VALUE(n, uint)

    unsigned int i = get_global_id(0);

    if (i >= n) return;

    REAL3 pos = vload3(i, position);

    vstore3(pos, i, originalpos);
}

////////// Kernels //////////

kernel void compute_rotations_and_strains (PSO_ARGS) {
    USE_GRID_PROPS

    USE_FIELD_FIRST_VALUE(n, uint)
    USE_FIELD_FIRST_VALUE(smoothingradius, REAL)
    USE_FIELD_FIRST_VALUE(mass, REAL)
    
    USE_FIELD(rotation, REAL) 
    USE_FIELD(strain, REAL) 
    USE_FIELD(position, REAL)
    USE_FIELD(originalpos, REAL) 
    USE_FIELD(density0, REAL) 
    USE_FIELD(density, REAL)

    uint i = get_global_id(0);

    if (i >= n) return;

    REAL3 ipos = vload3(i, position);
    REAL3 ipos0 = vload3(i, originalpos);

    REAL apq[9];

    for (uint d = 0; d < 9; ++d) apq[d] = 0;

    /*FOR_PARTICLES_IN_RANGE(i, j,
        if (j == i) continue;

        REAL3 p = vload3(j, position) - ipos;
        REAL3 q = vload3(j, originalpos) - ipos0;

        contributeApq(p, q, applyKernel(length(p), smoothingradius), apq);
    )*/  
    {
    int gx, gy, gz;
    for (gz = -1; gz <= 1; ++gz)
        for (gy = -1; gy <= 1; ++gy){
            for (gx = -1; gx <= 1; ++gx){
                int cell = get_cell_at_offset(gridres, gridcell[i], gx, gy, gz);
                if (cell == -1) continue;
                uint offset = celloffset[cell];
                for (uint jp = offset; jp < offset + gridcount[cell]; ++jp) {
                    uint j = cellparticles[jp];
                    j = cellparticles[jp];
                    if (j == i) continue;
                    REAL3 p = vload3(j, position) - ipos;
                    REAL3 q = vload3(j, originalpos) - ipos0;
                    contributeApq(p, q, applyKernel(length(p), smoothingradius), apq);
                }  
            }
        }
    }

    global REAL * r = rotation + i*3*3;

    getR(apq, r);

    REAL r_t[9];

    transpose(r, r_t);

    REAL deformation[9];

    for (uint d = 0; d < 9; ++d) deformation[d] = 0;

    /*FOR_PARTICLES_IN_RANGE(i, j,
        if (j == i) continue;

        REAL3 p = vload3(j, position) - ipos;
        p = multiplyMatrixVectorPrivate(r_t, p);

        REAL3 q = vload3(j, originalpos) - ipos0;

        REAL3 u = p-q;

        addOuterProduct(u*mass/density0[i], applyKernelGradient(q, smoothingradius), deformation);
    )*/ 
    {
    int gx, gy, gz;
    for (gz = -1; gz <= 1; ++gz)
        for (gy = -1; gy <= 1; ++gy){
            for (gx = -1; gx <= 1; ++gx){
                int cell = get_cell_at_offset(gridres, gridcell[i], gx, gy, gz);
                if (cell == -1) continue;
                uint offset = celloffset[cell];
                for (uint jp = offset; jp < offset + gridcount[cell]; ++jp) {
                    uint j = cellparticles[jp];
                    if (j == i) continue;
                    REAL3 p = vload3(j, position) - ipos;
                    p = multiplyMatrixVectorPrivate(r_t, p);
                    REAL3 q = vload3(j, originalpos) - ipos0;
                    REAL3 u = p-q;
                    addOuterProduct(u*mass/density0[i], applyKernelGradient(q, smoothingradius), deformation);
                }  
            }
        }
    }

    global REAL * s = strain + i*6;

    zeroArray(s, 6);

    computeCauchyStrain(deformation, s);
}
kernel void compute_original_density (PSO_ARGS) {
    USE_GRID_PROPS

    USE_FIELD_FIRST_VALUE(pnum, uint) USE_FIELD_FIRST_VALUE(n, uint) USE_FIELD_FIRST_VALUE(mass, REAL)
    USE_FIELD(position, REAL) USE_FIELD(density0, REAL) USE_FIELD_FIRST_VALUE(smoothingradius, REAL)

    unsigned int i = get_global_id(0);

    if (i >= n) return;

    REAL3 ipos = vload3(i, position);

    density0[i] = 0;

    /*FOR_PARTICLES_IN_RANGE(i, j,
        j = cellparticles[jp];

        REAL3 jpos = vload3(j, position);

        REAL3 diff = jpos - ipos;

        REAL dist = length(diff);

        density0[i] += applyKernel(dist, smoothingradius);
    )*/  
    {
    int gx, gy, gz;
    for (gz = -1; gz <= 1; ++gz)
        for (gy = -1; gy <= 1; ++gy){
            for (gx = -1; gx <= 1; ++gx){
                int cell = get_cell_at_offset(gridres, gridcell[i], gx, gy, gz);
                if (cell == -1) continue;
                uint offset = celloffset[cell];
                for (uint jp = offset; jp < offset + gridcount[cell]; ++jp) {
                    uint j = cellparticles[jp];
                    REAL3 jpos = vload3(j, position);
                    REAL3 diff = jpos - ipos;
                    REAL dist = length(diff);
                    density0[i] += applyKernel(dist, smoothingradius);
                }  
            }
        }
    }
    density0[i] *= mass;
}
kernel void compute_stresses (PSO_ARGS) {
    USE_FIELD(strain, REAL) 
    USE_FIELD(stress, REAL)
    
    USE_FIELD_FIRST_VALUE(bulk_modulus, REAL) 
    USE_FIELD_FIRST_VALUE(shear_modulus, REAL)
    USE_FIELD_FIRST_VALUE(n, uint)
    
    unsigned int i = get_global_id(0);
    if (i >= n) return;
    computeStress(strain + i*6, bulk_modulus, shear_modulus, stress + i*6);
}

kernel void compute_forces_solids (PSO_ARGS) {
    USE_GRID_PROPS
    USE_FIELD_FIRST_VALUE(n, uint)
    USE_FIELD_FIRST_VALUE(smoothingradius, REAL)
    USE_FIELD_FIRST_VALUE(mass, REAL)
    USE_FIELD_FIRST_VALUE(viscosity, REAL)
    
    USE_FIELD(originalpos, REAL) 
    USE_FIELD(stress, REAL)  
    USE_FIELD(density0, REAL)
    USE_FIELD(rotation, REAL) 
    USE_FIELD(force, REAL) 
    USE_FIELD(velocity, REAL)
    USE_FIELD(density, REAL) 
    USE_FIELD(position, REAL)
    
    unsigned int i = get_global_id(0);
    if (i >= n) return;
    REAL3 ipos = vload3(i, position);
    REAL3 ivel = vload3(i, velocity);
/*
// Repurpose for later batches
#define INIT_SOLIDS_FORCE_COMPUTATION \
    REAL3 ipos0 = vload3(i, originalpos);\
\
    global REAL * i_stress = stress + i*6;\
    global REAL * i_rotation = rotation + i*3*3;\
\
    REAL3 f_e = (REAL3)(0, 0, 0);\
    REAL3 f_v = (REAL3)(0, 0, 0);

    INIT_SOLIDS_FORCE_COMPUTATION
    */  
    REAL3 ipos0 = vload3(i, originalpos);
    global REAL * i_stress = stress + i*6;
    global REAL * i_rotation = rotation + i*3*3;
    REAL3 f_e = (REAL3)(0, 0, 0);
    REAL3 f_v = (REAL3)(0, 0, 0);

    /*FOR_PARTICLES_IN_RANGE(i, j,
        if (j == i) continue;

#define SOLIDS_FORCE_COMPUTATION \
        REAL3 jpos0 = vload3(j, originalpos);\
\
        global REAL * j_stress = stress + j*6;\
\
        REAL3 x = jpos0 - ipos0;\
        REAL vivj = mass*mass/density0[i]/density0[j];\
        REAL3 f_ji = -vivj * multiplySymMatrixVector(i_stress, applyKernelGradient(x, smoothingradius));\
        REAL3 f_ij = -vivj * multiplySymMatrixVector(j_stress, applyKernelGradient(-x, smoothingradius));\
\
        global REAL * j_rotation = rotation + j*3*3;\
\
        f_e += -multiplyMatrixVector(i_rotation, f_ji) + multiplyMatrixVector(j_rotation, f_ij);\
\
        REAL3 jpos = vload3(j, position);\
        REAL3 jvel = vload3(j, velocity);\
\
        f_v += mass/density[j] * (jvel - ivel) * applyLapKernel(length(jpos - ipos), smoothingradius);

        SOLIDS_FORCE_COMPUTATION
    )*/  
    
    {
    int gx, gy, gz;
    for (gz = -1; gz <= 1; ++gz)
        for (gy = -1; gy <= 1; ++gy){
            for (gx = -1; gx <= 1; ++gx){
                int cell = get_cell_at_offset(gridres, gridcell[i], gx, gy, gz);
                if (cell == -1) continue;
                uint offset = celloffset[cell];
                for (uint jp = offset; jp < offset + gridcount[cell]; ++jp) {
                    uint j = cellparticles[jp];
                    if (j == i) continue;
                    REAL3 jpos0 = vload3(j, originalpos);
                    global REAL * j_stress = stress + j*6;
                    REAL3 x = jpos0 - ipos0;
                    REAL vivj = mass*mass/density0[i]/density0[j];
                    REAL3 f_ji = -vivj * multiplySymMatrixVector(i_stress, applyKernelGradient(x, smoothingradius));
                    REAL3 f_ij = -vivj * multiplySymMatrixVector(j_stress, applyKernelGradient(-x, smoothingradius));
                    global REAL * j_rotation = rotation + j*3*3;
                    f_e += -multiplyMatrixVector(i_rotation, f_ji) + multiplyMatrixVector(j_rotation, f_ij);
                    REAL3 jpos = vload3(j, position);
                    REAL3 jvel = vload3(j, velocity);
                    f_v += mass/density[j] * (jvel - ivel) * applyLapKernel(length(jpos - ipos), smoothingradius);
                }  
            }
        }
    }
    

#define FINALISE_SOLIDS_FORCE_COMPUTATION \
    f_e *= 0.5f;\
    f_v *= viscosity;

    FINALISE_SOLIDS_FORCE_COMPUTATION

    REAL3 f = f_e + f_v;

    APPLY_CUBE_BOUNDS(ipos, f, -2.0f, 2.0f)

    vstore3(f, i, force);


}
