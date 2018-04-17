////////// Physics stuff //////////

inline void contributeApq (double3 p, double3 q, double weight, double apq[9]) {
    double3 wp = weight*p;

    apq[I_00] += wp.x*q.x; apq[I_01] += wp.x*q.y; apq[I_02] += wp.x*q.z;
    apq[I_10] += wp.y*q.x; apq[I_11] += wp.y*q.y; apq[I_12] += wp.y*q.z;
    apq[I_20] += wp.z*q.x; apq[I_21] += wp.z*q.y; apq[I_22] += wp.z*q.z;
}
inline void computeCauchyStrain (double deformation[9], global double strain[6]) {
    strain[S_00] = deformation[I_00];
    strain[S_01] = 0.5*(deformation[I_01] + deformation[I_10]);
    strain[S_02] = 0.5*(deformation[I_02] + deformation[I_20]);
    strain[S_11] = deformation[I_11];
    strain[S_12] = 0.5*(deformation[I_12] + deformation[I_21]);
    strain[S_22] = deformation[I_22];
}
inline void computeStress (global double strain[6], double bulk_modulus, double shear_modulus, global double stress[6]) {
    double lame_lambda = bulk_modulus - (2/3) * shear_modulus;
    double lame_mu = shear_modulus;

    double c[36] = {
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
    USE_FIELD(position, double) USE_FIELD(originalpos, double)
    USE_FIELD_FIRST_VALUE(n, uint)

    unsigned int i = get_global_id(0);

    if (i >= n) return;

    double3 pos = vload3(i, position);

    vstore3(pos, i, originalpos);
}

////////// Kernels //////////

kernel void compute_rotations_and_strains (PSO_ARGS) {
    USE_GRID_PROPS

    USE_FIELD_FIRST_VALUE(n, uint)

    USE_FIELD(rotation, double) USE_FIELD(strain, double) USE_FIELD(position, double)
    USE_FIELD(originalpos, double) USE_FIELD_FIRST_VALUE(smoothingradius, double)
    USE_FIELD(density0, double) USE_FIELD_FIRST_VALUE(mass, double)
    USE_FIELD(density, double)

    uint i = get_global_id(0);

    if (i >= n) return;

    double3 ipos = vload3(i, position);
    double3 ipos0 = vload3(i, originalpos);

    double apq[9];

    for (uint d = 0; d < 9; ++d) apq[d] = 0;

    FOR_PARTICLES_IN_RANGE(i, j,
        if (j == i) continue;

        double3 p = vload3(j, position) - ipos;
        double3 q = vload3(j, originalpos) - ipos0;

        contributeApq(p, q, applyKernel(length(p), smoothingradius), apq);
    )

    global double * r = rotation + i*3*3;

    getR(apq, r);

    double r_t[9];

    transpose(r, r_t);

    double deformation[9];

    for (uint d = 0; d < 9; ++d) deformation[d] = 0;

    FOR_PARTICLES_IN_RANGE(i, j,
        if (j == i) continue;

        double3 p = vload3(j, position) - ipos;
        p = multiplyMatrixVectorPrivate(r_t, p);

        double3 q = vload3(j, originalpos) - ipos0;

        double3 u = p-q;

        addOuterProduct(u*mass/density0[i], applyKernelGradient(q, smoothingradius), deformation);
    )

    global double * s = strain + i*6;

    zeroArray(s, 6);

    computeCauchyStrain(deformation, s);
}
kernel void compute_original_density (PSO_ARGS) {
    USE_GRID_PROPS

    USE_FIELD_FIRST_VALUE(pnum, uint) USE_FIELD_FIRST_VALUE(n, uint) USE_FIELD_FIRST_VALUE(mass, double)
    USE_FIELD(position, double) USE_FIELD(density0, double) USE_FIELD_FIRST_VALUE(smoothingradius, double)

    unsigned int i = get_global_id(0);

    if (i >= n) return;

    double3 ipos = vload3(i, position);

    density0[i] = 0;

    FOR_PARTICLES_IN_RANGE(i, j,
        j = cellparticles[jp];

        double3 jpos = vload3(j, position);

        double3 diff = jpos - ipos;

        double dist = length(diff);

        density0[i] += applyKernel(dist, smoothingradius);
    )

    density0[i] *= mass;
}
kernel void compute_stresses (PSO_ARGS) {
    USE_FIELD(strain, double) USE_FIELD(stress, double)
    USE_FIELD_FIRST_VALUE(bulk_modulus, double) USE_FIELD_FIRST_VALUE(shear_modulus, double)

    USE_FIELD_FIRST_VALUE(n, uint)

    unsigned int i = get_global_id(0);

    if (i >= n) return;

    computeStress(strain + i*6, bulk_modulus, shear_modulus, stress + i*6);
}

kernel void compute_forces_solids (PSO_ARGS) {
    USE_GRID_PROPS

    USE_FIELD_FIRST_VALUE(n, uint)

    USE_FIELD(originalpos, double) USE_FIELD_FIRST_VALUE(smoothingradius, double)
    USE_FIELD(stress, double) USE_FIELD_FIRST_VALUE(mass, double) USE_FIELD(density0, double)
    USE_FIELD(rotation, double) USE_FIELD(force, double) USE_FIELD(velocity, double)
    USE_FIELD_FIRST_VALUE(viscosity, double) USE_FIELD(density, double) USE_FIELD(position, double)

    unsigned int i = get_global_id(0);

    if (i >= n) return;

    double3 ipos = vload3(i, position);
    double3 ivel = vload3(i, velocity);

// Repurpose for later batches
#define INIT_SOLIDS_FORCE_COMPUTATION \
    double3 ipos0 = vload3(i, originalpos);\
\
    global double * i_stress = stress + i*6;\
    global double * i_rotation = rotation + i*3*3;\
\
    double3 f_e = (double3)(0, 0, 0);\
    double3 f_v = (double3)(0, 0, 0);

    INIT_SOLIDS_FORCE_COMPUTATION

    FOR_PARTICLES_IN_RANGE(i, j,
        if (j == i) continue;

#define SOLIDS_FORCE_COMPUTATION \
        double3 jpos0 = vload3(j, originalpos);\
\
        global double * j_stress = stress + j*6;\
\
        double3 x = jpos0 - ipos0;\
        double vivj = mass*mass/density0[i]/density0[j];\
        double3 f_ji = -vivj * multiplySymMatrixVector(i_stress, applyKernelGradient(x, smoothingradius));\
        double3 f_ij = -vivj * multiplySymMatrixVector(j_stress, applyKernelGradient(-x, smoothingradius));\
\
        global double * j_rotation = rotation + j*3*3;\
\
        f_e += -multiplyMatrixVector(i_rotation, f_ji) + multiplyMatrixVector(j_rotation, f_ij);\
\
        double3 jpos = vload3(j, position);\
        double3 jvel = vload3(j, velocity);\
\
        f_v += mass/density[j] * (jvel - ivel) * applyLapKernel(length(jpos - ipos), smoothingradius);

        SOLIDS_FORCE_COMPUTATION
    )

#define FINALISE_SOLIDS_FORCE_COMPUTATION \
    f_e *= 0.5;\
    f_v *= viscosity;

    FINALISE_SOLIDS_FORCE_COMPUTATION


    double3 f = f_e + f_v;

    APPLY_CUBE_BOUNDS(ipos, f, -2.0, 2.0)

    vstore3(f, i, force);


}
