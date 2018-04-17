#pragma OPENCL EXTENSION cl_khr_fp64 : enable

#define USE_FIELD(name, type) global type * name = (global type *) name##_m;
#define USE_FIELD_FIRST_VALUE(name, type) private type name = *((global type *) name##_m);

#define USE_GRID_PROPS \
    USE_FIELD(gridcell,      uint) USE_FIELD(gridcount, uint) USE_FIELD(celloffset, uint)\
    USE_FIELD(cellparticles, uint) USE_FIELD(gridres,   uint)

#define PSO_ARGS uint num_fields,\
                global char * names,         global uint * names_offsets,\
                global uint * dimensions,    global uint * num_dimensions, global uint * dimensions_offsets, global uint * entry_sizes,\
                global void * data,          global uint * data_sizes,     global uint * data_offsets

#define PI 3.1415926535

#define I_00 0
#define I_01 1
#define I_02 2
#define I_10 3
#define I_11 4
#define I_12 5
#define I_20 6
#define I_21 7
#define I_22 8

#define S_00 0
#define S_11 1
#define S_22 2
#define S_01 3
#define S_02 4
#define S_12 5
#define S_10 S_01
#define S_20 S_02
#define S_21 S_12

#define APPLY_CUBE_BOUNDS(position, force, minimum, maximum) \
    if (position.x < minimum+0.2) force.x += 1/(-(minimum)+position.x) - 5;\
    else if (position.x > maximum-0.2) force.x -= 1/(maximum-position.x) - 5;\
    if (position.y < minimum+0.2) force.y += 1/(-(minimum)+position.y) - 5;\
    else if (position.y > maximum-0.2) force.y -= 1/(maximum-position.y) - 5;\
    if (position.z < minimum+0.2) force.z += 1/(-(minimum)+position.z) - 5;\
    else if (position.z > maximum-0.2) force.z -= 1/(maximum-position.z) - 5;

// Can't use function pointers so silly macro hacks it is
// Needs gridres, gridcell, gridcount, celloffset, cellparticles from outer scope
// Whatever m_particlenum is changed to should already be defined, otherwise how do
// we know what particle you're asking about
#define FOR_PARTICLES_IN_RANGE(m_particlenum, m_otherparticlenum, m_block) {\
    int gx, gy, gz;\
\
    for (gz = -1; gz <= 1; ++gz)\
    for (gy = -1; gy <= 1; ++gy)\
    for (gx = -1; gx <= 1; ++gx)\
    {\
        int cell = get_cell_at_offset(gridres, gridcell[m_particlenum], gx, gy, gz);\
\
        if (cell == -1) continue;\
\
        uint offset = celloffset[cell];\
\
        for (uint jp = offset; jp < offset + gridcount[cell]; ++jp) {\
            uint m_otherparticlenum = cellparticles[jp];\
\
            m_block\
        }\
    }\
}

////////// Utility //////////

inline void zeroArray(global double * array, uint size) {
    for (uint i = 0; i < size; ++i) array[i] = 0;
}

////////// Kernel functions //////////

inline double applyKernel (double dist, double smoothingradius) {
    return dist > smoothingradius ? 0 : 315.0/64.0/PI/pow(smoothingradius, 9) * pow(pow(smoothingradius, 2) - pow(dist, 2), 3);
}
inline double3 applyKernelGradient (double3 dist, double smoothingradius) {
    double norm_dist = length(dist);
    double3 direction = dist/norm_dist;
    //return norm_dist > smoothingradius ? (double3)(0, 0, 0) : direction * -1890/64/PI/pow(smoothingradius, 9) * norm_dist * pow(pow(smoothingradius, 2) - pow(norm_dist, 2), 3);

    return norm_dist > smoothingradius ? (double3)(0, 0, 0) : direction * -45/PI/pow(smoothingradius, 6) * pow(smoothingradius - norm_dist, 2);
}
inline double applyLapKernel (double dist, double smoothingradius) {
    return dist > smoothingradius ? 0 : 45/PI/pow(smoothingradius, 6) * (smoothingradius - dist);
}

////////// Matrix functions //////////
// Many are for solids only but what harm does it do leaving them here

constant uint givens_i_xx_values[] = { 0, 0, 1 };
constant uint givens_i_yy_values[] = { 1, 2, 2 };
constant uint givens_i_ee_values[] = { 2, 1, 0 };
constant uint givens_i_ex_values[] = { 4, 3, 3 };
constant uint givens_i_ey_values[] = { 5, 5, 4 };

static inline void givensRotateToZero (double i[6], double pin[9], uint offdiag, double o[6], double pout[9]) {
    /* [ 0  3  4
         .  1  5
         .  .  2 ] */

    const uint ix = offdiag - 3;

/*
    const uint i_xx = xpy / 3;
    const uint i_yy = xpy - i_xx;
    const uint i_ee = 3-xpy;
    const uint i_xy = offdiag;
    const uint i_ex = 2+i_ee+i_xx;
    const uint i_ey = 2+i_ee+i_yy;
*/

    const uint i_xx = givens_i_xx_values[ix];
    const uint i_yy = givens_i_yy_values[ix];
    const uint i_ee = givens_i_ee_values[ix];
    const uint i_xy = offdiag;
    const uint i_ex = givens_i_ex_values[ix];
    const uint i_ey = givens_i_ey_values[ix];

    const double angle = 0.5*atan(2*i[i_xy]/(i[i_xx]-i[i_yy]));

    const double c = cos(angle);
    const double s = sin(angle);

    const double c2 = c*c;
    const double s2 = s*s;
    const double cs = c*s;

    o[i_xx] = i[i_xx]*c2 + i[i_yy]*s2 + 2*i[i_xy]*cs;
    o[i_xy] = -i[i_xx]*cs + i[i_yy]*cs + i[i_xy]*(c2-s2);
    o[i_yy] = i[i_xx]*s2 + i[i_yy]*c2 - 2*i[i_xy]*cs;
    o[i_ee] = i[i_ee];
    o[i_ex] = i[i_ex]*c + i[i_ey]*s;
    o[i_ey] = -i[i_ex]*s + i[i_ey]*c;

    const uint p_0x = i_xx, p_1x = 3+i_xx, p_2x = 6+i_xx,
               p_0y = i_yy, p_1y = 3+i_yy, p_2y = 6+i_yy,
               p_0e = i_ee, p_1e = 3+i_ee, p_2e = 6+i_ee;

    pout[p_0x] = pin[p_0x]*c + pin[p_0y]*s;
    pout[p_1x] = pin[p_1x]*c + pin[p_1y]*s;
    pout[p_2x] = pin[p_2x]*c + pin[p_2y]*s;

    pout[p_0y] = -pin[p_0x]*s + pin[p_0y]*c;
    pout[p_1y] = -pin[p_1x]*s + pin[p_1y]*c;
    pout[p_2y] = -pin[p_2x]*s + pin[p_2y]*c;

    pout[p_0e] = pin[p_0e];
    pout[p_1e] = pin[p_1e];
    pout[p_2e] = pin[p_2e];
}
inline uint maxAbsInArray (private double * array, uint size) {
    double m = 0;
    uint i_m = 0;
    for (uint i = 0; i < size; ++i) {
        double a = fabs(array[i]);
        if (a > m) {
            m = a;
            i_m = i;
        }
    }

    return i_m;
}
inline void multiplyBothSides (double d[3], double p[9], double out[6]) {
    // P*D*Pt
    const double d0p00 = d[0]*p[I_00], d1p01 = d[1]*p[I_01], d2p02 = d[2]*p[I_02],
                 d0p10 = d[0]*p[I_10], d1p11 = d[1]*p[I_11], d2p12 = d[2]*p[I_12],
                 d0p20 = d[0]*p[I_20], d1p21 = d[1]*p[I_21], d2p22 = d[2]*p[I_22];

    out[S_00] = p[I_00]*d0p00 + p[I_01]*d1p01 + p[I_02]*d2p02;
    out[S_01] = p[I_00]*d0p10 + p[I_01]*d1p11 + p[I_02]*d2p12;
    out[S_02] = p[I_00]*d0p20 + p[I_01]*d1p21 + p[I_02]*d2p22;

    //out[S_10] = p[I_10]*d0p00 + p[I_11]*d1p01 + p[I_12]*d2p02;
    out[S_11] = p[I_10]*d0p10 + p[I_11]*d1p11 + p[I_12]*d2p12;
    out[S_12] = p[I_10]*d0p20 + p[I_11]*d1p21 + p[I_12]*d2p22;

    //out[S_20] = p[I_20]*d0p00 + p[I_21]*d1p01 + p[I_22]*d2p02;
    //out[S_21] = p[I_20]*d0p10 + p[I_21]*d1p11 + p[I_22]*d2p12;
    out[S_22] = p[I_20]*d0p20 + p[I_21]*d1p21 + p[I_22]*d2p22;
}
inline void multiplyAS (double a[9], double s[6], global double out[9]) {
    out[I_00] = a[I_00]*s[S_00] + a[I_01]*s[S_10] + a[I_02]*s[S_20];
    out[I_01] = a[I_00]*s[S_01] + a[I_01]*s[S_11] + a[I_02]*s[S_21];
    out[I_02] = a[I_00]*s[S_02] + a[I_01]*s[S_12] + a[I_02]*s[S_22];

    out[I_10] = a[I_10]*s[S_00] + a[I_11]*s[S_10] + a[I_12]*s[S_20];
    out[I_11] = a[I_10]*s[S_01] + a[I_11]*s[S_11] + a[I_12]*s[S_21];
    out[I_12] = a[I_10]*s[S_02] + a[I_11]*s[S_12] + a[I_12]*s[S_22];

    out[I_20] = a[I_20]*s[S_00] + a[I_21]*s[S_10] + a[I_22]*s[S_20];
    out[I_21] = a[I_20]*s[S_01] + a[I_21]*s[S_11] + a[I_22]*s[S_21];
    out[I_22] = a[I_20]*s[S_02] + a[I_21]*s[S_12] + a[I_22]*s[S_22];
}
inline void multiplySymmetric (double a[6], double b[6], double out[6]) {
    out[S_00] = a[S_00]*b[S_00] + a[S_01]*b[S_10] + a[S_02]*b[S_20];
    out[S_01] = a[S_00]*b[S_01] + a[S_01]*b[S_11] + a[S_02]*b[S_21];
    out[S_02] = a[S_00]*b[S_02] + a[S_01]*b[S_12] + a[S_02]*b[S_22];
    out[S_11] = a[S_10]*b[S_01] + a[S_11]*b[S_11] + a[S_12]*b[S_21];
    out[S_12] = a[S_10]*b[S_02] + a[S_11]*b[S_12] + a[S_12]*b[S_22];
    out[S_22] = a[S_20]*b[S_02] + a[S_21]*b[S_12] + a[S_22]*b[S_22];
}
inline void multiplyMatrices (double a[9], double b[9], double out[9]) {
    out[I_00] = a[I_00]*b[I_00] + a[I_01]*b[I_10] + a[I_02]*b[I_20];
    out[I_01] = a[I_00]*b[I_01] + a[I_01]*b[I_11] + a[I_02]*b[I_21];
    out[I_02] = a[I_00]*b[I_02] + a[I_01]*b[I_12] + a[I_02]*b[I_22];

    out[I_10] = a[I_10]*b[I_00] + a[I_11]*b[I_10] + a[I_12]*b[I_20];
    out[I_11] = a[I_10]*b[I_01] + a[I_11]*b[I_11] + a[I_12]*b[I_21];
    out[I_12] = a[I_10]*b[I_02] + a[I_11]*b[I_12] + a[I_12]*b[I_22];

    out[I_20] = a[I_20]*b[I_00] + a[I_21]*b[I_10] + a[I_22]*b[I_20];
    out[I_21] = a[I_20]*b[I_01] + a[I_21]*b[I_11] + a[I_22]*b[I_21];
    out[I_22] = a[I_20]*b[I_02] + a[I_21]*b[I_12] + a[I_22]*b[I_22];
}
inline double3 multiplyMatrixVector (global double m[9], double3 v) {
    return (double3)(m[I_00]*v.x + m[I_01]*v.y + m[I_02]*v.z,
                     m[I_10]*v.x + m[I_11]*v.y + m[I_12]*v.z,
                     m[I_20]*v.x + m[I_21]*v.y + m[I_22]*v.z);
}
inline double3 multiplyMatrixVectorPrivate (double m[9], double3 v) {
    return (double3)(m[I_00]*v.x + m[I_01]*v.y + m[I_02]*v.z,
                     m[I_10]*v.x + m[I_11]*v.y + m[I_12]*v.z,
                     m[I_20]*v.x + m[I_21]*v.y + m[I_22]*v.z);
}
inline double3 multiplySymMatrixVector(global double s[6], double3 v) {
    return (double3)(s[S_00]*v.x + s[S_01]*v.y + s[S_02]*v.z,
                     s[S_10]*v.x + s[S_11]*v.y + s[S_12]*v.z,
                     s[S_20]*v.x + s[S_21]*v.y + s[S_22]*v.z);
}
inline void addOuterProduct (double3 a, double3 b, double out[9]) {
    out[I_00] += a.x*b.x; out[I_01] += a.x*b.y; out[I_02] += a.x*b.z;
    out[I_10] += a.y*b.x; out[I_11] += a.y*b.y; out[I_12] += a.y*b.z;
    out[I_20] += a.z*b.x; out[I_21] += a.z*b.z; out[I_22] += a.z*b.z;
}
inline void jacobi (private double i[6], double threshold, uint iterations, private double o[6], double p[9]) {
    uint particle = get_global_id(0);
    uint itix = particle*8;
    uint itnum = particle*8*3;

    uint _iterations = 0;

    p[I_00] = 1; p[I_01] = 0; p[I_02] = 0;
    p[I_10] = 0; p[I_11] = 1; p[I_12] = 0;
    p[I_20] = 0; p[I_21] = 0; p[I_22] = 1;

    double p2[9];

    for (uint ix = 0; ix < 9; ++ix) p2[ix] = p[ix];

    double temp[6];
    double temp2[6];

    private double * ioff = i + 3;
    double * tempoff = temp + 3;
    double * temp2off = temp2 + 3;

    for (uint ix = 0; ix < 6; ++ix) temp[ix] = i[ix];

    uint i_max = maxAbsInArray(tempoff, 3)+3;

    while ((_iterations < iterations) || (threshold > 0 && fabs(i[i_max]) > threshold)) {
        givensRotateToZero(temp, p, i_max, temp2, p2);

        i_max = maxAbsInArray(temp2off, 3)+3;

        givensRotateToZero(temp2, p2, i_max, temp, p);

        i_max = maxAbsInArray(tempoff, 3)+3;

        _iterations += 2;
    }

    for (uint ix = 0; ix < 6; ++ix) o[ix] = temp[ix];
}
inline void getR (double apq[9], global double r[9]) {
    double p[9];

    double s[6] = { apq[I_00]*apq[I_00] + apq[I_10]*apq[I_10] + apq[I_20]*apq[I_20],
                    apq[I_01]*apq[I_01] + apq[I_11]*apq[I_11] + apq[I_21]*apq[I_21],
                    apq[I_02]*apq[I_02] + apq[I_12]*apq[I_12] + apq[I_22]*apq[I_22],
                    apq[I_00]*apq[I_01] + apq[I_10]*apq[I_11] + apq[I_20]*apq[I_21],
                    apq[I_00]*apq[I_02] + apq[I_10]*apq[I_12] + apq[I_20]*apq[I_22],
                    apq[I_01]*apq[I_02] + apq[I_11]*apq[I_12] + apq[I_21]*apq[I_22] };

    double eigenvalues[6];
    double eigenvectors[9];

    jacobi(s, 0, 8, eigenvalues, eigenvectors);

    for (uint i = 0; i < 3; ++i) eigenvalues[i] = 1/sqrt(eigenvalues[i]);

    double temp[6];

    multiplyBothSides(eigenvalues, eigenvectors, temp);

    multiplyAS(apq, temp, r);
}
inline void transpose (global double m[9], double m_t[9]) {
    m_t[I_00] = m[I_00];
    m_t[I_11] = m[I_11];
    m_t[I_22] = m[I_22];

    m_t[I_01] = m[I_10];
    m_t[I_10] = m[I_01];

    m_t[I_02] = m[I_20];
    m_t[I_20] = m[I_02];

    m_t[I_12] = m[I_21];
    m_t[I_21] = m[I_12];
}

////////// Particle <-> grid mapping //////////

inline int get_cell_at_offset (global uint * gridres, uint original, int x, int y, int z) {
    uint rowlength = gridres[0];
    uint layersize = gridres[0] * gridres[1];

    uint zo = original / layersize;
    uint yo = (original - layersize * zo) / rowlength;
    uint xo = original - layersize * zo - rowlength * yo;

    int zn = zo + z;
    int yn = yo + y;
    int xn = xo + x;

    if (xn >= gridres[0] || xn < 0 ||
        yn >= gridres[1] || yn < 0 ||
        zn >= gridres[2] || zn < 0) return -1;

    return zn * layersize + yn * rowlength + xn;
}
kernel void find_particle_bins (PSO_ARGS) {
    USE_FIELD(position, double) USE_FIELD(gridbounds, double)
    USE_FIELD(gridres, uint)    USE_FIELD(gridcell, uint) USE_FIELD(gridcount, uint)
    USE_FIELD_FIRST_VALUE(pnum, uint) USE_FIELD_FIRST_VALUE(n, uint)
    USE_FIELD_FIRST_VALUE(smoothingradius, double)

    size_t i = get_global_id(0);

    if (i >= n) return;

    double3 pos = (double3) (position[i*3], position[i*3 + 1], position[i*3 + 2]);

    uint3 gc = clamp((uint3) ((pos.x - gridbounds[0]) / (gridbounds[1] - gridbounds[0]) * gridres[0],
                              (pos.y - gridbounds[2]) / (gridbounds[3] - gridbounds[2]) * gridres[1],
                              (pos.z - gridbounds[4]) / (gridbounds[5] - gridbounds[4]) * gridres[2]),
                     (uint3) (0, 0, 0),
                     (uint3) (gridres[0]-1, gridres[1]-1, gridres[2]-1));

    gridcell[i] = gc.z * gridres[0] * gridres[1] + gc.y * gridres[0] + gc.x;
}
kernel void zero_gridcount (PSO_ARGS) {
    USE_FIELD(gridcount, uint) USE_FIELD(gridres, uint)

    size_t i = get_global_id(0);

    if (i >= gridres[0]*gridres[1]*gridres[2]) return;

    gridcount[i] = 0;
}
kernel void count_particles_in_bins (PSO_ARGS) {
    USE_FIELD(gridcell, uint) USE_FIELD(gridcount, uint) USE_FIELD_FIRST_VALUE(n, uint)
    size_t i = get_global_id(0);

    if (i >= n) return;

    atomic_inc (&gridcount[gridcell[i]]);
}
kernel void prefix_sum (PSO_ARGS, local uint * temp, uint array_size, global uint * block_totals, uint num_blocks) {
    USE_FIELD(gridcount, uint) USE_FIELD(celloffset, uint) USE_FIELD(gridres, uint)

    size_t i_global = get_global_id(0);
    size_t i_local = get_local_id(0);

    /* Find bloody block we're in, binary search */
    /*
    if (num_blocks == 1) {
        s = 0;
    } else {
        s = num_blocks / 2;

        uint min_index = 0;
        uint max_index = num_blocks;

        int repeat = 1;
        while (repeat) {
            if (i_global < block_indices[s]) {
                max_index = s;
                s = (s + min_index) / 2;
            } else if (i_global >= block_indices[s+1]) {
                min_index = s+1;
                s = (s + max_index) / 2;
            } else {
                repeat = 0;
            }

            if (s == num_blocks - 1 || s == 0) repeat = 0;
        }
    }
    */

    size_t i = get_local_id(0);
    size_t n = get_local_size(0)*2;
    size_t g_id = get_group_id(0);
    size_t global_offset = get_group_id(0) * get_local_size(0) * 2;

    int in_array_1 = i + global_offset < array_size;
    int in_array_2 = i + n/2 + global_offset < array_size;

    /* Copy array to temp memory */
    if (in_array_1) temp[i] = gridcount[i + global_offset];
    else temp[i] = 0;
    if (in_array_2) temp[i + n/2] = gridcount[i + n/2 + global_offset];
    else temp[i + n/2] = 0;

    uint last_element_count;
    if (i == n/2 - 1) {
        last_element_count = temp[i + n/2];
    }

    /* Begin black magic */
    /* i.e. do the prefix sum on temp */
    uint stride = 1;
    for (uint d = n>>1; d > 0; d >>= 1) {
        barrier(CLK_LOCAL_MEM_FENCE);

        if (i < d) {
            uint ai = stride*(2*i+1) - 1;
            uint bi = ai + stride;

            temp[bi] += temp[ai];
        }

        stride <<= 1;
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    if (i == 0) temp[n - 1] = 0;

    for (uint d = 1; d < n; d <<= 1) {
        barrier(CLK_LOCAL_MEM_FENCE);

        stride >>= 1;

        if (i < d) {
            uint ai = stride*(2*i+1) - 1;
            uint bi = ai + stride;

            uint t = temp[bi];
            temp[bi] += temp[ai];
            temp[ai] = t;
        }
    }
    /* End black magic */

    barrier(CLK_LOCAL_MEM_FENCE);

    if (in_array_1) celloffset[i + global_offset] = temp[i];
    if (in_array_2) celloffset[i + n/2 + global_offset] = temp[i + n/2];

    if (i == n/2 - 1) {
        block_totals[g_id] = temp[i + n/2] + last_element_count;
        uint g_id_c;
        for (g_id_c = g_id+1; g_id_c < num_blocks; ++g_id_c) {
            block_totals[g_id_c] += block_totals[g_id];
        }
    }

    barrier(CLK_GLOBAL_MEM_FENCE);

    if (g_id > 0) {
        if (in_array_1) celloffset[i + global_offset] += block_totals[g_id-1];
        if (in_array_2) celloffset[i + n/2 + global_offset] += block_totals[g_id-1];
    }
}
kernel void copy_celloffset_to_backup (PSO_ARGS, global uint * backup_prefix_sum, uint num_grid_cells) {
    USE_FIELD(celloffset, uint)

    size_t i = get_global_id(0);

    if (i < num_grid_cells) backup_prefix_sum[i] = celloffset[i];
}
kernel void insert_particles_in_bin_array (PSO_ARGS, global uint * backup_prefix_sum) {
    USE_FIELD_FIRST_VALUE(n, uint) USE_FIELD(gridcell, uint) USE_FIELD(cellparticles, uint)

    size_t i = get_global_id(0);

    if (i < n) cellparticles[atomic_inc(backup_prefix_sum + gridcell[i])] = i;
}

////////// Particle creation & transformation //////////

kernel void populate_position_cuboid (PSO_ARGS, double3 corner1, double3 corner2, uint3 size) {
    if (get_work_dim() < 3) return;

    size_t i = get_global_id(0);
    size_t j = get_global_id(1);
    size_t k = get_global_id(2);

    if (i >= size.x || j >= size.y || k >= size.z) return;

    USE_FIELD(position, double) USE_FIELD_FIRST_VALUE(pnum, uint) USE_FIELD(n, uint)

    uint p_ix = i + size.x*j + size.x*size.y*k;

    if (p_ix >= pnum) return;
    else if (p_ix == pnum - 1) n[0] = pnum;
    else if (i == size.x - 1 && j == size.y - 1 && k == size.z - 1) n[0] = size.x*size.y*size.z;

    double3 p = (double3) (((double) i/ (double) (size.x-1)) * (corner2.x - corner1.x) + corner1.x,
                           ((double) j/ (double) (size.y-1)) * (corner2.y - corner1.y) + corner1.y,
                           ((double) k/ (double) (size.z-1)) * (corner2.z - corner1.z) + corner1.z);

    vstore3(p, p_ix, position);
}
kernel void rotate_particles (PSO_ARGS, double angle_x, double angle_y, double angle_z) {
    // ZYX
    USE_FIELD(position, double) USE_FIELD_FIRST_VALUE(n, uint)

    uint i = get_global_id(0);

    if (i >= n) return;

    double rotation_x[9];
    double rotation_y[9];
    double rotation_z[9];

    double c_x = cos(angle_x);
    double s_x = sin(angle_x);

    double c_y = cos(angle_y);
    double s_y = sin(angle_y);

    double c_z = cos(angle_z);
    double s_z = sin(angle_z);

    rotation_x[I_00] = 1; rotation_x[I_01] = 0;   rotation_x[I_02] = 0;
    rotation_x[I_10] = 0; rotation_x[I_11] = c_x; rotation_x[I_12] = -s_x;
    rotation_x[I_20] = 0; rotation_x[I_21] = s_x; rotation_x[I_22] = c_x;

    rotation_y[I_00] = c_y; rotation_y[I_01] = 0; rotation_y[I_02] = -s_y;
    rotation_y[I_10] = 0;   rotation_y[I_11] = 1;   rotation_y[I_12] = 0;
    rotation_y[I_20] = s_y; rotation_y[I_21] = 0; rotation_y[I_22] = c_y;

    rotation_z[I_00] = c_z; rotation_z[I_01] = -s_z; rotation_z[I_02] = 0;
    rotation_z[I_10] = s_z; rotation_z[I_11] = c_z;  rotation_z[I_12] = 0;
    rotation_z[I_20] = 0;   rotation_z[I_21] = 0;    rotation_z[I_22] = 1;

    double rotation_yz[9];
    double rotation[9];

    multiplyMatrices(rotation_y, rotation_z, rotation_yz);
    multiplyMatrices(rotation_x, rotation_yz, rotation);

    double3 pos = vload3(i, position);

    double3 rpos = (double3)(
        rotation[I_00]*pos.x + rotation[I_01]*pos.y + rotation[I_02]*pos.z,
        rotation[I_10]*pos.x + rotation[I_11]*pos.y + rotation[I_12]*pos.z,
        rotation[I_20]*pos.x + rotation[I_21]*pos.y + rotation[I_22]*pos.z
    );

    vstore3(rpos, i, position);
}

////////// Kernels //////////

kernel void compute_density (PSO_ARGS) {
    USE_GRID_PROPS

    USE_FIELD_FIRST_VALUE(pnum, uint) USE_FIELD_FIRST_VALUE(n, uint) USE_FIELD_FIRST_VALUE(mass, double)
    USE_FIELD(position, double) USE_FIELD(density, double) USE_FIELD_FIRST_VALUE(smoothingradius, double)

    unsigned int i = get_global_id(0);

    if (i >= n) return;

    double3 ipos = vload3(i, position);

    density[i] = 0;

    FOR_PARTICLES_IN_RANGE(i, j,
        j = cellparticles[jp];

        double3 jpos = vload3(j, position);

        double3 diff = jpos - ipos;

        double dist = length(diff);

        density[i] += applyKernel(dist, smoothingradius);
    )

    density[i] *= mass;
}
kernel void step_forward (PSO_ARGS) {
    USE_FIELD_FIRST_VALUE(mass, double) USE_FIELD_FIRST_VALUE(timestep, double)
    USE_FIELD_FIRST_VALUE(n, uint)

    USE_FIELD(acceleration, double) USE_FIELD(force, double) USE_FIELD(position, double)
    USE_FIELD(velocity, double) USE_FIELD(veleval, double) USE_FIELD(posnext, double)
    USE_FIELD(gravity, double)

    uint i = get_global_id(0);

    if (i >= n) return;

    double3 g = vload3(0, gravity);

    double3 f = vload3(i, force);
    double3 a = f / mass + g;

    vstore3(a, i, acceleration);

    double3 v = vload3(i, velocity);
    double3 vnext = v + timestep * a;
    double3 veval = (v + vnext) / 2;

    double3 p = vload3(i, position);

    double3 pnext = p + timestep * vnext;

    vstore3(pnext, i, position);
    vstore3(vnext, i, velocity);
}
