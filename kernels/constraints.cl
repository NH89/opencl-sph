// constraints.cl ///////////////////////////////////////////////////////////////////////////////////////
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

kernel void apply_plane_constraints(PSO_ARGS) {
    USE_FIELD(position, REAL) 
    USE_FIELD(velocity, REAL)
    USE_FIELD(plane_constraints, REAL) 
    USE_FIELD(plane_constraints_particles, uint)

    USE_FIELD_FIRST_VALUE(n, uint)

    uint i = get_global_id(0);

    if (i >= n || plane_constraints_particles[i] == 0) return;

    global REAL * plane_constraint = plane_constraints + (plane_constraints_particles[i]-1)*6;

    REAL3 point = vload3(0, plane_constraint);
    REAL3 normal = vload3(1, plane_constraint);

    REAL3 ipos = vload3(i, position);
    REAL3 ivel = vload3(i, velocity);

    ipos = ipos - dot(normal, ipos - point) * normal;
    ivel = ivel - dot(normal, ivel) * normal;

    vstore3(ipos, i, position);
    vstore3(ivel, i, velocity);
}
