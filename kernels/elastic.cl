#include "common.cl"

kernel void compute_forces_elastic (PSO_ARGS) {
    
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
    
    
    //  for master links
    for (      ){
        // get current distance
        
        
        // get type -> yield_len, rest_len, modulus, melt_pt,  etc 
        // get temp 
        
        
        if(cur_dist > yield){
         // release the link, nb remove from both particles 
            
        }else{
            // compute force
            force = cur_dist * modulus / rest_len;
            
            // add to both particles
            master_particle_forces += force ;//locally 
            
            follower_particle_forces -= force ;//atomic
            
        }
        // NB also consider new links: freezing/plasticity
        // If plastic/ductile/maleable will replace broken links with new local links
        // If melting break all links
        // If freezing form new local links
        // If locomoting automaton, select where to make/break links 
        //  &/or change rest_len & modulus (muscle)
        // If remodelling .. respond to forces & chem gradients
        
        
    }
    
     master_particle_forces;  // write to this particle
    
}


kernel void full_copy_elastic (PSO_ARGS) {
    
    // do general copy necessary values
    // NB could split data into swapped and non-swapped arrays.
    
    
    // for master links
    for (   ){
        // get new IDs from old
        
        
        // chk follower records link & not over connected.
        
    }
    

    
}

// Need buffer for "compute_forces_elastic(...)"

uchar master[8]; // list of IDs of master particles writing link forces to this one.
uchar follower[8];  // list of IDs of follower particles written to by this one.
uchar follower_idx[8]; // where in the follower's list of links.
uchar bond_type[8];  // ID of bond type in pallette 

struct bond_struct{float rest_len, modulus, yield_strain, melt_pt;}
const uchar pallette_size; // number of types of bonds in the simulation. 
bond_struct bond[pallette_size];


// Need a conf file


// Need a simlation dataset file + input/output.
    
    

