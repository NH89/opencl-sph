#!/bin/bash

FILE=$1

avg()
{
	if [ -z "$1" ]; then
		n=1
	else
		n=$1
	fi
	awk "function isnum(x){return(x==x+0)} { if(isnum(\$$n)) { sum+=\$$n; sumsq+=\$$n*\$$n ; n+=1;} } END { print sum/n, sum, n, sqrt(sumsq/n - sum*sum/n/n) }"
}

sum() {
	grep $1 $2 | cut -f 2 | avg
}

FUNCTIONS="benchmark-populate_position_cuboid_device_opencl
benchmark-init_original_position
benchmark-rotate_particles_device_opencl
benchmark-zero_gridcount
benchmark-bin_and_count_device_opencl
benchmark-prefix_sum_device_opencl
benchmark-copy_celloffset_to_backup_device_opencl
benchmark-insert_particles_in_bin_array_device_opencl
benchmark-compute_original_density
benchmark-compute_density
benchmark-compute_rotations_and_strains
benchmark-compute_stresses
benchmark-compute_forces_solids
"

# print table of:
# $function, avgTime, sumTime, count, stdevTime

for f in $FUNCTIONS ; do
	echo $f | tr '\n' '\t'
	sum $f $FILE | tr ' ' '\t'
done
