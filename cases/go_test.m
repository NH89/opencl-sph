clear all;
addpath('../bin/mex/');  
%init_wrap('/../../conf/solidstretch.conf');

% Set some initial conditions
[y, x, z] = meshgrid(-1:0.4:1);

xyz = cat(4, x, y, z);

position = reshape(xyz, [numel(xyz)/3 3])';
n = uint32(size(position, 2));
originalposition = position;


plane_constraints_particles = zeros(1, n);

plane_constraints_particles(position(2,:) == -1) = 1;
plane_constraints_particles(position(2,:) == 1) = 2;

plane_constraints_particles = uint32(plane_constraints_particles);

r = [1 0 0; 0 cos(-0.02) -sin(-0.02); 0 sin(-0.02) cos(-0.02)];
position = r*position;


setstate('position', position, 'n', n, 'originalpos', originalposition,...
         'plane_constraints_particles', plane_constraints_particles);

computebins();
callkernel('compute_original_density');

querystate('n')

num_steps = 2000;

[position, density, n] = querystate('position', 'density', 'n');

figure;
for i=1:num_steps
    scatter3(position(1,1:n), position(2,1:n), position(3,1:n), 10, density(1:n));
    axis([-2 2 -2 2 -2 2]);
    view(0, 0);
    
    % Do stuff related to changing boundary conditions
    callkernel('apply_plane_constraints')
    
    computebins();
    
    callkernel('compute_density');
    callkernel('compute_rotations_and_strains');
    callkernel('compute_stresses');
    callkernel('compute_forces_solids');
    
    callkernel('step_forward');
    
    [position, density] = querystate('position', 'density');
    
    drawnow;
end