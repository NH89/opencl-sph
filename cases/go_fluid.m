clear all;

num_steps = 2000;
addpath(genpath('..'));
testinit_fluid_wrap();

[position, density, n] = querystate('position', 'density', 'n');

figure;
for i = 1:num_steps
    scatter3(position(1,1:n), position(2,1:n), position(3,1:n), 10, density(1:n));
    axis([-2.0 2.0 -2.0 2.0 -2.0 2.0]);

    simstep_fluid();
    [position, density] = querystate('position', 'density');
    
    drawnow;
end
