clear all;

num_steps = 1;
addpath(genpath('..'));
testinit_wrap();

[position, density, n] = querystate('position', 'density', 'n');

figure;
for i = 1:num_steps
    scatter3(position(1,1:n), position(2,1:n), position(3,1:n), 10, density(1:n));
    axis([-2 2 -2 2 -2 2]);

    simstep();
    [position, density] = querystate('position', 'density');
    
    drawnow;
end
