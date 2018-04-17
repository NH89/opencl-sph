close all;
clear all;
clc;


addpath(genpath('..'));

path = '../positions/';

pause on;

for i =1:2000
    pause(0.01);
    M = csvread(strcat(path, 'position_',num2str(i),'.csv'));
    noOfParticles = sum(M(:,1)~=0 & M(:,2)~=0 & M(:,3)~=0);
    M = M(1:noOfParticles,1:3);
    scatter3(M(:,1), M(:,2), M(:,3), 'filled');
    axis([-2 2 -2 2 -2 2]);
    view(43,1)

    drawnow;
end