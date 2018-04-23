close all;
clear all;
clc;


addpath(genpath('..'));

%video = VideoWriter('fluid_3Dview.avi','Motion JPEG 2000');
%open(video)
path = '/media/aslab/data/hackthon_data/solid/positions/';
fig1 = figure('units','normalized','outerposition',[0 0 0.5, 1]);
pause on;

for i =1:3000
    %pause(0.01);
    M = csvread(strcat(path, 'position_',num2str(i-1),'.csv'));
    noOfParticles = sum(M(:,1)~=0 & M(:,2)~=0 & M(:,3)~=0);
    M = M(1:noOfParticles,1:3);
    size_m = size(M,1);
    size_v = ones(size_m,1)*100;
    depth = sqrt((-2 - M(:,1)).^2 +  (-2 - M(:,2)).^2 );
    scatter3(M(:,1), M(:,2), M(:,3),size_v, depth, 'filled');
    axis([-2 2 -2 2 -2 2]);
    %view(0,0);
    view(-114,2);
    grid off;
    set(gca,'XTickLabel',[])
    set(gca,'YTickLabel',[])
    set(gca,'ZTickLabel',[])
    drawnow;
    image = getframe(gcf);
    %writeVideo(video,image.cdata);  
    
    %im = imresize(image.cdata,[960,960]);
    %imwrite(im, strcat('/media/aslab/data/hackthon_data/solid/images/img',num2str(i,'%04.f'),'.png'));
end

%close(video);
