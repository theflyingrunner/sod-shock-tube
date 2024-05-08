% plot_sod_shock_tube.m
% MATLAB script to plot the solution output by the C program to solve
% the Sod Shock Tube problem using Steger-Warming scheme and MUSCL
% formalism

clc;
clear;
close all;

% Read solution files
pressureData = readmatrix('pressure.txt');
densityData = readmatrix('density.txt');
temperatureData = readmatrix('temperature.txt');
velocityData = readmatrix('velocity.txt');

% constants
nt = 12000;
nx = 1001;
dx = 10.0 / (nx - 1);
x = linspace(-5, 5, nx - 1);

% Plot
M(nt) = struct('cdata', [], 'colormap', []);

f = figure(1);
f.Visible = 'off';
myVideo = VideoWriter('solutionAnimationt12000'); %open video file
myVideo.FrameRate = 240;  %can adjust this, 5 - 10 works well for me
open(myVideo);

for i = 1:nt
    subplot(2, 2, 1);
    sgtitle('Sod Shock Tube w/ Steger-Warming and MUSCL');
    plot(x(1, :), pressureData(i, 1:1000), '-');
    hold on;
    xlabel('Location');
    ylabel('p');
    title('Pressure');
    axis([-5.0 5.0 0.0 1.2]);
    
    subplot(2, 2, 2);
    plot(x(1, :), densityData(i, 1:1000), '-');
    xlabel('Location');
    ylabel('\rho');
    title('Density');
    axis([-5.0 5.0 0.0 1.2]);
    
    subplot(2, 2, 3);
    plot(x(1, :), temperatureData(i, 1:1000), '-');
    xlabel('Location');
    ylabel('T');
    title('Temperature');
    axis([-5.0 5.0 0.7 1.2]);
    
    subplot(2, 2, 4);
    plot(x(1, :), velocityData(i, 1:1000), '-');
    hold on;
    xlabel('Location');
    ylabel('u');
    title('Velocity');
    axis([-5.0 5.0 0.0 1.0]);
    hold off;
    
    M(i) = getframe(f);
    writeVideo(myVideo, M(i));
    clf;
    %print(i);
end

f.Visible = 'on';
close(myVideo);

movie(M, 1, 240);


