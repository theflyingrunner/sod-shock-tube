% plot_sod_shock_tube.m
% MATLAB script to plot the solution output by the C program to solve
% the Sod Shock Tube problem using Steger-Warming scheme and MUSCL
% formalism

clc;
clear;
close all;

% Read solution file
data = readmatrix('solution_muscl.txt');
data_old = readmatrix('solution.txt');

% Extract parameters from data
x = data(:, 1);
p = data(:, 2);
rho = data(:, 3);
T = data(:, 4);
u = data(:, 5);

x_old = data_old(:, 1);
p_old = data_old(:, 2);
rho_old = data_old(:, 3);
T_old = data_old(:, 4);
u_old = data_old(:, 5);

% Plot
f = figure(1);
sgtitle({'Sod Shock Tube: Steger-Warming Scheme w/ and w/o MUSCL'});
subplot(2, 2, 1);
plot(x, p);
hold on;
plot(x_old, p_old);
title('Pressure');
xlabel('Location');
ylabel('p');
legend('MUSCL', 'No MUSCL', 'location', 'sw');
hold off;

subplot(2, 2, 2);
plot(x, T);
hold on;
plot(x_old, T_old);
title('Temperature');
xlabel('Location');
ylabel('T');
legend('MUSCL', 'No MUSCL', 'location', 'sw');
hold off;

subplot(2, 2, 3);
plot(x, rho);
hold on;
plot(x_old, rho_old);
title('Density');
xlabel('Location');
ylabel('\rho');
legend('MUSCL', 'No MUSCL', 'location', 'sw');
hold off;

subplot(2, 2, 4);
plot(x, u);
hold on;
plot(x_old, u_old);
title('Velocity');
xlabel('Location');
ylabel('u');
legend('MUSCL', 'No MUSCL', 'location', 'sw');
hold off;

f.Position = [100 100 1080 800];