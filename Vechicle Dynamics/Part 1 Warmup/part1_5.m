%% Part 1.5: Verifying answers

clear; clc; close all;
% Check using ode45 function (Given from Assignment)
odefun = @(x, y) y;

% Define the range of x values 
xspan = [0 5];

% Initial Condition
y0 = 1;

% Solve the ODE using ode45
[x, y] = ode45(odefun, xspan, y0);

%Plotting
figure(1);
plot(x, y);
xlabel('x');
ylabel('y');
title('Solution of the ODE using ode45');
grid on;
