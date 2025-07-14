%% Part 1.4: Forward Euler

clear; clc; close all;

% Define ODE & intialize parameters
f = @(x, y) y;

% Define step sizes 
step_sizes = [0.1, 0.05, 0.001];

% Intialize counter
i = 1; 

% Loop through each step size
while i <= length(step_sizes)
    % Define step size and create a range of x values from 0 to 5
    h = step_sizes(i);
    x = 0:h:5;

    y = zeros(size(x));
    y(1) = 1; % Initial condition

    % Euler method calcualtion
    for n = 1:length(x)-1
        % Update the value of the function
        y(n+1) = y(n) + h * f(x(n), y(n));
    end

    % Plotting Figure
    figure(1);
    plot(x, y, 'LineWidth', 1, 'DisplayName', ['Euler: \Deltax = ' num2str(h)]);
    hold on;
    
    % Increment Counter
    i = i + 1; 

end

% Legend
legend('\Deltax = 0.1', '\Deltax = 0.05', '\Deltax = 0.001');
xlabel('x'); ylabel('y');
title('Forward Euler Approximations');
grid on;