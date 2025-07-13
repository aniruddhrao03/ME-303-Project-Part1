%% Part 1.3: Plot solution from part 1 (y = e^x) and part 2 (power series)

clear; clc; close all;

% Define ODE & intialize parameters
f = @(x, y) y;

% Part 1.3: Plot solution from part 1 (y = e^x)
x = 0:0.1:5; % plotting 50 points between x = 0 and x = 5 used for both

% Define the exponential function
y_exact = exp(x);

% Plot 
figure(1);
plot(x, y_exact);
title('Plot of y = e^x');
xlabel('x');
ylabel('y');
grid on; 

% Part 1.3: Plot solution from part 2 (Power Series)
N_values = [1, 3, 5, 10, 100]; % Different numbers of terms in the series

figure(2);
hold on; % Hold on to plot multiple series on the same figure

% Intialize counter
i = 1; 

while i <= length(N_values)
    N = N_values(i); % Get the current number of terms
    a = zeros(1, N); % defining an array to hold all coefficients

    % Power series approximation
    y_series = zeros(size(x));
    for n = 1:N-1
        y_series = y_series + 1 / factorial(n) * (x.^(n)); % Sum the series
    end

    % Plot the power series approximation
    plot(x, y_series, 'DisplayName', ['N = ', num2str(N)]);

    % Increment counter
    i = i + 1;
end

title('Power Series Approximation of ODE Solution');
xlabel('x');
ylabel('y');
grid on;
legend show; % Show legend
hold off; % Release the hold on the figure