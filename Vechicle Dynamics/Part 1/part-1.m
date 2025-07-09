clear; clc; close all;

% Define ODE & intialize parameters
f = @(x, y) y;

% Part 1.3: Plot solution from part 1 (y = e^x)
    x = linspace(0, 5, 100); % plotting 100 points between x = 0 and x = 5    
    y_exact = exp(x);
    
    figure(1);
    plot(x, y_exact);
    title('Plot of y = e^x');
    xlabel('x');
    ylabel('y');
    grid on; % Adding a grid for clarity

% Part 1.3: Plot solution from part 2 (Power Series)
    N = 100; % Number of terms in the series
    a = zeros(1, N); % defining an array to hold all coefficients

    % Calculate the coefficients for the power series
    for n = 1:N-1
        a(n) = 1 / factorial(n); % Equation for coefficients, a(n), of power series
    end

    % Define the range for x
    x = linspace(0, 5, 100); % plotting 100 points between x=0 and 5

    % Calculate the power series approximation
    y_series = zeros(size(x));
    for n = 1:N-1
        y_series = y_series + a(n) * (x.^(n)); % Sum the series
    end
 
    figure(2);
    plot(x, y_series);
    title('Power Series Approximation of ODE Solution');
    xlabel('x');
    ylabel('y(x)');
    grid on;

% Part 1.4: Forward Euler
    % Define step sizes and colors
    step_sizes = [0.1, 0.05, 0.001];
    colours = ['r', 'g', 'b'];

    i = 1; 

    % Loop through each step size
    while i <= length(step_sizes)

        h = step_sizes(i);
        x = 0:h:5;
        y = zeros(size(x));

        y(1) = 1; % Initial condition

        for n = 1:length(x)-1
            y(n+1) = y(n) + h * f(x(n), y(n));
        end

        figure(3);
        plot(x, y, 'Color', colours(i), 'LineWidth', 1.5, 'DisplayName', ['Euler: \Deltax = ' num2str(h)]);
        hold on;
    
        i= i + 1; 

    end

    legend('\Deltax = 0.1', '\Deltax = 0.05', '\Deltax = 0.001');
    xlabel('x'); ylabel('y');
    title('Forward Euler Approximations for dy/dx = y');
    grid on;

%% Check using ode45 function (Given from Assignment)
    odefun = @(x, y) y;
    xspan = [0 5];
    [x, y] = ode45(odefun, xspan, y0);
    
    figure(4)
    plot(x, y);
    xlabel('x');
    ylabel('y');
    title('Solution of the ODE using ode45');
    grid on;

 
