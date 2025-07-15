clear; clc; close all;

% Default Vehicle Parameters
%m = 1400;            % Gross weight, kg
a = 1.14;            % cg to front axle, m
b = 1.33;            % cg to rear axle, m
Cf = 25000;          % front tyre cornering stiffness
Cr = 21000;          % front tyre cornering stiffness

% Constants
Iz = 2420;           % yaw inertia, kg*m^2
g = 9.81;            % acceleration due to gravity, m/s^2

% m values from 500 to 2500
m_values = 500:1:2500; % Increment by 1

% Preallocate arrays for results
v_critical = zeros(size(m_values));
Kus_values = zeros(size(m_values));


% Iterate through m values
for i = 1:length(m_values)
    m = m_values(i);

    % Calculate L, Wf, Wr
    L = a + b; % Wheelbase, m
    Wf = (b / L) * m * g; % Weight on front axle, N
    Wr = (a / L) * m * g; % Weight on rear axle, N

    % Calculate understeer coefficient
    Kus = Wf / Cf - Wr / Cr; % Understeer coefficient
    Kus_values(i) = Kus; % Store Kus for later use

    % Calculate highest stable speed only for negative Kus (oversteer)
    %if Kus < 0
        v_critical(i) = sqrt((g * L) / abs(Kus)) * 3.6; % Convert to km/h
    %end

    % Display results based on Kus
    if Kus == 0
        disp('Neutral Steer');
    elseif Kus > 0
        fprintf('Understeer, infinite stability, characteristic speed = %.2f km/h\n', sqrt((g * L) / Kus) * 3.6);
    elseif Kus < 0
        fprintf('Oversteer, highest stable speed = %.2f km/h\n', v_critical(i));
    end
end

% Plot v_critical against Cr values
figure;
plot(m_values, v_critical);
xlabel('Mass (m) [kg]');
ylabel('Highest Stable Speed (v_{critical}) [km/h]');
title('Highest Stable Speed vs Mass');
grid on;