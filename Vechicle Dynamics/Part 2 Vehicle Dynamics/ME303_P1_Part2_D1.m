clear; clc; close all;

% Default Vehicle Parameters
m = 1400;            % Gross weight, kg
a = 1.14;            % cg to front axle, m
b = 1.33;            % cg to rear axle, m
Cf = 25000;          % front tyre cornering stiffness

% Constants
Iz = 2420;           % yaw intertia, kg*m^2
g = 9.81;            % acceleration due to gravity, m/s^2

% Cr values from 21000 to 21800
Cr_values = 21000:1:21800; % Increment by 1

% Preallocate arrays for results
v_critical = zeros(size(Cr_values));
Kus_values = zeros(size(Cr_values));

% Calculate L, Wf, Wr
L = a + b; % Wheelbase, m
Wf = (b / L) * m * g; % Weight on front axle, N
Wr = (a / L) * m * g; % Weight on rear axle, N

% Iterate through Cr values
for i = 1:length(Cr_values)
    Cr = Cr_values(i);

    % Calculate understeer coefficient
    Kus = Wf / Cf - Wr / Cr; % Understeer coefficient
    Kus_values(i) = Kus; % Store Kus for later use

    % Calculate highest stable speed only for negative Kus (oversteer)
    %if Kus < 0
        v_critical(i) = sqrt((g * L) / -Kus) * 3.6; % Convert to km/h
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
plot(Cr_values, v_critical);
xlabel('Cornering Stiffness (Cr) [N/rad]');
ylabel('Highest Stable Speed (v_{critical}) [km/h]');
title('Highest Stable Speed vs Rear Cornering Stiffness');
grid on;