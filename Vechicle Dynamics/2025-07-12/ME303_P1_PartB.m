clear; clc; close all;

% Vehicle Parameters
m = 1400;
a = 1.14;
b = 1.33;
Cf = 25000;
Cr = 21000;
Iz = 2420;
delta = 0.1;

% Time setup
h = 0.01;
t = 0:h:20;
N = length(t);

f1 = figure;
f2 = figure;
f3 = figure;
f4 = figure;

One_g = -9.81* ones(size(t));

for u = [25, 50, 75, 100, 200, 300]* 1000 / 3600
    
    A_1 = [-(Cf + Cr)/(m*u), (-a*Cf + b*Cr)/(m*u) - u; -(a*Cf - b*Cr)/(Iz*u), -(a^2*Cf + b^2*Cr)/(Iz*u)];
    B_1 = [Cf/m; a*Cf/Iz];

    % Initialize vectors: [v_y; yaw rate]
    x_euler = zeros(2, N);
    x_rk4   = zeros(2, N);
    
    % Initialize vector: a_y (lateral acceleration)
    a_y_euler = zeros(1, N);
    a_y_rk4 = zeros(1, N);

    %Euler Method
    for n = 1:N-1
        x_euler(:,n+1) = x_euler(:,n) + h * (A_1 * x_euler(:,n) + B_1 * delta);
        a_y_euler (1, :) = (-(Cf + Cr)/(m*u))*x_euler(1,:) + ((-a*Cf + b*Cr)/(m*u) - u)*x_euler(2,:) + (Cf/m)*delta;
    end

    %RK4 Method 
    for n = 1:N-1
        k1 = A_1 * x_rk4(:,n) + B_1 * delta;
        k2 = A_1 * (x_rk4(:,n) + 0.5*h*k1) + B_1 * delta;
        k3 = A_1 * (x_rk4(:,n) + 0.5*h*k2) + B_1 * delta;
        k4 = A_1 * (x_rk4(:,n) + h*k3) + B_1 * delta;
        x_rk4(:,n+1) = x_rk4(:,n) + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
        a_y_rk4 (1, :) = (-(Cf + Cr)/(m*u))*x_rk4(1,:) + ((-a*Cf + b*Cr)/(m*u) - u)*x_rk4(2,:) + (Cf/m)*delta;
    end 

    %Plot of Lateral Accleration for Euler
    figure(f1);
    plot(t, a_y_euler(1,:), 'DisplayName', ['Euler a_y, u = ' num2str(u*3600/1000) 'km/h']); %Euler Plot
    hold on;
   
    %Plot of Lateral Velocity for RK4
    figure(f2);
    plot(t, a_y_rk4(1,:), 'DisplayName', ['RK4 a_y, u = ' num2str(u*3600/1000) 'km/h']); %RK4 Plot
    hold on;

    %Plot of Yaw Angle for Euler
    figure(f3);
    plot(t, x_euler(2,:), 'DisplayName', ['Euler yaw rate, u = '  num2str(u*3600/1000) 'km/h']);
    hold on;

    %Plot of Yaw Angle for RK4
    figure(f4);
    plot(t, x_rk4(2,:),'DisplayName', ['RK4 yaw rate, u = '  num2str(u*3600/1000) 'km/h']);
    hold on;
end

    figure(f1);
    plot(t,One_g, 'DisplayName', '1g');
    hold on;
    xlabel('Time (s)');
    ylabel('Lateral Acceleration (m/s^2)');
    title('Euler Lateral Accleration');
    legend; 
    grid on;

    figure(f2);
    plot(t,One_g, 'DisplayName', '1g');
    hold on;
    xlabel('Time (s)');
    ylabel('Lateral Acceleration (m/s^2)');
    title('RK4 Lateral Acceleration');
    legend; 
    grid on;

    figure(f3);
    xlabel('Time (s)');
    ylabel('Yaw Rate (rad/s)');
    title('Euler Yaw Rate');
    legend;
    grid on;
 
    figure(f4);
    xlabel('Time (s)');
    ylabel('Yaw Rate (rad/s)');
    title('RK4 Yaw Rate');
    legend;
    grid on;

speed = 200* 1000 / 3600:0.1:500* 1000 / 3600;
lambda_values = zeros(size(speed));

for i = 1:length(speed)
    u = speed(i);
       
    % Define vectors
    A_2 = [-(Cf + Cr)/(m*u), (-a*Cf + b*Cr)/(m*u) - u; -(a*Cf - b*Cr)/(Iz*u), -(a^2*Cf + b^2*Cr)/(Iz*u)];
    I = [1, 0; 0, 1];
  
    % Find A - lambda(I)
    syms lambda;
    C = A_2 - lambda*I;
    Determinant = C(1,1)*C(2,2) - C(1,2)*C(2,1);

    % Solve for lambda (eigenvalues) and store
    eigenvalues = solve(Determinant == 0,lambda);
    lambda_values(i) = max(double(eigenvalues));
end

% Find the highest stable speed (lambda < 0)
stable_speed = speed(lambda_values < 0);
if ~isempty(stable_speed)
    highest_stable_speed = max(stable_speed)*3600/1000;
    fprintf('The highest stable speed is: %d km/hr\n', highest_stable_speed);
else
    fprintf('No stable speeds found in the given range.\n');
end

% Car movement relative to ground frame
u = 75* 1000 / 3600;
A_3 = [-(Cf + Cr)/(m*u), (-a*Cf + b*Cr)/(m*u) - u; -(a*Cf - b*Cr)/(Iz*u), -(a^2*Cf + b^2*Cr)/(Iz*u)];
B_3 = [Cf/m; a*Cf/Iz];

% Initialize vectors: [y; yaw psi]
y_psi_values = zeros(2, N);
x_euler = zeros(2, N);

%Euler Method
for n = 1:N-1
    x_euler(:,n+1) = x_euler(:,n) + h * (A_3 * x_euler(:,n) + B_3 * delta);
end

for n = 1:N-1
    y_psi_values(:,n+1) = y_psi_values(:,n) + h * (x_euler(:,n));
end

psi_angle = y_psi_values(2,:);
y_dot = x_euler(1,:);
psi_dot = x_euler(2,:);
X_dot = zeros(size(psi_angle));
Y_dot = zeros(size(psi_angle));
X = zeros(size(psi_angle));
Y = zeros(size(psi_angle));

% 
for i = 1:length(psi_angle)
    psi = psi_angle(i);
    X_dot(i) = u * cos(psi) - (y_dot(i) + a * psi_dot(i)) * sin(psi);
    Y_dot(i) = (y_dot(i) + a * psi_dot(i)) * cos(psi) + u * sin(psi); 

    if i < length(psi_angle)
        X(i+1) = X(i) + h * X_dot(i);
        Y(i+1) = Y(i) + h * Y_dot(i);
    end
end

figure;
plot(X, Y, 'DisplayName', 'Path travelled by car');
hold on;
xlabel('Position in x (m)');
ylabel('Position in y (m)');
title('Path Travelled by Car');
legend;
grid on;
