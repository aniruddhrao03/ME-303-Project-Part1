clear; clc; close all;

% Vehicle Parameters
m = 1400;
a = 1.14;
b = 1.33;
Cf = 25000;
Cr = 21000;
Iz = 2420;
u = 75 * 1000 / 3600;   % m/s
delta = 0.1;

%A and B matrices 
A = [-(Cf + Cr)/(m*u), (-a*Cf + b*Cr)/(m*u) - u; -(a*Cf - b*Cr)/(Iz*u), -(a^2*Cf + b^2*Cr)/(Iz*u)];
B = [Cf/m; a*Cf/Iz];

% Time setup
h = 0.1;
t = 0:h:5;
N = length(t);

% Initialize vectors: [v_y; psi]
x_euler = zeros(2, N);
x_rk4   = zeros(2, N);

%Euler Method
for n = 1:N-1
    x_euler(:,n+1) = x_euler(:,n) + h * (A * x_euler(:,n) + B * delta);
end

%RK4 Method 
for n = 1:N-1
    k1 = A * x_rk4(:,n) + B * delta;
    k2 = A * (x_rk4(:,n) + 0.5*h*k1) + B * delta;
    k3 = A * (x_rk4(:,n) + 0.5*h*k2) + B * delta;
    k4 = A * (x_rk4(:,n) + h*k3) + B * delta;
    x_rk4(:,n+1) = x_rk4(:,n) + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
end

%Compare Plot of Lateral Velocity for Euler and RK4
figure;
plot(t, x_euler(1,:), 'r', 'DisplayName', ['Euler v_y, u = ' num2str(u*3600/1000) 'km/h']); %Euler Plot
hold on;
plot(t, x_rk4(1,:), 'b', 'DisplayName', ['RK4 v_y, u = ' num2str(u*3600/1000) 'km/h']); %RK4 Plot
xlabel('Time (s)');
ylabel('Lateral Velocity (m/s)');
title('Lateral Velocity Response: Euler vs RK4');
legend; 
grid on;

% Compare Plot of Yaw Angle for Euler and RK4
figure;
plot(t, x_euler(2,:), 'r', 'DisplayName', ['Euler yaw angle, u = '  num2str(u*3600/1000) 'km/h']);
hold on;
plot(t, x_rk4(2,:), 'b','DisplayName', ['RK4 yaw angle, u = '  num2str(u*3600/1000) 'km/h']);
xlabel('Time (s)');
ylabel('Yaw Rate (rad/s)');
title('Yaw Rate: Euler vs RK4');
legend;
grid on;

%Grid Independence Check

dt_vec = [0.1, 0.05, 0.025, 0.0125, 0.00625];  % Time intervals reduced by half each time 
% Initialize vectors  
euler_final = zeros(size(dt_vec));
rk4_final   = zeros(size(dt_vec));

% Loop through each time step and compute.
for i = 1:length(dt_vec)    
    % Set time step
    h = dt_vec(i);
    t = 0:h:5;
    N = length(t);
    x0 = [0; 0];

    % Euler
    x_euler = zeros(2, N);

    % Initial condition
    x_euler(:,1) = x0;
    for n = 1:N-1
        x_euler(:,n+1) = x_euler(:,n) + h * (A * x_euler(:,n) + B * delta);
    end

    % RK4
    x_rk4 = zeros(2, N);

    % Initial condition
    x_rk4(:,1) = x0;
    for n = 1:N-1 
        k1 = A * x_rk4(:,n) + B * delta;
        k2 = A * (x_rk4(:,n) + 0.5*h*k1) + B * delta;
        k3 = A * (x_rk4(:,n) + 0.5*h*k2) + B * delta;
        k4 = A * (x_rk4(:,n) + h*k3) + B * delta;
        x_rk4(:,n+1) = x_rk4(:,n) + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
    end

    % Store norm of final state (v_y and psi)
    euler_final(i) = norm(x_euler(:,end));
    rk4_final(i)   = norm(x_rk4(:,end));
end

% Error Calculation
euler_err = abs(euler_final - euler_final(end));
rk4_err   = abs(rk4_final   - rk4_final(end));

% Grid Independence Log-log plot
figure;
loglog(dt_vec, euler_err, 'r', 'DisplayName', 'Euler Error');
hold on;
loglog(dt_vec, rk4_err, 'b', 'DisplayName', 'RK4 Error');
xlabel('Time (s)');
ylabel('Error');
title('Grid Independence Check');
legend;
grid on;

%% Part B

f1 = figure;
f2 = figure;
f3 = figure;
f4 = figure;

One_g = -9.81* ones(size(t));

for u = [25, 50, 75, 100, 200, 300]* 1000 / 3600
    
    A = [-(Cf + Cr)/(m*u), (-a*Cf + b*Cr)/(m*u) - u; -(a*Cf - b*Cr)/(Iz*u), -(a^2*Cf + b^2*Cr)/(Iz*u)];
    B = [Cf/m; a*Cf/Iz];

    % Initialize vectors: [v_y; yaw rate]
    x_euler = zeros(2, N);
    x_rk4   = zeros(2, N);
    
    % Initialize vector: a_y (lateral acceleration)
    a_y_euler = zeros(1, N);
    a_y_rk4 = zeros(1, N);

    %Euler Method
    for n = 1:N-1
        x_euler(:,n+1) = x_euler(:,n) + h * (A * x_euler(:,n) + B * delta);
        a_y_euler (1, :) = (-(Cf + Cr)/(m*u))*x_euler(1,:) + ((-a*Cf + b*Cr)/(m*u) - u)*x_euler(2,:) + (Cf/m)*delta;
    end

    %RK4 Method 
    for n = 1:N-1
        k1 = A * x_rk4(:,n) + B * delta;
        k2 = A * (x_rk4(:,n) + 0.5*h*k1) + B * delta;
        k3 = A * (x_rk4(:,n) + 0.5*h*k2) + B * delta;
        k4 = A * (x_rk4(:,n) + h*k3) + B * delta;
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
    plot(t, x_euler(2,:), 'DisplayName', ['Euler yaw angle, u = '  num2str(u*3600/1000) 'km/h']);
    hold on;

    %Plot of Yaw Angle for RK4
    figure(f4);
    plot(t, x_rk4(2,:),'DisplayName', ['RK4 yaw angle, u = '  num2str(u*3600/1000) 'km/h']);
    hold on;
end

    figure(f1);
    plot(t,One_g, 'DisplayName', ['1g']);
    hold on;
    xlabel('Time (s)');
    ylabel('Lateral Acceleration (m/s^2)');
    title('Euler Lateral Accleration');
    legend; 
    grid on;

    figure(f2);
    plot(t,One_g, 'DisplayName', ['1g']);
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

speed = 220* 1000 / 3600:0.01:250* 1000 / 3600;
lambda_values = zeros(size(speed));

for i = 1:length(speed)
    u = speed(i);
       
    % Define vectors
    A = [-(Cf + Cr)/(m*u), (-a*Cf + b*Cr)/(m*u) - u; -(a*Cf - b*Cr)/(Iz*u), -(a^2*Cf + b^2*Cr)/(Iz*u)];
    I = [1, 0; 0, 1];
  
    % Find A - lambda(I)
    syms lambda;
    C = A - lambda*I;
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
  
%Car movement relative to ground frame

psi_angle = 0:1:360;
y_dot = x_rk4(1,:);
psi_dot = x_rk4(2,:);
X_dot = zeros(size(psi_angle));
Y_dot = zeros(size(psi_angle));
X = zeros(size(psi_angle));
Y = zeros(size(psi_angle));

% Needed to specific x_rk4 (1,1) think this is wrong
for i = 1:length(psi_angle)
    u = 100 * 1000/3600; 
    psi = psi_angle(i);
    X_dot(i) = u * cosd(psi) - (y_dot(i) + a * psi_dot(i)) * sind(psi);
    Y_dot(i) = (y_dot(i) + a * psi_dot(i)) * cosd(psi) + u * sind(psi); 

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


