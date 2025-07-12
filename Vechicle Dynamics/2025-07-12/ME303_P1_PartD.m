clear; clc; close all;

% Vehicle Parameters
m = 1400;
a = 1.14;
b = 1.33;
Iz = 2420;
delta = 0.1;
u = 400 * 1000/3600;

% Time setup
h = 0.01;
t = 0:h:20;
N = length(t);

%% Case 1: Cf = Cr

f1 = figure;

for Cf = [200, 400, 600, 800, 1000, 10000, 20000]
    Cr = Cf;

    A_1 = [-(Cf + Cr)/(m*u), (-a*Cf + b*Cr)/(m*u) - u; -(a*Cf - b*Cr)/(Iz*u), -(a^2*Cf + b^2*Cr)/(Iz*u)];
    B_1 = [Cf/m; a*Cf/Iz];

    % Initialize vectors
    y_psi_values = zeros(2, N);
    x_euler = zeros(2, N);

    for n = 1:N-1
        x_euler(:,n+1) = x_euler(:,n) + h * (A_1 * x_euler(:,n) + B_1 * delta);
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

    for i = 1:length(psi_angle)
        psi = psi_angle(i);
        X_dot(i) = u * cos(psi) - (y_dot(i) + a * psi_dot(i)) * sin(psi);
        Y_dot(i) = (y_dot(i) + a * psi_dot(i)) * cos(psi) + u * sin(psi); 

        if i < length(psi_angle)
            X(i+1) = X(i) + h * X_dot(i);
            Y(i+1) = Y(i) + h * Y_dot(i);
        end
    end

    figure(f1);
    plot(X, Y, 'DisplayName', ['Path travelled by car at constant speed with Cf = Cr = ' num2str(Cf)]);
    hold on;
    
    I = [1, 0; 0, 1];
  
    % Find A - lambda(I)
    syms lambda;
    C = A_1 - lambda*I;
    Determinant = C(1,1)*C(2,2) - C(1,2)*C(2,1);

    % Solve for lambda (eigenvalues) and store
    eigenvalues = solve(Determinant == 0,lambda);
    lambda_values = max(double(eigenvalues));

    % Find if stable (lambda < 0)
    if lambda_values < 0
        fprintf(['Stable at Cf = ' num2str(Cf) ' Cr = ' num2str(Cr) '\n']);
    else
        fprintf(['Not stable at Cf = ' num2str(Cf) ' Cr = ' num2str(Cr) '\n']);
    end
end

figure(f1);
xlabel('Position in x (m)');
ylabel('Position in y (m)');
title('Path Travelled by Car for Cf = Cr');
legend;
grid on;




%% Case 2: Cf > Cr

f2 = figure;

for Cf = [200, 400, 600, 800, 1000]
    Cr = Cf/2;

    A_2 = [-(Cf + Cr)/(m*u), (-a*Cf + b*Cr)/(m*u) - u; -(a*Cf - b*Cr)/(Iz*u), -(a^2*Cf + b^2*Cr)/(Iz*u)];
    B_2 = [Cf/m; a*Cf/Iz];

    % Initialize vectors
    y_psi_values = zeros(2, N);
    x_euler = zeros(2, N);

    for n = 1:N-1
        x_euler(:,n+1) = x_euler(:,n) + h * (A_2 * x_euler(:,n) + B_2 * delta);
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

    for i = 1:length(psi_angle)
        psi = psi_angle(i);
        X_dot(i) = u * cos(psi) - (y_dot(i) + a * psi_dot(i)) * sin(psi);
        Y_dot(i) = (y_dot(i) + a * psi_dot(i)) * cos(psi) + u * sin(psi); 

        if i < length(psi_angle)
            X(i+1) = X(i) + h * X_dot(i);
            Y(i+1) = Y(i) + h * Y_dot(i);
        end
    end

    figure(f2);
    plot(X, Y, 'DisplayName', ['Path travelled. Cf = ' num2str(Cf) ' Cr = ' num2str(Cr)]);
    hold on;

    I = [1, 0; 0, 1];
  
    % Find A - lambda(I)
    syms lambda;
    C = A_2 - lambda*I;
    Determinant = C(1,1)*C(2,2) - C(1,2)*C(2,1);

    % Solve for lambda (eigenvalues) and store
    eigenvalues = solve(Determinant == 0,lambda);
    lambda_values = max(double(eigenvalues));

    % Find if stable (lambda < 0)
    if lambda_values < 0
        fprintf(['Stable at Cf = ' num2str(Cf) ' Cr = ' num2str(Cr) '\n']);
    else
        fprintf(['Not stable at Cf = ' num2str(Cf) ' Cr = ' num2str(Cr) '\n']);
    end

end

figure (f2);
xlabel('Position in x (m)');
ylabel('Position in y (m)');
title('Path Travelled by Car at constant speed for Cf > Cr');
legend;
grid on;

%% Case 3: Cr > Cf

f3 = figure;

for Cf = [100, 200, 300, 400, 500]
    Cr = Cf*2;

    A_3 = [-(Cf + Cr)/(m*u), (-a*Cf + b*Cr)/(m*u) - u; -(a*Cf - b*Cr)/(Iz*u), -(a^2*Cf + b^2*Cr)/(Iz*u)];
    B_3 = [Cf/m; a*Cf/Iz];

    % Initialize vectors
    y_psi_values = zeros(2, N);
    x_euler = zeros(2, N);

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

    for i = 1:length(psi_angle)
        psi = psi_angle(i);
        X_dot(i) = u * cos(psi) - (y_dot(i) + a * psi_dot(i)) * sin(psi);
        Y_dot(i) = (y_dot(i) + a * psi_dot(i)) * cos(psi) + u * sin(psi); 

        if i < length(psi_angle)
            X(i+1) = X(i) + h * X_dot(i);
            Y(i+1) = Y(i) + h * Y_dot(i);
        end
    end

    figure(f3);
    plot(X, Y, 'DisplayName', ['Path travelled. Cf = ' num2str(Cf) ' Cr = ' num2str(Cr)]);
    hold on;

    I = [1, 0; 0, 1];
  
    % Find A - lambda(I)
    syms lambda;
    C = A_3 - lambda*I;
    Determinant = C(1,1)*C(2,2) - C(1,2)*C(2,1);

    % Solve for lambda (eigenvalues) and store
    eigenvalues = solve(Determinant == 0,lambda);
    lambda_values = max(double(eigenvalues));

    % Find if stable (lambda < 0)
    if lambda_values < 0
        fprintf(['Stable at Cf = ' num2str(Cf) ' Cr = ' num2str(Cr) '\n']);
    else
        fprintf(['Not stable at Cf = ' num2str(Cf) ' Cr = ' num2str(Cr) '\n']);
    end

end

figure (f3);
xlabel('Position in x (m)');
ylabel('Position in y (m)');
title('Path Travelled by Car at constant speed for Cf < Cr');
legend;
grid on;