clear; clc; close all;

% Part C Question 1

% Vehicle Parameters
Cf = 25000;
Cr = 21000;
Iz = 2420;
delta = 0.1;

% Time setup
h = 0.01;
t = 0:h:400;
N = length(t);

%% Center of gravity shifted forward

% Define speed range in meters per second
speed_cgforward = 200* 1000 / 3600:1:1000* 1000 / 3600;
lambda_values_cgforward = zeros(size(speed_cgforward));

% Iterate over the speed_cgforward array
for i = 1:length(speed_cgforward)

    % Get the current speed 
    u = speed_cgforward(i);

    % Define constants 
    a = 1.1;
    b = 1.37;
    m = 1450;
       
    % Define vectors
    A_1 = [-(Cf + Cr)/(m*u), (-a*Cf + b*Cr)/(m*u) - u; -(a*Cf - b*Cr)/(Iz*u), -(a^2*Cf + b^2*Cr)/(Iz*u)];
    I = [1, 0; 0, 1];
  
    % Find A - lambda(I)
    syms lambda_cgforward;
    C = A_1 - lambda_cgforward*I;
    Determinant_cgforward = C(1,1)*C(2,2) - C(1,2)*C(2,1);

    % Solve for lambda (eigenvalues) and store
    eigenvalues_cgforward = solve(Determinant_cgforward == 0,lambda_cgforward);
    lambda_values_cgforward(i) = max(double(eigenvalues_cgforward));
end

% Find the highest stable speed (lambda < 0)
stable_speed_cgforward = speed_cgforward(lambda_values_cgforward < 0);


% Check for stable speeds
if ~isempty(stable_speed_cgforward)
    highest_stable_speed_cgforward = max(stable_speed_cgforward)*3600/1000;
    fprintf('The highest stable speed with a forward shifted center of gravity: %d km/hr\n', highest_stable_speed_cgforward);
else
    fprintf('No stable speeds found in the given range.\n');
end


%% Center of gravity shifted rearward 

% Define speed range in meters per second
speed_cgback = 100* 1000 / 3600:0.1:500* 1000 / 3600;
lambda_values_cgback = zeros(size(speed_cgback));

% Iterate over the speed_cgforward array
for i = 1:length(speed_cgback)

    % Get the current speed 
    u = speed_cgback(i);

    %Define Constants
    a = 1.186;
    b = 1.284;
    m = 1450;

    % Define vectors
    A_2 = [-(Cf + Cr)/(m*u), (-a*Cf + b*Cr)/(m*u) - u; -(a*Cf - b*Cr)/(Iz*u), -(a^2*Cf + b^2*Cr)/(Iz*u)];
    I = [1, 0; 0, 1];
  
    % Find A - lambda(I)
    syms lambda_cgback;
    C_2 = A_2 - lambda_cgback*I;
    Determinant_cgback = C_2(1,1)*C_2(2,2) - C_2(1,2)*C_2(2,1);

    % Solve for lambda (eigenvalues) and store
    eigenvalues_cgback = solve(Determinant_cgback == 0,lambda_cgback);
    lambda_values_cgback(i) = max(double(eigenvalues_cgback));
end

% Find the highest stable speed (lambda < 0)
stable_speed_cgback = speed_cgback(lambda_values_cgback < 0);

% Check for stable speeds
if ~isempty(stable_speed_cgback)
    highest_stable_speed_cgback = max(stable_speed_cgback)*3600/1000;
    fprintf('The highest stable speed with a backward shifted center of gravity: %d km/hr\n', highest_stable_speed_cgback);
else
    fprintf('No stable speeds found in the given range.\n');
end

%% Why driving fast/turning hard on slippery roads is dangerous 

% With a delta of 0.1, 0.35, 0.06
for delta = [0.1, 0.35, 0.06]
    f1 = figure;

    %Define constants
    a = 1.14;
    b = 1.33;
    m = 1400;
    
    %Changing the front and rear cornering stiffness based on the steering
    %angle
    if delta == 0.1
        Cf = 100;
        Cr = 100;
    elseif delta == 0.35
        Cf = 0;
        Cr = 0;
    else 
        Cf = 20000;
        Cr = 20000;
    end 

     
    % Define speed on ice in meters per second
    speed_on_ice = 200* 1000 / 3600:0.1:500* 1000 / 3600;
    lambda_values_on_ice = zeros(size(speed_on_ice));

    for i = 1:length(speed_on_ice)
        u = speed_on_ice(i);

        % Define vectors
        A_3 = [-(Cf + Cr)/(m*u), (-a*Cf + b*Cr)/(m*u) - u; -(a*Cf - b*Cr)/(Iz*u), -(a^2*Cf + b^2*Cr)/(Iz*u)];
        I = [1, 0; 0, 1];
  
        % Find A - lambda(I)
        syms lambda_on_ice;
        C = A_3 - lambda_on_ice*I;
        Determinant_ice = C(1,1)*C(2,2) - C(1,2)*C(2,1);

        % Solve for lambda (eigenvalues) and store
        eigenvalues_ice = solve(Determinant_ice == 0,lambda_on_ice);
        lambda_values_on_ice(i) = max(double(eigenvalues_ice));
    end

    % Find the highest stable speed on ice (lambda < 0)
    stable_speed_ice = speed_on_ice(lambda_values_on_ice < 0);
    if ~isempty(stable_speed_ice)
        highest_stable_speed_ice = max(stable_speed_ice)*3600/1000;
        fprintf('The highest stable speed on ice with delta = %.2f: %d km/hr\n', delta, highest_stable_speed_ice);
    else
        fprintf('No stable speeds found for driving on an icy road (not winter tires) in the given range for delta = %.2f.\n', delta);
    end

    for u = [30, 50, 70, 80, 90, 100, 150]* 1000 / 3600;
    
        A_3 = [-(Cf + Cr)/(m*u), (-a*Cf + b*Cr)/(m*u) - u; -(a*Cf - b*Cr)/(Iz*u), -(a^2*Cf + b^2*Cr)/(Iz*u)];
        B_3 = [Cf/m; a*Cf/Iz];

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
 
        
        % Loop through the psi_angle array
        for i = 1:length(psi_angle)
            psi = psi_angle(i);

            
            % Calculate the derivatives of X and Y based on the current inputs
            X_dot(i) = u * cos(psi) - (y_dot(i) + a * psi_dot(i)) * sin(psi);
            Y_dot(i) = (y_dot(i) + a * psi_dot(i)) * cos(psi) + u * sin(psi); 

            %Update the next iteraiton
            if i < length(psi_angle)
                X(i+1) = X(i) + h * X_dot(i);
                Y(i+1) = Y(i) + h * Y_dot(i);
            end
        end
    
        figure(f1);
        plot(X, Y, 'DisplayName', ['Path travelled, u = ' num2str(u*3600/1000) 'km/h']);
        hold on;
    end

    figure(f1);
    xlabel('Position in x (m)');
    ylabel('Position in y (m)');
    title(['Path Travelled by Car on Icy Road (No Winter Tires)- Delta = ' num2str(delta)]);
    legend;
    grid on;

end 


%% Why winter tires help
% With a delta of 0.1, 0.35, 0.06
for delta = [0.1, 0.35, 0.06]
    f1 = figure;

    %Define constants
    a = 1.14;
    b = 1.33;
    m = 1400;
    
    %Changing the front and rear cornering stiffness based on the steering
    %angle
    if delta == 0.1
        Cf = 5000;
        Cr = 5000;
    elseif delta == 0.35
        Cf = 0;
        Cr = 0;
    else 
        Cf = 20000;
        Cr = 20000;
    end 

    % Define speed on ice in meters per second
    speed_on_ice = 200* 1000 / 3600:0.1:500* 1000 / 3600;
    lambda_values_on_ice = zeros(size(speed_on_ice));

    for i = 1:length(speed_on_ice)
        u = speed_on_ice(i);

        % Define vectors
        A_3 = [-(Cf + Cr)/(m*u), (-a*Cf + b*Cr)/(m*u) - u; -(a*Cf - b*Cr)/(Iz*u), -(a^2*Cf + b^2*Cr)/(Iz*u)];
        I = [1, 0; 0, 1];
  
        % Find A - lambda(I)
        syms lambda_on_ice;
        C = A_3 - lambda_on_ice*I;
        Determinant_ice = C(1,1)*C(2,2) - C(1,2)*C(2,1);

        % Solve for lambda (eigenvalues) and store
        eigenvalues_ice = solve(Determinant_ice == 0,lambda_on_ice);
        lambda_values_on_ice(i) = max(double(eigenvalues_ice));
    end

    % Find the highest stable speed on ice (lambda < 0)
    stable_speed_ice = speed_on_ice(lambda_values_on_ice < 0);
    if ~isempty(stable_speed_ice)
        highest_stable_speed_ice = max(stable_speed_ice)*3600/1000;
        fprintf('The highest stable speed on ice with winter tires. Delta = %.2f: %d km/hr\n', delta, highest_stable_speed_ice);
    else
        fprintf('No stable speeds found for icy conditions with winter tires in the given range for delta = %.2f.\n', delta);
    end

    for u = [30, 50, 70, 80, 90, 100, 150]* 1000 / 3600;
    
        %Define the matrices
        A_3 = [-(Cf + Cr)/(m*u), (-a*Cf + b*Cr)/(m*u) - u; -(a*Cf - b*Cr)/(Iz*u), -(a^2*Cf + b^2*Cr)/(Iz*u)];
        B_3 = [Cf/m; a*Cf/Iz];

        %Intialize 
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
 
        for i = 1:length(psi_angle)
            psi = psi_angle(i);

            % Calculate the derivatives of X and Y based on the current inputs
            X_dot(i) = u * cos(psi) - (y_dot(i) + a * psi_dot(i)) * sin(psi);
            Y_dot(i) = (y_dot(i) + a * psi_dot(i)) * cos(psi) + u * sin(psi); 

            % Update the next iteration
            if i < length(psi_angle)
                X(i+1) = X(i) + h * X_dot(i);
                Y(i+1) = Y(i) + h * Y_dot(i);
            end
        end
    
        figure(f1);
        plot(X, Y, 'DisplayName', ['Path travelled, u = ' num2str(u*3600/1000) 'km/h']);
        hold on;
    end

    figure(f1);
    xlabel('Position in x (m)');
    ylabel('Position in y (m)');
    title(['Path Travelled by Car with Winter Tires- Delta = ' num2str(delta)]);
    legend;
    grid on;

end 
