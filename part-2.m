m = 1400;               % mass (kg)
a = 1.14;               % front axle distance (m)
b = 1.33;               % rear axle distance (m)
Cf = 25000;             % front cornering stiffness (N/rad)
Cr = 21000;             % rear cornering stiffness (N/rad)
Iz = 2420;              % yaw inertia (kg.m^2)
u = 75 * 1000 / 3600;   % longitudinal velocity (m/s)
delta = 0.05;           % step steer input (rad)

% State-space matrices
A = [-(Cf+Cr)/(m*u), (-a*Cf + b*Cr)/(m*u) - u;
     (-a*Cf + b*Cr)/(Iz*u), -(a^2*Cf + b^2*Cr)/(Iz*u)];

B = [Cf/m;
     a*Cf/Iz];



