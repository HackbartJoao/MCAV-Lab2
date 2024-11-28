clear all; close all; clc

%% Vehicle parameters
m = 1050;  % kg
J = diag([1750, 700, 1750]);  % Inertia matrix

%% Thruster configuration
num_thrusters = 4;
thruster_positions = [...
    1,  2, 1;  % Thruster 1 position relative to CM
   -1,  2, 1;  % Thruster 2
   -1, -2, 1;  % Thruster 3
    1, -2, 1   % Thruster 4
];

%% Environment parameters
g_mars = 3.72;  % Mars gravity (m/s^2)
rho_mars = 0.020;  % Mars atmosphere density (kg/m^3)
Cd = 1.2;  % Drag coefficient
A = 8.0;  % Reference area (m^2)

%% Dynamics
% Initial conditions
X0 = zeros(12,1);  % Initial state: [p; v; euler_angles; omega]
X0(1) = 0;         % Initial X position (m)
X0(2) = 0;         % Initial Y position (m)
X0(3) = -2100;     % Initial Z position (m)
X0(4) = 0;         % Initial X velocity (m/s)
X0(5) = 0;         % Initial Y velocity (m/s)      
X0(6) = 0;         % Initial Z velocity (m/s)
X0(7) = 0;         % Initial roll (ϕ) angle (rad)
X0(8) = 0;         % Initial pitch (θ) angle (rad)
X0(9) = 0;         % Initial yaw (ψ) angle (rad)
X0(10) = 0;        % Initial roll angular velocity (rad/s)
X0(11) = 0;        % Initial pitch angular velocity (rad/s)
X0(12) = 0;        % Initial yaw angular velocity (rad/s)

Dt = 0.01;
t = 0:Dt:120;
Nsim = length(t);
X(:,1) = X0;

Uz = ones(12,1);   % Control input: [T; alpha; beta]
Uz(1:4) = 1000;       % Initial thrust force (N)
Uz(5) = 0;         % Initial T1 alpha angle (rad)
Uz(6) = 0;         % Initial T2 alpha angle (rad)
Uz(7) = 0;         % Initial T3 alpha angle (rad)
Uz(8) = 0;         % Initial T4 alpha angle (rad)
Uz(9) = 0;         % Initial T1 beta angle (rad)
Uz(10) = 0;        % Initial T2 beta angle (rad)
Uz(11) = 0;        % Initial T3 beta angle (rad)
Uz(12) = 0;        % Initial T4 beta angle (rad)

U = Uz * (t<30);

for k = 1:Nsim
    % State vector X = [p; v; euler_angles; omega]
    % Input vector U = [T1,T2,T3,T4, alpha1,alpha2,alpha3,alpha4, beta1,beta2,beta3,beta4]
    
    % Extract states
    p = X(1:3,k);       % Position in inertial frame
    v = X(4:6,k);       % Velocity in body frame
    euler = X(7:9,k);   % Euler angles
    omega = X(10:12,k); % Angular velocity in body frame
    
    % Extract inputs
    thrusts = U(1:4,k);
    alphas = U(5:8,k);   % x-axis rotation
    betas = U(9:12,k);   % y-axis rotation
    
    % Get rotation matrix (body to inertial)
    R = Euler2R(euler);
    
    % Calculate forces and moments   
    F_b = zeros(3,1);
    M_b = zeros(3,1);
    
    for i = 1:num_thrusters
        % Base thrust vector (aligned with z-axis)
        Ti = [0; 0; -thrusts(i)];
        % Apply rotations
        Ti = rot3D(betas(i),2) * rot3D(alphas(i),1) * Ti;
        % Add to total force
        F_b = F_b + Ti;
        % Calculate moment
        M_b = M_b + cross(thruster_positions(i,:)', Ti);
    end
    
    F_g = R' * [0; 0; m * g_mars];  % Gravity in body frame
    F_d = -0.5 * rho_mars * Cd * A * norm(v) * v; % Drag in body frame
    
    % Total forces and moments in body frame
    F_total = F_b + F_g + F_d;
    M_total = M_b;
    
    % State derivatives
    dp = R * v;  % Position derivative
    dv = F_total/m - cross(omega, v);  % Velocity derivative
    
    % Euler angle derivatives
    phi = euler(1); theta = euler(2);
    W = [1, sin(phi)*tan(theta), cos(phi)*tan(theta);
         0, cos(phi), -sin(phi);
         0, sin(phi)/cos(theta), cos(phi)/cos(theta)];
    deuler = W * omega;
    
    % Angular velocity derivative
    domega = J \ (M_total - cross(omega, J*omega));
    
    % Combine all derivatives
    dX = [dp; dv; deuler; domega];

    X(:,k+1) = X(:,k) + Dt*dX;

        % Termination condition
    if X(3,k) > 0
        X = X(:,1:k); % Trim unused states
        U = U(:,1:k); % Trim unused states
        t = 0:Dt:((k-1)*Dt);
        t_end = (k-1)*Dt;
        fprintf('Loop stopped at %.2fs because lander has crashed into terrain!\n',t_end);
        break;
    end
end
X = X(:,1:k);
U = U(:,1:k);

%% Plotting results
figure('Name', 'Mars Lander Simulation Results', 'Position', [100 100 1200 800]);

% Position plot
subplot(2,2,1);
plot(t, X(1:3,:)');
grid on;
legend('x', 'y', 'z');
title('Position');
xlabel('Time (s)');
ylabel('Position (m)');

% Velocity plot
subplot(2,2,2);
plot(t, X(4:6,:)');
grid on;
legend('Vx', 'Vy', 'Vz');
title('Velocity');
xlabel('Time (s)');
ylabel('Velocity (m/s)');

% Attitude plot
subplot(2,2,3);
plot(t, rad2deg(X(7:9,:)'));
grid on;
legend('\phi', '\theta', '\psi');
title('Euler Angles');
xlabel('Time (s)');
ylabel('Angle (deg)');

% Control inputs
subplot(2,2,4);
plot(t, U(1:4,:)');
grid on;
legend('T1', 'T2', 'T3', 'T4');
title('Thruster Forces');
xlabel('Time (s)');
ylabel('Thrust (N)');

% 3D trajectory visualization
figure('Name', 'Lander Trajectory', 'Position', [700 100 600 600]);
plot3(X(1,:), X(2,:), X(3,:), 'b-');
grid on;
hold on;
% Plot start and end points
plot3(X(1,1), X(2,1), X(3,1), 'go', 'MarkerSize', 10, 'LineWidth', 2);
% Flip Z Axis
set(gca, 'ZDir', 'reverse');
plot3(X(1,end), X(2,end), X(3,end), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
% Flip Z Axis
set(gca, 'ZDir', 'reverse');
legend('Trajectory', 'Start', 'End');
title('3D Trajectory');
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
view(45, 30);
axis equal;