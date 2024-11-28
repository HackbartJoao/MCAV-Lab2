clear all; close all; clc

% Vehicle parameters
m = 1050;  % kg
J = diag([1750, 700, 1750]);  % Inertia matrix

% Thruster configuration
num_thrusters = 4;
thruster_positions = [...
    1,  2, 0;  % Thruster 1 position relative to CM
   -1,  2, 0;  % Thruster 2
   -1, -2, 0;  % Thruster 3
    1, -2, 0   % Thruster 4
];

% Environment parameters
g_mars = 3.72;  % Mars gravity (m/s^2)
rho_mars = 0.020;  % Mars atmosphere density (kg/m^3)
Cd = 1.2;  % Drag coefficient
A = 8.0;  % Reference area (m^2)

%% OP.1(a) - Maximum vertical descent velocity
Vt = 89;
X_op1a = zeros(12,1);
X_op1a(3) = -2100;  % Altitude
X_op1a(6) = Vt;  % Vertical velocity
Teq_op1a = ((Euler2R(X_op1a(7:9)))' * [0; 0; m * g_mars])+(0.5 * rho_mars * norm(X_op1a(4:6)) * (-X_op1a(4:6)) * A * Cd);
U_op1a = zeros(12,1);
U_op1a(1:4) = Teq_op1a(3,1)/4;

%% OP.1(b) - Median vertical descent velocity
Vm = 45;
X_op1b = zeros(12,1);
X_op1b(3) = -2100;  % Altitude
X_op1b(6) = Vm;  % Vertical velocity
Teq_op1b = ((Euler2R(X_op1b(7:9)))' * [0; 0; m * g_mars])+(0.5 * rho_mars * norm(X_op1b(4:6)) * (-X_op1b(4:6)) * A * Cd);
U_op1b = zeros(12,1);
U_op1b(1:4) = Teq_op1b(3,1)/4;

%% OP.2(a) - Hover (zero velocity)
X_op2a = zeros(12,1);
X_op2a(3) = -21.3;  % Altitude
U_op2a = zeros(12,1);
U_op2a(1:4) = (m * g_mars) / 4;

%% OP.2(b) - Low forward velocity cruise
Vc = 10;  % m/s

X_op2b = zeros(12,1);
X_op2b(5) = Vc;  % Forward velocity

U_op2b = zeros(12,1);       % Control input: [T; alpha; beta]
U_op2b(1) = 975.5;          % Initial T1 force (N)
U_op2b(2) = 975.5;          % Initial T2 force (N)
U_op2b(3) = 976.5;          % Initial T3 force (N)
U_op2b(4) = 976.5;          % Initial T4 force (N)
U_op2b(5) = 0.0010235;      % Initial T1 alpha angle (rad)
U_op2b(6) = 0.0010235;      % Initial T2 alpha angle (rad)
U_op2b(7) = 0.0010235;      % Initial T3 alpha angle (rad)
U_op2b(8) = 0.0010235;      % Initial T4 alpha angle (rad)
U_op2b(9) = 0;              % Initial T1 beta angle (rad)
U_op2b(10) = 0;             % Initial T2 beta angle (rad)
U_op2b(11) = 0;             % Initial T3 beta angle (rad)
U_op2b(12) = 0;             % Initial T4 beta angle (rad)

% Extract inputs
thrusts = U_op2b(1:4);
alphas = U_op2b(5:8);   % x-axis rotation
betas = U_op2b(9:12);   % y-axis rotation

% Calculate forces and moments   
F_b = zeros(3,1);
M_b = zeros(3,1);

for i = 1:4
    % Base thrust vector (aligned with z-axis)
    T0 = [0; 0; -thrusts(i)];
    % Apply rotations
    T = rot3D(betas(i),2) * rot3D(alphas(i),1) * T0;
    % Add to total force
    F_b = F_b + T;
    % Calculate moment
    M_b = M_b + cross(thruster_positions(i,:)', T);
end