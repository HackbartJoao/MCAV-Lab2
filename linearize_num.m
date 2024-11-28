clear all; close all; clc

% Define symbolic variables
syms p_x p_y p_z v_x v_y v_z phi theta psi omega_x omega_y omega_z real
syms T alpha1 alpha2 alpha3 alpha4 beta1 beta2 beta3 beta4 real


m = 1050;           % kg
g_mars = 3.72;    % m/s^2
rho_mars = 0.02;  % kg/m^3
A = 8;          % m^2
Cd = 1.2;         % drag coefficient
Jxx = 1750; Jyy = 700; Jzz = 1750; % kg*m^2

% State vector
x_state = [p_x; p_y; p_z; v_x; v_y; v_z; phi; theta; psi; omega_x; omega_y; omega_z];

% Control input vector
u = [T; alpha1; alpha2; alpha3; alpha4; beta1; beta2; beta3; beta4];

% Output vector
y_t = [p_x; p_y; p_z; v_x; v_y; v_z; phi; theta; psi; omega_x; omega_y; omega_z];

% Inertia matrix
J = diag([Jxx, Jyy, Jzz]);

% Rotation matrix (Body to Inertial)
R = [cos(theta)*cos(psi), sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi), cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi);
     cos(theta)*sin(psi), sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi), cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi);
     -sin(theta),         sin(phi)*cos(theta),                         cos(phi)*cos(theta)];

% Base thrust direction vector
T_base = [0; 0; -1];

% Define thruster positions
thruster_positions = [...
         1,  2, 1;  % Thruster 1 position relative to CM
        -1,  2, 1;  % Thruster 2
        -1, -2, 1;  % Thruster 3
         1, -2, 1   % Thruster 4
];

% Initialize force and moment vectors
F_b = sym(zeros(3, 1));
M_b = sym(zeros(3, 1));

% Compute thrust forces and moments for each thruster
for i = 1:4
    % Thrust rotation (using corresponding alpha and beta for each thruster)
    T_rot = rot3D(u(1+i), 2) * rot3D(u(5+i), 1) * T_base;
    
    % Thrust force in body frame
    thrust_force = u(i) * T_rot;
    F_b = F_b + thrust_force;
    
    % Moment contribution from thrust
    M_b = M_b + cross(thruster_positions(i, :)', thrust_force);
end

% Gravity force in the body frame
F_g = R' * [0; 0; m * g_mars];

% Velocity vector
v = [v_x; v_y; v_z];

% Drag force in the body frame
F_d = 0.5 * rho_mars * A * Cd * norm(v) * (-v);

% Total force and moment
F_total = F_b + F_g + F_d;
M_total = M_b;

% Position derivative
dp = R * v;

% Velocity derivative
dv = F_total / m - cross([omega_x; omega_y; omega_z], v);

% Euler angle derivative matrix
W = [1, sin(phi)*tan(theta), cos(phi)*tan(theta);
     0, cos(phi), -sin(phi);
     0, sin(phi)/cos(theta), cos(phi)/cos(theta)];
deuler = W * [omega_x; omega_y; omega_z];

% Angular velocity derivative
domega = J \ (M_total - cross([omega_x; omega_y; omega_z], J * [omega_x; omega_y; omega_z]));

% Final state derivative
dx = [dp; dv; deuler; domega];


%% Jacobian Linearization
A_sym = jacobian(dx,x_state);
B_sym = jacobian(dx,u);
C_sym = jacobian(x_state,y_t);
D_sym = zeros(12,9);

%% OPERATING POINTS

% OP.1 - a)
x_op1 = [0; 0; -2100; 10; 10; 10; 0; 0; 0; 0; 0; 0];
u_op1 = [1000; 0; 0; 0; 0; 0; 0; 0; 0];
A1 = double(subs(A_sym, [x_state; u], [x_op1; u_op1]));
B1 = double(subs(B_sym, [x_state; u], [x_op1; u_op1]));
C1 = double(subs(C_sym, [x_state; u], [x_op1; u_op1]));
D1 = double(subs(D_sym, [x_state; u], [x_op1; u_op1]));

% OP.1 - b)
x_op2 = [0; 0; -2100; 10; 10; 10; 0; 0; 0; 0; 0; 0];
u_op2 = [1000; 0; 0; 0; 0; 0; 0; 0; 0];
A2 = double(subs(A_sym, [x_state; u], [x_op2; u_op2]));
B2 = double(subs(B_sym, [x_state; u], [x_op2; u_op2]));
C2 = double(subs(C_sym, [x_state; u], [x_op2; u_op2]));
D2 = double(subs(D_sym, [x_state; u], [x_op2; u_op2]));

% OP.2 - a)
x_op3 = [0; 0; -2100; 10; 10; 10; 0; 0; 0; 0; 0; 0];
u_op3 = [1000; 0; 0; 0; 0; 0; 0; 0; 0];
A3 = double(subs(A_sym, [x_state; u], [x_op3; u_op3]));
B3 = double(subs(B_sym, [x_state; u], [x_op3; u_op3]));
C3 = double(subs(C_sym, [x_state; u], [x_op3; u_op3]));
D3 = double(subs(D_sym, [x_state; u], [x_op3; u_op3]));

% OP.2 - b)
x_op4 = [0; 0; -2100; 10; 10; 10; 0; 0; 0; 0; 0; 0];
u_op4 = [1000; 0; 0; 0; 0; 0; 0; 0; 0];
A4 = double(subs(A_sym, [x_state; u], [x_op4; u_op4]));
B4 = double(subs(B_sym, [x_state; u], [x_op4; u_op4]));
C4 = double(subs(C_sym, [x_state; u], [x_op4; u_op4]));
D4 = double(subs(D_sym, [x_state; u], [x_op4; u_op4]));

%% Model Analysis

%% OP.1 - a)
A1_eig = eig(A1);
A1_jor = jordan(A1);
Mc1 = [B1,  A1*B1, (A1^2)*B1, (A1^3)*B1];
Mo1 = [C1;  C1*A1; C1*(A1^2); C1*(A1^3)];
R_mc1 = rank(Mc1);
R_mo1 = rank(Mo1);
sys1 = ss(A1,B1,C1,D1);
svd1 = svd(sys1(1,1));
G1 = zpk(tf(sys1));

%% OP.1 - b)
A2_eig = eig(A2);
A2_jor = jordan(A2);
Mc2 = [B2,  A2*B2, (A2^2)*B2, (A2^3)*B2];
Mo2 = [C2;  C2*A2; C2*(A2^2); C2*(A2^3)];
R_mc2 = rank(Mc2);
R_mo2 = rank(Mo2);
sys2 = ss(A2,B2,C2,D2);
svd2 = svd(sys2(1,1));
G2 = zpk(tf(sys2));

%% OP.2 - a)
A3_eig = eig(A3);
A3_jor = jordan(A3);
Mc3 = [B3,  A3*B3, (A3^2)*B3, (A3^3)*B3];
Mo3 = [C3;  C3*A3; C3*(A3^2); C3*(A3^3)];
R_mc3 = rank(Mc3);
R_mo3 = rank(Mo3);
sys3 = ss(A3,B3,C3,D3);
svd3 = svd(sys3(1,1));
G3 = zpk(tf(sys3));

%% OP.2 - b)
A4_eig = eig(A4);
A4_jor = jordan(A4);
Mc4 = [B4,  A4*B4, (A4^2)*B4, (A4^3)*B4];
Mo4 = [C4;  C4*A4; C4*(A4^2); C4*(A4^3)];
R_mc4 = rank(Mc4);
R_mo4 = rank(Mo4);
sys4 = ss(A4,B4,C4,D4);
svd4 = svd(sys4(1,1));
G4 = zpk(tf(sys4));

%% Helper function for 3D rotation matrices
function R = rot3D(angle, axis)
    % Create rotation matrix around specified axis
    % axis: 1 for x-rotation, 2 for y-rotation, 3 for z-rotation
    switch axis
        case 1 % x-axis rotation
            R = [1, 0, 0;
                 0, cos(angle), -sin(angle);
                 0, sin(angle), cos(angle)];
        case 2 % y-axis rotation
            R = [cos(angle), 0, sin(angle);
                 0, 1, 0;
                 -sin(angle), 0, cos(angle)];
        case 3 % z-axis rotation
            R = [cos(angle), -sin(angle), 0;
                 sin(angle), cos(angle), 0;
                 0, 0, 1];
        otherwise
            error('Invalid rotation axis');
    end
end
