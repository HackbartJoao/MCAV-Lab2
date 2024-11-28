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
D_sym = 0;

% Define a non-singular operating point
x_op = [0; 0; -2100; 10; 10; 10; 0; 0; 0; 0; 0; 0]; % Example
u_op = [1000; 0; 0; 0; 0; 0; 0; 0; 0]; % Example thrust and angles


% Evaluate symbolic matrices
A = double(subs(A_sym, [x_state; u], [x_op; u_op]));
B = double(subs(B_sym, [x_state; u], [x_op; u_op]));
C = double(subs(C_sym, [x_state; u], [x_op; u_op]));
D = double(subs(D_sym, [x_state; u], [x_op; u_op]));


%% Model Analysis
E = eig(A(1:4,1:4));
% J = jordan(E);

Mc = [B,  A*B, (A^2)*B, (A^3)*B];
Mo = [C;  C*A; C*(A^2); C*(A^3)];

fprintf('Mc rank:')
disp(rank(Mc))
fprintf('Mo rank:')
disp(rank(Mo))

sys = ss(A,B,C,D);
G = tf(sys);


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
