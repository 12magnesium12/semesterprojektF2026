clear;
clc;
syms th(t) x(t) F_in(t) ddx ddth x_s dx_s th_s dth_s
% syms l m_p m_r m_c I_p I_r g b_p b_c;

l = 0.35; m_p = 0.002; m_r = 0.082; m_c = 0.5; I_p = m_p*0.35^2; I_r = (1/3)*m_r*0.35^2; g = 9.82; b_p = 0.0012; b_c = 5;

% Energies and Lagrangian
E_kin_c = (1/2)*m_c*diff(x, t)^2;

dx_r = diff(x, t)+diff(th,t)*cos(th)*l/2;
dy_r = -diff(th, t)*sin(th)*l/2;
E_kin_r = (1/2)*(I_r*diff(th, t)^2 + m_r*dx_r^2 + m_r*dy_r^2);

dx_p = diff(x, t)+diff(th, t)*cos(th)*l;
dy_p = -diff(th, t)*sin(th)*l;
E_kin_p = (1/2)*(I_p*diff(th, t)^2 + m_p*dx_p^2 + m_p*dy_p^2);

E_kin = (E_kin_c + E_kin_r + E_kin_p);
% E_kin = (E_kin_c + E_kin_r);

E_pot_r = m_r*g*cos(th)*l/2;
E_pot_p = m_p*g*cos(th)*l;
E_pot = (E_pot_p + E_pot_r);
% E_pot = (E_pot_r);


L = simplify(E_kin - E_pot);

% Equations of Motion
EqTh = diff(diff(L, diff(th,t)), t) - diff(L, th) == -b_p*diff(th,t);
EqX  = diff(diff(L, diff(x,t)), t) - diff(L, x)  == F_in - b_c*diff(x,t);

% Solve for accelerations
Eqs_sub = subs([EqTh, EqX], [diff(x,t,2), diff(th,t,2), diff(x,t), diff(th,t)], [ddx, ddth, dx_s, dth_s]);

sol = solve(Eqs_sub, [ddx, ddth]);

% State-Space Formulation and Linearization
% Map functional states to simple symbols for Jacobian
f = [dx_s; 
     sol.ddx; 
     dth_s; 
     sol.ddth];
f = subs(f, [x(t), diff(x,t), th(t), diff(th,t)], [x_s, dx_s, th_s, dth_s]);

A_sym = jacobian(f, [x_s, dx_s, th_s, dth_s]);
B_sym = jacobian(f, F_in);

% Substitute equilibrium (position=0, velocities=0) and constants
A_num = double(subs(A_sym, [x_s, dx_s, th_s, dth_s], [0, 0, 0, 0]));
B_num = double(subs(B_sym, [x_s, dx_s, th_s, dth_s], [0, 0, 0, 0]));

C = [1 0 0 0
     0 0 1 0];

D = [0; 0];

sys = ss(A_num, B_num, C, D);

% LQR
Q = diag([10000 1 5 1]);
R = 5;

K = lqr(A_num,B_num,Q,R)


%% Reference tracking (From 0.5m to 0m)
C_ref = [1 0 0 0]; % Focus on cart position
Nbar = -inv(C_ref*inv(A_num - B_num*K)*B_num);

% Define the system
sys_cl = ss(A_num - B_num*K, B_num*Nbar, C, D);

% 1. Change the Reference (Target) to 0
t = 0:0.01:2; 
r = 0 * ones(size(t)); % Target position is 0m

% 2. Set Initial Conditions [x; dx; th; dth]
% We start at x = 0.5, all other states (velocity/angle) are 0
x0 = [0.5; 0; 0; 0]; 

% 3. Run simulation with initial condition x0
[y, t, x_states] = lsim(sys_cl, r, t, x0);

% Plotting the result
figure
plot(t, y(:,1), 'LineWidth', 2)
ylabel('Position (m)')
title('Cart Moving from 0.5m to 0m')
grid on
hold on;

plot(t, y(:,2), 'LineWidth', 2)
ylabel('Angle (rad)')
xlabel('Time (s)')
title('Pendulum Response')
grid on

%% Reference tracking WITH Controller Noise
C_ref = C(1,:);
Nbar = -inv(C_ref*inv(A_num - B_num*K)*B_num);

% Create noise signal
t = 0:0.01:50;
noise_magnitude = 0.5; % Tweak this to make the noise stronger/weaker
noise = noise_magnitude * randn(length(t), 1); % normally distributed random noise

% Create reference signal
r = 0.2 * ones(length(t), 1); % 20 cm step

% Augment the System to accept two inputs: [r, noise]
% The new B matrix handles B*Nbar*r and B*noise
B_aug = [B_num*Nbar, B_num]; 
D_aug = [D, D]; 

sys_cl_noisy = ss(A_num - B_num*K, B_aug, C, D_aug);

% Combine inputs into a two-column matrix
% Column 1: reference, Column 2: noise
inputs = [r, noise];

% Simulate
[y, t, x] = lsim(sys_cl_noisy, inputs, t);

% Compute noisy control force u(t) applied to the plant
u_actual = zeros(length(t), 1);
for i = 1:length(t)
    % u = -Kx + Nbar*r + noise
    u_actual(i) = -K * x(i,:)' + Nbar * r(i) + noise(i); 
end

% plot
figure
subplot(3,1,1)
plot(t, y(:,1), 'LineWidth', 2)
xlabel('Tid (s)')
ylabel('Position (m)')
title('Cart Position (with noisy controller)')
grid on

subplot(3,1,2)
plot(t, y(:,2), 'LineWidth', 2) 
xlabel('Tid (s)')
ylabel('Angle (rad)')
title('Pendulum Angle')
grid on

subplot(3,1,3)
plot(t, u_actual, 'r', 'LineWidth', 1)
xlabel('Tid (s)')
ylabel('Force (N)')
title('Actual Noisy Control Signal (u)')
grid on