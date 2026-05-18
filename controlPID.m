clear;
clc;
%%
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


G = tf(sys);

s = tf('s');

% Extract the functional transfer functions from your state-space
G_x = minreal(G(1,1));
G_th = minreal(G(2,1));



%%
% add integrator 1/s (pole at s=0) to remove zero in G_th
% using Rltool find suitable placement of two zeros from PID
% rltool(minreal(-G_th*(1/s)))
C_inner = -  19.788* (s+3.273) *(s+3.846)/s;
% C_inner = - 5.8725 *(s+5.443) *(s+4.991)/s;

T_inner = feedback(G_th * C_inner, 1);

C_inner_pid = pid(C_inner)
% step(T_inner)

%%
P_outer = minreal( (G_x * C_inner) / (1 + G_th * C_inner) );
% rltool(P_outer)
%%
C_outer = 0.0026139+ 0.018008*s + 0.0000364/s;
% C_outer = 0.00042394 + 3.635e-06/s + 0.012361*s;

T_cart_position = feedback(C_outer * P_outer, 1);
% step(T_cart_position);

C_outer_pid = pid(C_outer)

%% Discrete PLC Simulation with 5-Degree Initial Error

% 1. Discretize the Continuous Plant for 1ms PLC scan
Ts = 0.001; 
sys_c = ss(A_num, B_num, eye(4), zeros(4,1));
sys_d = c2d(sys_c, Ts, 'zoh'); % Zero-Order Hold discretization
Ad = sys_d.A; 
Bd = sys_d.B;

% Simulation Setup
t_sim = 50; 
N = t_sim / Ts;
t = 0:Ts:(t_sim-Ts);

% State vector: [x; dx; theta; dtheta]
x_state = [0; 0; 0; 0]; 

% Storage arrays for plotting
x_history = zeros(1, N);
th_history = zeros(1, N);
u_history = zeros(1, N);

% Initialize PLC Variables
e_x_old = 0; theta_ref_I_old = 0; theta_ref_D_old = 0;
e_theta_old = 0; u_I_old = 0; u_D_old = 0;
alpha = 0.01;

% Gains
outer_k_P = C_outer_pid.kp; outer_k_I = C_outer_pid.ki; outer_k_D = C_outer_pid.kd;
inner_k_P = C_inner_pid.kp; inner_k_I = C_inner_pid.ki; inner_k_D = C_inner_pid.kd;

% Execute 1ms Scan Cycle Loop
for k = 1:N
    if k == 5/Ts
        
    end
    % Read current state from the physical model
    x_pos = x_state(1);
    theta = x_state(3);
    
    % Log data
    x_history(k) = x_pos;
    th_history(k) = theta;
    
    % Outer Loop (Cart Position)
    e_x_new = 0 - x_pos;
    theta_ref_P = outer_k_P * e_x_new;
    theta_ref_I = theta_ref_I_old + (outer_k_I * e_x_new * Ts);
    theta_ref_D = (alpha * (outer_k_D * (e_x_new - e_x_old) / Ts)) + ((1.0 - alpha) * theta_ref_D_old);
    
    theta_ref = theta_ref_P + theta_ref_I + theta_ref_D;
    
    if theta_ref > 0.15
        theta_ref = 0.15;
    elseif theta_ref < -0.15
        theta_ref = -0.15;
    end
    
    % Inner Loop (Pendulum Angle)
    e_theta_new = theta_ref - theta; 
    
    u_P = inner_k_P * e_theta_new;
    u_I = u_I_old + (inner_k_I * e_theta_new * Ts);
    
    
    u_D = (alpha * (inner_k_D * (e_theta_new - e_theta_old) / Ts)) + ((1.0 - alpha) * u_D_old);
    
    u = (u_P + u_I + u_D);
    u = u + u*randn()*1.5;
    
    
    u_history(k) = u;
    
    % Variable Updates
    e_x_old = e_x_new;
    theta_ref_I_old = theta_ref_I;
    theta_ref_D_old = theta_ref_D;
    
    e_theta_old = e_theta_new;
    u_I_old = u_I;
    u_D_old = u_D;

    % add change in theta for disturbance
    
    x_state = Ad * x_state + Bd * u;
end

% Plot the Results
figure('Name', 'Inverted Pendulum 5-Degree Drop Simulation');

subplot(3,1,1);
plot(t, th_history * (180/pi), 'r', 'LineWidth', 1.5);
yline(0, 'k--');
title('Pendulum Angle (\theta)');
ylabel('Degrees');
grid on;

subplot(3,1,2);
plot(t, x_history, 'b', 'LineWidth', 1.5);
yline(0, 'k--');
title('Cart Position (x)');
ylabel('Meters');
grid on;

subplot(3,1,3);
plot(t, u_history, 'k', 'LineWidth', 1.5);
title('Motor Control Effort (Force)');
xlabel('Time (Seconds)');
ylabel('Force (Newtons)');
grid on;