clear;

syms I_p I_r m_p m_r m_c l g b_p b_c th(t) x(t) F_in(t) ddx ddth x_s dx_s th_s dth_s

% Energies and Lagrangian
E_kin_c = (1/2)*m_c*diff(x, t)^2;

dx_r = diff(x, t)+diff(th,t)*cos(th)*l/2;
dy_r = -diff(th, t)*sin(th)*l/2;
E_kin_r = (1/2)*(I_r*diff(th, t)^2 + m_r*dx_r^2 + m_r*dy_r^2);

dx_p = diff(x, t)+diff(th, t)*cos(th)*l;
dy_p = -diff(th, t)*sin(th)*l;
E_kin_p = (1/2)*(I_p*diff(th, t)^2 + m_p*dx_p^2 + m_p*dy_p^2);

E_kin = (E_kin_c + E_kin_r + E_kin_p);

E_pot_r = m_r*g*cos(th)*l/2;
E_pot_p = m_p*g*cos(th)*l;
E_pot = (E_pot_p + E_pot_r);

L = (E_kin - E_pot);

% Equations of Motion
EqTh = diff(diff(L, diff(th,t)), t) - diff(L, th) == -b_p*diff(th,t);
EqX  = diff(diff(L, diff(x,t)), t) - diff(L, x)  == -b_c*diff(x,t) + F_in;

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

% Numerical values
constVars = [l, m_p, m_r, m_c, I_p, I_r, g, b_p, b_c];
constVals = [0.35, 0.002, 0.082, 0.5, 0.002*0.35^2, (1/3)*0.082*0.35^2, 9.82, 0.0012, 5];

% Substitute equilibrium (position=0, velocities=0) and constants
A_num = double(subs(A_sym, [x_s, dx_s, th_s, dth_s, constVars], [0, 0, 0, 0, constVals]))
B_num = double(subs(B_sym, [x_s, dx_s, th_s, dth_s, constVars], [0, 0, 0, 0, constVals]))

C = [1 0 0 0
     0 0 1 0];
D = [0; 0]

sys = ss(A_num, B_num, C, D)
