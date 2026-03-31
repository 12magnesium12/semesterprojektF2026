%%matlabFunction
%%latex
clc;
clear;

%%change to real values
% I = 1;
% m = 1;
% M = 1;
% l = 1;
% g = 9.82;
% b_p = 1;
% b_c = 1;



syms I m_p m_c l g b_p b_c th x dth dx F_in
contVars = [I, m_p, m_c, l, g, b_p, b_c];
constValues = [1, 1, 1, 1, 9.82, 1, 1];

X_state = [th; x; dth; dx];

M_mat = [I+m_p*l^2 m_p*l*cos(th)
        m_p*l*cos(th) m_p+m_c ]

A_part = [m_p*g*l*sin(th) - b_p*dth;
          m_p*l*dth^2*sin(th) - b_c*dx];

f = [dth;
     dx;
     inv(M_mat) * A_part];

g_func = [0;
          0;
          inv(M_mat) * [0; 1]];

A_jacobian = jacobian(f, X_state);

% Equilibrium Point (Upright and stationary)
A_linearized = subs(A_jacobian, [th, x, dth, dx], [0, 0, 0, 0]);
B_linearized = subs(g_func,     [th, x, dth, dx], [0, 0, 0, 0]);

A = simplify(A_linearized)

B = simplify(B_linearized)

C = [1 0 0 0
     0 1 0 0]
D = [0; 0]

constVars = [I, m_p, m_c, l, g, b_p, b_c];
constValues = [1, 1, 1, 1, 9.82, 1, 1];

A = double(subs(A, constVars, constValues))
B = double(subs(B, constVars, constValues))

ss(A, B, C, D)