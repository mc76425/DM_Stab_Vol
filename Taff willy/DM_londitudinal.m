clear all; close all; clc

%% Data
% Mass
g = 9.81;
m = 6214;
Ixx = 30048.7;
Iyy = 26639.8;
Izz = 55299.8;
Ixz = 606;

% Flight conditions
rho = 0.302;
V = 94.8;
q_ = (1/2)*rho*V^2;
M = 0.6;
theta0 = 2.2 * pi/180; % Make it radian
U0 = 94.8;

% Wing
S = 21.81;
A = 6;
e = 0.89;
c_ = 2.05;

CL0 = 0.813;
CD0 = 0.135;
Cm0 = 0;

CL_alpha = 5.22;
CLq = 3.9;
CL_sige = 0.34;
CL_sigT = 0;

Cmu = 0;
Cm_alpha = -0.401;
Cm_alphadot = -5;
Cmq = -10;
Cm_sige = -0.89;
Cm_sigT = 0;

CDu = 0;
CD_sige = 0;
CD_sigT = 0;



%% Simplifcation term
Xu = -(2*CD0 + CDu)*q_*S/(m*U0);
Xw = (CL0 - 2*CL_alpha*CL0/(pi*e*A))*q_*S/(m*U0);

Zu = -(2*CL0 + (CL0*M^2)/(1-M^2))*q_*S/(m*U0);
Zw = -(CD0 + CL_alpha)*q_*S/(m*U0);
Zq = (CLq)*c_*q_*S/(2*m*U0);

Mu = (Cmu)*c_*q_*S/(Iyy*U0);
Mw = (Cm_alpha)*c_*q_*S/(Iyy*U0);
Mwdot = (Cm_alphadot)*(c_^2)*q_*S/(2*Iyy*(U0^2));
Mq = (Cmq)*(c_^2)*q_*S/(2*Iyy*U0);

%% Aircraft's matrix
matrix_A = [Xu              Xw          0               -g*cos(theta0);
            Zu              Zw          (U0)            -g*sin(theta0);
            Mu+(Mwdot*Zu) Mw+(Mwdot*Zw) Mq+(Mwdot*U0)   -Mwdot*g*sin(theta0);
            0               0           1               0];

%% Simplification term for control matrix B
X_sige = CD_sige*q_*S/(m*U0);
X_sigT = CD_sigT*q_*S/(m*U0);

Z_sige = CL_sige*q_*S/(m*U0);
Z_sigT = CL_sigT*q_*S/(m*U0);

M_sige = Cm_sige*q_*S*c_/(Iyy*U0);
M_sigT = Cm_sigT*q_*S*c_/(Iyy*U0);

%% Control matrix
matrix_B = [X_sige                  X_sigT;
            Z_sige                  Z_sigT;
            M_sige + Mwdot*Z_sige   M_sigT + Mwdot*Z_sigT;
            0                       0];

%% Characteristic equation
sympref('FloatingPointOutput',true);
syms s
p = tf('p')
coef = poly(matrix_A);
equation = coef(1)*s^4 + coef(2)*s^3 + coef(3)*s^2 + coef(4)*s + coef(5);

eig_vals = eig(matrix_A);

%% Simulation
control = [0.1 ; 0.1];
initial_state = [0;0;0;0];


den = det(s*eye(4)-matrix_A);
num = adjoint(s*eye(4)-matrix_A);

%% Short period mode
om_nsp = sqrt(eig_vals(1)*eig_vals(2));
Eps_sp = -(eig_vals(1)+eig_vals(2))/(2*om_nsp);

G_sp = 1/(s^2 + 2*Eps_sp*om_nsp*s + om_nsp^2);
G_sp_lti = 1/(p^2 + 2*Eps_sp*om_nsp*p + om_nsp^2);
x_sp = ilaplace(G_sp);
matlabFunctionBlock('DM_simulink/ShortPeriod_t', x_sp, 'outputs',{'Response'})
%% Phugoid mode
om_np = sqrt(eig_vals(3)*eig_vals(4));
Eps_p = -(eig_vals(3)+eig_vals(4))/(2*om_np);

G_p = 1 / (s^2 + 2*Eps_p*om_np*s + om_np^2);
G_p_lti = 1 / (p^2 + 2*Eps_p*om_np*p + om_np^2);
x_p = ilaplace(G_p);
matlabFunctionBlock('DM_simulink/Phugoid_t', x_p, 'outputs',{'Response'})

%% Transfer Functions
U_root = double(root(num(1,1)));
w_root = root(num(2,2));
q_root = root(num(3,3));
theta_root = root(num(4,4));

T_u = -1/U_root(1);
om_u = sqrt(U_root(2)*U_root(3));
eps_u = -(U_root(2)+U_root(3))/(2*om_u);
G_u = (s + 1/T_u)*(s^2 + 2*eps_u*om_u*s + om_u^2);
G_u_lti = (p + 1/T_u)*(p^2 + 2*eps_u*om_u*p + om_u^2);
U_s = G_u*(G_p*G_sp);
U_s_lti = G_u_lti*G_p_lti*G_sp_lti;

figure(1)
rlocus(U_s_lti)
%%
T_w = -1/w_root(3);
om_w = sqrt(w_root(1)*w_root(2));
eps_w = -(w_root(1)+w_root(2))/(2*om_w);
G_w = (s + 1/T_w)*(s^2 + 2*eps_w*om_w*s + om_w^2);

G_q = s*(s-q_root(2))*(s-q_root(3));
%inverse = num/den
%phi_t = ilaplace(num/den)