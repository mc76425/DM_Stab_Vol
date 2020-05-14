%Longitudinal Dynamic Stability of Airplane: X-15 première vitesse
%Aé421
%Magdalena CALKA
%Copyrigth Beer version: Si tu te sers de mon travail paie une bière

clear all; clc;
%% Let's set up airplane config:
%Nomainal Configuration
m = 7063; %[kg]
Xcg = 0.22; %[without unit]
Ixx = 4951.7; %[kg.m^2]
Iyy = 108511; %[kg.m^2]
Izz = 111244; %[kg.m^2]
g = 9.81; %[m.s^(-2)] or [N/kg]

%Flight conditions
H = 12192; %Altitude [m]
rho = 0.302; %Air Density [km/m^3]
M = 0.80; %Mach Number [without unit]
U0 = 235.9; %Velocity [m/s]
theta0 = 7.5 ; %Initial attitude theta0 [deg]

%Geometric Data
S = 18.58; %Wing Area [m^2]
b = 6.82;  %Wing Span [m]
C = 3.13; %Wing Chord [m]
AR = 2.5; %Aspect ratio [without unit]
e = 0.93; %Oswald factor [without unit]

%Steady state conditions [without unit]
CL0 = 0.42; CD0 = 0.13;
CTx0 = 0.025; Cm0 = 0; Cmt0 = 0;

%Aerodynamic derivatives [without unit]
Cmu = -0.1; Cmalpha = -1; Cmalphad = -1.6;
Cmq = -12.5; Cmtu = 0;Cmtalpha = 0;

CLu = 0; CLalpha = 4; CLalphad = 0; CLq = 0;

CDu = 0; CDalpha = 0.8; CDalphad = 0; CDq = 0;

CTu = 0; CLdeltae = 1.1; CDdeltae = 0; Cmdeltae = -1.5;

%And the calculation of q to reduce expressions of equations
q = (1/2)*(rho)*(U0)^2;

%% QUESTION 1
%Elements needed to construct the A matrix
Xu = -q*S/(m*U0)*(2*CD0+CDu);
Xw = q*S/(m*U0)*(CL0-2/(pi*e*AR)*CL0*CLalpha);

Zu = -q*S/(m*U0)*(2*CL0 + M^2/(1-M^2)*CL0);
Zw = -q*S/(m*U0)*(CD0+CLalpha);
Zq = q*S*C/(2*m*U0)*CLq;

Mu = q*S*C/(Iyy*U0)*M*Cmu;
Mw = q*S*C/(Iyy*U0)*Cmalpha;
Mwd = q*S*C^2/(2*Iyy*U0^2)*Cmalphad;
Mq = q*S*C^2/(2*Iyy*U0)*Cmq;

%Elements needed to construct the B matrix
CMdeltae = Cmdeltae/M;
CLdeltaT = 0; CDdeltaT = 0; CMdeltaT = 0;
%deltae: elevator deflection
Xdeltae = q*S/(m*U0)*CDdeltae;
Zdeltae= q*S/(m*U0)*CLdeltae;
Mdeltae = q*S*C/(Iyy*U0)*CMdeltae;
XdeltaT = q*S/(m*U0)*CDdeltaT;
ZdeltaT = q*S/(m*U0)*CLdeltaT;
MdeltaT = q*S*C/(Iyy*U0)*CMdeltaT;


%% QUESTION2
A = [Xu Xw 0 -g*cos(deg2rad(theta0));
    Zu Zw U0 -g*sin(deg2rad(theta0));
    Mu+Mwd*Zu Mw+Mwd*Zw Mq+U0*Mwd -Mwd*g*sin(deg2rad(theta0));
    0 0 1 0];

B = [Xdeltae XdeltaT;
    Zdeltae ZdeltaT;
    Mdeltae+Mwd*Zdeltae MdeltaT+Mwd*ZdeltaT;
    0 0];

%% QUESTION3
%CS: Characteristic equation
CS = poly(A);

%% QUESTION4
eig_vals = eig(A);

%% QUESTION5 
%Short period mode
wnSP = sqrt(eig_vals(1)*eig_vals(2)); %Natural Frequency
% or wnSP = sqrt(Zw*Mq-U0*Mw); 
zetaSP = -(eig_vals(1)+eig_vals(2))/(2*wnSP); %Damping Factor
% or -(Zw+Mq+U0*Mwd)/(2*wnSP); 

%Phugoid mode
wnPM = sqrt(eig_vals(3)*eig_vals(4));%Natural Frequency
zetaPM  = -(eig_vals(3)+eig_vals(4))/(2*wnPM);%Damping Factor
%wnPM = sqrt(-g/(U0*Zu)); %Natural Frequency
%zetaPM = -Xu/(2*wnPM); %Damping Factor

%% QUESTION6 Phugoid mode
h = 0.1;
t=(0:h:1000);
%Calculation of needed parameters according time:
%Pitch angle
%With initial conditions
a1pitch = deg2rad(theta0);
a2pitch = zetaPM*wnPM*10/(wnPM*sqrt(1 - zetaPM^2));
thetaPM = exp(-zetaPM*wnPM*t).*(a1pitch*cos(wnPM*sqrt(1-zetaPM^2)*t) + ...
    a2pitch*sin(wnPM*sqrt(1-wnPM^2)*t));
%Pitch rate, basicly the dervative of pitch angle
thetadPM = diff(thetaPM)/h;
%Angle of attack
alphaPM = thetaPM + deg2rad(theta0);
%Velocity (Mach number)
%with initial conditions
a1v = U0;
a2v = zetaPM*wnPM*U0/(wnPM*sqrt(1-zetaPM^2));
MPM = exp(-zetaPM*wnPM*t).*(a1v*cos(wnPM*sqrt(1-zetaPM^2)*t) + ...
    a2v*sin(wnPM*sqrt(1-zetaPM^2)*t));

figure(1)
plot(t,MPM, t, thetaPM, t(:,1:length(thetadPM)), thetadPM, t, alphaPM)
legend('axial velocity','Pitch angle (rad)', ...
    'Pitch rate (rad)','angle of attack (rad)');
grid on 
xlabel('Time(s)');
title('Phugoid Mode');

%% QUESTION6 Short period mode
h = 0.1;
tt=(0:h:200);
%Calculation of needed parameters according time:
a2SP = 0.1/(wnSP*sqrt(1-zetaSP^2)); 
%Pitch rate
thetadSP = exp(-zetaSP*wnSP*tt).*a2SP.*sin(wnSP*sqrt(1-zetaSP^2)*tt);
%Pitch angle
thetaSP = cumtrapz(tt,thetadSP); % The area under the curve
%Angle of attack
alphaSP = thetaSP + deg2rad(theta0);
%Axial velocity
a1vSP = U0;
a2vSP = zetaSP*wnSP*U0/(wnSP*sqrt(1-zetaSP^2));
MSP = exp(-zetaSP*wnSP*tt).*(a1vSP*cos(wnSP*sqrt(1-zetaSP^2)*tt) + ...
    a2vSP*sin(wnSP*sqrt(1-zetaSP^2)*tt));

figure(2)
plot(tt,thetadSP, tt, thetaSP, tt, alphaSP)
legend('Pitch rate (rad)', 'Pitch angle (rad)',...
    'angle of attack (rad)');
grid on 
xlabel('Time(s)');
title('Short Period mode');
commandwindow

%% Transfert functions
syms s ;

XPM=[s-Xu g; Zu/U0 s];
DetPM = det(vpa(XPM));

XSP = [s-Zw -U0; -(Mw+Mwd*Zw) s-(Mq+Mwd*U0)];
DetSP = det(vpa(XSP));

%U(s)/deltae(s)
numeratorU=s*Xdeltae-g*Zdeltae;
TF_U = vpa(numeratorU/DetPM)

%Theta(s)/deltae(s)
numeratorT=-(Zu/U0)*Xdeltae+(s-Xu)*Zdeltae;
TF_t = vpa(numeratorT/DetPM)

%W(s)/deltae(s)
numeratorW=(s-(Mq+Mwd*U0))*Zdeltae+U0*(Mdeltae+Mwd*Zdeltae);
TF_W = vpa(numeratorW/DetSP)

%q(s)/deltae(s)
numeratorQ = (Mw+Mwd*Zw)*Zdeltae+(s-Zw)*(Mdeltae+Mwd*Zdeltae);
TF_q = vpa(numeratorQ/DetSP)