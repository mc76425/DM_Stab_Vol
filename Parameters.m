%Let's set up airplane config:
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
theta0 = 7.5 ; %Initial attitude theta0 [deg] DON'T FORGET: NEED RAD

%Geometric Data
S = 18.58; %Wing Area [m^2]
b = 6.82;  %Wing Span [m]
C = 3.13; %Wing Chord [m]
AR = 2.5; %Aspect ratio [without unit]
e = 0.93; %Oswald factor [without unit]

%Steady state conditions [without unit]
CL0 = 0.42; CD0 = 0.13; CTx0 = 0.025; Cm0 = 0; Cmt0 = 0;

%Aerodynamic derivatives [without unit]
Cmu = -0.1; Cmalpha = -1; Cmalphad = -1.6; Cmq = -12.5; Cmtu = 0; Cmtalpha = 0;

CLu = 0; CLalpha = 4; CLalphad = 0; CLq = 0;

CDu = 0; CDalpha = 0.8; CDalphad = 0; CDq = 0;

CTu = 0; CLdeltae = 1.1; CDdeltae = 0; Cmdeltae = -1.5;

%And the calculation of q to reduce expressions of equations
q = (1/2)*(rho)*(U0)^2;