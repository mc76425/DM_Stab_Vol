%Tittle: Assignment : Lateral Dynamic Stability of Airplane
%Course: Mécanique du vol - Qualité de vol (Aé 421)
%Name: NGY Alexis
%Group: MS2

%Nominal configuration
m = 7392 ;
Xcg = 0.07 ;
Ixx = 4480.8 ;
Iyy = 73999.5 ;
Izz = 75390.8 ;
g = 9.81  

%Flight conditions
H = 10668 ; %Altitude [m]
rho = 0.379 ; %Air Density [km/m^3]
M = 0.90 ; %Mach Number [1]
U0 = 267 ; %Velocity [m/s]
theta = 2.86*(pi()/180) ; %Initial attitude theta0 [deg]
q = (1/2)*(rho)*(U0)^2 

%Geometric Data
S = 18.22 ; %Wing Area [m²]
b = 6.69 ;  %Wing Span [m]
c = 2.91 ; %Wing Chord [m]
AR = 2.45 ; %Aspect ratio [1]
e = 0.92 ; %Oswald factor [1]

%Steady state conditions
CL0 = 0.735 ;
CD0 = 0.263 ;
CTx0 = 0.025 ;
Cm0 = 0 ;
Cmt0 = 0 ;

%Aerodynamic derivatives
Clb = -0.175 ;
Clp = -0.285 ;
Clr = 0.265 ;
CldeltaA = 0.039 ;
CldetaR = 0.045 ;
Cnb = -0.5 ;
Cnp = -0.14 ;
Cnr = -0.75 ;
CndeltaA = 0.0042 ;
CndeltaR = -0.16 ;
Cyb = -1.17 ;
Cyp = 0 ;
Cyr = 0 ;
CydeltaA = 0 ;
CydeltaR = 0.208 ;

%Question 1: 
%Matrix A
Yv = ((q*S)/(m*U0))*Cyb ;
Yp = ((q*S*b)/(2*m*U0))*Cyp ;
Yr = ((q*S*b)/(2*m*U0))*Cyr ;

Lv = ((q*S*b)/(Ixx*U0))*Clb ;
Lp = ((q*S*b^2)/(2*Ixx*U0))*Clp ; 
Lr = ((q*S*b^2)/(2*Ixx*U0))*Clr ;

Nv = ((q*S*b)/(Izz*U0))*Cnb ;
Np = ((q*S*b^2)/(2*Izz*U0))*Cnp ;
Nr = ((q*S*b^2)/(2*Izz*U0))*Cnr ;

%Question 2:
disp('Matrix A')
A = [Yv Yp (Yr-U0) g*cos(theta); Lv Lp Lr 0; Nv Np Nr 0; 0 1 0 0] 

%Question 3
S1 = poly(A)
disp(S1)

%Question 4
[u,v] = eig(A)
lambda1 = v(1,1)
lambda2 = v(2,2)
lambda3 = v(3,3)
lambda4 = v(4,4)

