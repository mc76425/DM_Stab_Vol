clc
clear all
format short
 
%% DATA
%Flight condition
altitude=40000*0.3048; air_density=0.000588*515.383;%(kg/m3)
U0=677*0.3048 ;cg=0.32;
theta0=0; theta= 2.7*(pi()/180); alpha0=0; alpha=theta;
 
gamma = 1.4 ; R = 287 ; T0 = 15+ 273.15 ; T=216.5;
a=sqrt(gamma*R*T);M = U0/a ; g=9.81;
 
%Geometry and Inertias
wing_area=230*((0.3048)^2); % in m²
wing_span=34*0.3048; % in m
c=7*0.3048; % in m 
 
%Inertial Data
W =  564000* 0.453592*9.81 ;%Weight [N]
m = 13000*0.453592; % Mass [kg]
IxxB = 28000* 1.3558 ;% Inertia [kg.m?]
IyyB = 18800* 1.3558; % Inertia [kg.m?]
IzzB = 47000* 1.3558 ;% Inertia [kg.m?]
IxzB = 1300* 1.3558 ;% Inertia [kg.m?]
 
%Steady state Coefficients
cL0=0.41; cD=0.0335; CtX=0.0335; Cm=0;CmT=0;
 
%Longitudinal Derivatives
Cmu=0.05; Cmalpha=-0.64; Cmalphap=-6.7; Cmq=-15.5; CmTu=-0.003; CmTalpha=0;
CLu=0.40;  Clalpha=5.84;  CLalphap=2.2;  CLq=4.7;  CDalpha=0.30;
CDu=0.104;  CTXu=0;  CLdeltaE=0.46;  CDdeltaE=0;  CmdeltaE=-1.24;
 
CdM=Cmu/M
CmM=CDu/M
 
CLdeltaT = 0 ; CDdeltaT = 0 ; CMdeltaT = 0 ;
 
%Others
cL= cL0 + Clalpha*alpha;
aspect_ratio=((wing_span)^2)/wing_area;
%e=(2*cL *Clalpha)/(pi*aspect_ratio*CDalpha ); %oswald factor
e=0.85  %oswald factor
qbar= (1/2)*air_density*((U0)*U0);
 
 
%% Question 2:
%Matrix A's coefficient
Xu=-((qbar*wing_area)/(m*U0))*(2*cD+M*CdM);  
Xw = ((qbar*wing_area)/(m*U0))*(cL0 -((2/(pi*e*aspect_ratio))*(cL0*Clalpha)));
Zu = -((qbar*wing_area)/(m*U0))*((2*cL0)+(((M^2)/(1-(M^2)))*cL0));
Zw = -((qbar*wing_area)/(m*U0))*(cD+Clalpha);
Zw_dot = ((qbar*wing_area*c)/(2*m*(U0^2)))*(cD+CLalphap); 
Zq = ((qbar*wing_area*c)/(2*m*(U0)))*(CLq);
Mu = ((qbar*wing_area*c)/(IyyB*(U0)))*(CmM*M);
Mw = ((qbar*wing_area*c)/(IyyB*(U0)))*(Cmalpha);
Mw_dot = ((qbar*wing_area*(c^2))/(2*IyyB*(U0^2)))*(Cmalphap);
Mq = ((qbar*wing_area*(c^2))/(2*IyyB*(U0)))*(Cmq);
 
%Matrix A
A=[Xu Xw 0 -g*cos(theta0); Zu Zw U0 -g*sin(theta0); 
    Mu+(Mw_dot*Zu) Mw+(Mw_dot*Zw) Mq+(U0*Mw_dot) -Mw_dot*g*sin(theta);
    0 0 1 0];
disp('Matrice A')
disp('A=')
disp (A);
 
%%Matrix B's coeficient
Xdeltae=((qbar*wing_area)/(m*U0))*CDdeltaE;
Zdeltae=(qbar*wing_area)/(m*U0)*CLdeltaE;
Mdeltae=((qbar*wing_area*c)/(IyyB*(U0)))*CmdeltaE;
Xdeltat=((qbar*wing_area)/(m*U0))*CDdeltaT;
Zdeltat=((qbar*wing_area)/(m*U0))*CLdeltaT;
Mdeltat=((qbar*wing_area*c)/(IyyB*(U0)))*CMdeltaT;
 
%Matrix B
B=[Xdeltae Xdeltat;Zdeltae Zdeltat;
    Mdeltae+(Zdeltae*Mw_dot) Mdeltat+Zdeltat*Mw_dot;
    0 0 ];
disp('Matrice B')
disp('B=')
disp(B);
 
%% Question 3 : Charasteristic equation     
disp('Characteristic equation')
 
syms x ;
 
X=A-eye(4).*x
Det=vpa(det(X))
 
%% Question 4: Eigen values of A
disp('Eigenvalues of A')
 
eigenvalues=eig(A);
disp(eigenvalues)
 
%eigenvect1=l(1,1);
%eigenvect2=l(2,2);
%eigenvect3=l(3,3);
%eigenvect4=l(4,4);
 
%% Question 5: 
 
%Short-period mode
disp('Short Period mode');
wnSH=sqrt(Zw*Mq-U0*Mw)
zetaSH=(-(Zw+Mq+U0*Mw_dot)/(2*wnSH))
 
%Phugoid mode
display('Phugoid mode')
wnPH=sqrt((-g/U0)*Zu)
zetaPH=(-Xu)/(2*wnPH)
 
%% Curves of longitudinal motion
display('Curves of Longitudinal Mode')
 
T1=(0:1:1000);
T2=(0:0.1:10);
 
    %Phugoide Mode    
%Pitch angle
%theta1=exp(wnPH*(-zetaPH+1i*sqrt(1-(zetaPH)^2))*T1);
theta1=(exp(-T1*zetaPH*wnPH)).*(0*cos(wnPH*sqrt(1-(zetaPH^2))))+((zetaPH*wnPH*10)/(wnPH*sqrt(1-(zetaPH^2)))*sin(wnPH*sqrt(1-(zetaPH^2))*T1));
%Pitch rate
%thetapoint1=theta0*wnPH*(-zetaPH+1i*sqrt(1-(zetaPH)^2))*(theta1/theta0); 
thetapoint1=(exp(-T1*zetaPH*wnPH)).*(0*cos(wnPH*sqrt(1-(zetaPH^2))*T1))+((zetaPH*wnPH)/(wnPH*sqrt(1-(zetaPH^2)))*sin(wnPH*sqrt(1-(zetaPH^2))*T1));
%Angle of attack
%alpha1=alpha0*exp(wnPH*(-zetaPH+1i*sqrt(1-(zetaPH)^2))*T1);
alpha1=(exp(-T1*zetaPH*wnPH)).*(theta*cos(wnPH*sqrt(1-(zetaPH^2))*T1))+((zetaPH*wnPH*theta)/(wnPH*sqrt(1-(zetaPH^2)))*sin(wnPH*sqrt(1-(zetaPH^2))*T1));
%Axial velocity
M1=(exp(-T1*zetaPH*wnPH)).*((U0*cos(T1*wnPH*sqrt(1-(zetaPH^2))))+((zetaPH*wnPH*U0)/(wnPH*sqrt(1-(zetaPH^2)))*sin(T1*wnPH*sqrt(1-(zetaPH^2)))))
 
figure(1)
plot(T1,M1,T1,theta1,T1,thetapoint1,T1,alpha1)
legend('Axial Velocity','Angle of attack (rad)','Pitch rate','Pitch angle (rad)');
grid on 
xlabel('Time(s)');
title('Phugoid Mode');
 
figure(11)
plot(T1,M1,'blue');
grid on 
xlabel('Time(s)');
title('Axial Velocity with Phugoid')
 
figure(12)
plot(T1,alpha1,'red')
grid on 
xlabel('Time(s)')
title('Angle of attack with Phugoid')
 
figure(13)
plot(T1,thetapoint1,'yellow')
grid on 
legend('Pitch rate')
title('Pitch rate with Phugoid')
 
figure(14)
plot(T1,theta1,'cyan')
grid on 
legend('Pitch angle')
title('Pitch angle with Phugoid mode')
 
 
    % Short Period Mode
%Pitch angle
theta2=theta*exp(wnSH*(-zetaSH+1i*sqrt(1-(zetaSH)^2))*T2);
%Pitch rate
thetapoint2=theta*wnPH*(-zetaPH+1i*sqrt(1-(zetaPH)^2))*(theta2/theta); 
%Angle of attack
alpha2=alpha*exp(wnSH*(-zetaSH+1i*sqrt(1-(zetaSH)^2))*T2);
%Axial velocity
M2=0.1*exp(wnSH*(-zetaSH+1i*sqrt(1-(zetaSH)^2))*T2); 
%M2=(exp(-T2*zetaSH*wnSH)).*((U0*cos(T2*wnSH*sqrt(1-(zetaSH^2))))+((zetaSH*wnSH*U0)/(wnSH*sqrt(1-(zetaSH^2)))*sin(T2*wnSH*sqrt(1-(zetaSH^2)))))
 
 
figure(2)
plot(T2,M2,T2,theta2,T2,thetapoint2,T2,alpha2)
legend('Axial Velocity','Angle of attack (rad)','Pitch rate','Pitch angle (rad)');
grid on
xlabel('Time s');
title('Short Period Mode');
 
figure(21)
plot(T2,M2,'blue');
grid on 
xlabel('Time(s)');
title('Axial Velocity with Short Period')
 
figure(22)
plot(T2,alpha2,'red')
grid on 
xlabel('Time(s)')
title('Angle of attack with Short Period')
 
figure(23)
plot(T2,thetapoint2,'yellow')
grid on 
legend('Pitch rate')
title('Pitch rate with Short Period')
 
figure(24)
plot(T2,theta2,'cyan')
grid on 
legend('Pitch angle')
title('Pitch angle with Short Period')




%% Transfer function of each variable
 
syms s ;
format short e
 
%U(s)/deltae (s)
X1=[s-Xu g; Zu/U0 s]
Det1=det(vpa(X1))
num1=[s*Xdeltae-g*Zdeltae];
TF_U=vpa(num1/Det1)
disp(TF_U)
disp('TF_U/deltae=')
 
%Theta(s)/deltae (s)
X1=[s-Xu g; Zu/U0 s]
Det1=det(vpa(X1))
num2=[-(Zu/U0)*Xdeltae+(s-Xu)*Zdeltae];
TF_theta=vpa(num2/Det1)
disp(TF_theta)
disp('TF_U/deltae=')
 
%Short-period Transfer Functions Approximations
%W(s)/deltae(s)
X2=[s-Zw -U0; -(Mw+Mw_dot*Zw) s-(Mq+Mw_dot*U0)]
Det2=det(vpa(X2))
num3=[(s-(Mq+Mw_dot*U0))*Zdeltae+U0*(Mdeltae+Mw_dot*Zdeltae)];
TF_W=vpa(num3/Det2)
disp(TF_W)
disp('TF_W/deltae=')
 
%q(s)/deltae(s)
X2=[s-Zw -U0; -(Mw+Mw_dot*Zw) s-(Mq+Mw_dot*U0)]
Det2=det(vpa(X2))
num4=[(Mw+Mw_dot*Zw)*Zdeltae+(s-Zw)*(Mdeltae+Mw_dot*Zdeltae)];
TF_q=vpa(num4/Det2)
disp(TF_q)
disp('TF_q/deltae =')
