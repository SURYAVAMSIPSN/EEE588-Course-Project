% EEE588 - Robust Multivariable Control
% NAME: SIVANAGA SURYA VAMSI POPURI
% ASU ID: 1217319207
% Control of a Car-Like Robot with Differential Robot Model.

% 3 Controllers to be used
% 1. LQR
% 2. H_Infinity
% 3. LQG/LTRO

clc
close all
clear all

% User can choose the required controller.
choice_of_loop = input('Enter the loop you would like to try out\nAvailable loops: 1, 2, 2.1, 3, 3.1. Your loop? : ')
d = input ('Enter the distance d: ')
%% The Model of the robot. 
% d = distance from the mid-point on the line between the rear wheels to 
% the center of gravity.

% Note that as long as d > 0, the system remains stable as the center of
% gravity is in front. For d < 0, the system becomes unstable. Also, the
% system is non-linear for d != 0. 

%  --- Turtlebot Waffle's Params used here.
m = 1.8; mc = 1.6; I = 0.979 ; Iw =  7.2880e-04; L = 0.153; R = 0.015; 
w_o = 0.5; v_o = 0.5; % Equilibrium v ad w points for linearization.
p_1 = 1 / (m + 2*Iw/(R^2)); p_2 = 1 /(I + (2*Iw*L^2)/R^2);

% --- The plant's linearized state space without any motor torque considerations (Nominal Plant). 
A = [0 2*mc*d*w_o*p_1; -mc*d*w_o * p_2 -mc*d*v_o * p_2] ; B = [p_1/R p_1/R; L * p_2/R -L * p_2/R]; C = [1 0; 0 1]; D = [0 0;0 0];

nominal_sys = ss(A,B,C,D);
P_Nominal = tf(nominal_sys); % This gives us the nominal plant. 

%  --- Actuator Dynamics - conversion from voltage to the torque of the wheel's motor. 
% Based on METS Robots.
Kt = 9.5*10^(-3); Kg = 9.68; Kb = Kt; B = 3.29*10^(-6); La = 2300 * 10^-6; Ra = 13.7;
s = tf([1 0],[1]);
P_actuator = Kt*Kg*(s*Iw + B)/((s*Iw + B)*(s*La + Ra) + Kb*Kt); P_actuator = zpk(P_actuator);

% Plant including Actuator Dyanmics input voltage to V, W
Plant = P_Nominal*[P_actuator 0; 0 P_actuator]; Plant = zpk(Plant)
Plantss = ss(Plant,'minimal');

% --- Natural Modes - given by the Poles (Eigenvalues) and the Eigenvectors
 
[eig_vec,eig_val] = eig(Plantss.A)  % evec contains eigenvectors
                        % eval contains poles or eigenvalues
%  --- Modal Analysis: We want to examine the natural modes (tendencies) of the system.
%
t_start = 0;
t_final  = 1000;
t_step  = 0.001;

t_vec     = [t_start:t_step:t_final];   % Vector of uniformly spaced time points
u     = [0*t_vec' 0*t_vec'];         % Input was set to zero for zero-input response (initial conditions);
                             
x = lsim(ss(Plantss.A, Plantss.B, Plantss.C, Plantss.D), u, t_vec, abs(eig_vec(:,2))); 
figure; plot(t_vec,x)
grid
title(': Initial Condition States')
ylabel('States (deg, deg/sec)')
xlabel('Time (seconds)')
pause

% --- Calculation of Transmission Zeroes ---
%
P_Transmission_Zeroes = tzero(ss(Plantss.A, Plantss.B, Plantss.C, Plantss.D));       
%
% There are no finite transmission zeroes for the System. 

%  --- Controllability --- 
% 
Ap = Plantss.A; Bp = Plantss.B; Cp = Plantss.C; Dp = Plantss.D;
M_C  = [Bp Ap*Bp (Ap^2)*Bp  (Ap^3)*Bp (Ap^4)*Bp (Ap^5)*Bp];
rank_MC = rank(M_C); % Rank of Controllability Matrix
%
% system is controllable; that is, rank of controllability matrix is 4 

%  --- Observability ---
%
M_O = [Cp; Cp * Ap; Cp * (Ap^2); Cp * (Ap^3); Cp * (Ap^4); Cp * (Ap^5)];
rank_MO = rank(M_O); % Rank of Observability Matrix
%
% system is observable; that is, rank of observability matrix is 4 

% --- Singular Value Analysis ---
w_start  = -1;
w_final   =  2;
no_of_points  = 20000;  
w      = logspace(w_start,w_final,no_of_points);   % Form vector of logarithmically spaced freq points
singular_value = sigma(ss(Ap, Bp, Cp, Dp),w);
singular_value = 20*log10(singular_value);
figure; semilogx(w, singular_value)

title('Outputs: linear velocity (m/s), angular velocity (rad/s); Inputs: left actuator voltage, right actuator voltage')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause

% PLANT SVD Analysis at w = 4.5299 rad/sec
%
wsvd       = j*0.001;
plant_s    = Cp*inv(wsvd*eye(6)-Ap)*Bp + Dp;
[u, singular_value, v ] = svd(plant_s)
%% Final State Space Representation form of Nominal Plant and Actuator Dynamics:  

% Coefficients for simplification. 
a1 = R*(1/p_1)*(B*Ra + Kb*Kt); a2 = R*(1/p_1)*(La*B + Iw*Ra); a3 = (La*R*(1/p_1)*Iw); b1 = Kt*Kg*Iw; Bb = Kt*Kg*B;
c1 =  R*(1/p_2)*(B*Ra + Kb*Kt)/L; c2 = R*(1/p_2)*(La*B + Iw*Ra)/L; c3 = (La*R*(1/p_2)*Iw/L); 
A_LQR = [0 1 0 0 0 0; ...
      0 0 1/a3 0 0 0; ...     
      0 -a1 -a2/a3 0 0 0; ...
      0 0 0 0 1 0; ...
      0 0 0 0 0 1/c3; ...
      0 0 0 0 -c1 -c2/c3];
B_LQR = [0 0; b1/a3 b1/a3; Bb-b1*a2/a3 Bb-b1*a2/a3; 0 0; b1/c3 -b1/c3; Bb-b1*c2/c3 -Bb+b1*c2/c3];
C_LQR = [1 0 0 0 0 0; 0 0 0 1 0 0]; D_LQR = [0 0; 0 0];
P_SS_final = ss(A_LQR,B_LQR,C_LQR,D_LQR);
P_final = tf(P_SS_final); P_final = zpk(P_final);


As = [(-2*B*Kg^2)/((m*R^2)+Iw) 0 Kt*Kg/(La*R*(1/p_1)) Kt*Kg/(La*R*(1/p_1))
       0   -2*L*L*Kg*Kg*B/(2*R*R*(1/p_2)) Kt*Kg*L/(La*R*(1/p_2)) -Kt*Kg*L/(La*R*(1/p_2))
       -Kb*Kg/R -Kb*Kg*L/R -Ra/La 0
       -Kb*Kg/R Kb*Kg*L/R 0 -Ra/La];
Bs = [0 0; 0 0; 1 0; 0 1]; Cs = [1 0 0 0; 0 1 0 0]; Ds = [0 0;0 0];
P_SS_final_S = ss(As,Bs,Cs,Ds);
P_final_S = minreal(tf(P_SS_final_S)); P_SS_final_S = ss(P_final_S); P_final_S = zpk(P_final_S);

P_Transmission_Zeroes = tzero(ss(As, Bs, Cs, Ds)); 

[eig_vec,eig_val] = eig(As)

wsvd       = j*0.000;
plant_s    = Cs*inv(wsvd*eye(4)-As)*Bs + Ds;
[u, singular_value, v] = svd(plant_s)

w_start  = -1;
w_final   =  20;
no_of_points  = 20000;  
w      = logspace(w_start,w_final,no_of_points);   % Form vector of logarithmically spaced freq points
singular_value = sigma(ss(As, Bs, Cs, Ds),w);
singular_value = 20*log10(singular_value);
figure; semilogx(w, singular_value)
title('Plant Singular Values Plot')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause
%% LQR Controller Design.
if (choice_of_loop == 1)
    disp('Great. You have chosen LQR. Press Enter to continue')
    pause
t_vec = [0:0.1:6000];
x = step(ss(Ap, Bp, Cp, 0*ones(2,2)), t_vec); % Step response for each state and each input

A_LQR = As; B_LQR = Bs; C_LQR = Cs; D_LQR = Ds;
Aa = [ 0*ones(2,2) C_LQR
      0*ones(4,2) A_LQR ];

Bb = [0*ones(2,2)
     B_LQR ];
 
q = diag([0.4, 0.8, 10, 12, 0, 0]);
r = diag ([0.03, 0.03]);
[g, k, CLPOLES] = lqr(Aa, Bb, q, r);


% --- Open Loop Dynamical System

% State x = [ z' xp' ]' = [ z'  y'  xr' ]'
%
% where 
%       z is the integrator state
%       y is the output (theta and flight path angle)
%       xr contains the rest of the state (pitch rate and speed)
%
gz = g(:,1:2);
gy = g(:,3:4);
gr = g(:,5:6);

A_ol = [ 0*ones(2,2) 0*ones(2,4)
       -B_LQR*gz   A_LQR-B_LQR*[0*ones(2,2) gr] ];

B_ol = [ - eye(2,2)
          B_LQR*gy        ];

C_ol = [ 0*ones(2,2) C_LQR ];
       
D_ol = 0*ones(2,2);

% ---Dynamics for the Closed Loop System
A_cl = A_ol - B_ol*C_ol;
B_cl = B_ol;
C_cl = C_ol;
D_cl = D_ol;
CLS = ss(A_cl,B_cl,C_cl,D_cl);
% ---CLOSED LOOP RESPONSE TO STEP COMMAND
t_vec = [0:0.02:10];
[y, t_vec, x] = step(CLS,t_vec);

% --- Reference command: Linear Velocity (r = [1 0])

plot(t_vec,y(:,:,1))
grid
title('Linear and angular velocity response for r = [1  0] references')
xlabel('Time (seconds)')
ylabel('linear (m/s) & angular (rad/s)')
legend('v_linear','w_angular')
pause
    
% --- Control Law: u =  gz(r -Cm*x) - gr*Crx

u10 = [-g gy]*[x(:,:,1)'
             ones(1, size(x(:,:,1)')*[0 1]')
             0*ones(1, size(x(:,:,1)')*[0 1]')];
plot(t_vec,u10)
grid
title('Controller Response for r = [1  0]  reference');
ylabel('Voltage (V)')
xlabel('Time (seconds)')
legend('Right actuator I/P','Left Actuator I/P')
pause

% --- Angular Velocity Reference (r = [0 1])

plot(t_vec,y(:,:,2))
grid
title('Linear and angular velocity response for r = [0  1] reference')
xlabel('Time (seconds)')
ylabel('linear (m/s) & angular (rad/s)')
legend('Linear Vel','Angular Vel')
pause

% --- Control Law: u =  [ -g gy ] [x' r']'

u01 = [-g gy]*[x(:,:,2)'
             0*ones(1, size(x(:,:,2)')*[0 1]')
             ones(1, size(x(:,:,2)')*[0 1]') ];
plot(t_vec,u01)
grid
title('Response of controller for r = [0  1] reference')
xlabel('Time (seconds)')
ylabel('Voltage (V)')
legend('I/P - right actuator','I/P - left actuator')
pause

% --- Response to both input references (r = [1 1])

[y, t_vec, x] = lsim(CLS,[ones(size(t_vec)) ones(size(t_vec))],t_vec);

plot(t_vec,y)
grid
title('Linear and angular velocity response for r = [1  1] reference')
xlabel('Time (seconds)')
ylabel('linear (m/s) & angular (rad/s)')
pause
    
% --- Control Law: u =  [ -g gy ] [x' r']'
u11 = [-g gy]*[x'
               ones(1, size(x')*[0 1]')
               ones(1, size(x')*[0 1]') ];

plot(t_vec,u11)
grid
title('Right and left voltage response for r = [1  1] reference')
xlabel('Time (seconds)')
ylabel('Voltage (V) & Voltage (V)')
pause

% --- Open Loop frequency response of LQ Controller.  

w = logspace(-3,3,100);
singular_value = sigma(ss(Aa, Bb, g, 0*ones(2,2)),w);
singular_value = 20*log10(singular_value);
semilogx(w, singular_value)
%clear sv
title('Open Loop Singular Values: Plant Input')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause
%return

singular_value = sigma(ss(A_ol, B_ol, C_ol, D_ol),w);
singular_value = 20*log10(singular_value);
semilogx(w, singular_value)

title('Open Loop Singular Values: Error Signal')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause

% --- Closed Loop Frequency response of LQ controller.  

singular_value = sigma(ss(Aa-Bb*g, Bb, -g, eye(2,2)-0*ones(2,2)),w);
singular_value = 20*log10(singular_value);
semilogx(w, singular_value)
%clear sv
title('LQ  Sensitivity: Plant Input')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause

singular_value = sigma(ss(A_cl, B_cl, -C_cl, eye(2,2)-D_cl)*P_final_S,w);
singular_value = 20*log10(singular_value);
semilogx(w, singular_value)
title('Disturbance at Input to Output Singular Value')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause


singular_value = sigma(ss(A_cl, B_cl, C_cl, D_cl),w);
singular_value = 20*log10(singular_value);
semilogx(w, singular_value)

title('LQs T (Comp Sens): Plant Output')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause

singular_value = sigma(ss(Aa-Bb*g, Bb, g, 0*ones(2,2)),w);
singular_value = 20*log10(singular_value);
semilogx(w, singular_value)
title('LQ Complementary Sensitivity: Plant Input')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause
end


%% H-infinity Controller Design.
if (choice_of_loop == 2)
    disp('Damn brother. H infinity controller it is. Hit the Enter button!')
    pause
P=ss(P_SS_final_S.A,P_SS_final_S.B,P_SS_final_S.C,P_SS_final_S.D); 
 
 % Weights
% Standard first order weights are chosen
% For any non-standard/higher-order weights, form the transfer functions
% appropiately
% See AAR's book for details on higher-order weights
M11=1.3; w11=5; Eps11=0.01; M12=1.3; w12=10; Eps12=0.01;
W1 = [tf([1/M11 w11], [1 w11*Eps11]) 0; 0 tf([1/M12 w12], [1 w12*Eps12])];

M21=100; w21=1000; Eps21=0.1; M22=100; w22=1000; Eps22=0.1;
W2 = [tf([1 w21/M21], [Eps21 w21]) 0; 0 tf([1 w22/M22], [Eps22 w22])] ;

M31=1.1; w31=3; Eps31=0.001; M32=1.1; w32=18; Eps32=0.001;
W3 = [tf([1 w31/M31], [Eps31 w31]) 0; 0 tf([1 w32/M32], [Eps32 w32])] ;

G=augw(P,W1,W2,W3);

% Obtain controller using Matlab's hinfsyn command
% See help on hinfsyn for more options
K=hinfsyn(G);

% --- Closed Loop Maps
% helper.m is a matlab function for computing OL and CL maps in a standard
% output feedback structure. This file must be in the current Matlab folder
[Lo,Li,So,Si,To,Ti,KS,PS] = helper(P,K);


%% PLOTS
wvec=logspace(-3,3,10000);
%**************************************************************************
%
% --- Frequency response for Open Loop 
figure; 
sigma(Lo,wvec);
title('Open Loop SVs: Error')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause

figure;
sigma(Li,wvec);
title('Open Loop SVs: plant input')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause 


%  --- Frequency response for closed loop.

figure; 
sigma(So,wvec);
title('S (Sensitivity): Error')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause 

figure; 
sigma(Si,wvec);
title('S (Sensitivity): plant input')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause 

figure; 
sigma(To,wvec);
title('T (Comp Sens): plant output')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause 

figure; 
sigma(Ti,wvec);
title('T (Comp Sens): plant input')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')


% --- Time response for closed loop 

t_vec = [0:0.02:5];
[y, t_vec, x] = step(To,t_vec);

% --- Linear Velocity response: r = [ 1  0 ] 
figure;
plot(t_vec,y(:,:,1))
grid
title('Step response for Linear Velocity reference')
ylabel('lin vel (m/s)')
xlabel('Time (seconds)')

% --- Control Law:
[u, t_vec] = step(KS,t_vec);
u10=u(:,:,1);

% --- Control for r = [ 1  0 ] 
figure;
plot(t_vec,u10)
grid
title('Response to Step Reference Command for Velocity')
ylabel('Controls (V)')
xlabel('Time (seconds)')


%  r = [ 0 1 ]
figure;
plot(t_vec,y(:,:,2))
grid
title('Response to Step Reference Command for Ang Velocity')
ylabel('Outputs (rad/s)')
xlabel('Time (seconds)')

% Form control:
u01 = u(:,:,2);
figure;
plot(t_vec,u01)
grid
title('Response to Step Reference Command for Ang Velocity')
ylabel('Controls (V)')
xlabel('Time (seconds)')
end

%% H-infinity sub-part 1

if (choice_of_loop == 2.1)
    disp('Okay. H_infinity sub part - 1. Press Enter to continue')
    pause
M=1.8;
Mc=1.6;
I=0.979;
Iw=7.2880e-04;
L=0.153;
R=0.015;
Kt=9.5*10^-3;
Kg=9.68;
Kb=Kt;
Beta=3.29*10^(-6);
La=2300*10^(-6);
Ra=13.7;

p_1=M+(2*Iw/(R^2));
p_2=I+(2*(L^2)*Iw)/(R^2);
A_initial= [(-2*Beta*Kg^2)/((M*R^2)+Iw) 0 Kt*Kg/(La*R*p_1) Kt*Kg/(La*R*p_1)
       0   -2*L*L*Kg*Kg*Beta/(2*R*R*p_2) Kt*Kg*L/(La*R*p_2) -Kt*Kg*L/(La*R*p_2)
       -Kb*Kg/R -Kb*Kg*L/R -Ra/La 0
       -Kb*Kg/R Kb*Kg*L/R 0 -Ra/La]
B_initial= [0 0; 0 0; 1 0; 0 1]
C_initial= [1 0 0 0; 0 1 0 0]
D_initial= [0 0;0 0]

P_initial = ss(A_initial,B_initial,C_initial,D_initial); 

P_initial_tf =tf(P_initial);
% A=P_initial.A;    % 4*4
% B=P_initial.B;    % 4*2
% C=P_initial.C;    % 2*4
% D=P_initial.D; 

% --- Plant transfer function post minimal realization
P = minreal(tf(P_initial_tf)); 
P_zpk = zpk(P);


% --- State Space of Plant
P_SS=ss(P);
A=P_SS.A   % 4*4
B=P_SS.B    % 4*2
C=P_SS.C    % 2*4
D=P_SS.D;
    M11=1; w11=1.05; Eps11=0.001; M12=1.1; w12=1.05; Eps12=0.001;
W1 = [tf([1/M11 w11], [1 w11*Eps11]) 0; 0 tf([1/M12 w12], [1 w12*Eps12])]

M21=200; w21=1000; Eps21=0.1; M22=100; w22=1000; Eps22=0.1;
W2 = [tf([1 w21/M21], [Eps21 w21]) 0; 0 tf([1 w22/M22], [Eps22 w22])] 

M31=2; w31=2.1; Eps31=0.001; M32=2; w32=2.1; Eps32=0.001;
W3 = [tf([1 w31/M31], [Eps31 w31]) 0; 0 tf([1 w32/M32], [Eps32 w32])] 

G=augw(P,W1,W2,W3);

K_hinf=hinfsyn(G);

[Lo,Li,So,Si,To,Ti,KS,PS] = f_CLTFM(P,K_hinf)

wvec=logspace(-3,3,10000);

% OPEN LOOP FREQUENCY RESPONSE 
figure; 
sigma(Lo,wvec);
title('Open Loop Singular Values  with Hinf compensator: Error Signal')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

figure;
sigma(Li,wvec);
title('Open Loop Singular Values with Hinf compensator: Plant Input')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

% CLOSED LOOP FREQUENCY RESPONSE 

figure; 
sigma(So*P,wvec);
title('Input Disturbance to Output Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

figure; 
sigma(Si,wvec);
title('Sensitivity Hinf: Plant Input')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

figure; 
sigma(To,wvec);
title('Complementary Sensitivity Hinf: Plant Output')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

figure; 
sigma(Ti,wvec);
title('Complementary Sensitivity Hinf: Plant Input')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

figure; 
sigma(KS,wvec);
title('Reference to Control Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

% CLOSED LOOP TIME RESPONSE 
t_vec = [0:0.02:10];
[y1, t_vec, x] = step(To,t_vec);

clppoles=(eig(To))
tzero(ss(To))
damp(clppoles)

% r = [ 1  0 ] 
figure;
plot(t_vec,y1(:,:,1))
grid
title('Lin Velocity & Ang Velocity Response To r = [1  0]  Command')
ylabel('v (m/s) & \omega (rad/s)')
xlabel('Time (seconds)')
legend('Lin Vel','Ang Vel')
pause

% Form control:
[u, t_vec] = step(KS,t_vec);
u10=u(:,:,1);

% Controls: r = [ 1  0 ] 
figure;
plot(t_vec,u10)
grid
title('Controller Response To r = [1  0]  Command');
ylabel('Volts (V)')
xlabel('Time (seconds)')
legend('Input to Right actuator','Input to Left Actuator')
pause

%  r = [ 0 1 ]
figure;
plot(t_vec,y1(:,:,2))
grid
title('Lin Velocity & Ang Velocity Response To r = [0  1]  Command')
ylabel('v (m/s) & \omega (rad/s)')
xlabel('Time (seconds)')
legend('Lin Vel','Ang Vel')
pause 

% Form control:
u01 = u(:,:,2);
figure;
plot(t_vec,u01)
grid
title('Controller Response To r = [0  1]  Command')
ylabel('Volts (V)')
xlabel('Time (seconds)')
legend('Input to Right actuator','Input to Left Actuator')
pause
end

%% LQG / LTRO Design
if (choice_of_loop == 3)
    
    disp('Hmm. LQG / LTRO. Interesting choice.')
Ap = As; Bp = Bs; Cp = Cs; Dp = Ds;

% Augment Plant with Integrators at Plant Input and Plot Singular Values
%
[NI NC] = size(Bp);                      % ns = number of inputs;  nc = number of controls;   
Aa = [ Ap             Bp
      0*ones(NC,NI)    0*ones(NC,NC) ]

Bb = [ 0*ones(NI,NC)
      eye(NC)      ]

c = [ Cp  0*ones(NC,NC) ]

d = 0*ones(NC,NC)

singular_value = sigma(ss(Aa, Bb, c, d),w);
singular_value = 20*log10(singular_value);
semilogx(w, singular_value)
%clear sv
title('Design Plant Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause

% Design of Target Loop Singular Values Using Kalman Filter
%
ll =  inv(Cp*inv(-Ap)*Bp + Dp);     % Choose ll and lh to match singular values at all frequencies
lh = -inv(Ap)*Bp*ll;
l = [lh 
     ll];                           % ll, lh - for low and high frequency loop shaping

singular_value = sigma(ss(Aa, l, c, d),w);
singular_value = 20*log10(singular_value);
semilogx(w, singular_value)
%clear sv
title('Filter Open Loop (G_{FOL}) Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause

pnint = eye(NC)                                    % Process Noise Intensity Matrix
mu = 0.01;                                         % Measurement Noise Intesity; Used to adjust Kalman Filter Bandwidth
                                                   % Small mu - expensive sensor   - large bandwidth
                                                   % Large mu - inexpensive sensor - small bandwidth
mnint = mu*eye(NC)                                 % Measurement Noise Intensity Matrix 
[kest, h, sig]= kalman(ss(Aa, [Bb l], c, [d 0*ones(NC,NC)]),pnint, mnint);  % Compute Filter Gain Matrix h
%[sig, poles, g1, rr] = care(a',c',l*l', mnint);                          % Alternate Method for Computing h
%h = g1';
singular_value = sigma(ss(Aa, h, c, d),w);
tsv = 20*log10(singular_value);
semilogx(w, tsv)
%clear sv
title('Target Loop (G_{KF}) Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause

singular_value = sigma(ss(Aa-h*c, h, c, 0*eye(NC)),w);
singular_value = 20*log10(singular_value);
semilogx(w, singular_value, w, 20*log10(10./w))
%clear sv
title('Target Complementary (T_{KF}) Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause

%
% Recover Target Loop By Solving Cheap LQR Problem
%
q = c'*c;                                            % State Weighting Matrix
rho = 1e-9;                                          % Cheap control recovery parameter;
                                                     % The smaller the parameter, the better the recovery.
r = rho*eye(NC)                                      % Control Weigthing Matrix
[k, poles, g, rr] = care(Aa,Bb,q,r);                   % Compute Control Gain Matrix G

%
% Compensator Analysis
%
A_k = [ Aa-Bb*g-h*c  0*ones(NI+NC,NC)
       g          0*ones(NC,NC) ]

B_K = [ h
       0*ones(NC,NC) ]

C_K = [0*ones(NC, NI+NC) eye(NC,NC) ]

Comp_poles = eig(A_k)                               % Compensator Poles
Comp_zeroes = tzero(Aa, h, g, 0*ones(NC,NC))         % Compensator Zeros
zerocheck = tzero(A_k, B_K, C_K, 0*ones(NC,NC))   % Check Compensator Zeros

singular_value = sigma(ss(A_k, B_K, C_K, 0*eye(NC)),w);
singular_value = 20*log10(singular_value);
semilogx(w, singular_value)
%clear sv
title('Compensator Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause

%
% Open Loop Analysis
%
A_L = [ Ap                     Bp*C_K
       0*ones(NI+NC+NC,NI)    A_k    ]

B_L = [ 0*ones(NI,NC)
       B_K ]
    
C_L = [ Cp  0*ones(NC,NI+NC+NC) ]
    
OL_poles = eig(A_L)                          % Open Loop Poles
OL_zeroes = tzero(A_L,B_L,C_L,0*ones(NC,NC))    % Open Loop Zeros
    
singular_value = sigma(ss(A_L, B_L, C_L, 0*eye(NC)),w);
singular_value = 20*log10(singular_value);
semilogx(w, singular_value, w, tsv)
%clear sv
title('Open Loop Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause   

% Closed Loop Analysis
%
CLPOLES = eig(A_L-B_L*C_L)           % Closed Loop Poles
CLP_Kalman = eig(Aa - h*c)              % Closed Loop Poles Due to Kalman Filter
CLP_Reg = eig(Aa - Bb*g)             % Closed Loop Poles Due to Regulator


singular_value = sigma(ss(A_L-B_L*C_L, B_L, -C_L, eye(NC)),w);
singular_value = 20*log10(singular_value);
semilogx(w, singular_value)
%clear sv
title('Sensitivity Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause 

end

%% LQG / LTRO Design - subpart (a). 
if (choice_of_loop ==3.1)
    disp('Had a tough day, huh? Sub part of LQR / LTRO it is.')
M=1.8;
Mc=1.6;
I=0.979;
Iw=7.2880e-04;
L=0.153;
R=0.0153;
Kt=9.5*10^-3;
Kg=9.68;
Kb=Kt;
Beta=3.29*10^(-6);
La=2300*10^(-6);
Ra=13.7;

p_1=M+(2*Iw/(R^2));
p_2=I+(2*(L^2)*Iw)/(R^2);
A_initial= [(-2*Beta*Kg^2)/((M*R^2)+Iw) 0 Kt*Kg/(La*R*p_1) Kt*Kg/(La*R*p_1)
       0   -2*L*L*Kg*Kg*Beta/(2*R*R*p_2) Kt*Kg*L/(La*R*p_2) -Kt*Kg*L/(La*R*p_2)
       -Kb*Kg/R -Kb*Kg*L/R -Ra/La 0
       -Kb*Kg/R Kb*Kg*L/R 0 -Ra/La]
B_initial= [0 0; 0 0; 1 0; 0 1]
C_initial= [1 0 0 0; 0 1 0 0]
D_initial= [0 0;0 0]

P_initial = ss(A_initial,B_initial,C_initial,D_initial); 

P_initial_tf =tf(P_initial);
% A=P_initial.A;    % 4*4
% B=P_initial.B;    % 4*2
% C=P_initial.C;    % 2*4
% D=P_initial.D; 

% Plant transferfunction after minimal realization
P = minreal(tf(P_initial_tf)); 
P_zpk = zpk(P);


%Plant StateSpace
P_SS=ss(P);
A=P_SS.A   % 4*4
B=P_SS.B    % 4*2
C=P_SS.C    % 2*4
D=P_SS.D;    % 2*2
NC=2;
no=2;
NI=4;
% Dynamic plant augmentation with integrator
%
                    % NI = number of inputs;  NC = number of controls;   
Aa = [ A            B
      0*ones(NC,NI)    0*ones(NC,NC) ]

Bb = [ 0*ones(NI,NC)
      eye(NC)      ]

c = [ C  0*ones(NC,NC) ]

d = 0*ones(NC,NC)


w = logspace(-2,3,100);
singular_value = sigma(ss(Aa, Bb, c, d),w);
singular_value = 20*log10(singular_value);
figure;
semilogx(w, singular_value)
title('Design Plant Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
% pause

% --- Design of Target Loop Singular Values Using Kalman Filter

ll =  inv(C*inv(-A)*B + D);     % Choose ll and lh to match singular values at all frequencies
lh = -inv(A)*B*ll;
l = [lh 
     ll];                           % ll, lh - for low and high frequency loop shaping

w_start  = -1;
w_final   =  2;
no_of_points  = 200;  
w      = logspace(w_start,w_final,no_of_points); 
singular_value = sigma(ss(Aa, l, c, d),w);
singular_value = 20*log10(singular_value);
figure;semilogx(w, singular_value)
%clear sv
title('Filter Open Loop (G_{FOL}) Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')




pnint = eye(NC)                                    % Process Noise Intensity Matrix
mu = 0.01; 
                                      % Measurement Noise Intesity; Used to adjust Kalman Filter Bandwidth
                                                   % Small mu - expensive sensor   - large bandwidth
                                                   % Large mu - inexpensive sensor - small bandwidth
mnint = mu*eye(NC)                                 % Measurement Noise Intensity Matrix 
[kest, h, sig]= kalman(ss(Aa, [Bb l], c, [d 0*ones(NC,NC)]),pnint, mnint);  % Compute Filter Gain Matrix h

singular_value = sigma(ss(Aa, h, c, d),w);
tsv = 20*log10(singular_value);
figure;semilogx(w, tsv)
%clear sv
title('Target Loop (G_{KF}) Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

tolpoles = eig(Aa)                           % Target Open Loop Poles
targzeros = tzero(Aa,h,c,0*ones(NC,NC))      % Target Open Loop Zeros
tclpoles = eig(Aa-h*c)                       % Target Closed Loop Poles


z  = 2.7; 

% Form Pre-Filter
%
fil    = ss(tf({z 0; 0 z}, {[1 z] 1; 1 [1 z]}));


singular_value = sigma(ss(Aa-h*c, h, -c, eye(NC)),w);
singular_value = 20*log10(singular_value);
figure;semilogx(w, singular_value)
%clear sv
title('Target Sensitivity (S_{KF}) Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')



singular_value = sigma(ss(Aa-h*c, h, c, 0*eye(NC)),w);
singular_value = 20*log10(singular_value);
figure;semilogx(w, singular_value, w, 20*log10(10./w))
%clear sv
title('Target Complementary (T_{KF}) Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')


% --- Target Closed Loop Singular Values
w_start  = -1;
w_final   =  2;
no_of_points  = 200;  
w = logspace(w_start,w_final,no_of_points);    % vector of logarithmically spaced freq points
[NTarA_CL, NTarB_CL, NTarC_CL, NTarD_CL ] = series(fil.a, fil.b, fil.c, fil.d, Aa-h*c,h,c,d);
singular_value     = sigma(ss(NTarA_CL, NTarB_CL, NTarC_CL, NTarD_CL),w);
singular_value     = 20*log10(singular_value);
figure;semilogx(w, singular_value)
title('Target Closed Loop Singular Values (r to y)')
grid
xlabel('Frequency (rad/sec)')
ylabel('T_{ry}, Singular Values (dB)')

% --- Target Step Responses
t_start = 0;
t_step  = 0.05;
t_final  = 10.0;
t_vec     = [t_start:t_step:t_final]'; % Vector of uniformly spaced time points
r1    = [ones(size(t_vec)) zeros(size(t_vec))];
r2    = [zeros(size(t_vec)) ones(size(t_vec))];
ty1   = lsim(ss(NTarA_CL, NTarB_CL, NTarC_CL, NTarD_CL), r1,t_vec);
ty2   = lsim(ss(NTarA_CL, NTarB_CL, NTarC_CL, NTarD_CL), r2,t_vec);

figure;plot(t_vec,ty1(:,1),'b', t_vec,ty1(:,2),'r')
grid
title('Target Responses to lin_vel Step Reference Command')
ylabel('lin_vel, ang_vel')
xlabel('Time (seconds)')



plot(t_vec,ty2(:,1),'b', t_vec,ty2(:,2),'r')
grid
figure;title('Target Responses to ang_vel Step Reference Command')
ylabel('lin_vel, ang_vel (degrees)')
xlabel('Time (seconds)')





% --- Recover Target Loop By Solving Cheap LQR Problem

q = c'*c;                                            % State Weighting Matrix
rho = 0.5e-04;                                       % Cheap control recovery parameter;
                                                     % The smaller the parameter, the better the recovery.
r = rho*eye(NC)                                      % Control Weigthing Matrix
[k, poles, g, rr] = care(Aa,Bb,q,r);                   % Compute Control Gain Matrix G



% --- Analysis of Compensator.

A_k = [ Aa-Bb*g-h*c  0*ones(NI+NC,NC)
       g          0*ones(NC,NC) ]

B_K = [ h
       0*ones(NC,NC) ]

C_K = [0*ones(NC, NI+NC) eye(NC,NC) ]

D_K=0*ones(NC,NC)
Comp_poles = eig(A_k)                               % Compensator Poles
Comp_zeroes = tzero(Aa, h, g, 0*ones(NC,NC))         % Compensator Zeros
zerocheck = tzero(A_k, B_K, C_K, 0*ones(NC,NC))   % Check Compensator Zeros

singular_value = sigma(ss(A_k, B_K, C_K, 0*eye(NC)),w);
singular_value = 20*log10(singular_value);
figure;semilogx(w, singular_value)
%clear sv
title('Compensator Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')



% --- Open Loop Analysis

A_L = [ A                     B*C_K
       0*ones(NI+NC+NC,NI)    A_k    ]

B_L = [ 0*ones(NI,NC)
       B_K ]
C_L = [ C  0*ones(NC,NI+NC+NC) ]
 
D_L=0*ones(NC,NC);
OL_poles = eig(A_L)                          % Open Loop Poles
OL_zeroes = tzero(A_L,B_L,C_L,0*ones(NC,NC))    % Open Loop Zeros
    
singular_value = sigma(ss(A_L, B_L, C_L, 0*eye(NC)),w);
singular_value = 20*log10(singular_value);
figure;semilogx(w, singular_value, w, tsv)

title('Open Loop Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')  

% --- Closed Loop System
%
A_cl     = A_L-B_L*C_L;
B_cl     = B_L;
C_cl     = C_L;
D_cl     = D_L;

% --- Closed Loop Analysis

CLPOLES = eig(A_L-B_L*C_L)           % Closed Loop Poles
CLP_Kalman = eig(Aa - h*c)              % Closed Loop Poles Due to Kalman Filter
CLP_Reg = eig(Aa - Bb*g)             % Closed Loop Poles Due to Regulator


singular_value = sigma(ss(A_L-B_L*C_L, B_L, -C_L, eye(NC)),w);
singular_value = 20*log10(singular_value);
figure;semilogx(w, singular_value)
%clear sv
title('Sensitivity: Error Signal')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

singular_value = sigma(ss(A_L-B_L*C_L, B_L, C_L, 0*eye(NC)),w);
singular_value = 20*log10(singular_value);
figure;semilogx(w, singular_value)
%clear sv
title('T (Comp Sens) : Plant Output')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

% --- Pre-filter to Closed Loop Try


[NA_CL, NB_CL, NC_CL, ND_CL ] = series(fil.a, fil.b, fil.c, fil.d, A_cl,B_cl,C_cl,D_cl);
singular_value = sigma(ss(NA_CL, NB_CL, NC_CL, ND_CL),w);
singular_value = 20*log10(singular_value);
figure; semilogx(w, singular_value)

title('T (Comp Sens) : Plant Output')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)') 

% --- Reference to control (r to u)

[A_ru, B_ru, C_ru, D_ru ] = series(A_cl, B_cl, -C_cl, eye(2), A_k, B_K, C_K, D_K);
[A_ru, B_ru, C_ru, D_ru ] = series(fil.a, fil.b, fil.c, fil.d, A_ru, B_ru, C_ru, D_ru );    
w_start                 = -1;
w_final                  = 2;
no_of_points                 = 200;
w                     = logspace(w_start,w_final,no_of_points);
singular_value  = sigma(ss(A_ru, B_ru, C_ru, D_ru),w);
singular_value  = 20*log10(singular_value);
figure; semilogx(w, singular_value)

title('Reference to Control Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')



% --- Input Disturbance to Output

[A_diy, B_diy, C_diy, D_diy ] = series(A,B,C,D, A_cl, B_cl, -C_cl, eye(2));
w_start  = -1;
w_final   =  2;
no_of_points  = 200;  
w      = logspace(w_start,w_final,no_of_points); 
singular_value  = sigma(ss(A_diy, B_diy, C_diy, D_diy),w);
singular_value  = 20*log10(singular_value);
figure; semilogx(w, singular_value)
%clear sv
title('Input Disturbance to Output Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')


% --- DC analysis of input disturbance.
%
s = j*7
tdy              = C_diy*inv(s*eye(size(A_diy))-A_diy)*B_diy +  D_diy;
[utdy stdy vtdy] = svd(tdy)


% --- Finally, the closed Loop Step Responses
y1    = lsim(ss(NA_CL,NB_CL, NC_CL, ND_CL), r1, t_vec); 
figure;plot(t_vec,y1)
grid
title('Lin Velocity & Ang Velocity Response To r = [1  0]  Command')
ylabel('v (m/s) & \omega (rad/s)')
xlabel('Time (seconds)')
legend('Lin Vel','Ang Vel')

figure; plot(t_vec,y1, t_vec, ty1)
legend('y','ty')
title('Target & Recovery target Comaprison for step reference Lin Velocity')
ylabel('v (m/s) & \omega (rad/s)')
xlabel('Time (seconds)')
legend('Lin Vel','Ang Vel','Target Lin Vel')

y2    = lsim(ss(NA_CL,NB_CL, NC_CL, ND_CL), r2, t_vec); 
figure;plot(t_vec,y2)
grid
title('Lin Velocity & Ang Velocity Response To r = [0  1]  Command')
ylabel('v (m/s) & \omega (rad/s)')
xlabel('Time (seconds)')
legend('Lin Vel','Ang Vel')

figure; plot(t_vec,y2, t_vec, ty2)
legend('target ang vel','ang vel','lin vel')
title('Target & Recovery target Comparison for step reference ang vel')
ylabel('lin_vel (m/s) & ang_vel (rad/s)')
xlabel('Time (seconds)')

u1    = lsim(ss(A_ru, B_ru, C_ru, D_ru), r1, t_vec); 
figure; plot(t_vec,u1)
grid
title('Controller Response for r = [1  0]  reference');
ylabel('Voltage (V)')
xlabel('Time (seconds)')
legend('I/P to Right actuator','I/P to Left Actuator')

u2    = lsim(ss(A_ru, B_ru, C_ru, D_ru), r2, t_vec); 
figure;plot(t_vec,u2)
grid
title('Controller Response for r = [0  1] reference');
ylabel('Voltage (V)')
xlabel('Time (seconds)')
legend('I/p to right actuator','i/p to left actuator')   
end  
    