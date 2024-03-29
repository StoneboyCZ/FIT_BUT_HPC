% Problem: y' = (sin(t)cos(t)/e^t)

%% init

clc
close all

tmax = 2*pi;
dt = pi/10;

%% SOLUTION OF THE ORIGINAL PROBLEM
y0 = 1;
tspan = [0 tmax];
options = odeset('RelTol',1e-11,'AbsTol',1e-13);

tic;
[T_ODE45,Y_ODE45] = ode45(@(t,y) ode(t,y),tspan,y0,options);
TIME_ODE45 = toc;

figure
plot(T_ODE45,Y_ODE45(:,1));
grid on;

%% VARIANT 1: FK
% y_1' &= y_5y_2y_3 & y_1(0) &= 5 \\
% y_2' &= y_3  & y_2(0) &= \sin(0)\\
% y_3' &= -y_2 & y_3(0) &= \cos(0) \\
% y_4' &= y_4 & y_4(0) &= e^0 \\
% y_5' & = -y_5y_5y_4 & y_5(0) &= \frac{1}{y_4(0)}

% ode solvers
y0 = [1;sin(0);cos(0);1;1]; 
tspan = [0 tmax];
options = odeset('RelTol',1e-11,'AbsTol',1e-13);

tic;
[T_ODE45_FK,Y_ODE45_FK] = ode45(@(t,y) odeFK(t,y),tspan,y0,options);
TIME_FK_ODE45 = toc;

figure
plot(T_ODE45_FK,Y_ODE45_FK(:,1));
grid on;

ne = 5;
A = zeros(ne);
A(2,3) = 1;
A(3,2) = -1;
A(4,4) = 1;

ij = [];
B2 = [];

ijk = [
    5,2,3;
    5,5,4;
];
B3 = zeros(ne, size(ijk,1));
B3(1,1) = 1;
B3(5,2) = -1;

ijkl = [];
B4 = [];

ijklm = [];
B5 = [];

b = zeros(ne,1);


maxORD = 63;
minORD = 10;
hScaleFactor = 1;
eps = 1e-7;

tic;
[T_MTSM_FK,Y_MTSM_FK,ORD] = explicitTaylorMult_GNPV_ver14(dt,tspan,y0,A,B2,B3,B4,B5,b,ij,ijk,ijkl,ijklm,eps,maxORD,minORD,hScaleFactor);
T_MTSM=toc;

% figure
% plot(T_MTSM_FK,Y_MTSM_FK(1,:));
% grid on;

% error in last step
err = abs(Y_MTSM_FK(:,end)-Y_ODE45_FK(end,:)');


%% VARIANT 2: OH
%     K1 = 1/y(4);
%     K2 = (y(2)*y(3))/y(4);
% 
%     dy = zeros(5,1);
%     
%     dy(1) = y(5);
%     dy(2) = y(3);
%     dy(3) = -y(2);
%     dy(4) = y(4);
%     dy(5) = K1*y(3)*y(3) - K1*y(2)*y(2) - K2*K1*y(4);


% ode solvers
y10 = 1;
y20 = sin(0);
y30 = cos(0);
y40 = 1;
y50 = (sin(0)*cos(0))/1;

y0 = [y10;y20;y30;y40;y50]; 
tspan = [0 tmax];
options = odeset('RelTol',1e-11,'AbsTol',1e-13);

TIMES_OH_ODE45 = [];
for i=1:100
    tic;
    [T_ODE45_OH,Y_ODE45_OH] = ode45(@(t,y) odeOH(t,y),tspan,y0,options);
    TIMES_OH_ODE45(i) = toc;
end
median(TIMES_OH_ODE45)

% comparison between systems -- equal
err = abs(Y_ODE45_OH(end,1)-Y_ODE45_FK(end,1))

figure
plot(T_ODE45_OH,Y_ODE45_OH(:,1));
grid on;

ne = 5;
A = zeros(ne);
A(1,5) = 1;
A(2,3) = 1;
A(3,2) = -1;
A(4,4) = 1;

ij = [];
B2 = [];

ijk = [];
B3 = [];

ijkl = [];
B4 = [];

ijklm = [];
B5 = [];

b = zeros(ne,1);

maxORD = 63;
minORD = 10;
hScaleFactor = 1;
eps = 1e-7;

tic;
[T_MTSM_OH,Y_MTSM_OH,ORD] = explicitTaylorMult_GNPV_ver14_div2(dt,tspan,y0,A,B2,B3,B4,B5,b,ij,ijk,ijkl,ijklm,eps,maxORD,minORD,hScaleFactor);
T_MTSM=toc;

figure
plot(T_MTSM_OH,Y_MTSM_OH(1,:));
grid on;





function dy = ode(t, y)
    dy = zeros(1,1);
    dy(1) = (sin(t)*cos(t))/exp(t);
end

function dy =  odeFK(~, y)
    dy = zeros(5,1);
    dy(1) = y(5)*y(2)*y(3);
    dy(2) = y(3);
    dy(3) = -y(2);
    dy(4) = y(4);
    dy(5) = -y(5)*y(5)*y(4);
end

function dy = odeOH(t, y)
    K1 = 1/y(4);
    K2 = (y(2)*y(3))/y(4);

    dy = zeros(5,1);
    
    dy(1) = y(5);
    dy(2) = y(3);
    dy(3) = -y(2);
    dy(4) = y(4);
    dy(5) = K1*y(3)*y(3) - K1*y(2)*y(2) - K2*K1*y(4);
end