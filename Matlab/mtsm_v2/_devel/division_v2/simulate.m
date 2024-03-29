%% PROBLEM
% y1' = y1/y2 y1(0) = 1; 
% y2' = a*y2; y2(0) = 1;
%%

clc
close all



a = 2;
analfun = @(x) exp(-exp(-a*x)/a)*exp(1/a); 
tmax = 10;
dt = 0.1;
tspan = [0 tmax];
options = odeset('RelTol',1e-11,'AbsTol',1e-13);

%% FK division
% y0 = [1;1];
% [T_ODE451,Y_ODE451] = ode45(@(t,y) orig(t,y,a),tspan,y0,options);
% figure
% plot(T_ODE451,Y_ODE451(:,1));
% grid on;
y0 = [1;1;1];

tic;
[T_ODE45,Y_ODE45] = ode45(@(t,y) FK(t,y,a),tspan,y0,options);
TIME_ODE45 = toc;

figure
plot(T_ODE45, abs(analfun(T_ODE45)-Y_ODE45(:,1)),'*')
grid on;
title('ERR - FK ODE45')

%    dy(1) = y(3)*y(1);
%    dy(2) = a*y(2);
%    dy(3) = -a*y(3)*y(3)*y(2)

ne = 3;

A = zeros(ne);
A(2,2) = a;

ij = [
    3,1
];
B2 = zeros(ne, size(ij,1));
B2(1,1) = 1;

ijk = [
    3,3,2;
];
B3 = zeros(ne, size(ijk,1));
B3(3,1) = -a;


ijkl = []; B4 = [];

ijklm = []; B5 = [];

b = zeros(ne,1);

maxORD = 63;
minORD = 10;
hScaleFactor = 1;
eps = 1e-7;

tic;
[T_MTSM_FK,Y_MTSM_FK,ORD] = explicitTaylorMult_GNPV_ver14(dt,tspan,y0,A,B2,B3,B4,B5,b,ij,ijk,ijkl,ijklm,eps,maxORD,minORD,hScaleFactor);
TIME_MTSM_FK=toc;

figure
plot(T_MTSM_FK, abs(analfun(T_MTSM_FK)-Y_MTSM_FK(1,:)),'*')
grid on;
title('ERR - FK MTSM')

%% OH
y0 = [1;1;1];

tic;
[T_ODE45_OH,Y_ODE45_OH] = ode45(@(t,y) OH(t,y,a),tspan,y0,options);
TIME_ODE45 = toc;

figure
plot(T_ODE45_OH, abs(analfun(T_ODE45_OH)-Y_ODE45_OH(:,1)),'*')
grid on;
title('ERR - OH ODE45')

% K1 = 1/y(2);
% K2 = y(3);
% 
% dy = zeros(3,1);
% dy(1) = y(1)/y(2);
% dy(2) = a*y(2);

y0 = [1;1];
ne = 2;

A = zeros(ne);
A(2,2) = a;

ij = []; B2 = [];

ijk = []; B3 = [];

ijkl = []; B4 = [];

ijklm = []; B5 = [];

b = zeros(ne,1);

maxORD = 63;
minORD = 10;
hScaleFactor = 1;
eps = 1e-7;

tic;
[T_MTSM_OH,Y_MTSM_OH,ORD] = explicitTaylorMult_GNPV_ver14_div3(dt,tspan,y0,A,B2,B3,B4,B5,b,ij,ijk,ijkl,ijklm,eps,maxORD,minORD,hScaleFactor);
TIME_MTSM_OH=toc;

figure
plot(T_MTSM_OH, abs(analfun(T_MTSM_OH)-Y_MTSM_OH(1,:)),'*')
grid on;
title('ERR - OH MTSM')

function dy = FK(~,y,a)
    dy = zeros(3,1);
    dy(1) = y(3)*y(1);
    dy(2) = a*y(2);
    dy(3) = -y(3)*y(3)*a*y(2);
end

function dy = OH(~,y,a)
    K1 = 1/y(2);
    K2 = y(3);

    dy = zeros(3,1);
    dy(1) = y(3);
    dy(2) = a*y(2);
    dy(3) = K1*y(3)-K2*K1*a*y(2);
end