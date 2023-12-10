%% PROBLEM
% y1' = cos(t)/e^t y1(0) = 1; 
%%

clc;
clear all;
close all;


%% original solution using ODE45
y0 = [1; 1; 0; 1];
dt = 0.1;
tmax = 1;
tspan = [0:dt:tmax];
options = odeset('RelTol',1e-13,'AbsTol',1e-15);

[T_ODE451,Y_ODE451] = ode45(@(t,y) orig(t,y),tspan,y0,options);

figure
plot(T_ODE451, Y_ODE451(:,1),'*');
grid on

%% FK
y0 = [1; 1; 0; 1; 1];
dt = 0.1;
tmax = 1;
tspan = [0 tmax];

% dy(1) = y(5)*y(2);
% dy(2) = -y(3);
% dy(3) = y(2);
% dy(4) = y(4);
% dy(5) = -y(5)*y(5)*y(4);

ne = 5;

A = zeros(ne);
A(2,3) = -1;
A(3,2) = 1;
A(4,4) = 1;
ij = [
    5,2
];
B2 = zeros(ne, size(ij,1));
B2(1,1) = 1;

ijk = [
    5,5,4;
];
B3 = zeros(ne, size(ijk,1));
B3(5,1) = -1;

ijkl = []; B4 = [];

ijklm = []; B5 = [];

b = zeros(ne,1);

maxORD = 63;
minORD = 10;
hScaleFactor = 1;
eps = 1e-9;

tic;
[T_MTSM_FK,Y_MTSM_FK,ORD] = explicitTaylorMult_GNPV_ver14(dt,tspan,y0,A,B2,B3,B4,B5,b,ij,ijk,ijkl,ijklm,eps,maxORD,minORD,hScaleFactor);
TIME_MTSM_FK=toc;

figure
plot(T_MTSM_FK, Y_MTSM_FK(1,:),'*')
grid on;
title('ERR - FK MTSM')

%% OH
y0 = [1; 1; 0; 1];
dt = 0.1;
tmax = 1;
tspan = [0 tmax];

% dy(1) = y(2)/y(4);
% dy(2) = -y(3);
% dy(3) = y(2);
% dy(4) = y(4);

ne = 4;

A = zeros(ne);
A(2,3) = -1;
A(3,2) = 1;
A(4,4) = 1;
ij = []; B2 = [];

ijk = []; B3 = [];

ijkl = []; B4 = [];

ijklm = []; B5 = [];

b = zeros(ne,1);

maxORD = 63;
minORD = 10;
hScaleFactor = 1;
eps = 1e-9;

tic;
[T_MTSM_OH,Y_MTSM_OH,ORD] = explicitTaylorMult_GNPV_ver14_div3(dt,tspan,y0,A,B2,B3,B4,B5,b,ij,ijk,ijkl,ijklm,eps,maxORD,minORD,hScaleFactor);
TIME_MTSM_OH=toc;

figure
plot(T_MTSM_OH, Y_MTSM_OH(1,:),'*')
grid on;
title('ERR - OH MTSM')

Y_MTSM_OH(1,:)



function dy = orig(~,y)
    dy = zeros(4,1);
    dy(1) = y(2)/y(4);
    dy(2) = -y(3);
    dy(3) = y(2);
    dy(4) = y(4);
end

function dy = FK(~,y)
    dy = zeros(5,1);
    dy(1) = y(5)*y(2);
    dy(2) = -y(3);
    dy(3) = y(2);
    dy(4) = y(4);
    dy(5) = -y(5)*y(5)*y(4);
end
