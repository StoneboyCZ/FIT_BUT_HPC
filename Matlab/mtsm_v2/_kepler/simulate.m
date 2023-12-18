%%% KEPLER PROBLEM
clc;
close all;

analfun = @(y1,y2,e) (y1+e).*(y1+e) + (y2.*y2)/(1-e*e) - 1; 

%% parameters
e = 0.75;
RUNS = 10;

tmax = 2*pi;
tspan = [0 tmax];
options = odeset('RelTol',1e-13,'AbsTol',1e-15);
dt = 0.025;

maxORD = 63;
minORD = 10;
hScaleFactor = 1;
eps = 1e-9;

ind=load('DY_indexes_maxORD_GN_ordered_all_60','DY_ij','DY_ijk','DY_ijkl','DY_ijklm');

%% TEST THE ORIGINAL FUNCTION

y10 = 1-e;
y20 = 0;
y30 = 0;
y40 = sqrt((1+e)/(1-e));
y0 = [y10; y20; y30; y40];


TIMES_ODE451 = zeros(1,RUNS);
for i=1:RUNS
    tic;
    [T_ODE451,Y_ODE451] = ode45(@(t,y) ode1(t,y),tspan,y0,options);
    TIMES_ODE451(1,i) = toc;
end

mean(TIMES_ODE451)
% 
% figure
% plot(Y_ODE451(:,1),Y_ODE451(:,2));
% grid on

% err = analfun(Y_ODE451(:,1),Y_ODE451(:,2),e)';

%%% AUX with square root removed
y10 = 1-e;
y20 = 0;
y30 = 0;
y40 = sqrt((1+e)/(1-e));
y50 = sqrt((y10*y10 + y20*y20)*(y10*y10 + y20*y20)*(y10*y10 + y20*y20));
y60 = sqrt((y10*y10 + y20*y20));
y0 = [y10; y20; y30; y40; y50; y60];

TIMES_ODE452 = zeros(1,RUNS);
for i=1:RUNS
    tic;
    [T_ODE452,Y_ODE452] = ode45(@(t,y) ode2(t,y),tspan,y0,options);
    TIMES_ODE452(1,i) = toc;
end

mean(TIMES_ODE452)

% figure
% plot(Y_ODE452(:,1),Y_ODE452(:,2));
% grid on
% 
% err = analfun(Y_ODE452(:,1),Y_ODE452(:,2),e)';

%%% AUX with division removed
y10 = 1-e;
y20 = 0;
y30 = 0;
y40 = sqrt((1+e)/(1-e));
y50 = sqrt((y10*y10 + y20*y20)*(y10*y10 + y20*y20)*(y10*y10 + y20*y20));
y60 = sqrt((y10*y10 + y20*y20));
y70 = 1/y50;
y80 = 1/y60;
y0 = [y10; y20; y30; y40; y50; y60; y70; y80];

ne = 8;

A = zeros(ne);
A(1,3) = 1;
A(2,4) = 1;

ij = [
    1,7
    2,7
];
B2 = zeros(ne, size(ij,1));
B2(3,1) = -1;
B2(4,2) = -1;

ijk = [
    6,1,3;
    6,2,4;
    1,3,8;
    2,4,8;    
];
B3 = zeros(ne, size(ijk,1));
B3(5,1) = 3;
B3(5,2) = 3;
B3(6,3) = 1;
B3(6,4) = 1;

ijkl = []; B4 = [];

ijklm = [
    6,1,3,7,7;
    6,2,4,7,7;
    1,3,8,8,8;
    2,4,8,8,8;
]; 
B5 = zeros(ne, size(ijk,1));
B5(7,1) = -3;
B5(7,2) = -3;
B5(8,3) = -1;
B5(8,4) = -1;

b = zeros(ne,1);

% OH method
ne = 14;

y10 = 1-e;
y20 = 0;
y30 = 0;
y40 = sqrt((1+e)/(1-e));
y50 = sqrt((y10*y10 + y20*y20)*(y10*y10 + y20*y20)*(y10*y10 + y20*y20));
y60 = sqrt((y10*y10 + y20*y20));
y70 = y10*y30;
y80 = y20*y40;
y90 = y60*y70; 
y100 = y60*y80; 
y110 = y10/y50;
y120 = y20/y50;
y130 = y70/y60;
y140 = y80/y90;
y01 = [y10; y20; y30; y40; y50; y60; y70; y80; y90; y100; y110; y120; y130; y140];

% dy(1) = y(3);
% dy(2) = y(4);
% dy(3) = -y(11);
% dy(4) = -y(12);    
% dy(5) = 3*y(9) + 3*y(10);
% dy(6) = y(13) + y(14);
% 
% y(7) = y(1)*y(3)
% y(8) = y(2)*y(4)
% y(9) = y(6)*y(7)
% y(10) = y(6)*y(8)
% y(11) = y(1)/y(5)
% y(12) = y(2)/y(5)
% y(13) = y(7)/y(6)
% y(14) = y(8)/y(6)

index_l = 1:6;
index_m = 7:8;
index_d = 11:14;

A1 = zeros(6,ne);
A1(1,3) = 1;
A1(2,4) = 1;
A1(3,11) = -1;
A1(4,12) = -1;
A1(5,9) = 3;
A2(5,10) = 3;
A2(6,13) = 1;
A2(6,14) = 1;

b1 = zeros(6,1);

m = [
    1,3;
    2,4;    
];

d = [
    1,5;
    2,5;
    7,6;
    8,6;
];


TIMES_ODE453 = zeros(1,RUNS);
TIMES_MTSM_VS = zeros(1,RUNS);
TIMES_MTSM_OPT = zeros(1,RUNS);
TIMES_MTSM_OH = zeros(1,RUNS);
for i=1:RUNS
    tic;
    [T_ODE453,Y_ODE453] = ode45(@(t,y) ode3(t,y),tspan,y0,options);
    TIMES_ODE453(1,i) = toc;

    tic
    [T_MTSM_VS,Y_MTSM_VS,ORD_VS] = explicitTaylorMult(dt,tspan,y0,A,B2,B3,B4,B5,b,ij,ijk,ijkl,ijklm,eps,ind,maxORD);
    TIMES_MTSM_VS(1,i) = toc;
    
    tic
    [T_MTSM,Y_MTSM,ORD] = explicitTaylorMult_GNPV_ver14(dt,tspan,y0,A,B2,B3,B4,B5,b,ij,ijk,ijkl,ijklm,eps,maxORD,minORD,hScaleFactor);
    TIMES_MTSM_OPT(1,i) = toc;

    tic;
    [T_MTSM_OH,Y_MTSM_OH,ORD] = taylor_v4(dt,tspan,y01,eps,A1,b1,m,d,index_l,index_m,index_d,maxORD,minORD,hScaleFactor);
    TIMES_MTSM_OH(1,i)=toc;
end

mean(TIMES_ODE453)
mean(TIMES_MTSM_VS)
mean(TIMES_MTSM_OPT)
mean(TIMES_MTSM_OH)

err45 = norm(analfun(Y_ODE453(:,1),Y_ODE453(:,2),e))
errvs = norm(analfun(Y_MTSM_VS(1,:),Y_MTSM_VS(2,:),e))
erropt = norm(analfun(Y_MTSM(1,:),Y_MTSM(2,:),e))
erroh = norm(analfun(Y_MTSM_OH(1,:),Y_MTSM_OH(2,:),e))

function dy = ode1(t,y)
    dy = zeros(4,1);

    r = sqrt(y(1)*y(1) + y(2)*y(2));
    r3 = r*r*r;
    
    dy(1) = y(3);
    dy(2) = y(4);
    dy(3) = -(y(1)/r3);
    dy(4) = -(y(2)/r3);
end

function dy = ode2(t,y)
    dy = zeros(6,1);
    dy(1) = y(3);
    dy(2) = y(4);
    dy(3) = -(y(1)/y(5));
    dy(4) = -(y(2)/y(5));    
    dy(5) = 3*y(6)*y(1)*y(3) + 3*y(6)*y(2)*y(4);
    dy(6) = (y(1)*y(3) + y(2)*y(4))/y(6);
end

function dy = ode3(t,y)
    dy = zeros(8,1);
    dy(1) = y(3);
    dy(2) = y(4);
    dy(3) = -y(1)*y(7);
    dy(4) = -y(2)*y(7);    
    dy(5) = 3*y(6)*y(1)*y(3) + 3*y(6)*y(2)*y(4);
    dy(6) = y(1)*y(3)*y(8) + y(2)*y(4)*y(8);
    dy(7) = -3*y(6)*y(1)*y(3)*y(7)*y(7) - 3*y(6)*y(2)*y(4)*y(7)*y(7);
    dy(8) = -y(1)*y(3)*y(8)*y(8)*y(8) - y(2)*y(4)*y(8)*y(8)*y(8);
end

function dy = ode4(t,y)
    dy = zeros(14,1);
    dy(1) = y(3);
    dy(2) = y(4);
    dy(3) = -y(1)*y(7);
    dy(4) = -y(2)*y(7);    
    dy(5) = 3*y(6)*y(9) + 3*y(6)*y(10);
    dy(6) = y(9)*y(8) + y(10)*y(8);
    dy(7) = -3*y(6)*y(9)*y(7)*y(7) - 3*y(6)*y(10)*y(7)*y(7);
    dy(8) = -y(9)*y(8)*y(8)*y(8) - y(10)*y(8)*y(8)*y(8);
    dy(9) = y(11) - y(12)*y(7);
    dy(10) = y(13) - y(14)*y(7);
    dy(11) = -2*y(9)*y(7);
    dy(12) = 2*y(9);
    dy(13) = -2*y(10)*y(7);
    dy(14) = 2*y(10);

end