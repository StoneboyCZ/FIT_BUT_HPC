%%% KEPLER PROBLEM
clc;
close all;

analfun = @(y1,y2,e) (y1+e).*(y1+e) + (y2.*y2)/(1-e*e) - 1; 

%% parameters
e = 0.75;
cycles=2; % number of rotations

RUNS = 100;


% number of steps for one rotation
% nSteps=20; % e=0.25
% nSteps=50; % e=0.5
nSteps=100; % e=0.75

tmax = cycles*2*pi;
tspan = [0 tmax];
options = odeset('RelTol',1e-13,'AbsTol',1e-15);

dt=(2*pi)/nSteps; % integration stepsize

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

% mean(TIMES_ODE451);
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

% mean(TIMES_ODE452);

% figure
% plot(Y_ODE452(:,1),Y_ODE452(:,2));
% grid on
% 
% err = analfun(Y_ODE452(:,1),Y_ODE452(:,2),e)';

%%% AUX from INFORMATICS
y10 = 1-e;
y20 = 0;
y30 = 0;
y40 = sqrt((1+e)/(1-e));
y50 = sqrt((y10*y10 + y20*y20)*(y10*y10 + y20*y20)*(y10*y10 + y20*y20));
y60 = sqrt((y10*y10 + y20*y20));
y70 = 1/y50;
y80 = 1/y60;
y90 = y10*y30 + y20*y40;
y100 = y60*y90;
y110 = y80*y90;
y120 = y10*y10+y20*y20;

%% FK definition
y0_orig = [y10; y20; y30; y40; y50; y60; y70; y80; y90; y100; y110; y120];

ne = 12;

A1 = zeros(ne);
A1(1,3) = 1;
A1(2,4) = 1;
A1(5,10) = 3;
A1(6,11) = 1;

ij = [
    1,7
    2,7
    3,3
    4,4
    7,12
    9,11
    1,3
    2,4
];
B2 = zeros(ne, size(ij,1));
B2(3,1) = -1;
B2(4,2) = -1;
B2(9,3) = 1;
B2(9,4) = 1;
B2(9,5) = -1;
B2(10,6) = 1;
B2(12,7) = 2;
B2(12,8) = 2;

% dy(1) = y(3);
% dy(2) = y(4);
% dy(3) = -y(7)*y(1);
% dy(4) = -y(7)*y(2);
% dy(5) = 3*y(10);
% dy(6) = y(11);





% dy(7) = -3*y(7)*y(7)*y(10);
% dy(8) = -y(8)*y(8)*y(11);
% dy(9) = y(3)*y(3) + y(4)*y(4) - y(7)*y(12);
% dy(10) = y(9)*y(11) + y(6)*y(3)*y(3) + y(6)*y(4)*y(4) - y(6)*y(7)*y(12);
% dy(11) = -y(8)*y(11)*y(11) + y(8)*y(3)*y(3) + y(8)*y(4)*y(4) - y(8)*y(7)*y(12);
% dy(12) = 2*y(1)*y(3) + 2*y(2)*y(4);

ijk = [
    7,7,10
    8,8,11
    6,3,3
    6,4,4
    6,7,12
    8,11,11
    8,3,3
    8,4,4
    8,7,12
];
B3 = zeros(ne, size(ijk,1));
B3(7,1) = -3;
B3(8,2) = -1;
B3(10,3) = 1;
B3(10,4) = 1;
B3(10,5) = -1;
B3(11,6) = -1;
B3(11,7) = 1;
B3(11,8) = 1;
B3(11,9) = -1;

ijkl = []; B4 = [];

ijklm = []; B5 = [];

b1 = zeros(ne,1);


%% OH definition
ne = 12;
rows = 31;

A = zeros(ne,rows);
A(1,3) = 1;
A(2,4) = 1;
A(3,13) = -1;
A(4,14) = -1;
A(5,10) = 3;
A(6,11) = 1;
A(7,23) = -3;
A(8,24) = -1;

A(9,15) = 1;
A(9,16) = 1;
A(9,17) = -1;

A(10,18) = 1;
A(10,25) = 1;
A(10,26) = 1;
A(10,27) = -1;

A(11,28) = -1;
A(11,29) = 1;
A(11,30) = 1;
A(11,31) = -1;

A(12,19) = 2;
A(12,20) = 2;

index_l = 1:12;
index_m1 = 13:22;
index_m2 = 23:31;
index_d = [];

m1 = [
    7,1;
    7,2;
    3,3;
    4,4;
    7,12;
    9,11;
    1,3;
    2,4;
    7,10;
    8,11
];

m2 = [
    7,21;
    8,22;
    6,15;
    6,16;
    6,17;
    11,22;
    8,15;
    8,16;
    8,17;
];

d = [];

% init
y10 = 1-e;
y20 = 0;
y30 = 0;
y40 = sqrt((1+e)/(1-e));
y50 = sqrt((y10*y10 + y20*y20)*(y10*y10 + y20*y20)*(y10*y10 + y20*y20));
y60 = sqrt((y10*y10 + y20*y20));
y70 = 1/y50;
y80 = 1/y60;
y90 = y10*y30 + y20*y40;
y100 = y60*y90;
y110 = y80*y90;
y120 = y10*y10+y20*y20;

y130 = y70*y10;
y140 = y70*y20;
y150 = y30*y30;
y160 = y40*y40;
y170 = y70*y120;
y180 = y90*y110;
y190 = y10*y30;
y200 = y20*y40;
y210 = y70*y100;
y220 = y80*y110;

y230 = y70*y210;
y240 = y80*y220;
y250 = y60*y150;
y260 = y60*y160;
y270 = y60*y170;
y280 = y110*y220;
y290 = y80*y150;
y300 = y80*y160;
y310 = y80*y170;

b = zeros(ne,1);

y0 = [y10; y20; y30; y40; y50; y60; y70; y80; y90; y100; y110; y120; y130; y140; y150; y160; y170; y180; y190; y200; y210;
    y220; y230; y240; y250; y260; y270; y280; y290; y300; y310];

TIMES_MTSM_OH = zeros(1, RUNS);
TIMES_ODE453 = zeros(1,RUNS);
TIMES_MTSM_FK = zeros(1, RUNS);
TIMES_MTSM_VS = zeros(1, RUNS);
for i = 1:RUNS
    tic;
    [T_ODE453,Y_ODE453] = ode45(@(t,y) ode4(t,y),tspan,y0_orig,options);
    TIMES_ODE453(1,i) = toc;
    
    tic;
    [T_MTSM_FK,Y_MTSM_FK,ORD_FK] = explicitTaylorMult_GNPV_ver14(dt,tspan,y0_orig,A1,B2,B3,B4,B5,b1,ij,ijk,ijkl,ijklm,eps,maxORD,minORD,hScaleFactor);
    TIMES_MTSM_FK(1,i)=toc;
    
    tic;
    [T_MTSM_VS,Y_MTSM_VS,ORD_VS] = explicitTaylorMult(dt,tspan,y0_orig,A1,B2,B3,B4,B5,b1,ij,ijk,ijkl,ijklm,eps,ind,maxORD);
    TIMES_MTSM_VS=toc;

    tic;
    [T_MTSM_OH,Y_MTSM_OH,ORD_OH] = taylor_v52(dt,tspan,y0,eps,A,b,m1,m2,d,index_l,index_m1,index_m2,index_d,maxORD,minORD,hScaleFactor);
    TIMES_MTSM_OH(1,i)=toc;
end

err = analfun(Y_ODE453(:,1),Y_ODE453(:,2),e)';
err = norm(analfun(Y_MTSM_OH(1,:),Y_MTSM_OH(2,:),e));
err = norm(analfun(Y_MTSM_FK(1,:),Y_MTSM_FK(2,:),e));

mean(TIMES_ODE453)
mean(TIMES_MTSM_OH)
mean(TIMES_MTSM_FK)
mean(TIMES_MTSM_VS)

ratio 


% ne = 12;
% 
% A = zeros(ne);
% A(1,3) = 1;
% A(2,4) = 1;
% 
% ij = [
%     1,7
%     2,7
% ];
% B2 = zeros(ne, size(ij,1));
% B2(3,1) = -1;
% B2(4,2) = -1;
% 
% ijk = [
%     6,1,3;
%     6,2,4;
%     1,3,8;
%     2,4,8;    
% ];
% B3 = zeros(ne, size(ijk,1));
% B3(5,1) = 3;
% B3(5,2) = 3;
% B3(6,3) = 1;
% B3(6,4) = 1;
% 
% ijkl = []; B4 = [];
% 
% ijklm = [
%     6,1,3,7,7;
%     6,2,4,7,7;
%     1,3,8,8,8;
%     2,4,8,8,8;
% ]; 
% B5 = zeros(ne, size(ijk,1));
% B5(7,1) = -3;
% B5(7,2) = -3;
% B5(8,3) = -1;
% B5(8,4) = -1;
% 
% b = zeros(ne,1);



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
    dy = zeros(12,1);
    dy(1) = y(3);
    dy(2) = y(4);
    dy(3) = -y(7)*y(1);
    dy(4) = -y(7)*y(2);
    dy(5) = 3*y(10);
    dy(6) = y(11);
    dy(7) = -3*y(7)*y(7)*y(10);
    dy(8) = -y(8)*y(8)*y(11);
    dy(9) = y(3)*y(3) + y(4)*y(4) - y(7)*y(12);
    dy(10) = y(9)*y(11) + y(6)*y(3)*y(3) + y(6)*y(4)*y(4) - y(6)*y(7)*y(12);
    dy(11) = -y(8)*y(11)*y(11) + y(8)*y(3)*y(3) + y(8)*y(4)*y(4) - y(8)*y(7)*y(12);
    dy(12) = 2*y(1)*y(3) + 2*y(2)*y(4);
end