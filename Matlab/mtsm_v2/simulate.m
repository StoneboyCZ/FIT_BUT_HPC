%% PROBLEM
% y' = cos(t)/exp(t); y(0) = 1; TMAX = 1; h = 0.1
% 
% ODE
% y1' = y5 y1(0) = 1
% y2' = -y3 y2(0) = 1 cos(t)
% y3' = y2 y3(0) = 0 sin(t)
% y4' = y4 y4(0) = 1 e^t
% 
% DAE
% y5 = y2/y4;

%%
clc;
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
y0 = [1; 1; 1; 0; 1; 1];
dt = 0.1;
tmax = 1;
tspan = [0 tmax];

% dy(1) = y(5)*y(2);
% dy(2) = -y(3);
% dy(3) = y(2);
% dy(4) = y(4);
% dy(5) = -y(5)*y(5)*y(4);

ne = 6;

A = zeros(ne);
A(3,4) = -1;
A(4,3) = 1;
A(5,5) = 1;
ij = [
    3,6
    3,5
];
B2 = zeros(ne, size(ij,1));
B2(1,1) = 1;
B2(2,2) = 1;

ijk = [
    6,6,5;
];
B3 = zeros(ne, size(ijk,1));
B3(6,1) = -1;

ijkl = []; B4 = [];

ijklm = []; B5 = [];

b = zeros(ne,1);

maxORD = 63;
minORD = 10;
hScaleFactor = 1;
eps = 1e-9;

RUNS = 1000;

TIMES_MTSM_FK = zeros(1, RUNS);
for i=1:100
    tic;
    [T_MTSM_FK,Y_MTSM_FK,ORD] = explicitTaylorMult_GNPV_ver14(dt,tspan,y0,A,B2,B3,B4,B5,b,ij,ijk,ijkl,ijklm,eps,maxORD,minORD,hScaleFactor);
    TIMES_MTSM_FK(1,i)=toc;
end

% mean(TIMES_MTSM_FK)
% Y_MTSM_FK(1,:)

figure
plot(T_MTSM_FK, Y_MTSM_FK(1,:),'*')
grid on;
title('ERR - FK MTSM')


%% OH - one division, one multiplication
% dy(1) = y(7);
% dy(2) = y(6);
% dy(3) = -y(4);
% dy(4) = y(3);
% dy(5) = y(5);
% y(6) = y(3)*y(5)
% y(7) = y(3)/y(5)


y10 = 1; y20 = 1; y30 = 1; y40 = 0; y50 = 1;    
y0 = [y10; y20; y30; y40; y50;y30*y50;y30/y50];
dt = 0.1;
tmax = 1;
tspan = [0 tmax];

A = zeros(5, 7);
A(1,7) = 1;
A(2,6) = 1;
A(3,4) = -1;
A(4,3) = 1;
A(5,5) = 1;

b = zeros(5,1);

m = [
    3,5;
];

d = [
    3,5;
];


% indexes of linear part
index_l = 1:5;
index_m = 6;
index_d = 7;

maxORD = 63;
minORD = 10;
hScaleFactor = 1;
eps = 1e-9;

TIMES_MTSM_OH = zeros(1, RUNS);
for i=1:100
    tic;
    [T_MTSM_OH,Y_MTSM_OH,ORD] = taylor_v5(dt,tspan,y0,eps,A,b,m,d,index_l,index_m,index_d,maxORD,minORD,hScaleFactor);
    TIMES_MTSM_OH(1,i)=toc;
end
err = norm(Y_MTSM_FK(1,:) - Y_MTSM_OH(1,:))

%% OH - multiple equations


% y10 = 1; y20 = 1; y30 = 1; y40 = 0; y50 = 1; y60 = 1; y70 = 1;   
% y0 = [y10; y20; y30; y40; y50;y60;y70;y30*y50;y30*y50;y30/y50;y30/y50];
% dt = 0.1;
% tmax = 1;
% tspan = [0 tmax];
% 
% % dy(1) = y(10);
% % dy(2) = y(8);
% % dy(3) = -y(4);
% % dy(4) = y(3);
% % dy(5) = y(5);
% % dy(6) = y(11);
% % dy(7) = y(9);
% 
% % y(8) = y(3)*y(5)
% % y(9) = y(3)*y(5)
% % y(10) = y(3)/y(5)
% % y(11) = y(3)/y(5)
% 
% A = zeros(7, 11);
% A(1,10) = 1;
% A(2,8) = 1;
% A(3,4) = -1;
% A(4,3) = 1;
% A(5,5) = 1;
% A(6,11) = 1;
% A(7,9) = 1;
% 
% b = zeros(7,1);
% 
% m = [
%     3,5;
%     3,5;
% ];
% 
% d = [
%     3,5;
%     3,5
% ];
% 
% 
% % indexes of linear part
% index_l = 1:7;
% index_m = 8:9;
% index_d = 10:11;
% 
% maxORD = 63;
% minORD = 10;
% hScaleFactor = 1;
% eps = 1e-9;
%  
% TIMES_MTSM_OH = zeros(1, RUNS);
% for i=1:100
%     tic;
%     [T_MTSM_OH,Y_MTSM_OH,ORD] = taylor_v5(dt,tspan,y0,eps,A,b,m,d,index_l,index_m,index_d,maxORD,minORD,hScaleFactor);
%     TIMES_MTSM_OH(1,i)=toc;
% end
% mean(TIMES_MTSM_FK)/mean(TIMES_MTSM_OH)
% err = abs(Y_MTSM_FK(1,:) - Y_MTSM_OH(1,:))

%  Y_MTSM_OH(1,:)
% 
% figure
% plot(T_MTSM_OH, Y_MTSM_OH(2,:),'*')
% %plot(T_MTSM_OH, abs(Y_MTSM_FK(1,:)-Y_MTSM_OH(1,:)),'*')
% grid on;
% title('OH MTSM')

% Y_MTSM_OH(1,:)
% Y_MTSM_OH(2,:)

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
