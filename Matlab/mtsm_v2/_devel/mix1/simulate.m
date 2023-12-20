% y' = sin(t)*cos(t)*exp(t); y(0) = 1; TMAX = 1; h = 0.1
% 
% ODE
% y1' = y6 y1(0) = 1
% y2' = -y3 y2(0) = 1 cos(t)
% y3' = y2 y3(0) = 0 sin(t)
% y4' = y4 y4(0) = 1 e^t
% 
% DAE
% y5 = y2*y4;
% y6 = y3*y5;


clc;
close all;

tmax = 1;
dt = 0.1;
tspan = 0:dt:tmax;
options = odeset('RelTol',1e-13,'AbsTol',1e-15);


y0 = 1;

tic;
[T_ODE451,Y_ODE451] = ode45(@(t,y) ode(t,y),tspan,y0,options);
TIMES_ODE451 = toc;

Y_ODE451

% y1' = y6 y1(0) = 1
% y2' = -y3 y2(0) = 1 cos(t)
% y3' = y2 y3(0) = 0 sin(t)
% y4' = y4 y4(0) = 1 e^t
% 
% DAE
% y5 = y2*y4;
% y6 = y3*y5;

% ne = 4;
% A = zeros(ne, 6);
% A(1,6) = 1;
% A(2,3) = -1;
% A(3,2) = 1;
% A(4,4) = 1;
% 
% m = [
%     2,4;
% ];
% 
% d = [];
% 
% index_l = 1:4;
% index_m = 5;
% index_d = [];
% 
% b = zeros(4,1);
% 
% y10 = 1;
% y20 = 1;
% y30 = 0;
% y40 = 1;
% y50 = y20*y40;
% y60 = y30*y50;
% y0 = [y10;y20;y30;y40;y50;y60];
% 
% dt = 0.1;
% tspan = [0 tmax];
% minORD = 10;
% hScaleFactor = 1;
% maxORD = 64;
% 
% tic;
% [T_MTSM_OH,Y_MTSM_OH,ORD] = taylor_v51(dt,tspan,y0,eps,A,b,m,d,index_l,index_m,index_d,maxORD,minORD,hScaleFactor);
% TIMES_MTSM_OH=toc;
% 
% Y_MTSM_OH



function dy = ode(t,y)
    dy = zeros(1,1);
    dy(1) = (sin(t)*cos(t))/exp(t);
end
