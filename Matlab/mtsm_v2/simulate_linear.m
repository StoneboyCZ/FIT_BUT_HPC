clc
close all

ne = 2;
A = zeros(ne);
A(1,2) = 1;
A(2,1) = -1;

b = zeros(ne,1);

y0 = [0; 1];
dt = 0.1;
tmax = 2*pi;
options = odeset('RelTol',1e-13,'AbsTol',1e-15);
maxORD = 63;
minORD = 10;
hScaleFactor = 1;
eps = 1e-9;
tspan = [0 tmax];



d = [];
m = [];
index_l = [];
index_d = [];
index_m = [];

RUNS = 1000;

TIMES_MTSM_OH = zeros(1, RUNS);
TIMES_ODE45 = zeros(1,RUNS);
TIMES_MTSM_LIN = zeros(1,RUNS);
for i=1:100
    tic
    [T_ODE451,Y_ODE451] = ode45(@(t,y) ode(t,y),tspan,y0,options);
    TIMES_ODE45(1,i)=toc;

    tic;
    [T_MTSM_LIN,Y_MTSM_LIN,ORD] = explicitTaylorLinear(dt,tspan,y0,A,b,eps);
    TIMES_MTSM_LIN(1,i)=toc;

    tic;
    [T_MTSM_OH,Y_MTSM_OH,ORD] = taylor_v3(dt,tspan,y0,A,b,d,m,index_l, index_d, index_m, eps,maxORD,minORD,hScaleFactor);
    TIMES_MTSM_OH(1,i)=toc;
end

figure
plot(T_MTSM_OH, Y_MTSM_OH)
grid on

mean(TIMES_ODE45)
mean(TIMES_MTSM_LIN)
mean(TIMES_MTSM_OH)

function dy = ode(t,y)
    dy = zeros(2,1);
    dy(1) = y(2);
    dy(2) = -y(1);
end