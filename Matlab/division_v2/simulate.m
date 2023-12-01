% Problem: y' = (sin(t)cos(t)/e^t)

%% init

clc
close all

tmax = 2*pi;
eps = 1e-8;
dt = pi/10;
%% VARIANT 1: FK
% y_1' &= y_5y_2y_3 & y_1(0) &= \frac{\sin(0)\cos(0)}{e^0} \\
% y_2' &= y_3  & y_2(0) &= \sin(0)\\
% y_3' &= -y_2 & y_3(0) &= \cos(0) \\
% y_4' &= y_4 & y_4(0) &= e^0 \\
% y_5' & = -y_5y_5y_4 & y_5(0) &= \frac{1}{y_4(0)}

% ode solvers
y0 = [(sin(0)*cos(0)/1);sin(0);cos(0);1;1]; 
tspan = [0 tmax];
options = odeset('RelTol',1e-11,'AbsTol',1e-13);

tic;
[T,Y] = ode45(@(t,y) odeFK(t,y),tspan,y0,options);
TIME_FK_ODE45 = toc;

figure
plot(T,Y(:,1));
grid on;

TIME_FK_ODE45

ne = 5;
A = zeros(ne);
A(2,3) = 1;
A(3,2) = -1;
A(4,4) = 1;

ij = [];
B2 = zeros(ne, size(ij,1));

ijk = [
    5,2,3;
    5,5,4;
];
B3 = zeros(ne, size(ijk,1));

ijkl = [];
B4 = zeros(ne, size(ijkl,1));

ijklm = [];
B5 = zeros(ne, size(ijklm,1));







%% VARIANT 2: OH







function dy =  odeFK(~, y)
    dy = zeros(5,1);
    dy(1) = y(5)*y(2)*y(3);
    dy(2) = y(3);
    dy(3) = -y(2);
    dy(4) = y(4);
    dy(5) = -y(5)*y(5)*y(4);
end