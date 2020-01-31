%%% Double pendulum - auxillary and non-auxillary system
function [T_OUT,T_AUX_OUT,T_AUX_CONST_OUT,T_AUX_CONST_DIV_OUT,T_TAYLOR_OUT,T_TAYLOR_GN_V1_OUT,T_TAYLOR_GN_V2_OUT] = simulate(display,m1,m2,L1,L2,g,h,h_mtsm,eps,eps_mtsm,tmax)
    % clear
    %clc
    close all
    
    %%  parameters of the system
    data.m1 = m1; %0.1;
    data.m2 = m2; %0.5;
    data.L1 = L1; %1;
    data.L2 = L2; %2;
    data.g = g; %9.81;
       
    data.A = data.L1*data.m1;
    data.B = data.L1*data.m2;
    data.C = data.L2*data.m1;
    data.D = data.L2*data.m2;
    data.E = data.m1*data.g;
    data.F = data.m2*data.g;
    
    phi1_0 = 0;
    phi2_0 = 0;
    phi1dot_0 = 1;
    phi2dot_0 = 5;
    
    y0 = [
        phi1_0;
        phi2_0;
        phi1dot_0;
        phi2dot_0
        ];
    
     y0_aux = [
        y0;
        sin(phi1_0-phi2_0);
        cos(phi1_0-phi2_0);
        sin(phi1_0);
        cos(phi1_0);
        sin(phi2_0);
        cos(phi2_0);
     ];
 
     y0_aux_div = [
         y0_aux;
         1/(data.A+data.B - data.B*cos(phi1_0-phi2_0)^2);
         1/(data.C+data.D - data.D*cos(phi1_0-phi2_0)^2);
     ];
    
    ne=length(y0); % number of equations
    ne_aux=length(y0_aux); % number of equations
    ne_aux_div=length(y0_aux_div); % number of equations




%%% Solvers options
tspan = [0 tmax];
tspan_ode = 0:h:tmax;



% % MATLAB solvers options
ABSTOL=eps*ones(1,ne);
ABSTOL_aux=eps*ones(1,ne_aux);
ABSTOL_aux_div=eps*ones(1,ne_aux_div);
% options = odeset('RelTol',tol,'AbsTol',ABSTOL, 'Vectorized','on');
options = odeset('RelTol',eps,'AbsTol',ABSTOL);
options_aux = odeset('RelTol',eps,'AbsTol',ABSTOL_aux);
options_aux_div = odeset('RelTol',eps,'AbsTol',ABSTOL_aux_div);

% options = odeset('RelTol',eps,'AbsTol',ABSTOL,'Refine',1,'MaxStep',h);
% options_aux = odeset('RelTol',eps,'AbsTol',ABSTOL_aux,'Refine',1,'MaxStep',h);
% options_aux_div = odeset('RelTol',eps,'AbsTol',ABSTOL_aux_div,'Refine',1,'MaxStep',h);
    
    %% analytical solution for phi_1
    
    %% call an ODE solver - ode45
    tic
    [T,Y] = ode45(@(t,y) pendulum(t,y,data),tspan_ode,y0,options);
    T_OUT = toc; 
    
%     figure;
%     plot(T,Y, '-x');
%     grid on;
%     legend('pos_1', 'pos_2', 'speed_1', 'speed_2');
%     title('Original approach')

    tic
    [T_aux,Y_aux] = ode45(@(t,y) pendulum_aux(t,y,data),tspan_ode,y0_aux,options_aux);
    T_AUX_OUT = toc; 
%     
%     figure;
%     %'-rx', '-mx', '-bx', '-gx'
%     plot(T,Y(:,1:4), '-x');
%     grid on;
%     legend('pos_1', 'pos_2', 'speed_1', 'speed_2');
%     title('Generating eqs')
    
     tic
    [T_aux_const,Y_aux_const] = ode45(@(t,y) pendulum_aux_const(t,y,data),tspan_ode,y0_aux,options_aux);
    T_AUX_CONST_OUT = toc; 
    
%     figure;
%     %'-rx', '-mx', '-bx', '-gx'
%     plot(T,Y(:,1:4), '-x');
%     grid on;
%     legend('pos_1', 'pos_2', 'speed_1', 'speed_2');
%     title('Generating eqs (susbstitued constants')
    
        
    tic
    [T_aux_const_div,Y_aux_const_div] = ode45(@(t,y) pendulum_aux_const_div(t,y,data),tspan_ode,y0_aux_div,options_aux_div);
    T_AUX_CONST_DIV_OUT = toc; 
    
%      figure;
%      %'-rx', '-mx', '-bx', '-gx'
%      plot(T,Y(:,1:4), '-x');
%      grid on;
%      legend('pos_1', 'pos_2', 'speed_1', 'speed_2');
%      title('Generating eqs (susbstitued constants) and removed division')

 %% call an ODE solver - ode113
    tic
    [T_ODE113,Y_ODE113] = ode113(@(t,y) pendulum(t,y,data),tspan_ode,y0,options);
    T_OUT_ODE113 = toc; 

    tic
    [T_aux_ODE113,Y_aux_ODE113] = ode113(@(t,y) pendulum_aux(t,y,data),tspan_ode,y0_aux,options_aux);
    T_AUX_OUT_ODE113 = toc; 
    
     tic
    [T_aux_const_ODE113,Y_aux_const_ODE113] = ode113(@(t,y) pendulum_aux_const(t,y,data),tspan_ode,y0_aux,options_aux);
    T_AUX_CONST_OUT_ODE113 = toc; 
    
        
    tic
    [T_aux_const_div_ODE113,Y_aux_const_div_ODE113] = ode113(@(t,y) pendulum_aux_const_div(t,y,data),tspan_ode,y0_aux_div,options_aux_div);
    T_AUX_CONST_DIV_OUT_ODE113 = toc; 
    

%% Taylor %%%%%
%     system of ODES
%     dz(1) = z(3);
%     dz(2) = z(4);
%     dz(3) = z(11)*(-B*(z(3)^2)*z(5)*z(6) + F*z(9)*z(6) - D*(z(4)^2)*z(5) - (E+F)*z(7));
%     dz(4) = z(12)*( D*(z(4)^2)*z(5)*z(6) + (E+F)*z(7)*z(6) + (A+B)*(z(3)^2)*z(5) - (E+F)*z(9));
%     dz(5) = z(6)*(z(3)-z(4)); % sin(z1-z2)
%     dz(6) = z(5)*(z(4)-z(3)); % cos(z1-z2)
%     dz(7) = z(8)*z(3); % sin(z1)
%     dz(8) = -z(7)*z(3); % cos(z1)
%     dz(9) = z(10)*z(4); % sin(z2)
%     dz(10) = -z(9)*z(4); %cos(z2)
%     dz(11) = -(-2*B*dz(6)*z(6))*z(11)^2;
%     dz(12) = -(-2*D*dz(6)*z(6))*z(12)^2;
%

% Matrices and indexes
% y' = A*y + A2*y_ij + A3*y_ijk + A4*y_ijkl + A5*y_ijklm
A=zeros(12,12);
A(1,3)=1; A(2,4)=1;

ij=[11,7;
    12,9;
    6,3;
    6,4;
    5,4;
    5,3;
    8,3;
    7,3;
    10,4;
    9,4];
A2=zeros(12,size(ij,1));
A2(3,1)=-(data.E+data.F);
A2(4,2)=-(data.E+data.F);
A2(5,3)=1;
A2(5,4)=-1;
A2(6,5)=1;
A2(6,6)=-1;
A2(7,7)=1;
A2(8,8)=-1;
A2(9,9)=1;
A2(10,10)=-1;


ijk=[11 6 9;
    12 6 7];
A3=zeros(12,size(ijk,1));
A3(3,1)=data.F;
A3(4,2)=data.E+data.F;


ijkl=[11 4 4 5;
      12 3 3 5];
A4=zeros(12,size(ijkl,1));
A4(3,1)=-data.D;
A4(4,2)=data.A+data.B;

ijklm=[11 3 3 5 6;
       12 4 4 5 6;
       11 11 6 5 4;
       11 11 6 5 3;
       12 12 6 5 4;
       12 12 6 5 3];
A5=zeros(12,size(ijklm,1));
A5(3,1)=-data.B;
A5(4,2)=data.D;
A5(11,3)=2*data.B;
A5(11,4)=-2*data.B;
A5(12,5)=2*data.D;
A5(12,6)=-2*data.D;

b=zeros(12,1);

%%%%%% PRINT %%%%%
% disp('MATRICES');
% A
% A2
% A3
% A4
% A5
% figure;
% spy(A);
% grid on;
% title('A');
% 
% figure;
% spy(A2);
% grid on;
% title('A2');
% 
% figure;
% spy(A3);
% grid on;
% title('A3');
% 
% figure;
% spy(A4);
% grid on;
% title('A4');
% 
% figure;
% spy(A5);
% grid on;

% sizes of the matrices, vectors
% disp('sizes');
% s_A = size(A)
% s_A2 = size(A2)
% s_A3 = size(A3)
% s_A4 = size(A4)
% s_A5 = size(A5)
% 
% s_b = size(b)
% 
% s_y0_aux_div = size(y0_aux_div)

%%%%%% TEST %%%%%%
% odefun_test = @(t,y) A*y+A2*(y(ij(:,1)).*y(ij(:,2)))+A3*(y(ijk(:,1)).*y(ijk(:,2)).*y(ijk(:,3)))+...
%     A4*(y(ijkl(:,1)).*y(ijkl(:,2)).*y(ijkl(:,3)).*y(ijkl(:,4)))+b;
% tic;
% [T_ode45test,Y_ode45test] = ode45(odefun_test,tspan,y0_aux_div, options_aux_div);
% T_ode45=toc;
% 
% % A3=zeros(12,size(ijk,1));
% % A4=zeros(12,size(ijkl,1));
%  A5=zeros(12,size(ijklm,1));
% 
% [T_Taylor,Y_Taylor,ORD,DY_all] = explicitTaylorMult(h,tspan,y0_aux_div,A,A2,A3,A4,A5,b,ij,ijk,ijkl,ijklm,eps);
% 
% fprintf('||Y(tmax)-Y_Taylor(tmax): %g \n',norm(Y_ode45test(end,:)-Y_Taylor(:,end)'));

%%% end TEST %%%%
maxORD=60;
%ind=load('DY_indexes_maxORD_25.mat','DY_ij','DY_ijk','DY_ijkl','DY_ijklm');
ind=load('DY_indexes_maxORD_GN_ordered_all_60','DY_ij','DY_ijk','DY_ijkl','DY_ijklm');


tic
[T_Taylor,Y_Taylor,ORD_VS,~] = explicitTaylorMult(h_mtsm,tspan,y0_aux_div,A,A2,A3,A4,A5,b,ij,ijk,ijkl,ijklm,eps_mtsm,ind,maxORD);
% [T_Taylor,Y_Taylor,ORD,DY_all] = explicitTaylorMult_ver1(h,tspan,y0_aux_div,A,A2,A3,A4,A5,b,ij,ijk,ijkl,ijklm,eps);
T_TAYLOR_OUT=toc;

tic
[T_Taylor_GN_v1,Y_Taylor_GN_v1,ORD_GN_v1,~] = explicitTaylorMult_GNPV_ver1(h_mtsm,tspan,y0_aux_div,A,A2,A3,A4,A5,b,ij,ijk,ijkl,ijklm,eps_mtsm,ind,maxORD);
% [T_Taylor,Y_Taylor,ORD,DY_all] = explicitTaylorMult_ver1(h,tspan,y0_aux_div,A,A2,A3,A4,A5,b,ij,ijk,ijkl,ijklm,eps);
T_TAYLOR_GN_V1_OUT=toc;

tic
[T_Taylor_GN_v2,Y_Taylor_GN_v2,ORD_GN_v2,~] = explicitTaylorMult_GNPV_ver2_full(h_mtsm,tspan,y0_aux_div,A,A2,A3,A4,A5,b,ij,ijk,ijkl,ijklm,eps_mtsm,ind,maxORD);
% [T_Taylor,Y_Taylor,ORD,DY_all] = explicitTaylorMult_ver1(h,tspan,y0_aux_div,A,A2,A3,A4,A5,b,ij,ijk,ijkl,ijklm,eps);
T_TAYLOR_GN_V2_OUT=toc;

%%% ERRORS %%%%
if display
    fprintf('\n****** ERRORS t=TMAX ******\n');
    fprintf('||Y(tmax)-Y_aux_ode45(tmax)||: %g \n',cnorm(Y(end,:),Y_aux(end,1:4)));
    fprintf('||Y(tmax)-Y_aux_const_ode45(tmax)||: %g \n',cnorm(Y(end,:),Y_aux_const(end,1:4)));
    fprintf('||Y(tmax)-Y_aux_const_div_ode45(tmax)||: %g \n',cnorm(Y(end,:),Y_aux_const_div(end,1:4)));
    
    fprintf('||Y(tmax)-Y_aux_ode113(tmax)||: %g \n',cnorm(Y_ODE113(end,:),Y_aux_ODE113(end,1:4)));
    fprintf('||Y(tmax)-Y_aux_const_ode113(tmax)||: %g \n',cnorm(Y_ODE113(end,:),Y_aux_const_ODE113(end,1:4)));
    fprintf('||Y(tmax)-Y_aux_const_div_ode113(tmax)||: %g \n',cnorm(Y_ODE113(end,:),Y_aux_const_div_ODE113(end,1:4)));
    
    fprintf('||Y(tmax)-Y_Taylor(tmax)||: %g \n',cnorm(Y(end,:),Y_Taylor(1:4,end)'));
    fprintf('||Y(tmax)-Y_Taylor_GN v1(tmax)||: %g \n',cnorm(Y(end,:),Y_Taylor_GN_v1(1:4,end)'));
    fprintf('||Y(tmax)-Y_Taylor_GN v2(tmax)||: %g \n',cnorm(Y(end,:),Y_Taylor_GN_v2(1:4,end)'));


    %%% TIME OF COMPUTATIONS %%%%
    fprintf('\n****** TIME OF COMPUTATIONS [s] ******\n');
    fprintf('ode45: %g \n',T_OUT);
    fprintf('ode45_aux: %g \n',T_AUX_OUT);
    fprintf('ode45_div: %g \n',T_AUX_CONST_DIV_OUT);
    
    fprintf('ode113: %g \n',T_OUT_ODE113);
    fprintf('ode113_aux: %g \n',T_AUX_OUT_ODE113);
    fprintf('ode113_div: %g \n',T_AUX_CONST_DIV_OUT_ODE113);
    
    fprintf('Taylor: %g \n',T_TAYLOR_OUT);
    fprintf('Taylor_GN v1: %g \n',T_TAYLOR_GN_V1_OUT);
    fprintf('Taylor_GN v2: %g \n',T_TAYLOR_GN_V2_OUT);
end


%%% FIGURES %%%
if display
    figure
    subplot(4,1,1);
    plot(T_Taylor,Y_Taylor(1,:));
    grid on;
    xlabel('t')

    subplot(4,1,2);
    plot(T_Taylor,Y_Taylor(2,:));
    grid on;
    xlabel('t')

    subplot(4,1,3);
    plot(T_Taylor,Y_Taylor(3,:));
    grid on;
    xlabel('t')

    subplot(4,1,4);
    plot(T_Taylor,Y_Taylor(4,:));
    grid on;
    xlabel('t')

    figure
    plot(T_Taylor(2:end),ORD_VS(2:end),'+');
    grid on;
    legend('ORD')
    xlabel('t')

    figure
    plot(T_Taylor_GN_v1(2:end),ORD_GN_v1(2:end),'+');
    grid on;
    legend('ORD GN v1')
    xlabel('t')
    
    figure
    plot(T_Taylor_GN_v2(2:end),ORD_GN_v2(2:end),'+');
    grid on;
    legend('ORD GN v2')
    xlabel('t')
end

    
end