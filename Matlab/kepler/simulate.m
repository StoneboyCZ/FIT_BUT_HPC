% function [T_ODE45_sub, T_ODE45_div,T_ODE45_basic,T_ODE23_sub, T_ODE15s_sub, FULL_GN_v2_TIME, DIV_GN_v2_TIME, DIV_VS_TIME] =  simulate(display,tol,tol_ode,tmax,dt,e)
function simulate(display,tol,tol_ode,tmax,dt,e) 
    % simulate(true, 1e-9,1e-12,4*pi,4*pi/180,0.75)
    close all;
    
    global results
    
    Y_ANAL = 1;
    tspan = [0, tmax];
    
    %% MTSM nonlinear solver
    maxORD = 60;
    
    %% initial conditions
    y1_init = 1-e;
    y2_init = 0;
    y3_init = 0;
    y4_init = sqrt((1+e)/(1-e));

    y0 = [y1_init;y2_init;y3_init;y4_init];

    y5_init = sqrt(y1_init*y1_init+y2_init*y2_init)^3;
    y6_init = sqrt(y1_init*y1_init+y2_init*y2_init);
    y0_aux_sqrt = [y0;
        y5_init;
        y6_init
        ];

    y7_init = 1/y5_init;
    y8_init = 1/y6_init;
    y0_aux_div = [y0_aux_sqrt;
         y7_init;
         y8_init
    ];
    
    y9_init = y1_init*y3_init;
    y10_init = y2_init*y4_init;
    y11_init = y1_init^2;
    y12_init = y2_init^2;
    
    
    y0_aux_div_brackets = [y0_aux_div;
         y9_init;
         y10_init;
         y11_init;
         y12_init
    ];
    
    y13_init = y9_init + y10_init;
    y14_init = y11_init + y12_init;
    y15_init = y6_init*y13_init;
    y16_init = y8_init*y13_init;
    
    y0_aux_div_full = [y0_aux_div_brackets;
        y13_init;
        y14_init;
        y15_init;
        y16_init;
    ];
 

    ne=length(y0);
    ne_aux_sqrt=length(y0_aux_sqrt);
    ne_aux_div = length(y0_aux_div);
    ne_aux_div_brackets = length(y0_aux_div_brackets);
    ne_aux_div_full = length(y0_aux_div_full);
    
    % % MATLAB solvers options
    ABSTOL=tol_ode*ones(1,ne);
    ABSTOL_aux_sqrt=tol_ode*ones(1,ne_aux_sqrt);
    ABSTOL_aux_div=tol_ode*ones(1,ne_aux_div);
    ABSTOL_aux_div_brackets=tol_ode*ones(1,ne_aux_div_brackets);
    ABSTOL_aux_div_full=tol_ode*ones(1,ne_aux_div_full);
    % options = odeset('RelTol',tol,'AbsTol',ABSTOL, 'Vectorized','on');
    options = odeset('RelTol',tol_ode,'AbsTol',ABSTOL);
    options_aux_sqrt = odeset('RelTol',tol_ode,'AbsTol',ABSTOL_aux_sqrt);
    options_aux_div = odeset('RelTol',tol_ode,'AbsTol',ABSTOL_aux_div);
    options_aux_div_brackets = odeset('RelTol',tol_ode,'AbsTol',ABSTOL_aux_div_brackets);
    options_aux_div_full = odeset('RelTol',tol_ode,'AbsTol',ABSTOL_aux_div_full);

    %% Results are saved in the variable results
    % results.<type>.<solver>.time
    % type: kepler, kepler_aux_sqrt,kepler_aux_div, kepler_aux_div_brackets, kepler_aux_div_full 
    
        
    %% basic systems
    %     r = sqrt(y(1)*y(1) + y(2)*y(2));
    %     r3 = r*r*r;
    % 
    %     dy = zeros(4,1);
    %     dy(1) = y(3);
    %     dy(2) = y(4);
    %     dy(3) = -(y(1)/r3); 
    %     dy(4) = -(y(2)/r3);
    
    % ODE45 
    tic;
    [T_ODE45,Y_ODE45] = ode45(@(t,y) kepler(t,y),tspan,y0,options);
    T_ODE45_basic = toc;

    
    results.kepler.n = length(y0);
    results.kepler.ode45.time = T_ODE45_basic;
    results.kepler.ode45.T = T_ODE45;
    results.kepler.ode45.Y = Y_ODE45;
    results.kepler.ode45.norm = cnorm(Y_ANAL,(results.kepler.ode45.Y(:,1)+e).^2 + results.kepler.ode45.Y(:,2).^2/(1-e^2));

    if display
        figure
        plot(T_ODE45,Y_ODE45(:,2));
        grid on;
        xlabel('t');
        TITLE=sprintf("Kepler problem (e=%f) - ODE45",e);
        title(TITLE)
    end
    
    %% ODE solver - aux sqrt
    %     dy(1) = y(3);
    %     dy(2) = y(4);
    %     dy(3) = -y(1)/y(5); 
    %     dy(4) = -y(2)/y(5);
    %     dy(5) = 3*y(6)*(y(1)*y(3)+y(2)*y(4));
    %     dy(6) = (y(1)*y(3)+y(2)*y(4))/y(6);
    
    tic;
    [T_ODE45_AUX_SQRT,Y_ODE45_AUX_SQRT] = ode45(@(t,y) kepler_aux_sqrt(t,y),tspan,y0_aux_sqrt,options_aux_sqrt);
    T_ODE45_sqrt = toc;

    results.kepler_sqrt.n = length(y0_aux_sqrt);
    results.kepler_sqrt.ode45.time = T_ODE45_sqrt;
    results.kepler_sqrt.ode45.T = T_ODE45_AUX_SQRT;
    results.kepler_sqrt.ode45.Y = Y_ODE45_AUX_SQRT;
    results.kepler_sqrt.ode45.norm = cnorm(Y_ANAL,(results.kepler_sqrt.ode45.Y(:,1)+e).^2 + results.kepler_sqrt.ode45.Y(:,2).^2/(1-e^2));
    
    if display
        figure
        plot(T_ODE45_AUX_SQRT,Y_ODE45_AUX_SQRT(:,2));
        grid on;
        xlabel('t');
        TITLE=sprintf("Kepler problem (e=%f) - ODE45 (AUX SQRT)",e);
        title(TITLE)
    end
    
    %%   ODE solver - div
    %     dy(1) = y(3);
    %     dy(2) = y(4);
    %     dy(3) = -y(1)*y(7); 
    %     dy(4) = -y(2)*y(7); 
    %     dy(5) = 3*y(6)*y(1)*y(3)+3*y(6)*y(2)*y(4);
    %     dy(6) = y(1)*y(3)*y(8)+y(2)*y(4)*y(8); 
    %     dy(7) = -3*y(6)*y(1)*y(3)*y(7)*y(7)-3*y(6)*y(2)*y(4)*y(7)*y(7); 
    %     dy(8) = -y(1)*y(3)*y(8)*y(8)*y(8)-y(2)*y(4)*y(8)*y(8)*y(8);
    
    results.kepler_div.n = length(y0_aux_div);
    
    tic;
    [T_ODE45_AUX_DIV,Y_ODE45_AUX_DIV] = ode45(@(t,y) kepler_aux_div(t,y),tspan,y0_aux_div,options_aux_div);
    T_ODE45_div = toc;

    results.kepler_div.ode45.time = T_ODE45_div;
    results.kepler_div.ode45.T = T_ODE45_AUX_DIV;
    results.kepler_div.ode45.Y = Y_ODE45_AUX_DIV;
    results.kepler_div.ode45.norm = cnorm(Y_ANAL,(results.kepler_div.ode45.Y(:,1)+e).^2 + results.kepler_div.ode45.Y(:,2).^2/(1-e^2));
    
    tic;
    [T_ODE23_AUX_DIV,Y_ODE23_AUX_DIV] = ode23(@(t,y) kepler_aux_div(t,y),tspan,y0_aux_div,options_aux_div);
    T_ODE23_div = toc;

    results.kepler_div.ode23.time = T_ODE23_div;
    results.kepler_div.ode23.T = T_ODE23_AUX_DIV;
    results.kepler_div.ode23.Y = Y_ODE23_AUX_DIV;
    results.kepler_div.ode23.norm = cnorm(Y_ANAL,(results.kepler_div.ode23.Y(:,1)+e).^2 + results.kepler_div.ode23.Y(:,2).^2/(1-e^2));
    
    tic;
    [T_ODE15s_AUX_DIV,Y_ODE15s_AUX_DIV] = ode15s(@(t,y) kepler_aux_div(t,y),tspan,y0_aux_div,options_aux_div);
    T_ODE15s_div = toc;

    results.kepler_div.ode15s.time = T_ODE15s_div;
    results.kepler_div.ode15s.T = T_ODE15s_AUX_DIV;
    results.kepler_div.ode15s.Y = Y_ODE15s_AUX_DIV;
    results.kepler_div.ode15s.norm = cnorm(Y_ANAL,(results.kepler_div.ode15s.Y(:,1)+e).^2 + results.kepler_div.ode15s.Y(:,2).^2/(1-e^2));
    
    tic;
    [T_ODE113_AUX_DIV,Y_ODE113_AUX_DIV] = ode113(@(t,y) kepler_aux_div(t,y),tspan,y0_aux_div,options_aux_div);
    T_ODE113_div = toc;

    results.kepler_div.ode113.time = T_ODE113_div;
    results.kepler_div.ode113.T = T_ODE113_AUX_DIV;
    results.kepler_div.ode113.Y = Y_ODE113_AUX_DIV;
    results.kepler_div.ode113.norm = cnorm(Y_ANAL,(results.kepler_div.ode113.Y(:,1)+e).^2 + results.kepler_div.ode113.Y(:,2).^2/(1-e^2)); 
    
    if display
        figure
        plot(T_ODE45_AUX_DIV,Y_ODE45_AUX_DIV(:,2));
        grid on;
        xlabel('t');
        TITLE=sprintf("Kepler problem (e=%f) - ODE45 (AUX DIV)",e);
        title(TITLE)
    end
    
    init = y0_aux_div;
%     [DIV_VS_T,DIV_VS_Y,DIV_VS_TIME,DIV_VS_ORD,DIV_VS_ANAL,DIV_GN_v1_T,DIV_GN_v1_Y,DIV_GN_v1_TIME,DIV_GN_v1_ORD,DIV_GN_v1_ANAL,DIV_GN_v2_T,DIV_GN_v2_Y,DIV_GN_v2_TIME,DIV_GN_v2_ORD,DIV_GN_v2_ANAL] = taylor_divAux(dt,tspan,init,tol,maxORD,e,display);
    taylor_div(dt,tspan,init,tol,maxORD,e,display);
    
    %% Division with brackets (12 equations)
    %     dy(1) = y(3);
    %     dy(2) = y(4);
    %     dy(3) = -y(1)*y(7); 
    %     dy(4) = -y(2)*y(7); 
    %     dy(5) = 3*y(6)*(y(9)+y(10));
    %     dy(6) = y(8)*(y(9)+y(10));
    %     dy(7) = -3*y(6)*y(7)*y(7)*(y(9)+y(10));
    %     dy(8) = (-y(8))*(-y(8))*(-y(8))*(y(9)+y(10));
    %     dy(9) = y(3)*y(3) - y(7)*y(11);
    %     dy(10) = y(4)*y(4) -y(7)*y(12);
    %     dy(11) = 2*y(9);
    %     dy(12) = 2*y(10);
    
    results.kepler_div_brackets.n = length(y0_aux_div_brackets);
    
    tic;
    [T_ODE45_DIV_BRACKETS,Y_ODE45_DIV_BRACKETS] = ode45(@(t,y) kepler_aux_div_brackets(t,y),tspan,y0_aux_div_brackets,options_aux_div_brackets);
    T_ODE45_div_brakets = toc;

    results.kepler_div_brackets.ode45.time = T_ODE45_div_brakets;
    results.kepler_div_brackets.ode45.T = T_ODE45_DIV_BRACKETS;
    results.kepler_div_brackets.ode45.Y = Y_ODE45_DIV_BRACKETS;
    results.kepler_div_brackets.ode45.norm = cnorm(Y_ANAL,(results.kepler_div_brackets.ode45.Y(:,1)+e).^2 + results.kepler_div_brackets.ode45.Y(:,2).^2/(1-e^2));
    
    tic;
    [T_ODE23_AUX_DIV_BRACKETS,Y_ODE23_AUX_DIV_BRACKETS] = ode23(@(t,y) kepler_aux_div_brackets(t,y),tspan,y0_aux_div_brackets,options_aux_div_brackets);
    T_ODE23_div_brackets = toc;

    results.kepler_div_brackets.ode23.time = T_ODE23_div_brackets;
    results.kepler_div_brackets.ode23.T = T_ODE23_AUX_DIV_BRACKETS;
    results.kepler_div_brackets.ode23.Y = Y_ODE23_AUX_DIV_BRACKETS;
    results.kepler_div_brackets.ode23.norm = cnorm(Y_ANAL,(results.kepler_div_brackets.ode23.Y(:,1)+e).^2 + results.kepler_div_brackets.ode23.Y(:,2).^2/(1-e^2));
    
    tic;
    [T_ODE15s_AUX_DIV_BRACKETS,Y_ODE15s_AUX_DIV_BRACKETS] = ode15s(@(t,y) kepler_aux_div_brackets(t,y),tspan,y0_aux_div_brackets,options_aux_div_brackets);
    T_ODE15s_div_brackets = toc;

    results.kepler_div_brackets.ode15s.time = T_ODE15s_div_brackets;
    results.kepler_div_brackets.ode15s.T = T_ODE15s_AUX_DIV_BRACKETS;
    results.kepler_div_brackets.ode15s.Y = Y_ODE15s_AUX_DIV_BRACKETS;
    results.kepler_div_brackets.ode15s.norm = cnorm(Y_ANAL,(results.kepler_div_brackets.ode15s.Y(:,1)+e).^2 + results.kepler_div_brackets.ode15s.Y(:,2).^2/(1-e^2));
    
    tic;
    [T_ODE113_AUX_DIV_BRACKETS,Y_ODE113_AUX_DIV_BRACKETS] = ode113(@(t,y) kepler_aux_div_brackets(t,y),tspan,y0_aux_div_brackets,options_aux_div_brackets);
    T_ODE113_div_brackets = toc;

    results.kepler_div_brackets.ode113.time = T_ODE113_div_brackets;
    results.kepler_div_brackets.ode113.T = T_ODE113_AUX_DIV_BRACKETS;
    results.kepler_div_brackets.ode113.Y = Y_ODE113_AUX_DIV_BRACKETS;
    results.kepler_div_brackets.ode113.norm = cnorm(Y_ANAL,(results.kepler_div_brackets.ode113.Y(:,1)+e).^2 + results.kepler_div_brackets.ode113.Y(:,2).^2/(1-e^2)); 
    
    if display
        figure
        plot(T_ODE45_DIV_BRACKETS,Y_ODE45_DIV_BRACKETS(:,2));
        grid on;
        xlabel('t');
        TITLE=sprintf("Kepler problem (e=%f) - ODE45 (DIV BRACKETS)",e);
        title(TITLE)
    end
    
    init = y0_aux_div_brackets;
    taylor_div_brackets(dt,tspan,init,tol,maxORD,e,display);
    %% Full substitution
    %     dy(1) = y(3);
    %     dy(2) = y(4);
    %     dy(3) = -y(1)*y(7); 
    %     dy(4) = -y(2)*y(7); 
    %     dy(5) = 3*y(15);
    %     dy(6) = y(16);
    %     dy(7) = -3*y(15)*y(7)*y(7);
    %     dy(8) = -y(8)*y(8)*y(16);
    %     dy(9) = y(3)*y(3) - y(7)*y(11);
    %     dy(10) = y(4)*y(4) - y(7)*y(12);
    %     dy(11) = 2*y(9);
    %     dy(12) = 2*y(10);
    %     dy(13) = y(3)*y(3) + y(4)*y(4) - y(7)*y(14);
    %     dy(14) = 2*y(9) + 2*y(10);
    %     dy(15) = y(8)*y(13)*y(13) + y(6)*y(3)*y(3) + y(6)*y(4)*y(4) - y(6)*y(7)*y(14);
    %     dy(16) = -y(8)*y(16)*y(16) + y(8)*y(3)*y(3) + y(8)*y(4)*y(4) - y(8)*y(7)*y(14);    
    
    results.kepler_div_full.n = length(y0_aux_div_full);
    
    tic;
    [T_ODE45_DIV_FULL,Y_ODE45_DIV_FULL] = ode45(@(t,y) kepler_aux_div_full(t,y),tspan,y0_aux_div_full,options_aux_div_full);
    T_ODE45_div_full = toc;

    results.kepler_div_full.ode45.time = T_ODE45_div_full;
    results.kepler_div_full.ode45.T = T_ODE45_DIV_FULL;
    results.kepler_div_full.ode45.Y = Y_ODE45_DIV_FULL;
    results.kepler_div_full.ode45.norm = cnorm(Y_ANAL,(results.kepler_div_full.ode45.Y(:,1)+e).^2 + results.kepler_div_full.ode45.Y(:,2).^2/(1-e^2));    
    
    % ODE23
    tic;
    [T_ODE23_DIV_FULL,Y_ODE23_DIV_FULL] = ode23(@(t,y) kepler_aux_div_full(t,y),tspan,y0_aux_div_full,options_aux_div_full);
    T_ODE23_div_full = toc;
    
    results.kepler_div_full.ode23.time = T_ODE23_div_full;
    results.kepler_div_full.ode23.T = T_ODE23_DIV_FULL;
    results.kepler_div_full.ode23.Y = Y_ODE23_DIV_FULL;
    results.kepler_div_full.ode23.norm = cnorm(Y_ANAL,(results.kepler_div_full.ode23.Y(:,1)+e).^2 + results.kepler_div_full.ode23.Y(:,2).^2/(1-e^2)); 
    
    % ODE15s
    tic;
    [T_ODE15s_DIV_FULL,Y_ODE15s_DIV_FULL] = ode15s(@(t,y) kepler_aux_div_full(t,y),tspan,y0_aux_div_full,options_aux_div_full);
    T_ODE15s_div_full = toc;
    
    results.kepler_div_full.ode15s.time = T_ODE15s_div_full;
    results.kepler_div_full.ode15s.T = T_ODE15s_DIV_FULL;
    results.kepler_div_full.ode15s.Y = Y_ODE15s_DIV_FULL;
    results.kepler_div_full.ode15s.norm = cnorm(Y_ANAL,(results.kepler_div_full.ode23.Y(:,1)+e).^2 + results.kepler_div_full.ode23.Y(:,2).^2/(1-e^2)); 
    
    % ODE15s
    tic;
    [T_ODE113_DIV_FULL,Y_ODE113_DIV_FULL] = ode113(@(t,y) kepler_aux_div_full(t,y),tspan,y0_aux_div_full,options_aux_div_full);
    T_ODE113_div_full = toc;
    
    results.kepler_div_full.ode113.time = T_ODE113_div_full;
    results.kepler_div_full.ode113.T = T_ODE113_DIV_FULL;
    results.kepler_div_full.ode113.Y = Y_ODE113_DIV_FULL;
    results.kepler_div_full.ode113.norm = cnorm(Y_ANAL,(results.kepler_div_full.ode113.Y(:,1)+e).^2 + results.kepler_div_full.ode113.Y(:,2).^2/(1-e^2));
    
    if display
        figure
        plot(T_ODE45_DIV_FULL,Y_ODE45_DIV_FULL(:,2));
        grid on;
        xlabel('t');
        TITLE=sprintf("Kepler problem (e=%f) - ODE45 (DIV FULL)",e);
        title(TITLE)
    end
   
    %MTSM 
    init = y0_aux_div_full;
%     [FULL_VS_T,FULL_VS_Y,FULL_VS_TIME,FULL_VS_ORD,FULL_VS_ANAL,FULL_GN_v2_T,FULL_GN_v2_Y,FULL_GN_v2_TIME,FULL_GN_v2_ORD,FULL_GN_v2_ANAL] = taylor_fullAux(dt,tspan,init,tol,maxORD,e,display);
    taylor_div_full(dt,tspan,init,tol,maxORD,e,display);
   
    
    if display
        % Output 
        fprintf('==== SYSTEM OF %d EQUATIONS (BASIC SYSTEM) ====\n',results.kepler.n); 
        fprintf('ode45: %g steps: %d ||norm||: %g \n', results.kepler.ode45.time, length(results.kepler.ode45.T)-1, results.kepler.ode45.norm);
        
        fprintf('\n\n==== SYSTEM OF %d EQUATIONS (WITHOUT SQUARE ROOT)====\n',results.kepler_sqrt.n); 
        fprintf('ode45: %g steps: %d ||norm||: %g \n',results.kepler_sqrt.ode45.time, length(results.kepler_sqrt.ode45.T)-1,results.kepler_sqrt.ode45.norm);

        fprintf('\n\n==== SYSTEM OF %d EQUATIONS (WITHOUT SQUARE ROOT AND DIVISION)====\n',results.kepler_div.n); 
        fprintf('ode45:  %g  steps: %d  ||norm||: %g \n',results.kepler_div.ode45.time, length(results.kepler_div.ode45.T)-1,results.kepler_div.ode45.norm);
        fprintf('ode23:  %g  steps: %d  ||norm||: %g \n',results.kepler_div.ode23.time, length(results.kepler_div.ode23.T)-1,results.kepler_div.ode23.norm);
        fprintf('ode15s: %g  steps: %d  ||norm||: %g \n',results.kepler_div.ode15s.time, length(results.kepler_div.ode15s.T)-1,results.kepler_div.ode15s.norm);
        fprintf('ode113: %g  steps: %d  ||norm||: %g \n',results.kepler_div.ode113.time, length(results.kepler_div.ode113.T)-1,results.kepler_div.ode113.norm);
        fprintf('mtsm_basic\ntime: %g  steps: %d ORD: %g ||norm||: %g \n',results.kepler_div.mtsm_basic.time, length(results.kepler_div.mtsm_basic.T)-1,mean(results.kepler_div.mtsm_basic.ORD),results.kepler_div.mtsm_basic.norm);       
        fprintf('ode45/mtsm_basic: %g\n',results.kepler_div.ode45.time/results.kepler_div.mtsm_basic.time);
        fprintf('ode23/mtsm_basic: %g\n',results.kepler_div.ode23.time/results.kepler_div.mtsm_basic.time);
        fprintf('ode15s/mtsm_basic: %g\n',results.kepler_div.ode15s.time/results.kepler_div.mtsm_basic.time);
        fprintf('ode113/mtsm_basic: %g\n',results.kepler_div.ode113.time/results.kepler_div.mtsm_basic.time);
        
        fprintf('mtsm_v2\ntime: %g  steps: %d ORD: %g ||norm||: %g \n',results.kepler_div.mtsm_v2.time, length(results.kepler_div.mtsm_v2.T)-1,mean(results.kepler_div.mtsm_v2.ORD),results.kepler_div.mtsm_v2.norm);       
        fprintf('ode45/mtsm_v2: %g\n',results.kepler_div.ode45.time/results.kepler_div.mtsm_v2.time);
        fprintf('ode23/mtsm_v2: %g\n',results.kepler_div.ode23.time/results.kepler_div.mtsm_v2.time);
        fprintf('ode15s/mtsm_v2: %g\n',results.kepler_div.ode15s.time/results.kepler_div.mtsm_v2.time);
        fprintf('ode113/mtsm_v2: %g\n',results.kepler_div.ode113.time/results.kepler_div.mtsm_v2.time);
        
        fprintf('\n\n==== SYSTEM OF %d EQUATIONS (WITHOUT SQUARE ROOT AND DIVISION, NOT FULLY SUBSTITUTED)====\n',results.kepler_div_brackets.n); 
        fprintf('ode45:  %g  steps: %d  ||norm||: %g \n',results.kepler_div_brackets.ode45.time, length(results.kepler_div_brackets.ode45.T)-1,results.kepler_div_brackets.ode45.norm);
        fprintf('ode23:  %g  steps: %d  ||norm||: %g \n',results.kepler_div_brackets.ode23.time, length(results.kepler_div_brackets.ode23.T)-1,results.kepler_div_brackets.ode23.norm);
        fprintf('ode15s: %g  steps: %d  ||norm||: %g \n',results.kepler_div_brackets.ode15s.time, length(results.kepler_div_brackets.ode15s.T)-1,results.kepler_div_brackets.ode15s.norm);
        fprintf('ode113: %g  steps: %d  ||norm||: %g \n',results.kepler_div_brackets.ode113.time, length(results.kepler_div_brackets.ode113.T)-1,results.kepler_div_brackets.ode113.norm);
        fprintf('mtsm_basic\ntime: %g  steps: %d ORD: %g ||norm||: %g \n',results.kepler_div_brackets.mtsm_basic.time, length(results.kepler_div_brackets.mtsm_basic.T)-1,mean(results.kepler_div_brackets.mtsm_basic.ORD),results.kepler_div_brackets.mtsm_basic.norm);       
        fprintf('ode45/mtsm_basic: %g\n',results.kepler_div_brackets.ode45.time/results.kepler_div_brackets.mtsm_basic.time);
        fprintf('ode23/mtsm_basic: %g\n',results.kepler_div_brackets.ode23.time/results.kepler_div_brackets.mtsm_basic.time);
        fprintf('ode15s/mtsm_basic: %g\n',results.kepler_div_brackets.ode15s.time/results.kepler_div_brackets.mtsm_basic.time);
        fprintf('ode113/mtsm_basic: %g\n',results.kepler_div_brackets.ode113.time/results.kepler_div_brackets.mtsm_basic.time);
        
        fprintf('mtsm_v2\ntime: %g  steps: %d ORD: %g ||norm||: %g \n',results.kepler_div_brackets.mtsm_v2.time, length(results.kepler_div_brackets.mtsm_v2.T)-1,mean(results.kepler_div_brackets.mtsm_v2.ORD),results.kepler_div_brackets.mtsm_v2.norm);       
        fprintf('ode45/mtsm_v2: %g\n',results.kepler_div_brackets.ode45.time/results.kepler_div_brackets.mtsm_v2.time);
        fprintf('ode23/mtsm_v2: %g\n',results.kepler_div_brackets.ode23.time/results.kepler_div_brackets.mtsm_v2.time);
        fprintf('ode15s/mtsm_v2: %g\n',results.kepler_div_brackets.ode15s.time/results.kepler_div_brackets.mtsm_v2.time);
        fprintf('ode113/mtsm_v2: %g\n',results.kepler_div_brackets.ode113.time/results.kepler_div_brackets.mtsm_v2.time);
        
        
        
        fprintf('\n\n==== SYSTEM OF %d EQUATIONS (WITHOUT SQUARE ROOT AND DIVISION, FULLY SUBSTITUTED)====\n',results.kepler_div_full.n); 
        fprintf('ode45:  %g  steps: %d  ||norm||: %g \n',results.kepler_div_full.ode45.time, length(results.kepler_div_full.ode45.T)-1,results.kepler_div_full.ode45.norm);
        fprintf('ode23:  %g  steps: %d  ||norm||: %g \n',results.kepler_div_full.ode23.time, length(results.kepler_div_full.ode23.T)-1,results.kepler_div_full.ode23.norm);
        fprintf('ode15s: %g  steps: %d  ||norm||: %g \n',results.kepler_div_full.ode15s.time, length(results.kepler_div_full.ode15s.T)-1,results.kepler_div_full.ode15s.norm);
        fprintf('ode113: %g  steps: %d  ||norm||: %g \n',results.kepler_div_full.ode113.time, length(results.kepler_div_full.ode113.T)-1,results.kepler_div_full.ode113.norm);
        fprintf('mtsm_basic\ntime: %g  steps: %d ORD: %g ||norm||: %g \n',results.kepler_div_full.mtsm_basic.time, length(results.kepler_div_full.mtsm_basic.T)-1,mean(results.kepler_div_full.mtsm_basic.ORD),results.kepler_div_full.mtsm_basic.norm);       
        fprintf('ode45/mtsm_basic: %g\n',results.kepler_div_full.ode45.time/results.kepler_div_full.mtsm_basic.time);
        fprintf('ode23/mtsm_basic: %g\n',results.kepler_div_full.ode23.time/results.kepler_div_full.mtsm_basic.time);
        fprintf('ode15s/mtsm_basic: %g\n',results.kepler_div_full.ode15s.time/results.kepler_div_full.mtsm_basic.time);
        fprintf('ode113/mtsm_basic: %g\n',results.kepler_div_full.ode113.time/results.kepler_div_full.mtsm_basic.time);
        
        fprintf('mtsm_v2\ntime: %g  steps: %d ORD: %g ||norm||: %g \n',results.kepler_div_full.mtsm_v2.time, length(results.kepler_div_full.mtsm_v2.T)-1,mean(results.kepler_div_full.mtsm_v2.ORD),results.kepler_div_full.mtsm_v2.norm);       
        fprintf('ode45/mtsm_v2: %g\n',results.kepler_div_full.ode45.time/results.kepler_div_full.mtsm_v2.time);
        fprintf('ode23/mtsm_v2: %g\n',results.kepler_div_full.ode23.time/results.kepler_div_full.mtsm_v2.time);
        fprintf('ode15s/mtsm_v2: %g\n',results.kepler_div_full.ode15s.time/results.kepler_div_full.mtsm_v2.time);
        fprintf('ode113/mtsm_v2: %g\n',results.kepler_div_full.ode113.time/results.kepler_div_full.mtsm_v2.time);
    end
end

function taylor_div(dt,tspan,init,tol,maxORD,e,display)
    global results
% function [VS_T,VS_Y,VS_TIME,VS_ORD,VS_ANAL,GN_v1_T,GN_v1_Y,GN_v1_TIME,GN_v1_ORD,GN_v1_ANAL,GN_v2_T,GN_v2_Y,GN_v2_TIME,GN_v2_ORD,GN_v2_ANAL] = taylor_divAux(dt,tspan,init,tol,maxORD,e,display)
%     dy(1) = y(3);
%     dy(2) = y(4);
%     dy(3) = -y(1)*y(7); 
%     dy(4) = -y(2)*y(7); 
%     dy(5) = 3*y(6)*y(1)*y(3)+3*y(6)*y(2)*y(4);
%     dy(6) = y(1)*y(3)*y(8)+y(2)*y(4)*y(8); 
%     dy(7) = 3*y(6)*y(1)*y(3)*y(7)*y(7)+3*y(6)*y(2)*y(4)*y(7)*y(7); 
%     dy(8) = y(1)*y(3)*y(8)*y(8)*y(8)+y(2)*y(4)*y(8)*y(8)*y(8);

    ne=8;
    A = zeros(ne,ne);
    A(1,3) = 1;
    A(2,4) = 1;
    
     % 2 multiplications -- ij
    % rhs - indeces of y(i) * y(j)
    ij = [
      1,7; %1
      2,7; %2
    ];
    
    % A(equation_index, ij_index) = multiplication constant 
    A2=zeros(ne,size(ij,1));
    A2(3,1) = -1;
    A2(4,2) = -1;

    % 3 multiplications -- ijk
    ijk = [
        6,1,3; % 1
        6,2,4; % 2
        1,3,8; % 3
        2,4,8; % 4
    ];

    A3=zeros(ne,size(ijk,1));
    A3(5,1) = 3;
    A3(5,2) = 3;
    A3(6,3) = 1;
    A3(6,4) = 1;
    
    ijkl = [];
    A4=zeros(ne,size(ijkl,1));
    
    % 5 multiplications -- ijkl
    ijklm = [
        6,1,3,7,7; % 1
        6,2,4,7,7; % 2
        1,3,8,8,8; % 3
        2,4,8,8,8; % 4
    ];

    A5=zeros(ne,size(ijklm,1));
    A5(7,1) = -3;
    A5(7,2) = -3;
    A5(8,3) = -1;
    A5(8,4) = -1;

    b=zeros(ne,1);
    
    ind=load('DY_indexes_maxORD_GN_ordered_all_60','DY_ij','DY_ijk','DY_ijklm');
    
    % standard implementation
    tic
    [VS_T,VS_Y,VS_ORD,~] = explicitTaylorMult(dt,tspan,init,A,A2,A3,A4,A5,b,ij,ijk,ijkl,ijklm,tol,ind,maxORD);
    VS_TIME=toc;
    
    results.kepler_div.mtsm_basic.time = VS_TIME;
    results.kepler_div.mtsm_basic.T = VS_T;
    results.kepler_div.mtsm_basic.Y = VS_Y;
    results.kepler_div.mtsm_basic.ORD = VS_ORD;
    results.kepler_div.mtsm_basic.norm = cnorm(1,(results.kepler_div.mtsm_basic.Y(1,:)+e).^2 + results.kepler_div.mtsm_basic.Y(2,:).^2/(1-e^2));
    
    % GNPV implemetation
    tic
    [GN_v2_T,GN_v2_Y,GN_v2_ORD,~] = explicitTaylorMult_GNPV_ver2_full(dt,tspan,init,A,A2,A3,A4,A5,b,ij,ijk,ijkl,ijklm,tol,ind,maxORD);
    GN_v2_TIME=toc;
    
    results.kepler_div.mtsm_v2.time = GN_v2_TIME;
    results.kepler_div.mtsm_v2.T = GN_v2_T;
    results.kepler_div.mtsm_v2.Y = GN_v2_Y;
    results.kepler_div.mtsm_v2.ORD = GN_v2_ORD;
    results.kepler_div.mtsm_v2.norm = cnorm(1,(GN_v2_Y(1,:)+e).^2 + GN_v2_Y(2,:).^2/(1-e^2));
 
    if display
        figure
        plot(GN_v2_Y(1,:),GN_v2_Y(2,:));
        grid on;
        xlabel('x');
        ylabel('y');
        hold on;
        TITLE=sprintf("MTSM (div) - Orbit e=%f",e);
        title(TITLE)
           
        figure
        plot(GN_v2_T,GN_v2_ORD,'*');
        grid on;
        TITLE=sprintf("Kepler problem (e=%f) - TAYLOR GNPV (div ORD)",e);
        title(TITLE)

        figure
        plot(VS_T,VS_ORD,'*');
        grid on;
        TITLE=sprintf("Kepler problem (e=%f) - TAYLOR VS (div ORD)",e);
        title(TITLE)
    end
end

function taylor_div_brackets(dt,tspan,init,tol,maxORD,e,display)
    %     dy(1) = y(3);
    %     dy(2) = y(4);
    %     dy(3) = -y(1)*y(7); 
    %     dy(4) = -y(2)*y(7); 
    %     dy(5) = 3*y(6)*y(9)+3*y(6)*y(10);
    %     dy(6) = y(8)*y(9)+y(8)*y(10);
    %     dy(7) = -3*y(6)*y(7)*y(7)*y(9)-3*y(6)*y(7)*y(7)*y(10);
    %     dy(8) = -y(8)*y(8)*y(8)*y(9) - y(8)*y(8)*y(8)*y(10);
    %     dy(9) = y(3)*y(3) - y(7)*y(11);
    %     dy(10) = y(4)*y(4)- y(7)*y(12);
    %     dy(11) = 2*y(9);
    %     dy(12) = 2*y(10);
    
    global results;
    ne = 12;

    A = zeros(ne,ne);
    A(1,3) = 1;
    A(2,4) = 1;
    A(11,9) = 2;
    A(12,10) = 2;
    
    ij = [
        1,7;
        2,7;
        6,9;
        6,10;
        8,9;
        8,10;
        3,3;
        7,11;
        4,4;
        7,12
        ];
    A2=zeros(ne,size(ij,1));
    A2(3,1) = -1;
    A2(4,2) = -1;
    A2(5,3) = 3;
    A2(5,4) = 3;
    A2(6,5) = 1;
    A2(6,6) = 1;
    A2(9,7) = 1;
    A2(9,8) = -1;
    A2(10,9) = 1;
    A2(10,10) = -1;
    
    ijk = [];
    A3=zeros(ne,size(ijk,1)); 
    
    ijkl = [
        6,7,7,9;
        6,7,7,10;
        8,8,8,9;
        8,8,8,10;
        ];
    A4=zeros(ne,size(ijkl,1));
    A4(7,1) = -3;
    A4(7,2) = -3;
    A4(8,3) = -1;
    A4(8,4) = -1;
    
    ijklm = [];
    A5=zeros(ne,size(ijklm,1));
    
    b=zeros(ne,1);
    
    ind=load('DY_indexes_maxORD_GN_ordered_all_60','DY_ij','DY_ijkl');
    
    % VS implementation
    tic
    [VS_T,VS_Y,VS_ORD,~] = explicitTaylorMult(dt,tspan,init,A,A2,A3,A4,A5,b,ij,ijk,ijkl,ijklm,tol,ind,maxORD);
    VS_TIME=toc;
    
    results.kepler_div_brackets.mtsm_basic.time = VS_TIME;
    results.kepler_div_brackets.mtsm_basic.T = VS_T;
    results.kepler_div_brackets.mtsm_basic.Y = VS_Y;
    results.kepler_div_brackets.mtsm_basic.ORD = VS_ORD;
    results.kepler_div_brackets.mtsm_basic.norm = cnorm(1,(results.kepler_div_brackets.mtsm_basic.Y(1,:)+e).^2 + results.kepler_div_brackets.mtsm_basic.Y(2,:).^2/(1-e^2));    
    
    % GNPV implemetation
    tic
    [GN_T,GN_Y,GN_ORD,~] = explicitTaylorMult_GNPV_ver2_full(dt,tspan,init,A,A2,A3,A4,A5,b,ij,ijk,ijkl,ijklm,tol,ind,maxORD);
    GN_TIME=toc;
    
    results.kepler_div_brackets.mtsm_v2.time = GN_TIME;
    results.kepler_div_brackets.mtsm_v2.T = GN_T;
    results.kepler_div_brackets.mtsm_v2.Y = GN_Y;
    results.kepler_div_brackets.mtsm_v2.ORD = GN_ORD;
    results.kepler_div_brackets.mtsm_v2.norm = cnorm(1,(results.kepler_div_brackets.mtsm_v2.Y(1,:)+e).^2 + results.kepler_div_brackets.mtsm_v2.Y(2,:).^2/(1-e^2));
    
end

function taylor_div_full(dt,tspan,init,tol,maxORD,e,display)
    global results
    % System of ODEs
    %     dy = zeros(16,1);
    %     dy(1) = y(3);
    %     dy(2) = y(4);
    %     dy(3) = -y(1)*y(7); 
    %     dy(4) = -y(2)*y(7); 
    %     dy(5) = 3*y(15);
    %     dy(6) = y(16);
    %     dy(7) = -3*y(15)*y(7)*y(7);
    %     dy(8) = -y(8)*y(8)*y(16);
    %     dy(9) = y(3)*y(3) - y(7)*y(11);
    %     dy(10) = y(4)*y(4) - y(7)*y(12);
    %     dy(11) = 2*y(9);
    %     dy(12) = 2*y(10);
    %     dy(13) = y(3)*y(3) + y(4)*y(4) - y(7)*y(14);
    %     dy(14) = 2*y(9) + 2*y(10);
    %     dy(15) = y(8)*y(13)*y(13) + y(6)*y(3)*y(3) + y(6)*y(4)*y(4) - y(6)*y(7)*y(14);
    %     dy(16) = -y(8)*y(16)*y(16) + y(8)*y(3)*y(3) + y(8)*y(4)*y(4) - y(8)*y(7)*y(14);


    % matrices and indexes y' = A*y + A2*y_ij + A3*y_ijk + A4*y_ijkl + A5*y_ijklm
    % A(equation_index, rhs index) = multiplication constant 
    ne = 16;
    A = zeros(ne,ne);
    A(1,3) = 1;
    A(2,4) = 1;
    A(5,15) = 3;
    A(6,16) = 1;
    A(11,9) = 2;
    A(12,10) = 2;
    A(14,9) = 2;
    A(14,10) = 2;

    % 2 multiplications -- ij
    % rhs - indeces of y(i) * y(j)
    ij = [
      1,7; %1
      2,7; %2
      3,3; %3
      7,11; %4
      4,4; %5
      7,12; %6
      7,14; %7
    ];
    
    % A(equation_index, ij_index) = multiplication constant 
    A2=zeros(ne,size(ij,1));
    A2(3,1) = -1;
    A2(4,2) = -1;
    A2(9,3) = 1;
    A2(9,4) = -1;
    A2(10,5) = 1;
    A2(10,6) = -1;
    A2(13,3) = 1;
    A2(13,5) = 1;
    A2(13,7) = -1;

    % 3 multiplications -- ijk
    ijk = [
        15,7,7; % 1
        8,8,16; % 2
        8,13,13; % 3
        6,3,3; % 4
        6,4,4; % 5 
        6,7,14; % 6 
        8,16,16; % 7
        8,3,3; %8 
        8,4,4; %9
        8,7,14 %10
    ];

    A3=zeros(ne,size(ijk,1));
    A3(7,1) = -3;
    A3(8,2) = -1;
    A3(15,3) = 1;
    A3(15,4) = 1;
    A3(15,5) = 1;
    A3(15,6) = -1;
    A3(16,7) = -1;
    A3(16,8) = 1;
    A3(16,9) = 1;
    A3(16,10) = -1;
    
    ijkl = [];
    A4=zeros(ne,size(ijkl,1));
    
    ijklm = [];
    A5=zeros(ne,size(ijklm,1));

    b=zeros(ne,1);

    ind=load('DY_indexes_maxORD_GN_ordered_all_60','DY_ij','DY_ijk');
    
    % VS implementation
    tic
    [VS_T,VS_Y,VS_ORD,~] = explicitTaylorMult(dt,tspan,init,A,A2,A3,A4,A5,b,ij,ijk,ijkl,ijklm,tol,ind,maxORD);
    VS_TIME=toc;
    
    results.kepler_div_full.mtsm_basic.time = VS_TIME;
    results.kepler_div_full.mtsm_basic.T = VS_T;
    results.kepler_div_full.mtsm_basic.Y = VS_Y;
    results.kepler_div_full.mtsm_basic.ORD = VS_ORD;
    results.kepler_div_full.mtsm_basic.norm = cnorm(1,(results.kepler_div_full.mtsm_basic.Y(1,:)+e).^2 + results.kepler_div_full.mtsm_basic.Y(2,:).^2/(1-e^2));
    
    % GNPV implemetation
    tic
    [GN_T,GN_Y,GN_ORD,~] = explicitTaylorMult_GNPV_ver2_full(dt,tspan,init,A,A2,A3,A4,A5,b,ij,ijk,ijkl,ijklm,tol,ind,maxORD);
    GN_TIME=toc;
    
    results.kepler_div_full.mtsm_v2.time = GN_TIME;
    results.kepler_div_full.mtsm_v2.T = GN_T;
    results.kepler_div_full.mtsm_v2.Y = GN_Y;
    results.kepler_div_full.mtsm_v2.ORD = GN_ORD;
    results.kepler_div_full.mtsm_v2.norm = cnorm(1,(results.kepler_div_full.mtsm_v2.Y(1,:)+e).^2 + results.kepler_div_full.mtsm_v2.Y(2,:).^2/(1-e^2));
    
    if display
        figure
        plot(GN_Y(1,:),GN_Y(2,:));
        grid on;
        xlabel('x');
        ylabel('y');
        hold on;
        TITLE=sprintf("MTSM (full substitution) - Orbit e=%f",e);
        title(TITLE)
           
        figure
        plot(GN_T,GN_ORD,'*');
        grid on;
        TITLE=sprintf("Kepler problem (e=%f) - TAYLOR GNPV (ORD)",e);
        title(TITLE)

        figure
        plot(VS_T,VS_ORD,'*');
        grid on;
        TITLE=sprintf("Kepler problem (e=%f) - TAYLOR VS (ORD)",e);
        title(TITLE)
    end
end


