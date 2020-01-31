function TEST(runs)
    display = false;
    
    %% double pendulum
    m1 = 1;
    m2 = 5;
    L1 = 1;
    L2 = 2;
    g = 9.81;
    
    % solver options
    h = 0.01;
    h_mtsm = 0.01;
    
    eps=1e-13;
    eps_mtsm = 1e-9; 
    
    tmax = 3;
    
    T = [];
    T_AUX = [];
    T_AUX_CONST = [];
    T_AUX_CONST_DIV = [];
    T_TAYLOR = [];
    T_TAYLOR_GN_V1 = [];
    T_TAYLOR_GN_V2 = [];

    for n=1:runs
        n
        [T_OUT,T_AUX_OUT,T_AUX_CONST_OUT,T_AUX_CONST_DIV_OUT,T_TAYLOR_OUT,T_TAYLOR_GN_V1_OUT,T_TAYLOR_GN_V2_OUT] = simulate(display,m1,m2,L1,L2,g,h,h_mtsm,eps,eps_mtsm,tmax);
        T = [T,T_OUT];
        T_AUX = [T_AUX,T_AUX_OUT];
        T_AUX_CONST = [T_AUX_CONST,T_AUX_CONST_OUT];
        T_AUX_CONST_DIV = [T_AUX_CONST_DIV, T_AUX_CONST_DIV_OUT];
        T_TAYLOR = [T_TAYLOR, T_TAYLOR_OUT];
        T_TAYLOR_GN_V1 = [T_TAYLOR_GN_V1, T_TAYLOR_GN_V1_OUT];
        T_TAYLOR_GN_V2 = [T_TAYLOR_GN_V2, T_TAYLOR_GN_V2_OUT];
    end

fn = sprintf('runs_%d_m1_%g_m2_%g_L1_%g_L2_%g_g_%g_h_%g_hMTSM_%g_eps_%g_epsMTSM_%g_tmax_%g.txt',runs,m1,m2,L1,L2,g,h,h_mtsm,eps,eps_mtsm,tmax);
f = fopen(fn,'w');    
fprintf(f,'runs: %d\n m1: %g\n m2: %g\n L1: %g\n L2: %g\n g: %g\n h: %g\n h_mtsm: %g\n eps: %g\n eps_mtsm: %g\n tmax: %g\n',runs,m1,m2,L1,L2,g,h,h_mtsm,eps,eps_mtsm,tmax);
fprintf(f,'===================================================================\n');
fprintf(f,'time ode45 solver (original): %g\n', median(T));
fprintf(f,'time ode45 solver (auxillary system): %g\n', median(T_AUX));
fprintf(f,'time ode45 solver (auxillary system + constants): %g\n', median(T_AUX_CONST));
fprintf(f,'time ode45 solver (auxillary system + constants + no division): %g\n', median(T_AUX_CONST_DIV_OUT));    
fprintf(f,'time Taylor solver: %g\n', median(T_TAYLOR));    
fprintf(f,'time Taylor solver (GN v1): %g\n', median(T_TAYLOR_GN_V1));    
fprintf(f,'time Taylor solver (GN v2): %g\n', median(T_TAYLOR_GN_V2));    
fclose(f);
end