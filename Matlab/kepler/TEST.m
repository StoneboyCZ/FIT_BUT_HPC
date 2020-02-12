function TEST(numruns)
    display = false;
    
    %% kepler problem
    e = 0.75; % eccentricity
    
    % parameters of calculation
    tol = 1e-12;
    tol_ode = 3e-14;

    tmax = 14; %2*pi; does not work for ODE solvers
    
    dt = 0.07;
    dt_ode = 0.01;
    
    clc
    close all
    
    global data;
    
    T_ODE45_sub_ALL = [];
    T_ODE45_div_ALL = [];
    T_ODE45_basic_ALL = [];
    T_ODE23_sub_ALL = [];
    T_ODE15s_sub_ALL = [];
    T_TAYLOR_GN_ordered_v2_ALL = [];
    T_TAYLOR_GN_ordered_v2_div_ALL = [];
    T_TAYLOR_VS_ALL = [];
    for i=1:numruns
        i
        [T_ODE45_sub, T_ODE45_div,T_ODE45_basic,T_ODE23_sub, T_ODE15s_sub, T_TAYLOR_GN_ordered_v2, T_TAYLOR_GN_ordered_div, T_TAYLOR_VS]  =  simulate(display,tol,tol_ode,tmax,dt,dt_ode,e);
        T_ODE45_sub_ALL = [T_ODE45_sub_ALL,T_ODE45_sub];
        T_ODE45_div_ALL = [T_ODE45_div_ALL,T_ODE45_div];
        T_ODE45_basic_ALL = [T_ODE45_basic_ALL,T_ODE45_basic];
        T_ODE23_sub_ALL = [T_ODE23_sub_ALL, T_ODE23_sub];
        T_ODE15s_sub_ALL = [T_ODE15s_sub_ALL, T_ODE15s_sub];
        T_TAYLOR_GN_ordered_v2_ALL = [T_TAYLOR_GN_ordered_v2_ALL, T_TAYLOR_GN_ordered_v2];
        T_TAYLOR_GN_ordered_v2_div_ALL = [T_TAYLOR_GN_ordered_v2_div_ALL, T_TAYLOR_GN_ordered_div];
        T_TAYLOR_VS_ALL = [T_TAYLOR_VS_ALL, T_TAYLOR_VS];
    end
    fn = sprintf('runs_%d_e_%g_tolMTSM_%g_tolODE_%g_dtMTSM_%g_dtODE_%g_tmax_%g.txt',numruns,e,tol,tol_ode,dt,dt_ode,tmax);
    f = fopen(fn,'w');
    fprintf(f,'runs: %d\n e: %g\n tol_MTSM: %g\n tol_ODE: %g\n dt_MTSM: %g\n dt_ODE: %g\n tmax: %g\n',numruns,e,tol,tol_ode,dt,dt_ode,tmax);
    fprintf(f,'=========================================\n');
    fprintf(f,'ode45 sub: %g\n',median(T_ODE45_sub_ALL));
    fprintf(f,'ode45 div: %g\n',median(T_ODE45_div_ALL));
    fprintf(f,'ode45 basic: %g\n',median(T_ODE45_basic_ALL));
    fprintf(f,'ode23 sub: %g\n',median(T_ODE23_sub_ALL));
    fprintf(f,'ode15s sub: %g\n',median(T_ODE15s_sub_ALL));
    fprintf(f,'taylor GN ordered v2 (division not substituted): %g\n',median(T_TAYLOR_GN_ordered_v2_div_ALL));
    fprintf(f,'taylor GN ordered v2: %g\n',median(T_TAYLOR_GN_ordered_v2_ALL));
    fprintf(f,'taylor VS: %g\n',median(T_TAYLOR_VS_ALL));
    fclose(f);
end