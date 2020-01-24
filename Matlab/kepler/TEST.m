function TEST(numruns)
    display = false;
    
    clc
    close all
    
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
        [T_ODE45_sub, T_ODE45_div,T_ODE45_basic,T_ODE23_sub, T_ODE15s_sub, T_TAYLOR_GN_ordered_v2, T_TAYLOR_GN_ordered_div, T_TAYLOR_VS]  =  simulate(false);
        T_ODE45_sub_ALL = [T_ODE45_sub_ALL,T_ODE45_sub];
        T_ODE45_div_ALL = [T_ODE45_div_ALL,T_ODE45_div];
        T_ODE45_basic_ALL = [T_ODE45_basic_ALL,T_ODE45_basic];
        T_ODE23_sub_ALL = [T_ODE23_sub_ALL, T_ODE23_sub];
        T_ODE15s_sub_ALL = [T_ODE15s_sub_ALL, T_ODE15s_sub];
        T_TAYLOR_GN_ordered_v2_ALL = [T_TAYLOR_GN_ordered_v2_ALL, T_TAYLOR_GN_ordered_v2];
        T_TAYLOR_GN_ordered_v2_div_ALL = [T_TAYLOR_GN_ordered_v2_div_ALL, T_TAYLOR_GN_ordered_div];
        T_TAYLOR_VS_ALL = [T_TAYLOR_VS_ALL, T_TAYLOR_VS];
    end
    
    fprintf('runs: %d\n',numruns);
    fprintf('=========================================\n');
    fprintf('ode45 sub: %g\n',median(T_ODE45_sub_ALL));
    fprintf('ode45 div: %g\n',median(T_ODE45_div_ALL));
    fprintf('ode45 basic: %g\n',median(T_ODE45_basic_ALL));
    fprintf('ode23 sub: %g\n',median(T_ODE23_sub_ALL));
    fprintf('ode15s sub: %g\n',median(T_ODE15s_sub_ALL));
    fprintf('taylor GN ordered v2: %g\n',median(T_TAYLOR_GN_ordered_v2_ALL));
    fprintf('taylor VS: %g\n',median(T_TAYLOR_VS_ALL));
end