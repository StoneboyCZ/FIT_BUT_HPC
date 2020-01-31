function TEST(runs)
    display = false;
    T = [];
    T_AUX = [];
    T_AUX_CONST = [];
    T_AUX_CONST_DIV = [];
    T_TAYLOR = [];
    T_TAYLOR_GN_V1 = [];
    T_TAYLOR_GN_V2 = [];

    for n=1:runs
        n
        [T_OUT,T_AUX_OUT,T_AUX_CONST_OUT,T_AUX_CONST_DIV_OUT,T_TAYLOR_OUT,T_TAYLOR_GN_V1_OUT,T_TAYLOR_GN_V2_OUT] = simulate(display);
        T = [T,T_OUT];
        T_AUX = [T_AUX,T_AUX_OUT];
        T_AUX_CONST = [T_AUX_CONST,T_AUX_CONST_OUT];
        T_AUX_CONST_DIV = [T_AUX_CONST_DIV, T_AUX_CONST_DIV_OUT];
        T_TAYLOR = [T_TAYLOR, T_TAYLOR_OUT];
        T_TAYLOR_GN_V1 = [T_TAYLOR_GN_V1, T_TAYLOR_GN_V1_OUT];
        T_TAYLOR_GN_V2 = [T_TAYLOR_GN_V2, T_TAYLOR_GN_V2_OUT];
    end

fprintf('runs: %d\n',runs);
fprintf('time ode45 solver (original): %g\n', median(T));
fprintf('time ode45 solver (auxillary system): %g\n', median(T_AUX));
fprintf('time ode45 solver (auxillary system + constants): %g\n', median(T_AUX_CONST));
fprintf('time ode45 solver (auxillary system + constants + no division): %g\n', median(T_AUX_CONST_DIV_OUT));    
fprintf('time Taylor solver: %g\n', median(T_TAYLOR));    
fprintf('time Taylor solver (GN v1): %g\n', median(T_TAYLOR_GN_V1));    
fprintf('time Taylor solver (GN v2): %g\n', median(T_TAYLOR_GN_V2));    

end