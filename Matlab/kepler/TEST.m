function TEST(numruns)
    global results

    display = false;
    % parameters of calculation
    tol = 1e-9;
    tol_ode = 1e-12;

    tmax = 4*pi; %2*pi; does not work for ODE solvers
    
    dt = 4*pi/180;
    
    e = 0.75; % eccentricity
    
    
    
    clc
    close all

    test.kepler.ode45.time = [];
    
    test.kepler_sqrt.ode45.time = [];
    
    test.kepler_div.ode45.time = [];
    test.kepler_div.ode23.time = [];
    test.kepler_div.ode15s.time = [];
    test.kepler_div.ode113.time = [];
    test.kepler_div.mtsm_basic.time = [];
    test.kepler_div.mtsm_v2.time = [];
    
    test.kepler_div_brackets.ode45.time = [];
    test.kepler_div_brackets.ode23.time = [];
    test.kepler_div_brackets.ode15s.time = [];
    test.kepler_div_brackets.ode113.time = [];
    test.kepler_div_brackets.mtsm_basic.time = [];
    test.kepler_div_brackets.mtsm_v2.time = [];
    
    test.kepler_div_full.ode45.time = [];
    test.kepler_div_full.ode23.time = [];
    test.kepler_div_full.ode15s.time = [];
    test.kepler_div_full.ode113.time = [];
    test.kepler_div_full.mtsm_basic.time = [];
    test.kepler_div_full.mtsm_v2.time = [];
    
    for i=1:numruns
        i
%         [T_ODE45_sub, T_ODE45_div,T_ODE45_basic,T_ODE23_sub, T_ODE15s_sub, T_TAYLOR_GN_ordered_v2, T_TAYLOR_GN_ordered_div, T_TAYLOR_VS]  =  simulate(display,tol,tol_ode,tmax,dt,dt_ode,e);
        simulate(display,tol,tol_ode,tmax,dt,e);
        test.kepler.ode45.time = [test.kepler.ode45.time,results.kepler.ode45.time];
        
        test.kepler_sqrt.ode45.time = [test.kepler_sqrt.ode45.time,results.kepler_sqrt.ode45.time];
        
        test.kepler_div.ode45.time = [test.kepler_div.ode45.time,results.kepler_div.ode45.time];
        test.kepler_div.ode23.time = [test.kepler_div.ode23.time,results.kepler_div.ode23.time];
        test.kepler_div.ode15s.time = [test.kepler_div.ode15s.time,results.kepler_div.ode15s.time];
        test.kepler_div.ode113.time = [test.kepler_div.ode113.time,results.kepler_div.ode113.time];
        test.kepler_div.mtsm_basic.time = [test.kepler_div.mtsm_basic.time,results.kepler_div.mtsm_basic.time];
        test.kepler_div.mtsm_v2.time = [test.kepler_div.mtsm_v2.time,results.kepler_div.mtsm_v2.time];
        
        test.kepler_div_brackets.ode45.time = [test.kepler_div_brackets.ode45.time,results.kepler_div_brackets.ode45.time];
        test.kepler_div_brackets.ode23.time = [test.kepler_div_brackets.ode23.time,results.kepler_div_brackets.ode23.time];
        test.kepler_div_brackets.ode15s.time = [test.kepler_div_brackets.ode15s.time,results.kepler_div_brackets.ode15s.time];
        test.kepler_div_brackets.ode113.time = [test.kepler_div_brackets.ode113.time,results.kepler_div_brackets.ode113.time];
        test.kepler_div_brackets.mtsm_basic.time = [test.kepler_div_brackets.mtsm_basic.time,results.kepler_div_brackets.mtsm_basic.time];
        test.kepler_div_brackets.mtsm_v2.time = [test.kepler_div_brackets.mtsm_v2.time,results.kepler_div_brackets.mtsm_v2.time];
        
        test.kepler_div_full.ode45.time = [test.kepler_div_full.ode45.time,results.kepler_div_full.ode45.time];
        test.kepler_div_full.ode23.time = [test.kepler_div_full.ode23.time,results.kepler_div_full.ode23.time];
        test.kepler_div_full.ode15s.time = [test.kepler_div_full.ode15s.time,results.kepler_div_full.ode15s.time];
        test.kepler_div_full.ode113.time = [test.kepler_div_full.ode113.time,results.kepler_div_full.ode113.time];
        test.kepler_div_full.mtsm_basic.time = [test.kepler_div_full.mtsm_basic.time,results.kepler_div_full.mtsm_basic.time];
        test.kepler_div_full.mtsm_v2.time = [test.kepler_div_full.mtsm_v2.time,results.kepler_div_full.mtsm_v2.time];
    end
    fn = sprintf('runs_%d_e_%g_tolMTSM_%g_tolODE_%g_dtMTSM_%g_tmax_%g.txt',numruns,e,tol,tol_ode,dt,tmax);
    f = fopen(fn,'w');
    fprintf(f,'runs: %d\n e: %g\n tol_MTSM: %g\n tol_ODE: %g\n dt_MTSM: %g\n tmax: %g\n',numruns,e,tol,tol_ode,dt,tmax);
    fprintf(f,'=========================================\n');
    fprintf(f,'BASIC SYSTEM (%d equations): ode45 %g\n',results.kepler.n,median(test.kepler.ode45.time));
    fprintf(f,'SYSTEM WITHOUT SQUARE ROOT (%d equations): ode45 %g\n',results.kepler_sqrt.n,median(test.kepler_sqrt.ode45.time));
    fprintf(f,'SYSTEM WITHOUT DIVISION (%d equations)\n======\n',results.kepler_div.n);
    fprintf(f,'ode45: %g\n',median(test.kepler_div.ode45.time));
    fprintf(f,'ode23: %g\n',median(test.kepler_div.ode23.time));
    fprintf(f,'ode15s: %g\n',median(test.kepler_div.ode15s.time));
    fprintf(f,'ode113: %g\n',median(test.kepler_div.ode113.time));
    fprintf(f,'mtsm_basic: %g\n',median(test.kepler_div.mtsm_basic.time));
    
    fprintf(f,'ode45/mtsm_basic: %g\n',median(test.kepler_div.ode45.time)/median(test.kepler_div.mtsm_basic.time));
    fprintf(f,'ode23/mtsm_basic: %g\n',median(test.kepler_div.ode23.time)/median(test.kepler_div.mtsm_basic.time));
    fprintf(f,'ode15s/mtsm_basic: %g\n',median(test.kepler_div.ode15s.time)/median(test.kepler_div.mtsm_basic.time));
    fprintf(f,'ode113/mtsm_basic: %g\n',median(test.kepler_div.ode113.time)/median(test.kepler_div.mtsm_basic.time));
    
    fprintf(f,'mtsm_v2: %g\n',median(test.kepler_div.mtsm_v2.time));
    
    fprintf(f,'ode45/mtsm_v2: %g\n',median(test.kepler_div.ode45.time)/median(test.kepler_div.mtsm_v2.time));
    fprintf(f,'ode23/mtsm_v2: %g\n',median(test.kepler_div.ode23.time)/median(test.kepler_div.mtsm_v2.time));
    fprintf(f,'ode15s/mtsm_v2: %g\n',median(test.kepler_div.ode15s.time)/median(test.kepler_div.mtsm_v2.time));
    fprintf(f,'ode113/mtsm_v2: %g\n',median(test.kepler_div.ode113.time)/median(test.kepler_div.mtsm_v2.time));
    
    
    fprintf(f,'SYSTEM WITHOUT DIVISION, NOT FULLY SUBSTITUTED (%d equations)\n======\n',results.kepler_div_brackets.n);
    fprintf(f,'ode45: %g\n',median(test.kepler_div_brackets.ode45.time));
    fprintf(f,'ode23: %g\n',median(test.kepler_div_brackets.ode23.time));
    fprintf(f,'ode15s: %g\n',median(test.kepler_div_brackets.ode15s.time));
    fprintf(f,'ode113: %g\n',median(test.kepler_div_brackets.ode113.time));
    fprintf(f,'mtsm_basic: %g\n',median(test.kepler_div_brackets.mtsm_basic.time));
    
    fprintf(f,'ode45/mtsm_basic: %g\n',median(test.kepler_div_brackets.ode45.time)/median(test.kepler_div_brackets.mtsm_basic.time));
    fprintf(f,'ode23/mtsm_basic: %g\n',median(test.kepler_div_brackets.ode23.time)/median(test.kepler_div_brackets.mtsm_basic.time));
    fprintf(f,'ode15s/mtsm_basic: %g\n',median(test.kepler_div_brackets.ode15s.time)/median(test.kepler_div_brackets.mtsm_basic.time));
    fprintf(f,'ode113/mtsm_basic: %g\n',median(test.kepler_div_brackets.ode113.time)/median(test.kepler_div_brackets.mtsm_basic.time));
    
    fprintf(f,'mtsm_v2: %g\n',median(test.kepler_div_brackets.mtsm_v2.time));
    
    fprintf(f,'ode45/mtsm_v2: %g\n',median(test.kepler_div_brackets.ode45.time)/median(test.kepler_div_brackets.mtsm_v2.time));
    fprintf(f,'ode23/mtsm_v2: %g\n',median(test.kepler_div_brackets.ode23.time)/median(test.kepler_div_brackets.mtsm_v2.time));
    fprintf(f,'ode15s/mtsm_v2: %g\n',median(test.kepler_div_brackets.ode15s.time)/median(test.kepler_div_brackets.mtsm_v2.time));
    fprintf(f,'ode113/mtsm_v2: %g\n',median(test.kepler_div_brackets.ode113.time)/median(test.kepler_div_brackets.mtsm_v2.time));
    
        

    fprintf(f,'COMPLETELY SUBSTITUTED SYSTEM (%d equations)\n======\n',results.kepler_div_full.n);
    fprintf(f,'ode45: %g\n',median(test.kepler_div_full.ode45.time));
    fprintf(f,'ode23: %g\n',median(test.kepler_div_full.ode23.time));
    fprintf(f,'ode15s: %g\n',median(test.kepler_div_full.ode15s.time));
    fprintf(f,'ode113: %g\n',median(test.kepler_div_full.ode113.time));
    fprintf(f,'mtsm_basic: %g\n',median(test.kepler_div_full.mtsm_basic.time));

    fprintf(f,'ode45/mtsm_basic: %g\n',median(test.kepler_div_full.ode45.time)/median(test.kepler_div_full.mtsm_basic.time));
    fprintf(f,'ode23/mtsm_basic: %g\n',median(test.kepler_div_full.ode23.time)/median(test.kepler_div_full.mtsm_basic.time));
    fprintf(f,'ode15s/mtsm_basic: %g\n',median(test.kepler_div_full.ode15s.time)/median(test.kepler_div_full.mtsm_basic.time));
    fprintf(f,'ode113/mtsm_basic: %g\n',median(test.kepler_div_full.ode113.time)/median(test.kepler_div_full.mtsm_basic.time));
    
    fprintf(f,'mtsm_v2: %g\n',median(test.kepler_div_full.mtsm_v2.time));
    
    fprintf(f,'ode45/mtsm_v2: %g\n',median(test.kepler_div_full.ode45.time)/median(test.kepler_div_full.mtsm_v2.time));
    fprintf(f,'ode23/mtsm_v2: %g\n',median(test.kepler_div_full.ode23.time)/median(test.kepler_div_full.mtsm_v2.time));
    fprintf(f,'ode15s/mtsm_v2: %g\n',median(test.kepler_div_full.ode15s.time)/median(test.kepler_div_full.mtsm_v2.time));
    fprintf(f,'ode113/mtsm_v2: %g\n',median(test.kepler_div_full.ode113.time)/median(test.kepler_div_full.mtsm_v2.time));    

    fclose(f);
end