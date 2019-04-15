function [t, y, ORD] = MTSM_explicit_linear_precalc_one_mtx_parallel(h, tspan, init, A, B, eps, maxORD,S)
%% MTSM_explicit_linear_precalc_one_mtx_parallel
% Numerical solution of ODEs with explicit Taylor series method
% problem y' = A*y + b
% where A is Jacobian matrix of constants and b is a vector of constants
% ---------------------------------------------------------------
% PARAMETERS: 
% h         integration step size
% tspan 	interval of numerical solution [0,tmax]
% init      vector with initial conditions
% A         Jacobian matrix A 
% b         vector b 
% eps       precision of the calculation 
% maxORD    maximm number of Taylor series terms
% 
% RETURN:
% ---------------
% t         vector with simulation times
% y         vector with solutions 
% ORD       vector with orders for each integration step

    % number of equations
    ne=size(A,2); % number of equations
    %AT=A';
    % number of integration steps
    steps = round((tspan(2)-tspan(1))/h); % number of integration steps
    % vector for simulation times
    t = zeros(steps+1,1);
    % vector for orders
    ORD = zeros(steps+1,1);

    % stoping rule: maximum number DY terms are smaller then eps
    stopping = 3;
    % max ORD for double precision
    %maxORD=200;
    % last step tolerance
    ls_tol = h/10;%10^-10;

    % timestep index
    i = 1; 
    t(i) = tspan(1,1);
    % initial values
    y(1,:) = init'; 

    i = i+1;
    t(i) = t(i-1)+h;

    spmd
        AA = A(S(labindex,1):S(labindex,2),:);
        YY = y(:,:);
        while t(i-1)+ls_tol<tspan(1,2) %t(i-1)+ls_tol<tspan(1,2)
            % y_{i+1}=y_{i}+...
            %YY(i,:) = y
            YY(i,:)=YY(i-1,:);
            %y(i,:)=y(i-1,:);
            YYC = YY(i,:);
            AY = YYC';
            size(AA)
            size(AY)
            size(AA*AY)
            size(B)
            
            
            % by one multiplication (and addition) we obtain results for current
            % integration step 
            YYC(i,S(labindex,1):S(labindex,2)) = AA*AY + B(S(labindex,1):S(labindex,2),1)

            % shift simulation time
            i = i+1;
            t(i) = t(i-1)+h;

        end
    end
    %t(end)=[];
end



