function [t, y, ORD] = MTSM_explicit_linear_serial(h, tspan, init, A, b, eps, maxORD)
%% MTSM_explicit_linear_serial
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
    ne = size(A,1); 
    AT = A';
    % number of integration steps
    steps = round((tspan(2)-tspan(1))/h); 
    % vector for simulation times
    t = zeros(steps+1,1);
    % vector for orders
    ORD = zeros(steps+1,1);

    % stoping rule: maximum number DY terms are smaller then eps
    stopping = 3;
    % max ORD for double precision
    %maxORD = 200;
    % last step tolerance
    ls_tol = h/10; %10^-10;

    % timestep index
    i = 1; 
    t(i) = tspan(1,1);
    % initial values
    y(1,:) = init'; 

    i = i+1;
    t(i) = t(i-1)+h;

    while t(i-1)+ls_tol<tspan(1,2) %t(i-1)+ls_tol<tspan(1,2)
        % matrix for Taylor series terms for each equation
        DY = zeros(maxORD,ne); 
        % y_{i+1}=y_{i}+...
        y(i,:) = y(i-1,:); 
        Ay = A*y(i,:)';
        % first Taylor series term
        DY(1,:) = h*(Ay+b);

        % add first Taylor series term
        y(i,:)=y(i,:)+DY(1,:);

        k=1; % TODO: zbytecne??
        %maxDY = ones(1,stopping)*10^10;
        maxDY = [ones(1,stopping)*10^10 zeros(1,maxORD+stopping)];

        k = 2;
        % other Taylor series terms
        while norm(maxDY(k-1:stopping+k-2))>eps
            ADy = DY(k-1,:)*AT;
            DY(k,:) = (h/k)*ADy;

            y(i,:) = y(i,:) + DY(k,:);

            maxDY(k+stopping-1) = abs(max(DY(k,:)));
            k = k+1;
            % max ORD was reached
            if k > maxORD 
                break;
            end    
        end

        % max ORD was reached
        if k > maxORD 
            h = h/2; % we use half-size integration step
        else
            % end of current integration step 
            % save order and shift the simulation time
            ORD(i) = k;
            i = i+1;
            t(i) = t(i-1) + h;

        end
    end

    t(end)=[];
end



