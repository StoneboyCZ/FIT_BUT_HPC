function [MTSM_serial, MTSM_one_serial, MTSM_one_parallel] = simulateExportTimes(nRLC, nWorkers)
    %% simulate
    % Solves the Telegraph Line Equation (2nd order PDE)
    % PDE for voltage and current - simplified model
    % -> resistance of the wire (R(x)) and conductance between wires (G(x)) are
    % ommited
    % 
    % L*C \partial^2 u(x,t) / \partial t^2 - \partial^2 u(x,t) / \partial x^2 = 0
    % L*C \partial^2 i(x,t) / \partial t^2  - \partial^2 i(x,t) / \partial x^2} = 0
    % 
    % PDEs are transformed into the system of ODEs 
    % y' = A*y + b 
    % 
    % where 
    % A = [A11 A12; A21 A22]	parse block matrices for currents and voltages
    %   note: matrix A also includes system of 2 coupled linear ODEs for the 
    %   generation of sine function 
    %   U_0 = U_0*sin(omega*t)
    %   u_0' = omega*x      u_0(0) = 0;
    %   x' = \omega*u_0     x_0(0) = U_0;
    %
    % y0 = [0...0 1]            vector with initial conditions 
    % b = [0...0]
    % 
    % PARAMETERS:
    % ------------
    % nRLC  number of RLC segments
    % TODO: add more parameters - R1, R2, type of input, eps, dt?
    %
    % TODO: For more information, see: Model of the Telegraph Line and its Numerical
    % Solution, journal OpenCS, 2018.
    %close all
    %clc

    %% for now, just adjusted line
    tmax = 1;
    eps = 1e-9; % presicion of the calculation
    dt = 0.1; % size of the integration step

    % simulation time
    tspan = [0, tmax];  
    maxORD = 64;

    [A,b,init,nRLC] = loadFromFile();
    
    %% generate different matrix representations for parallel solution
    % Approach 1: 3D matrix - precalculation of the matrix A for the higher derivatives
    % each Taylor series term [(h^i/i!) * A^i] for i=1..maxORD is in special
    % matrix, so we have an array of these matrices -> 3D matrix 
%   [A_3D, B_3D] = generate_3D_mtx(nRLC, A, b, dt, maxORD); 

    % Approach 2: one big matrix A_hat
    % all Taylor series terms [(h^i/i!) * A^i] for i=1..maxORD are added
    % together, so we obtain one big matrix 

    % serial version
    [A_one_serial, b_one_serial,T_prealloc_serial] = generate_one_mtx(nRLC, A, b, dt, maxORD);
    % parallel version
    [A_one_parallel, b_one_parallel,T_prealloc_parallel] = generate_one_mtx_parallel(nRLC, A, b, dt, maxORD,nWorkers);

       
    %% comparison 
    %disp('Comparison')
    %tf = isequal(A_one,A_one_serial)
    
    %% MTSM solution
    % MTSM - serial
    tic;
    [T,Y,ORD] = MTSM_explicit_linear_serial(dt,tspan,init,A,b,eps,maxORD);
    MTSM_t_serial = toc;
    
%     % MTSM - serial - 3D matrix 
%     tic;
%     [T_3D,Y_3D,ORD_3D] = MTSM_explicit_linear_precalc_3D_mtx(dt,tspan,init,A_3D,B_3D,eps,maxORD);
%     MTSM_t_serial_3D = toc;

    % MTSM - serial - one matrix A_hat 
    tic;
    [T_2D_serial,Y_2D_serial,ORD_2D_serial] = MTSM_explicit_linear_precalc_one_mtx(dt,tspan,init,A_one_serial,b_one_serial,eps,maxORD);
    MTSM_t_serial_one = toc;

    % MTSM - parallel - one matrix A_hat 
    tic;
    [T_2D_parallel,Y_2D_parallel,ORD_2D_parallel] = MTSM_explicit_linear_precalc_one_mtx(dt,tspan,init,A_one_parallel,b_one_parallel,eps,maxORD);
    MTSM_t_parallel_one = toc;

    u1=Y(end,1:(nRLC-3)/2);
    u2=Y(end,(nRLC-3)/2+1:end-3);
    
    %     
    % MTSM - parallel - one matrix A_hat 
%     tic;
%     [T_2D_parallel,Y_2D_parallel,ORD_2D_parallel] = MTSM_explicit_linear_precalc_one_mtx_parallel(dt,tspan,init,A_one_parallel,b_one_parallel,eps,maxORD,S);
%     MTSM_t_parallel_one = toc;
    
   

%     %% Plot: MTSM serial
%     % MTSM - UC1
%     figure;
%     %subplot(2,1,1);
%     plot(T, Y(:,1), '-r', 'LineWidth',1.5);
%     grid on;
%     legend('UC1')
% 
%     hold on; 
% 
%     % MTSM - UC100
%     plot(T, Y(:,nRLC), '--b','LineWidth',1.5); % uc last
%     grid on;
%     legend('UC1' ,['UC',num2str(nRLC,'%u')])
%     TITLE=sprintf('MTSM Serial Original Approach');
%     title(TITLE);
% 
%     %% Plot: MTSM serial - 3D matrix
% %     % MTSM - UC1
% %     figure;
% %     %subplot(2,1,1);
% %     plot(T_3D, Y_3D(:,1), '-r', 'LineWidth',1.5);
% %     grid on;
% %     legend('UC1')
% 
% %     hold on; 
% % 
% %     % MTSM - UC100
% %     plot(T_3D, Y_3D(:,nRLC), '--b','LineWidth',1.5); % uc last
% %     grid on;
% %     legend('UC1' ,['UC',num2str(nRLC,'%u')])
% %     TITLE=sprintf('MTSM Serial 3D matrix');
% %     title(TITLE);
%  
%     %% Plot: MTSM serial - one matrix
%     % MTSM - UC1
%     figure;
%     %subplot(2,1,1);
%     plot(T_2D_serial, Y_2D_serial(:,1), '-r', 'LineWidth',1.5);
%     grid on;
%     legend('UC1')
% 
%     hold on; 
% 
%     % MTSM - UC100
%     plot(T_2D_serial, Y_2D_serial(:,nRLC), '--b','LineWidth',1.5); % uc last
%     grid on;
%     legend('UC1' ,['UC',num2str(nRLC,'%u')])
%     TITLE=sprintf('MTSM Serial One Matrix');
%     title(TITLE);    
%     
%      %% Plot: MTSM parallel - one matrix
%     % MTSM - UC1
%     figure;
%     %subplot(2,1,1);
%     plot(T_2D_parallel, Y_2D_parallel(:,1), '-r', 'LineWidth',1.5);
%     grid on;
%     legend('UC1')
% 
%     hold on; 

    % MTSM - UC100
%     plot(T_2D_parallel, Y_2D_parallel(:,nRLC), '--b','LineWidth',1.5); % uc last
%     grid on;
%     legend('UC1' ,['UC',num2str(nRLC,'%u')])
%     TITLE=sprintf('MTSM Parallel One Matrix');
%     title(TITLE);    
% 
%     %% Plots: sparse matrices
%     % Original A matrix
%     figure;
%     spy(A);
%     title('MTSM Original Matrix');
%     nz = nnz(A);
%     pct = 100 / numel(A);
%     xlabel(sprintf('nonzeros = %d (%.3f%%)',nz,nz*pct));
%     %set(gca,'XTick',1:1:2*(S+N-1));
%     %set(gca,'YTick',1:1:2*(S+N-1));
%     grid on;
%     set(gca, 'GridLineStyle', '-');
%     
% %     % 3D matrix 
% %     figure;
% %     spy(A_3D);
% %     title('MTSM Serial 3D Matrix');
% %     nz = nnz(A_3D);
% %     pct = 100 / numel(A_3D);
% %     xlabel(sprintf('nonzeros = %d (%.3f%%)',nz,nz*pct));
% %     %set(gca,'XTick',1:1:2*(S+N-1));
% %     %set(gca,'YTick',1:1:2*(S+N-1));
% %     grid on;
% %     set(gca, 'GridLineStyle', '-');
%     
%     % One big matrix (Ahat) serial
%     figure;
%     spy(A_one_serial);
%     title('MTSM Serial One Matrix');
%     nz = nnz(A_one_serial);
%     pct = 100 / numel(A_one_serial);
%     xlabel(sprintf('nonzeros = %d (%.3f%%)',nz,nz*pct));
%     %set(gca,'XTick',1:1:2*(S+N-1));
%     %set(gca,'YTick',1:1:2*(S+N-1));
%     grid on;
%     set(gca, 'GridLineStyle', '-');
%     
%     % One big matrix (Ahat) parallel
%     figure;
%     spy(A_one_parallel);
%     title('MTSM Parallel One Matrix');
%     nz = nnz(A_one_parallel);
%     pct = 100 / numel(A_one_parallel);
%     xlabel(sprintf('nonzeros = %d (%.3f%%)',nz,nz*pct));
%     %set(gca,'XTick',1:1:2*(S+N-1));
%     %set(gca,'YTick',1:1:2*(S+N-1));
%     grid on;
%     set(gca, 'GridLineStyle', '-');
    
    %% return measured times
    MTSM_serial = MTSM_t_serial;
    MTSM_one_serial = MTSM_t_serial_one++T_prealloc_serial;
    MTSM_one_parallel = MTSM_t_parallel_one+T_prealloc_parallel;
    
    %% Console print
    % Time of solutions
%     fprintf('\n*** TELEGRAPH LINE PROBLEM : linear system [%u, %u] *** \n', size(A,1), size(A,2));
%     fprintf('Time of solution Serial MTSM solver: %d seconds \n', MTSM_t_serial);
% %     fprintf('Time of solution Serial MTSM solver - 3D mtx: %d seconds \n', MTSM_t_serial_3D);
%     fprintf('Time of solution Serial MTSM solver - one mtx: %d seconds \n', MTSM_t_serial_one+T_prealloc_serial);
%     fprintf('Time of solution Parallel MTSM solver - one mtx: %d seconds \n', MTSM_t_parallel_one+T_prealloc_parallel);
% 
%     % errors 
%      fprintf('\n||MTSM serial - MTSM serial one mtx||: %g\n', norm(Y(end,:)-Y_2D_serial(end,:))/(norm(Y_2D_serial(end,:))+1));
%      fprintf('\n||MTSM serial - MTSM parallel one mtx||: %g\n', norm(Y(end,:)-Y_2D_parallel(end,:))/(norm(Y_2D_parallel(end,:))+1));
% %     fprintf('||MTSM serial - MTSM serial 3D mtx||: %g\n', norm(Y(end,:)-Y_3D(end,:))/(norm(Y_3D(end,:))+1));
% 
end