function simulate2(nRLC)
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
    close all
    clc

    L = 10^-8;     % inductance (H, Henry)
    C = 10^-12;    % capacitance (F, Farad)
    om = 3*10^9;   % omega (rad/s)

    % total delay of input signal
    tdelay = nRLC*sqrt(L*C);

    %% for now, just adjusted line
    tmax = 2*tdelay;
    R1 = 100; % input load (Ohm)
    R2 = 100; % output load (Ohm)
    eps = 1e-8; % presicion of the calculation
    dt = sqrt(L*C); % size of the integration step

    % simulation time
    tspan = [0, tmax];  
    maxORD = 2;

    %% initial conditions
    % set initial condition for the input function
    init = zeros(2*nRLC+2,1);
    init(2*nRLC+2) = 1; 

    %% assembly of matrix A
    % assembly_A_b;
    [A, b] = assembly_A_b_sparse(nRLC, L, C, om, R1, R2);
    
    size(A)
    
    %% generate different matrix representations for parallel solution
    % Approach 1: 3D matrix - precalculation of the matrix A for the higher derivatives
    % each Taylor series term [(h^i/i!) * A^i] for i=1..maxORD is in special
    % matrix, so we have an array of these matrices -> 3D matrix 
    [A_3D, B_3D] = generate_3D_mtx(nRLC, A, b, dt, maxORD); 

    % Approach 2: one big matrix A_hat
    % all Taylor series terms [(h^i/i!) * A^i] for i=1..maxORD are added
    % together, so we obtain one big matrix 
    [A_one, b_one] = generate_one_mtx_parallel(nRLC, A, b, dt, maxORD);

    A_one
    % serial version
    [A_one_serial, b_one_serial] = generate_one_mtx(nRLC, A, b, dt, maxORD);
    
    A_one_serial
    
    figure;
    spy(A_one_serial);
    title('MTSM One Matrix (Ahat) serial');
    nz = nnz(A_one_serial);
    pct = 100 / numel(A_one_serial);
    xlabel(sprintf('nonzeros = %d (%.3f%%)',nz,nz*pct));
    %set(gca,'XTick',1:1:2*(S+N-1));
    %set(gca,'YTick',1:1:2*(S+N-1));
    grid on;
    set(gca, 'GridLineStyle', '-');
    
    
    %% comparison 
    disp('Comparison')
    tf = isequal(A_one,A_one_serial)
    
    %% MTSM solution
    % MTSM - parallel - 3D matrix 
%     tic;
%     [T_3D,Y_3D,ORD_3D] = MTSM_explicit_linear_precalc_3D_mtx(dt,tspan,init,A_3D,B_3D,eps,maxORD);
%     MTSM_t_parallel_3D = toc;

%     % MTSM - parallel - ona matrix A_hat 
%     tic;
%     [T_2D,Y_2D,ORD_2D] = MTSM_explicit_linear_precalc_one_mtx(dt,tspan,init,A_one,b_one,eps,maxORD);
%     MTSM_t_parallel_One = toc;
%     
%     % MTSM - serial - ona matrix A_hat 
%     tic;
%     [T_2D_serial,Y_2D_serial,ORD_2D_serial] = MTSM_explicit_linear_precalc_one_mtx(dt,tspan,init,A_one_serial,b_one_serial,eps,maxORD);
%     MTSM_t_parallel_One = toc;
% 
%     % MTSM - serial
%     tic;
%     [T_serial,Y_serial,ORD_serial] = MTSM_explicit_linear_serial(dt,tspan,init,A,b,eps,maxORD);
%     MTSM_t_serial = toc;
% 
%     %% Plot: MTSM serial
%     % MTSM - UC1
%     figure;
%     %subplot(2,1,1);
%     plot(T_serial, Y_serial(:,1), '-r', 'LineWidth',1.5);
%     grid on;
%     legend('UC1')
% 
%     hold on; 
% 
%     % MTSM - UC100
%     plot(T_serial, Y_serial(:,nRLC), '--b','LineWidth',1.5); % uc last
%     grid on;
%     legend('UC1' ,['UC',num2str(nRLC,'%u')])
%     TITLE=sprintf('MTSM Serial');
%     title(TITLE);
% 
% %     %% Plot: MTSM parallel - 3D matrix
% %     % MTSM - UC1
% %     figure;
% %     %subplot(2,1,1);
% %     plot(T_3D, Y_3D(:,1), '-r', 'LineWidth',1.5);
% %     grid on;
% %     legend('UC1')
% % 
% %     hold on; 
% % 
% %     % MTSM - UC100
% %     plot(T_3D, Y_3D(:,nRLC), '--b','LineWidth',1.5); % uc last
% %     grid on;
% %     legend('UC1' ,['UC',num2str(nRLC,'%u')])
% %     TITLE=sprintf('MTSM Parallel 3D matrix');
% %     title(TITLE);
%  
%     %% Plot: MTSM parallel - one matrix
%     % MTSM - UC1
%     figure;
%     %subplot(2,1,1);
%     plot(T_2D, Y_2D(:,1), '-r', 'LineWidth',1.5);
%     grid on;
%     legend('UC1')
% 
%     hold on; 
% 
%     % MTSM - UC100
%     plot(T_2D, Y_2D(:,nRLC), '--b','LineWidth',1.5); % uc last
%     grid on;
%     legend('UC1' ,['UC',num2str(nRLC,'%u')])
%     TITLE=sprintf('MTSM Parallel One Matrix');
%     title(TITLE);    
%     
%      %% Plot: MTSM parallel - one matrix
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

    %% Plots: sparse matrices
    % Original A matrix
    figure;
    spy(A);
    title('MTSM - Original Matrix');
    nz = nnz(A);
    pct = 100 / numel(A);
    xlabel(sprintf('nonzeros = %d (%.3f%%)',nz,nz*pct));
    %set(gca,'XTick',1:1:2*(S+N-1));
    %set(gca,'YTick',1:1:2*(S+N-1));
    grid on;
    set(gca, 'GridLineStyle', '-');
    
    % 3D matrix 
    figure;
    spy(A_3D);
    title('MTSM Parallel 3D Matrix');
    nz = nnz(A_3D);
    pct = 100 / numel(A_3D);
    xlabel(sprintf('nonzeros = %d (%.3f%%)',nz,nz*pct));
    %set(gca,'XTick',1:1:2*(S+N-1));
    %set(gca,'YTick',1:1:2*(S+N-1));
    grid on;
    set(gca, 'GridLineStyle', '-');
    
    % One big matrix (Ahat)
    figure;
    spy(A_one);
    title('MTSM One Matrix (Ahat)');
    nz = nnz(A_one);
    pct = 100 / numel(A_one);
    xlabel(sprintf('nonzeros = %d (%.3f%%)',nz,nz*pct));
    %set(gca,'XTick',1:1:2*(S+N-1));
    %set(gca,'YTick',1:1:2*(S+N-1));
    grid on;
    set(gca, 'GridLineStyle', '-');
    
%     %% Console print
%     % Time of solutions
%     fprintf('\n*** TELEGRAPH LINE PROBLEM : linear system [%u, %u] *** \n', size(A,1), size(A,2));
%     fprintf('Time of solution Serial MTSM solver: %d seconds \n', MTSM_t_serial);
% %     fprintf('Time of solution Parallel MTSM solver - 3D mtx: %d seconds \n', MTSM_t_parallel_3D);
%     fprintf('Time of solution Parallel MTSM solver - One mtx: %d seconds \n', MTSM_t_parallel_One);
% 
%     % errors 
%     fprintf('\n||MTSM serial - MTSM parallel One mtx||: %g\n', norm(Y_serial(end,:)-Y_2D(end,:))/(norm(Y_2D(end,:))+1));
% %     fprintf('||MTSM serial - MTSM parallel 3D mtx||: %g\n', norm(Y_serial(end,:)-Y_3D(end,:))/(norm(Y_3D(end,:))+1));

end