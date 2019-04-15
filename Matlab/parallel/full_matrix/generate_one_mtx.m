function [A_one, b_one,T] = generate_one_mtx(nRLC, A, b, dt, maxORD)
%% generate_one_mtx 
% PARAMETERS:
% ------------
% nRLC  number of RLC segments
% A         original matrix A 
% b         original vector b 
% dt        integration step size
% maxORD    maximum ORD for MTSM
%
% RETURN:
% ------------
% A_one     all Taylor series terms [(h^i/i!) * A^i] for i=1..maxORD in one
%           matrix 
% b_one     vector b
    A_one = zeros(nRLC,nRLC);
    b_one = zeros(nRLC,1);
    I = eye(nRLC,nRLC);

    Ynominator = A*dt;
    Bnominator = b*dt;
    denominator = 1;
   
    tic
    for i=1:1:maxORD 
        % we sum all Taylor series terms together to create one big matrix
        % A_one
        A_one = A_one + Ynominator/denominator;
        b_one = Bnominator/denominator;

        Ynominator = Ynominator*(A*dt);
        Bnominator = A*dt*Bnominator;

        denominator = (i+1)*denominator;
    end
    T = toc;
    % we add identity matrix I, because the Taylor series has the form: 
    % y_i+1 = y_i + dt*A + ... 
    % we have to add "y_i"
    A_one = A_one + I;
    
end