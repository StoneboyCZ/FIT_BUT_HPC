function [A_one, b_one] = generate_one_mtx_eig(nRLC, A, b, dt, maxORD)
%% generate_one_mtx_eig 
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
    A_one = zeros(2*nRLC+2,2*nRLC+2);
    b_one = zeros(2*nRLC+2,1);
    I = eye(2*nRLC+2,2*nRLC+2);

    Ynominator = A*dt;
    Bnominator = b*dt;
    denominator = 1;
    
    [V,D] = eigs(A);
   
    for i=1:1:maxORD 
        % we sum all Taylor series terms together to create one big matrix
        % A_one
        A_one = A_one + Ynominator/denominator;
        b_one = Bnominator/denominator;

        Ynominator = V*D^i*V' * dt; %Ynominator*(A*dt);
        Bnominator = A*dt*Bnominator;

        denominator = (i+1)*denominator;
    end
    % we add identity matrix I, because the Taylor series has the form: 
    % y_i+1 = y_i + dt*A + ... 
    % we have to add "y_i"
    A_one = A_one + I;
    
end