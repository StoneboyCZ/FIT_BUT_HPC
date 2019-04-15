function [A_3D, B_3D] = generate_3D_mtx(nRLC, A, b, dt, maxORD)
%% generate_3D_mtx 
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
% A_3D      3D matrix with all Taylor series terms [(h^i/i!) * A^i] for i=1..maxORD 
% B_3D      3D vector

    % 3D matrix for the result
    A_3D = zeros(2*nRLC+2,2*nRLC+2,200); % 200 - maxORD 
    B_3D = zeros(2*nRLC+2,1,200);

    Ynominator = A*dt;
    Bnominator = b*dt;
    denominator = 1;
    for i=1:1:maxORD  
        % first derivative will have the b vector
        A_3D(:,:,i) = Ynominator/denominator;
        B_3D(:,:,i) = Bnominator/denominator;

        Ynominator = Ynominator*(A*dt);
        Bnominator = A*dt*Bnominator;

        denominator = (i+1)*denominator;
    end

end