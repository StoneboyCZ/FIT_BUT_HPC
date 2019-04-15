function [A, b] = assembly_A_b_sparse(nRLC, L, C, om, R1, R2)
%% assembly_A_b_sparse
% Assembles sparse block matrix A for currents and voltages and 
% also creates vector b (zeros in our case) 
% PARAMETERS:
% ---------------
% nRLC      number of segments of the line 
% L         inductance (H, Henry)
% C         capacitance (C, Farad)
% om        omega (rad/s)
% R1        input load (Ohm)
% R2        output load (Ohm)
%
% RETURN:
% ------------
% A         sparse block matrix A = [A11 A12; A21 A22]
%           matrices for currents and voltages 
% b         vector b 

    % sparse matrix A
    A=sparse(nRLC*2+2,nRLC*2+2);  
    b=zeros(nRLC*2+2,1);      

    % A = ( A11  A12 )
    %     ( A21  A22 )
    A11=spalloc(nRLC,nRLC,1);
    A12=spalloc(nRLC,nRLC+2,1);
    A21=spalloc(nRLC+2,nRLC,1);
    A22=spalloc(nRLC+2,nRLC+2,1);

    %% load coefficients
    e = ones(nRLC, 1);
    A11(nRLC,nRLC)=-1/(R2*C);
    A12(1:nRLC, 1:nRLC) = spdiags([e*1/C e*-1/C], [0 1], nRLC, nRLC);
    A21(1:nRLC, 1:nRLC) = spdiags([e*1/L e*-1/L], [-1 0], nRLC, nRLC);
    A22(1,1) = -R1/L;

    % i = 1;
    % for i=1:nRLC
    %     if (i~=1) && (i~=nRLC)
    %         A12(i,i:i+1)=[1/C,-1/C];
    %         A21(i,i-1:i)=[1/L,-1/L];
    %     elseif i==1 % first row
    %         A12(i,i:i+1)=[1/C,-1/C];
    %         A21(i,i)=-1/L;
    %         A22(i,i)=-R1/L;
    %     elseif i==nRLC % last row
    %         A11(i,i)=-1/(R2*C); 
    %         A12(i,i)=1/C;
    %         A21(i,i-1:i)=[1/L,-1/L];
    %     end
    % end

    %% load the function for u0
    A22(1,nRLC+1) = 1/L;

    % equations for sine
    A22(nRLC+1,nRLC+2) = om;
    A22(nRLC+2,nRLC+1) = -om;

    % final matrix
    A=[A11, A12;
      A21, A22];

    %full(A)
    %spy(A)
    %grid on;
    %title('vysledna matice');
end