function [A_one, b_one,T,S] = generate_one_mtx_parallel(nRLC, A, b, dt, maxORD,nWorkers)
%% generate_one_mtx_parallel 
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

    %A_one = zeros(2*nRLC+2,2*nRLC+2);
    %b_one = zeros(2*nRLC+2,1);
    I = eye(2*nRLC+2,2*nRLC+2);
    
   
    % precalculate factorials 1! to maxORD!
    denominators = zeros(maxORD,1);
    denominators(1) = 1;
    for i=2:maxORD
        denominators(i) = i*denominators(i-1);    
    end
   
 
    m = size(A,1); % number of rows of the input matrix
    n = size(A,2); % number of cols of the input matrix

    %% start parallel pool and set number of workers
    tic
    spmd(nWorkers)
      
        ilo = floor ( ( ( labindex ( ) - 1 ) * n ) / numlabs ( ) ) + 1;
        ihi = floor ( (   labindex ( )       * n ) / numlabs ( ) );  
        %fprintf('%d %d\n', ilo, ihi);
        % get block of the matrix A
        A_one_block = get_block(A, ilo, ihi);
        A_one_block_tmp = sparse(size(A_one_block,1), size(A_one_block,2));
        % save how the matrix was split
%         split = [ilo,ihi];
        
        
        % compute power of the matrix block 
        Ynominator = A_one_block * dt;
        for i=1:1:maxORD
             A_one_block_tmp = A_one_block_tmp + Ynominator/denominators(i);
             Ynominator = Ynominator*(A*dt);
        end
    % spmd end
     
    end
    T = toc;
    %% assembly powered blocks together
    A_one = vertcat(A_one_block_tmp{:});
    A_one = A_one + I;
    
    %% vector b
    b_one = A_one*b;
    
    % assembly the split (to know how the matrix was spit)
%     S = vertcat(split{:});
end

%============================================
% get some rows of the matrix A  (removed col indexing (j))
function Ablock = get_block (A, ilo, ihi)
  Ablock = A(ilo:ihi,:);
  
  %Ablock = sparse( ihi + 1 - ilo, jhi + 1 - jlo);


  %i = ilo:ihi;

  %for j = jlo : jhi
    %Ablock(i+1-ilo,j+1-jlo) = A(i,j);
    
  %end

end

 