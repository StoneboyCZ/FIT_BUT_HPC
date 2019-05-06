function [A_one, b_one] = generate_one_mtx_parallel(nRLC, A, b, dt, maxORD)
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

    A_one = zeros(2*nRLC+2,2*nRLC+2);
    b_one = zeros(2*nRLC+2,1);
    I = eye(2*nRLC+2,2*nRLC+2);
    dt_orig = dt;
    
    % precalculate factorials 1! to maxORD!
    denominators = zeros(maxORD,1);
    denominators(1) = 1;
    for i=2:maxORD
        denominators(i) = i*denominators(i-1);    
    end
   
 
   m = size(A,1); % number of rows of the input matrix
   n = size(A,2); % number of cols of the input matrix

  workers = 2;
  %parpool ( 'open', 'local', workers )
  poolobj = parpool('local', workers);


  
  spmd
     ilo = floor ( ( ( labindex ( ) - 1 ) * n ) / numlabs ( ) ) + 1;
     ihi = floor ( (   labindex ( )       * n ) / numlabs ( ) );  
     jlo = 1;
     jhi = m; 
%     ilo = 1;
%     ihi = m;
%     jlo = floor ( ( ( labindex ( ) - 1 ) * n ) / numlabs ( ) ) + 1;
%     jhi = floor ( (   labindex ( )       * n ) / numlabs ( ) );
     fprintf('%d %d %d %d\n', ilo, ihi, jlo, jhi);
     A_one_block = get_block(A, ilo, jlo, ihi, jhi);
    % A_one_block
     
   
%      figure;
%      spy(A_one_block)
%      title('block')
     %A_one_block = zeros(size(A_block,1), size(A_block,2));
  
     % compute power of the matrix block 
    
    
    
%     A_block_T = A_block';
%     Ynominator = A_block*dt;
%     
%     size(A_block)
%     size(A_block_T)
%    size(A_one_block)
%    
%     
%     denominator = 1;
%     size(A_one_block)
%     size(Ynominator/denominator)
%     
    %%  for i=1:1:maxORD
        %A_one_block_3D(:,:,i) = A_block^i;%*dt/denominators(i);
        
        %dt = dt*dt;
        %Ynominator = Ynominator*(A_block*dt);
     %%   A_one_block = A_one_block*dt/denominators(i) * A;
     %%   dt = dt*dt;
     %% end
  end
  
  
  
  %  A_one = [ A_one_block{:} ];
    A_one = [ A_one_block{1};  A_one_block{2}];
    %A_one = reshape(A_one,[2*nRLC+2,2*nRLC+2]);
    figure;
    spy(A_one)
    title('A one');
    
    figure;
    spy(A_one_block{1})
    title('A one block 1');
    
    figure;
    spy(A_one_block{2})
    title('A one block 2');
    
  output = A_one_block{:}
  
  %dt = dt_orig;
  %A_one = [ A_one_block{1}; A_one_block{2}; A_one_block{3}; A_one_block{4}]
  
  %spy(A_one_block{2});
  %A_one_block{2}
  
  %A_one = [ A_one_block{ : } ]';
  %Ares = [ A_block{:} ];
  %A_one_res = [ A_one_block{:} ];
%   figure;
%   spy(A_one_res)
%   grid on;
%   title('Joined matrix');
%   
%   figure;
%   spy(Ablock{1})
%   grid on;
%   
%   figure;
%   spy(Ablock{2})
%   grid on;

   % return
   
   delete(poolobj);
end

%============================================
% ziskani prislusneho bloku z matice A 
function Ablock = get_block (A, ilo, jlo, ihi, jhi)
  Ablock = zeros ( ihi + 1 - ilo, jhi + 1 - jlo );

  i = ilo:ihi;

  for j = jlo : jhi
    Ablock(i+1-ilo,j+1-jlo) = A(i,j);
  end

  %return
end

    %===========================================
    % all workers will calculate
%     spmd 
%           A_d = codistributed(A,codistributor('1d',1));
%           A_local = getLocalPart(A_d);
%           A_one_d = codistributed(A_one,codistributor('1d',1));  
%           A_one_local = getLocalPart(A_one_d);
%           
% %         A_d_local = getLocalPart(A_d);
% %         A_one_local = getLocalPart(A_one_d);
%           size(A_local)
%           size(A_one_local)
%           
%           A_one_local = A_one_local';
%          
%          Ynominator = A_local*dt;    
%          %size(A_d_local)
% %         size(Ynominator/denominators(1))
% % 
%          %for i=1:1:maxORD
% % 
%          %    A_one_local = A_one_local + Ynominator/denominators(i);
%          %    disp('A_one_local_size')
%          %    size(A_d_local)
%          %    disp('Ynominator')
%          %    size(Ynominator)
% 
%           %   Ynominator = Ynominator*(A_d_local*dt);
%          %end
% 
% %         %%%
% %         A_one = A_one + Ynominator/denominator;
% %         b_one = Bnominator/denominator;
% % 
% %         Ynominator = Ynominator*(A*dt);
% %         Bnominator = A*dt*Bnominator;
% % 
% %         denominator = (i+1)*denominator;
%     end

    %A_one = gather(A_one_local);
    