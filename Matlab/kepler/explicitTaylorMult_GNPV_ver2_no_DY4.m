function [t,y,ORD,DY_all]=explicitTaylorMult_GNPV_ver2_no_DY4(h,tspan,init,A,A2,A3,A5,b,ij,ijk,ijklm,eps,ind,maxORD)
%  [t,y,ORD]=explicitTaylor(h,tspan,init,A,B3,C,d,eps)
%  Numerical solution of ODEs with explicit Taylor series method
%  problem y' = A*y^2 + B3*y_iy_jy_k + C*y + d
%  problem y' = A*y + A2*y_ij + A3*y_ijk + A4*y_ijkl + A5*y_ijklm
%  where A is linear term matrix,
%        A2 quadrature term matrix,
%        A3 cubature term combinations matrix,
%        ....
%        ij,ijk,... indexes on y_i*y_j... multiplications
%  and   b is a vector of constants (right-hand side)
%  ---------------------------------------------------------------
%  h ... integration step size
%  tmax ... interval of numerical solution [0,tmax]
%  ORD ... order of Taylor series - number of Taylor series Terms
%  DY ... cell of all derivatives in each time step DY [y; dy; ddy; d3y; ...]

% VERSION 2: recurrent calculation of Taylor series terms among the same
% number of multiplication members - specially, DY3_terms, DY4_terms,
% DY5_terms (further optimalization for DY4 and DY5

ne=length(init); % number of equations

steps = round((tspan(2)-tspan(1))/h); % number of integration steps

% preallocation ORD, DY_all
ORD=zeros(steps*30,1); % 30 possibilities of halfing integration stepsize
DY_all=[];

% stoping rule: maximum number DY terms are smaller then eps
stopping=3;

% max ORD for double precision
% maxORD=30;

% last step tolerance
ls_tol = 10^-10;

%% precalculation_DY_mult (we do it in simulate.m script)
% --> precalculation_DY_mult(maxORD)
% ind=load('DY_indexes_maxORD_25.mat','DY_ij','DY_ijk','DY_ijkl','DY_ijklm');
% ind=load('DY_indexes_maxORD_65.mat','DY_ij','DY_ijk','DY_ijkl','DY_ijklm');

%% indexes on multiplications terms
i2=ij(:,1);
j2=ij(:,2);

i3=ijk(:,1);
j3=ijk(:,2);
k3=ijk(:,3);

i5=ijklm(:,1);
j5=ijklm(:,2);
k5=ijklm(:,3);
l5=ijklm(:,4);
m5=ijklm(:,5);


%% MTSM calculation
i=1; % timestep index
t(i)=tspan(1,1);
y(:,1)=init'; % initial values

i=i+1;


while t(i-1)+ls_tol<tspan(1,2)
    DY=zeros(ne,maxORD+1);
    
    % already calculated submatrices
    DY_3terms = zeros(size(ijk,1), maxORD+stopping);
    
    DY_5terms = zeros(size(ijklm,1), maxORD+stopping);
    
    % TODO: arrays for already calculated submatrices of the currently
    % calculated term
    DY_5terms_smaller = zeros(size(ijklm,1), maxORD+stopping);
    
    y(:,i)=y(:,i-1); % first term of Taylor series y_{i+1}=y_{i}+...
    DY(:,1)=y(:,i);
%     for n=1:ne
%        DY(1,n)=h*(sum(A(n,:).*y(i,:))+b(n));
%     end
       
    % linear term
    Ay=A*DY(:,1)+b;
    
    % TODO: for DY4 and DY5, save two / three multiplicications into the
    % new array (1 1 / 1 1 1)
    DY_5terms_smaller(:,1) = DY(ijklm(:,2),1).*DY(ijklm(:,3),1).*DY(ijklm(:,4),1);
    
    % computation of DY1 (indeces 1 1 1 1 / 1 1 1 1 1) 
    DY_3terms(:,1) = DY(ijk(:,2),1).*DY(ijk(:,3),1);
    DY_5terms(:,1) =  DY_5terms_smaller(:,1).*DY(ijklm(:,5),1);
    
    A2y=A2*(DY(ij(:,1),1).*DY(ij(:,2),1));
    A3y=A3*(DY(ijk(:,1),1).*DY_3terms(:,1));
    A5y=A5*(DY(ijklm(:,1),1).*DY_5terms(:,1));
    
    DY(:,2)=h*(Ay+A2y+A3y+A5y); % first derivative  
    y(:,i)=y(:,i)+DY(:,2); % first term (first derivative)

    maxDY = [ones(1,stopping)*10^10 zeros(1,maxORD+stopping)];

    k=2;
    
    % first columns indeces of submatrices with new values
    t3_new_start = 2;
    t5_new_start = 2;
    
    % TODO: indeces that indicate the start of the completely new coeficients
    % can be calculated using the number of newly added terms in every
    % h
    t5_new_start_smaller = 3;

    
    while norm(maxDY(k-1:stopping+k-2))>eps
        
        % linear term
        Ay=A*DY(:,k); % opraveno Vasek 23.10.2017
        % -------------------------------------------

        % 2 terms
        A2y=A2*sum(DY(i2,ind.DY_ij{1,k-1}).*DY(j2,ind.DY_ij{2,k-1}),2);
        
        % -------------------------------------------
        % unique 1st row of TS terms indeces (equivalent with unique(DY_ijk..{1,k-1})
        k_array = k:-1:2;
        % 3 terms - original solution
        %A3y=A3*sum(DY(i3,ind.DY_ijk{1,k-1}).*DY(j3,ind.DY_ijk{2,k-1}).*DY(k3,ind.DY_ijk{3,k-1}),2);
        
        % 3 terms - PVGN solution
        % already calculated terms - we use the submatrices and multiply them by the unique 1st DY_ijk row  
        A3y_calculated = sum(DY(i3,k_array).*DY_3terms(:,1:k-1),2);
        % new terms - new submatrix multiplied only by the corresponding
        % rows from DY, we do not consider DY_ijk indeces, because there
        % are always ones 
        DY_3terms(:,k) = sum(DY(j3,ind.DY_ijk{2,k-1}(t3_new_start:end)).*DY(k3,ind.DY_ijk{3,k-1}(t3_new_start:end)),2);
%         k
%         r2 = ind.DY_ijk{2,k-1}(i3_new_start:end)
%         r3 = ind.DY_ijk{3,k-1}(i3_new_start:end)
        
        % new  calculated terms have to be multiplied with the 1st row 
        A3y_new = DY(i3,1).*DY_3terms(:,k);
        A3y=A3*sum([A3y_calculated,A3y_new],2);
        % new column index - number of columns of the previous one
        t3_new_start = size(ind.DY_ijk{1,k-1},2)+1;             
        
        % -------------------------------------------
        % 5 terms  - PVGN solution
        %A5y=A5*sum(DY(i5,ind.DY_ijklm{1,k-1}).*DY(j5,ind.DY_ijklm{2,k-1}).*DY(k5,ind.DY_ijklm{3,k-1}).*DY(l5,ind.DY_ijklm{4,k-1}).*DY(m5,ind.DY_ijklm{5,k-1}),2);
        % already calculated terms - we use the submatrices and multiply them by the 1st DY_ijklm row  
        A5y_calculated = sum(DY(i5,k_array).*DY_5terms(:,1:k-1),2);
        % new terms - new submatrix multiplied only by the corresponding
        % rows from DY, we do not consider DY_ijk indeces, because there
        % are always ones 
        A5y_calculated_smaller = sum(DY(m5,k_array).*DY_5terms_smaller(:,1:k-1),2);
        
        DY_5terms_smaller(:,k) = sum(DY(j5,ind.DY_ijklm{2,k-1}(t5_new_start_smaller:end)).*DY(k5,ind.DY_ijklm{3,k-1}(t5_new_start_smaller:end)).*DY(l5,ind.DY_ijklm{4,k-1}(t5_new_start_smaller:end)),2);    
        
        DY_5terms_new_part = DY_5terms_smaller(:,k).*DY(m5,1);
        
        DY_5terms(:,k) = sum([A5y_calculated_smaller,DY_5terms_new_part],2);
        A5y_new = DY(i5,1).*DY_5terms(:,k);
        A5y=A5*sum([A5y_calculated,A5y_new],2);
        
        t5_new_start_smaller = t5_new_start;
        t5_new_start = size(ind.DY_ijklm{1,k-1},2)+1;
        t5_new_start_smaller = (t5_new_start - t5_new_start_smaller) + t5_new_start;
         
        % -------------------------------------------
%         fprintf('DEBUG: k = %d\n', k);
%         DY(:,k+1)=(h/k)*(Ay+A2y+A3y+A4y+A5y);
%         VS = DY(:,1:k+1)
%         
         DY(:,k+1)=(h/k)*(Ay+A2y+A3y+A5y);
         %GNPV = DY(:,1:k+1)
         %return
%         GNPV = DY(:,1:k+1)
%         
%         res = VS - GNPV
%         
%         A3_diff = (A3y - A3y_GN)'
%         A4_diff = (A4y - A4y_GN)'
%         A5_diff = (A5y - A5y_GN)'
        

        y(:,i)=y(:,i)+DY(:,k+1);
        
        maxDY(k+stopping-1)=max(abs(DY(:,k+1)));
        k=k+1;
        
        if k > maxORD % max ORD was reached
            break;
        end
        
    end

    if k > maxORD % max ORD was reached
        h=h/2; % we use half-size integration step
        fprintf('Halving integration step size.\n')
    else
%         DY_all{i-1}=DY;
%         if k<(stopping+4)
%             h=h*2;
%             Hk=(h.^iPascal_mixed)./fact_k;
%             fprintf('Increasing integration step size.\n')
%         end
        %DY(:,1:k+1)
        ORD(i)=k;
        t(i)=t(i-1)+h;
        i=i+1;
        %return
    end
end

% shrink the rest (empty) array
ORD(i:end)=[];


