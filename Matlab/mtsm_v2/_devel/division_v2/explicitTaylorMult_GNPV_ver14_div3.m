function [t,y,ORD,DY_all]=explicitTaylorMult_GNPV_ver14_div3( ... 
    h, ... size of the integration step
    tspan, ... times of calculation
    init, ... initial condition
    A, ... Jacobian for the linear part of the system
    A2, ... quadrature term matrix,
    A3, ... three term multiplications combination matrix
    A4, ... four term multiplications combination matrix
    A5, ... five term multiplications combination matrix
    b, ...  vector of constants (right-hand side)
    ij, ...  indexes on y_i*y_j multiplications
    ijk, ... indexes on y_i*y_j*y_k multiplications
    ijkl, ... indexes on y_i*y_j*y_k*y_l multiplications
    ijklm, ... indexes on y_i*y_j*y_k*y_l*y_m multiplications
    eps, ...  accuracy
    maxORD, ... maximum order of the Taylor series
    minORD, ... minimum ord needed for step increase 
    hScaleFactor ... how h is scaled (hnew = h*hscalefactor) when ORD is stable
    )
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

% VERSION 4: last step fix, while optimization

% VERSION 5: removal of cell array

% VERSION 6: removal of ind variable, misc comments

% VERSION 7: increase of ORD if it´s low enough long enough

% VERSION 8: settings for the ORD increase 

% VERSION 9: placeholder

% VERSION 10: optimized ORD increase

% VERSION 11: fixed ORD increase, optimized ORD handling, fixed accuracy

% VERSION 14: updated handling of the step size scaling, optimizations

ne=length(init); % number of equations

steps = round((tspan(2)-tspan(1))/h); % number of integration steps

% preallocation ORD, DY_all
ORD=zeros(steps*30,1); % 30 possibilities of halfing integration stepsize
ORD_tmp = zeros(3,1);
DY_all=[];

% stoping rule: maximum number DY terms are smaller then eps
stopping=3;

% last step tolerance
ls_tol = 10^-10;

scaleIndex = 1;
mult2 = true;
mult5 = false;

%% precalculation_DY_mult (we do it in simulate.m script)
% --> precalculation_DY_mult(maxORD)
% ind=load('DY_indexes_maxORD_25.mat','DY_ij','DY_ijk','DY_ijkl','DY_ijklm');
% ind=load('DY_indexes_maxORD_65.mat','DY_ij','DY_ijk','DY_ijkl','DY_ijklm');

%% indexes on multiplications terms
if isempty(ij)
    A2y=zeros(ne,1);
else
    i2=ij(:,1);
    j2=ij(:,2);
end

if isempty(ijk)
    A3y=zeros(ne,1);
else
    i3=ijk(:,1);
    j3=ijk(:,2);
    k3=ijk(:,3);
    mult2 = false;
end

if isempty(ijkl)
    A4y=zeros(ne,1);
else
    i4=ijkl(:,1);
    j4=ijkl(:,2);
    k4=ijkl(:,3);
    l4=ijkl(:,4);
    mult2 = false;
end

if isempty(ijklm)
    A5y=zeros(ne,1);
else
    i5=ijklm(:,1);
    j5=ijklm(:,2);
    k5=ijklm(:,3);
    l5=ijklm(:,4);
    m5=ijklm(:,5);
    mult2 = false;
    mult5 = true;
end

%% MTSM calculation
i=1; % timestep index
t(i)=tspan(1,1);
y(:,1)=init'; % initial values
i=i+1;

% while t(i-1)+ls_tol<tspan(1,2)
  %while t(i-1) <= tspan(1,2)
 while true
  %while t(i-1)+ls_tol<=tspan(1,2) 
%     scaleCount = 0;
    DY=zeros(ne,maxORD+1);
    
    y(:,i)=y(:,i-1); % first term of Taylor series y_{i+1}=y_{i}+...
    DY(:,1)=y(:,i);
    
    K1 = 1/DY(2,1);
    K2 = DY(1,1);
    
    % linear term
    Ay=A*DY(:,1)+b;

    ij1 = 2:-1:1;
    ij2 = 1:2;
    
    ijk2 = ij1;
    ijk3 = ij2;

    ijkl2 = ij1;
    ijkl3 = ij2;

    ijklm2 = [2 1 1];
    ijklm3 = [1 2 1];
    ijklm4 = [1 1 2];


    if ~isempty(ij)
        A2y=A2*(DY(ij(:,1),1).*DY(ij(:,2),1));
    end
    
    if ~isempty(ijk)
        % already calculated submatrices
        DY_3terms = zeros(size(ijk,1), maxORD+stopping);
        
        % computation of DY1 (indeces 1 1 1 1 / 1 1 1 1 1) 
        DY_3terms(:,1) = DY(j3,1).*DY(k3,1);
        
        A3y=A3*(DY(i3,1).*DY_3terms(:,1));
    end
    
    if ~isempty(ijkl)
        % already calculated submatrices
        DY_4terms = zeros(size(ijkl,1), maxORD+stopping);
        
        % calculated term
        DY_4terms_smaller = zeros(size(ijkl,1), maxORD+stopping);
        
        % new array (1 1 / 1 1 1)
        DY_4terms_smaller(:,1) = DY(ijkl(:,2),1).*DY(ijkl(:,3),1);
        
        % new array (1 1 / 1 1 1)
        DY_4terms(:,1) = DY_4terms_smaller(:,1).*DY(ijkl(:,4),1);
        
        A4y=A4*(DY(ijkl(:,1),1).*DY_4terms(:,1));
    end
    
    if ~isempty(ijklm) 
        
        DY_5terms = zeros(size(ijklm,1), maxORD+stopping);

        % TODO: arrays for already calculated submatrices of the currently
        DY_5terms_smaller = zeros(size(ijklm,1), maxORD+stopping);
        
        % TODO: for DY4 and DY5, save two / three multiplicications into the
        % new array (1 1 / 1 1 1)
        DY_5terms_smaller(:,1) = DY(ijklm(:,2),1).*DY(ijklm(:,3),1).*DY(ijklm(:,4),1);

        DY_5terms(:,1) =  DY_5terms_smaller(:,1).*DY(ijklm(:,5),1);
        
        A5y=A5*(DY(ijklm(:,1),1).*DY_5terms(:,1));
    end
    
    D = zeros(2,1);
    D(1) = DY(1,1);    
    
    DY(:,2)=h*(Ay+A2y+A3y+A4y+A5y+D); % first derivative  
    y(:,i)=y(:,i)+DY(:,2); % first term (first derivative)

    maxDY = ones(1,stopping)*10^10;

    k=2;
    mi = 0;
    while norm(maxDY)>eps
        % linear term
        Ay=A*DY(:,k); % opraveno Vasek 23.10.2017
        % -------------------------------------------
        if ~isempty(ij)
            % 2 terms
%           A2y=A2*sum(DY(i2,ind.DY_ij{1,k-1}).*DY(j2,ind.DY_ij{2,k-1}),2);
%             ij1 = ind.DY_ij{1,k-1};
%             ij2 = ind.DY_ij{2,k-1};
            A2y=A2*sum(DY(i2,ij1).*DY(j2,ij2),2);
        end
        
        
        % -------------------------------------------
        % unique 1st row of TS terms indeces (equivalent with unique(DY_ijk..{1,k-1})
         
        if ~mult2
            k_array = k:-1:2;
        end
        % 3 terms - original solution
        %A3y=A3*sum(DY(i3,ind.DY_ijk{1,k-1}).*DY(j3,ind.DY_ijk{2,k-1}).*DY(k3,ind.DY_ijk{3,k-1}),2);
        
        if ~isempty(ijk)
            % 3 terms - PVGN solution
            % already calculated terms - we use the submatrices and multiply them by the unique 1st DY_ijk row  
            A3y_calculated = sum(DY(i3,k_array).*DY_3terms(:,1:k-1),2);
            
            % new terms - new submatrix multiplied only by the corresponding
            % rows from DY, we do not consider DY_ijk indeces, because there
            % are always ones 
%             DY_3terms(:,k) = sum(DY(j3,ind.DY_ijk{2,k-1}(t3_new_start:end)).*DY(k3,ind.DY_ijk{3,k-1}(t3_new_start:end)),2);
%             DY_3terms(:,k) = sum(DY(j3,ijk2(t3_new_start:end)).*DY(k3,ijk3(t3_new_start:end)),2);
            DY_3terms(:,k) = sum(DY(j3,ijk2).*DY(k3,ijk3),2);

            % new  calculated terms have to be multiplied with the 1st row 
            A3y_new = DY(i3,1).*DY_3terms(:,k);
            A3y=A3*sum([A3y_calculated,A3y_new],2);
            % new column index - number of columns of the previous one
%             t3_new_start = size(ijk2,2)+1;
        end

        if ~isempty(ijkl)
            % 4 terms  - PVGN solution
            %A4y=A4*sum(DY(i4,ind.DY_ijkl{1,k-1}).*DY(j4,ind.DY_ijkl{2,k-1}).*DY(k4,ind.DY_ijkl{3,k-1}).*DY(l4,ind.DY_ijkl{4,k-1}),2);
            % already calculated terms - we use the submatrices and multiply them by the 1st DY_ijkl row  
            % DY4: 2 1 1 1
            A4y_calculated = sum(DY(i4,k_array).*DY_4terms(:,1:k-1),2);

            % TODO: calculate the multiplications for already calculated
            % submatrices
            % DY4: 1 1 2
            % DY_4terms_smaller = 1 1 (k-1)
            % A4y_calculated_smaller = 1 1 2
            A4y_calculated_smaller = sum(DY(l4,k_array).*DY_4terms_smaller(:,1:k-1),2);

            % TODO: calculate new columns, save result
            % DY_4terms_smaller(:,k) = 2 1 + 1 2 (green) (ind. k)
%             DY_4terms_smaller(:,k) = sum(DY(j4,ind.DY_ijkl{2,k-1}(t4_new_start_smaller:end)).*DY(k4,ind.DY_ijkl{3,k-1}(t4_new_start_smaller:end)),2);    
            DY_4terms_smaller(:,k) = sum(DY(j4,ijkl2).*DY(k4,ijkl3),2);   

            %DY_4tmp = 2 1 1 + 1 2 1
            DY_4terms_new_part = DY_4terms_smaller(:,k).*DY(l4,1);
            % new terms - new submatrix multiplied only by the corresponding
            % rows from DY, we do not consider DY_ijkl indeces, because there
            % are always ones 
            DY_4terms(:,k) = sum([A4y_calculated_smaller,DY_4terms_new_part],2);

            % new part, multiply by 1 in the 1st row (1112 + 1211 + 1121)
            A4y_new = DY(i4,1).*DY_4terms(:,k);
            A4y=A4*sum([A4y_calculated,A4y_new],2);
        end
              
        if ~isempty(ijklm)
            % 5 terms  - PVGN solution
            %A5y=A5*sum(DY(i5,ind.DY_ijklm{1,k-1}).*DY(j5,ind.DY_ijklm{2,k-1}).*DY(k5,ind.DY_ijklm{3,k-1}).*DY(l5,ind.DY_ijklm{4,k-1}).*DY(m5,ind.DY_ijklm{5,k-1}),2);
            % already calculated terms - we use the submatrices and multiply them by the 1st DY_ijklm row  
            A5y_calculated = sum(DY(i5,k_array).*DY_5terms(:,1:k-1),2);
            % new terms - new submatrix multiplied only by the corresponding
            % rows from DY, we do not consider DY_ijk indeces, because there
            % are always ones 
            A5y_calculated_smaller = sum(DY(m5,k_array).*DY_5terms_smaller(:,1:k-1),2);

%             DY_5terms_smaller(:,k) = sum(DY(j5,ind.DY_ijklm{2,k-1}(t5_new_start_smaller:end)).*DY(k5,ind.DY_ijklm{3,k-1}(t5_new_start_smaller:end)).*DY(l5,ind.DY_ijklm{4,k-1}(t5_new_start_smaller:end)),2); 
            DY_5terms_smaller(:,k) = sum(DY(j5,ijklm2).*DY(k5,ijklm3).*DY(l5,ijklm4),2);
   

            DY_5terms_new_part = DY_5terms_smaller(:,k).*DY(m5,1);

            DY_5terms(:,k) = sum([A5y_calculated_smaller,DY_5terms_new_part],2);
            A5y_new = DY(i5,1).*DY_5terms(:,k);
            A5y=A5*sum([A5y_calculated,A5y_new],2);
        end
        y_index = k-1:-1:1;
        g_index = 2:1:k; 
        D(1) = K1*(DY(1,k) - sum(DY(1,y_index).*DY(2,g_index)));
        
        DY(:,k+1)=(h/k)*(Ay+A2y+A3y+A4y+A5y+D);
        y(:,i)=y(:,i)+DY(:,k+1);
        
        mi = mi+1;
        maxDY(mi)=max(abs(DY(:,k+1)));
        k=k+1;
        
        if mi == stopping
            mi = 0;
        end

        if k > maxORD % max ORD was reached
            break;
        end
        
        ij1 = k:-1:1; 
        ij2 = 1:k;
        

        if ~mult2
            ijk2 = ij1;
            ijk3 = ij2;

            ijkl2 = ij1;
            ijkl3 = ij2;
            if mult5
%                 neg = k:-1:1;
%                 pos = 1:k;
                ijklm2 = repelem(ij1,ij2);
                ijklm3 = [ijklm3,ij1];
                ijklm4 = [ijklm4,ij2];
            end
        end
    end

    if k > maxORD % max ORD was reached
        h=h/2; % we use half-size integration step
        fprintf('GN: Halving integration step size.\n')
    else
%         DY_all{i-1}=DY;
%         if k<(stopping+4)
%             h=h*2;
%             Hk=(h.^iPascal_mixed)./fact_k;
%             fprintf('Increasing integration step size.\n')
%         end
        %DY(:,1:k+1)
        ORD(i)=k;
        ORD_tmp(scaleIndex) = k;
        scaleIndex = scaleIndex + 1;
        if scaleIndex == 4
            scaleIndex = 1;
        end
        %scaleCount = scaleCount+k;
        %if scaleCount > 3*minORD  
        %    scaleCount = 0;
        %end
        t(i)=t(i-1)+h;
        
        if i > stopping 
            if hScaleFactor > 1
        %       if sum(ORD(i-stopping+1:i) < minORD) == 3
                scaleCount = ORD_tmp(1)+ORD_tmp(2)+ORD_tmp(3); 
                if scaleCount == 3*minORD
                    h = h*hScaleFactor;
    %                 scaleCount = 0;
                end
            end
        end
        
        if t(i)+h > tspan(1,2)
            h = tspan(1,2) - t(i);
        end
        
        if t(i)+ls_tol<tspan(1,2)
            i=i+1;
        else
            break
        end
    end
end

% shrink the rest (empty) array
% ORD(i:end)=[];
ORD(i+1:end)=[];


