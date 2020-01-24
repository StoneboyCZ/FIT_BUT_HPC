function [t,y,ORD,DY_all]=explicitTaylorMult(h,tspan,init,A,A2,A3,A4,A5,b,ij,ijk,ijkl,ijklm,eps,ind,maxORD)
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


%%% indexes on multiplications DY terms
% DY_ij=cell(2,maxORD+stopping);
% DY_ijk=cell(3,maxORD+stopping);
% DY_ijkl=cell(4,maxORD+stopping);
% DY_ijklm=cell(5,maxORD+stopping);
% 
% for k=2:maxORD+stopping-1
%     
%     mm1=k:-1:1;
%     mm2=1:k;
%     DY_ij{1,k-1}=mm1;
%     DY_ij{2,k-1}=mm2;
%     
%     for jj=1:k
%         mm1=k+1-jj:-1:1;
%         mm2=(1:k-jj+1);
%         DY_ijk{1,k-1}=[DY_ijk{1,k-1},mm1];
%         DY_ijk{2,k-1}=[DY_ijk{2,k-1},mm2];
%         DY_ijk{3,k-1}=[DY_ijk{3,k-1},jj*ones(1,length(mm1))];
%     end
%     
%     for ii=1:k
%         for jj=1:k-ii+1
%             mm1=k+2-jj-ii:-1:1;
%             mm2=(1:k-jj-ii+2);
%             DY_ijkl{1,k-1}=[DY_ijkl{1,k-1},mm1];
%             DY_ijkl{2,k-1}=[DY_ijkl{2,k-1},mm2];
%             DY_ijkl{3,k-1}=[DY_ijkl{3,k-1},jj*ones(1,length(mm1))];
%             DY_ijkl{4,k-1}=[DY_ijkl{4,k-1},ii*ones(1,length(mm1))];
%         end
%     end
%     
%     for iii=1:k
%         for ii=1:k-iii+1
%             for jj=1:k-ii-iii+2
%                 mm1=k+3-jj-ii-iii:-1:1;
%                 mm2=(1:k-jj-ii-iii+3);
%                 DY_ijklm{1,k-1}=[DY_ijklm{1,k-1},mm1];
%                 DY_ijklm{2,k-1}=[DY_ijklm{2,k-1},mm2];
%                 DY_ijklm{3,k-1}=[DY_ijklm{3,k-1},jj*ones(1,length(mm1))];
%                 DY_ijklm{4,k-1}=[DY_ijklm{4,k-1},ii*ones(1,length(mm1))];
%                 DY_ijklm{5,k-1}=[DY_ijklm{5,k-1},iii*ones(1,length(mm1))];
%             end
%         end
%     end
% end % k (ORD)
% 
% filename=['DY_indexes_maxORD_',int2str(maxORD)];
% save(filename,'DY_ij','DY_ijk','DY_ijkl','DY_ijklm');

% ind=load('DY_indexes_maxORD_25.mat','DY_ij','DY_ijk','DY_ijkl','DY_ijklm');
% ind=load('DY_indexes_maxORD_65.mat','DY_ij','DY_ijk','DY_ijkl','DY_ijklm');

%%% indexes on multiplications terms

i2=ij(:,1);
j2=ij(:,2);

i3=ijk(:,1);
j3=ijk(:,2);
k3=ijk(:,3);

if isempty(ijkl)
    A4y=zeros(ne,1);
else
    i4=ijkl(:,1);
    j4=ijkl(:,2);
    k4=ijkl(:,3);
    l4=ijkl(:,4);
end

if isempty(ijklm)
   A5y=zeros(ne,1);
else 
    i5=ijklm(:,1);
    j5=ijklm(:,2);
    k5=ijklm(:,3);
    l5=ijklm(:,4);
    m5=ijklm(:,5);
end


i=1; % timestep index
t(i)=tspan(1,1);
y(:,1)=init'; % initial values

i=i+1;

while t(i-1)+ls_tol<tspan(1,2)
    DY=zeros(ne,maxORD+1);
    y(:,i)=y(:,i-1); % first term of Taylor series y_{i+1}=y_{i}+...
    DY(:,1)=y(:,i);
%     for n=1:ne
%        DY(1,n)=h*(sum(A(n,:).*y(i,:))+b(n));
%     end
    
    % linear term
    Ay=A*DY(:,1)+b;
    

    % multiplications terms
    A2y=A2*(DY(ij(:,1),1).*DY(ij(:,2),1));
    A3y=A3*(DY(ijk(:,1),1).*(DY(ijk(:,2),1).*DY(ijk(:,3),1)));
    if ~isempty(ijkl)
    A4y=A4*(DY(ijkl(:,1),1).*DY(ijkl(:,2),1).*DY(ijkl(:,3),1).*DY(ijkl(:,4),1));
    end
    if ~isempty(ijklm)
    A5y=A5*(DY(ijklm(:,1),1).*DY(ijklm(:,2),1).*DY(ijklm(:,3),1).*DY(ijklm(:,4),1).*DY(ijklm(:,5),1));
    end
    
    DY(:,2)=h*(Ay+A2y+A3y+A4y+A5y); % first derivative
    
    y(:,i)=y(:,i)+DY(:,2); % first term (first derivative)
    
%     k=1;
%     maxDY = ones(1,stopping)*10^10;
    maxDY = [ones(1,stopping)*10^10 zeros(1,maxORD+stopping)];

    k=2;
    while norm(maxDY(k-1:stopping+k-2))>eps
        
        % linear term
        Ay=A*DY(:,k); % opraveno Vasek 23.10.2017
        
        A2y=A2*sum(DY(i2,ind.DY_ij{1,k-1}).*DY(j2,ind.DY_ij{2,k-1}),2);
        A3y=A3*sum(DY(i3,ind.DY_ijk{1,k-1}).*DY(j3,ind.DY_ijk{2,k-1}).*DY(k3,ind.DY_ijk{3,k-1}),2);
        if ~isempty(ijkl)
            A4y=A4*sum(DY(i4,ind.DY_ijkl{1,k-1}).*DY(j4,ind.DY_ijkl{2,k-1}).*DY(k4,ind.DY_ijkl{3,k-1}).*DY(l4,ind.DY_ijkl{4,k-1}),2);
        end
        if ~isempty(ijklm)
            A5y=A5*sum(DY(i5,ind.DY_ijklm{1,k-1}).*DY(j5,ind.DY_ijklm{2,k-1}).*DY(k5,ind.DY_ijklm{3,k-1}).*DY(l5,ind.DY_ijklm{4,k-1}).*DY(m5,ind.DY_ijklm{5,k-1}),2);
        end
        DY(:,k+1)=(h/k)*(Ay+A2y+A3y+A4y+A5y);
        
        y(:,i)=y(:,i)+DY(:,k+1);
        
        maxDY(k+stopping-1)=max(abs(DY(:,k+1)));
        k=k+1;
        
        if k > maxORD % max ORD was reached
            break;
        end
        
    end
    
    if k > maxORD % max ORD was reached
        h=h/2; % we use half-size integration step
        fprintf('VS: Halving integration step size.\n')
    else
%         DY_all{i-1}=DY;
%         if k<(stopping+4)
%             h=h*2;
%             Hk=(h.^iPascal_mixed)./fact_k;
%             fprintf('Increasing integration step size.\n')
%         end
        ORD(i)=k;
        t(i)=t(i-1)+h;
        i=i+1;
    end
end

% shrink the rest (empty) array
ORD(i:end)=[];


