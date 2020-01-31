function [t,y,ORD,DY_all]=explicitTaylorMult_ver2(h,tspan,init,A,A2,A3,A4,A5,b,ij,ijk,ijkl,ijklm,eps)
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
% AT=A';
steps = round((tspan(2)-tspan(1))/h); % number of integration steps
% t=tspan(1):h:steps*h; % prealocation only for fixed step
%y=zeros([steps+1,ne]);
%DY=zeros(ORD,ne); % Taylor series terms - for linear system of n ODEs should be size [ORDxne]

% preallocation ORD, DY_all
ORD=zeros(steps*30,1); % 30 possibilities of halfing integration stepsize
DY_all=[];


% stoping rule: maximum number DY terms are smaller then eps
stopping=3;
% max ORD for double precision
maxORD=60;
% last step tolerance
ls_tol = 10^-10;


% indexes for mixed quadratic term B3*yi.*yj.*yk
% ijk=nmultichoosek(1:ne,3);
% ii=ijk(:,1);
% jj=ijk(:,2);
% kk=ijk(:,3);

%%% indexes on multiplications terms
% i2=ij(:,1);
% j2=ij(:,2);
% 
% i3=ijk(:,1);
% j3=ijk(:,2);
% k3=ijk(:,3);


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
    A4y=A4*(DY(ijkl(:,1),1).*DY(ijkl(:,2),1).*DY(ijkl(:,3),1).*DY(ijkl(:,4),1));
    A5y=A5*(DY(ijklm(:,1),1).*DY(ijklm(:,2),1).*DY(ijklm(:,3),1).*DY(ijklm(:,4),1).*DY(ijklm(:,5),1));
    
    DY(:,2)=h*(Ay+A2y+A3y+A4y+A5y); % first derivative
    
    y(:,i)=y(:,i)+DY(:,2); % first term (first derivative)
    
%     k=1;
%     maxDY = ones(1,stopping)*10^10;
    maxDY = [ones(1,stopping)*10^10 zeros(1,maxORD+stopping)];

    k=2;
    while norm(maxDY(k-1:stopping+k-2))>eps
        
        % linear term
        Ay=A*DY(:,k); % opraveno Vasek 23.10.2017
        

        ijkl_length=length(ijkl(:,1));
        ijklm_length=length(ijklm(:,1));
        A4y=zeros(ijkl_length,1);
        A5y=zeros(ijklm_length,1);
        for iii=1:k
            A5ysum_iii=zeros(ijklm_length,1);
            for ii=1:k-iii+1
                if iii==1
                    A4ysum_ii=zeros(ijkl_length,1);
                end
                A5ysum_ii=zeros(ijklm_length,1);
                for jj=1:k-ii-iii+2
                    if  (iii==1) && (ii==1) && (jj==1)   % iii==ii==jj==1 % iii=2, ii= 2, jj=1 -> TRUE
                        mm1=k:-1:1;
                        mm2=1:k;
                        A2y=sum(DY(ij(:,1),mm1).*DY(ij(:,2),mm2),2);
                        A3y=DY(ijk(:,1),jj).*sum(DY(ijk(:,2),mm1).*DY(ijk(:,3),mm2),2);
                        A4ysum_ii=A4ysum_ii+DY(ijkl(:,1),jj).*sum(DY(ijkl(:,2),mm1).*DY(ijkl(:,3),mm2),2);
                    elseif (iii==1) && (ii==1)
                        mm1=k+1-jj:-1:1;
                        mm2=(1:k-jj+1);
                        A3y=A3y+DY(ijk(:,1),jj).*sum(DY(ijk(:,2),mm1).*DY(ijk(:,3),mm2),2);
                        A4ysum_ii=A4ysum_ii+DY(ijkl(:,1),jj).*sum(DY(ijkl(:,2),mm1).*DY(ijkl(:,3),mm2),2);
                    elseif iii==1
                        mm1=k+2-jj-ii:-1:1;
                        mm2=(1:k-jj-ii+2);
                        A4ysum_ii=A4ysum_ii+DY(ijkl(:,1),jj).*sum(DY(ijkl(:,2),mm1).*DY(ijkl(:,3),mm2),2);
                    else
                        mm1=k+3-jj-ii-iii:-1:1;
                        mm2=(1:k-jj-ii-iii+3);
                    end % end if
                    A5ysum_ii=A5ysum_ii+DY(ijklm(:,1),jj).*sum(DY(ijklm(:,2),mm1).*DY(ijklm(:,3),mm2),2);
                end % end for jj
                if iii==1
                    A4y=A4y+A4ysum_ii.*DY(ijkl(:,4),ii);
                end
                A5ysum_iii=A5ysum_iii+A5ysum_ii.*DY(ijklm(:,4),ii);
            end %end for ii
            A5y=A5y+A5ysum_iii.*DY(ijklm(:,5),iii);
        end % end for iii
        A2y=A2*A2y;
        A3y=A3*A3y;
        A4y=A4*A4y;
        A5y=A5*A5y;
        
        
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
        fprintf('Halving integration step size.\n')
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


