function [t,y,ORD,DY_all]=taylor( ... 
    h, ... size of the integration step
    tspan, ... times of calculation
    init, ... initial condition
    A, ... Jacobian for the linear part of the system
    b, ...  vector of constants (right-hand side)
    eps, ...  accuracy
    maxORD, ... maximum order of the Taylor series
    minORD, ... minimum ord needed for step increase 
    hScaleFactor ... how h is scaled (hnew = h*hscalefactor) when ORD is stable
    )
%  [t,y,ORD]=taylor(h,tspan,init,A,B3,C,d,eps)
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
ORD_tmp = zeros(3,1);
DY_all=[];

% stoping rule: maximum number DY terms are smaller then eps
stopping=3;

% last step tolerance
ls_tol = 10^-10;

scaleIndex = 1;

%% MTSM calculation
i=1; % timestep index
t(i)=tspan(1,1);
y(:,1)=init'; % initial values

i=i+1;

while true
    DY=zeros(ne,maxORD+1);
    
    y(:,i)=y(:,i-1); % first term of Taylor series y_{i+1}=y_{i}+...
    DY(:,1)=y(:,i);
    
    % linear term
    Ay=A*DY(:,1)+b;


%     if ~isempty(ij)
%         A2y=A2*(DY(ij(:,1),1).*DY(ij(:,2),1));
%     end
    
    
    DY(:,2)=h*(Ay+A2y+A3y+A4y+A5y); % first derivative  
    y(:,i)=y(:,i)+DY(:,2); % first term (first derivative)

    maxDY = ones(1,stopping)*10^10;

    k=2;
    mi = 0;
    while norm(maxDY)>eps
        % linear term
        Ay=A*DY(:,k); % opraveno Vasek 23.10.2017
        % -------------------------------------------
%         if ~isempty(ij)
%             A2y=A2*sum(DY(i2,ij1).*DY(j2,ij2),2);
%         end

        DY(:,k+1)=(h/k)*(Ay+A2y+A3y+A4y+A5y);

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