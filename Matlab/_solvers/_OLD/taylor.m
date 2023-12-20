function [t,y,ORD,DY_all]=taylor( ... 
    h, ... size of the integration step
    tspan, ... times of calculation
    init, ... initial conditions for the system
    A, ... Jacobian for the linear part of the system
    b, ...  vector of constants (right-hand side)
    d, ... terms in the divisions
    m, ... terms in the multiplications
    index_l, ... indexes of linear part,
    index_d, ... indexes for division,
    index_m, ... indexes for multiplication,
    eps, ...  accuracy
    maxORD, ... maximum order of the Taylor series
    minORD, ... minimum ord needed for step increase 
    hScaleFactor ... how h is scaled (hnew = h*hscalefactor) when ORD is stable
    )

ne_ode = size(A, 1); % number of ODEs
ne_div = size(d, 1); % number of divisions
ne_mult = size(m,1); % number of multiplications

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

if ~isempty(d)
    d1=d(:,1);
    d2=d(:,2);
end

if ~isempty(m)
    m1 = m(:,1);
    m2 = m(:,2);
end

%% MTSM calculation
i=1; % timestep index
t(i)=tspan(1,1);
y(:,1)=init; % initial values

i=i+1;


while true
    DY=zeros(ne_ode+ne_div+ne_mult,maxORD+1);
    DY_lin=zeros(ne_ode,maxORD+1);
    DY_div=zeros(ne_div,maxORD+1);
    DY_mult=zeros(ne_mult,maxORD+1);

    y(:,i)=y(:,i-1); % first term of Taylor series y_{i+1}=y_{i}+...
    y(index_d,i) = y(d1,i)./y(d2,i);
    y(index_m,i) = y(m1,i).*y(m2,i);

    DY_lin(:,1) = y(index_l,i);
    DY_div(:,1) = y(index_d,i);
    DY_mult(:,1) = y(index_m,i);
    DY(:,1)=[DY_lin(:,1);DY_div(:,1);DY_mult(:,1)];

    % division constant
    
    B1 = eye(ne_div)*(1./DY(d2,1));


    %% first  derivative
    Ay=A*DY(:,1)+b;
    DY_lin(:,2)=h*(Ay); % first derivative
    DY_div(:,2)=B1*(DY_lin(d1,2)-DY_div(:,1).*DY_lin(d2,2));
    DY_mult(:,2)=DY_lin(m1,2).*DY_lin(m2,1) + DY_lin(m1,1).*DY_lin(m2,2);
    DY(:,2) = [DY_lin(:,2);DY_div(:,2);DY_mult(:,2)];

    y(:,i)=y(:,i)+DY(:,2); % first term (first derivative)

    maxDY = ones(1,stopping)*10^10;

    %%% second and additional derivatives
    k=2;
    mi = 0;
    while norm(maxDY)>eps
        i1 = k:-1:1;
        i2 = 2:1:k+1;
        j1 = k+1:-1:1;
        j2 = 1:k+1;
        Ay=A*DY(:,k); % opraveno Vasek 23.10.2017
        
        DY_lin(:,k+1)= (h/k)*(Ay);
        DY_div(:,k+1)= B1*(DY_lin(d1,k+1) - sum(DY_div(:,i1).*DY_lin(d2,i2)));
        DY_mult(:,k+1) = sum(DY_lin(m1,j1).*DY_lin(m2,j2));
        DY(:,k+1) = [DY_lin(:,k+1);DY_div(:,k+1);DY_mult(:,k+1)];
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
        
%         ij1 = k:-1:1; 
%         ij2 = 1:k;
%         
    end

    if k > maxORD % max ORD was reached
        h=h/2; % we use half-size integration step
        fprintf('GN: Halving integration step size.\n')
    else
        ORD(i)=k;
        ORD_tmp(scaleIndex) = k;
        scaleIndex = scaleIndex + 1;
        if scaleIndex == 4
            scaleIndex = 1;
        end
        t(i)=t(i-1)+h;
        
        if i > stopping 
            if hScaleFactor > 1
                scaleCount = ORD_tmp(1)+ORD_tmp(2)+ORD_tmp(3); 
                if scaleCount == 3*minORD
                    h = h*hScaleFactor;
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
ORD(i+1:end)=[];


