function [t,y,ORD,DY_all]=taylor_v6( ... 
    h, ... size of the integration step
    tspan, ... times of calculation
    init, ... initial conditions for the system
    eps, ...  accuracy
    A, ... Jacobian for the linear part of the system
    b, ...  vector of constants (right-hand side)
    m, ... terms in the multiplications
    m_extra, ... terms in the multiplication 2
    d, ... terms in the divisions
    maxORD, ... maximum order of the Taylor series
    minORD, ... minimum ord needed for step increase 
    hScaleFactor ... how h is scaled (hnew = h*hscalefactor) when ORD is stable
)

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

%% number of equations
ne = size(A,2); % number of equations
ne_ode = size(A, 1); % number of ODEs

index_l = 1:ne_ode;

linear = false;
if ne == ne_ode
    linear = true;
end

multiplication = true;
if ~isempty(m)
    ne_mult1 = size(m,1); 
   
    
    index_m = ne_ode+1:ne_ode+ne_mult1;
    
    
    m1 = m(:,1);
    m2 = m(:,2);

    mult1 = init(m1,1).*init(m2,1);
    init = [init;mult1];
else
    multiplication = false;
end

multiplication2 = true;
if ~isempty(m_extra)
    ne_mult2 = size(m_extra,1);

    index_m2 = ne_ode+ne_mult1+1:ne_ode+ne_mult1+ne_mult2;

    m11 = m_extra(:,1);
    m12 = m_extra(:,2);
    
    mult2 = init(m11,1).*init(m12,1);
    init = [init;mult2];
else
    multiplication2 = false;
end

division = true;
if ~isempty(d)
    ne_div = size(d, 1); % number of divisions
    
    index_d = ne_ode + ne_mult1 + ne_mult2 +1:ne_ode+ne_mult1+ne_mult2+ne_div;
    
    d1=d(:,1);
    d2=d(:,2);
    
    div = init(d1,1)./init(d2,1);
    init = [init;div];
else
    division = false;
end

%% MTSM calculation
i=1; % timestep index
t(i)=tspan(1,1);
y(:,1)=init; % initial values

i=i+1;

while true
    DY=zeros(ne,maxORD+1);

    y(:,i)=y(:,i-1); % first term of Taylor series y_{i+1}=y_{i}+...
    
    if i > 2
        if multiplication
            y(index_m,i) = y(m1,i).*y(m2,i);
        end


        if multiplication2
            y(index_m2,i) = y(m11,i).*y(m12,i);
        end

        if division
            y(index_d,i) = y(d1,i)./y(d2,i);
        end
    end
    DY(:,1) = y(:,i);

    % division constant
    
    if division
        B1 = eye(ne_div).*((1./DY(d2,1))');
    end
    
    %% first  derivative
    Ay=A*DY(:,1)+b;
    
    if linear
        DY(:,2)=h*(Ay); % first derivative
    else
        DY(index_l,2)=h*(Ay); % first derivative
    end
    
    if multiplication
        DY(index_m,2)=DY(m1,2).*DY(m2,1) + DY(m1,1).*DY(m2,2);
    end
    
    
    if multiplication2
        DY(index_m2,2)=DY(m11,2).*DY(m12,1) + DY(m11,1).*DY(m12,2);
    end

    if division
        DY(index_d,2)=B1*(DY(d1,2)-DY(index_d,1).*DY(d2,2));
    end
    
    maxDY = ones(1,stopping)*10^10;

    %%% second and additional derivatives
    k=2;
    mi = 0;
    while norm(maxDY)>eps
        if division
            i1 = k:-1:1;
            i2 = 2:1:k+1;
        end
        
        if multiplication
            j1 = k+1:-1:1;
            j2 = 1:k+1;
        end
        Ay=A*DY(:,k);
        
        if linear
            DY(:,k+1) = (h/k)*(Ay);
        else
            DY(index_l,k+1) = (h/k)*(Ay);
        end
        
        if multiplication
            DY(index_m,k+1) = sum(DY(m1,j1).*DY(m2,j2),2);
            
        end
        
        if multiplication2
            DY(index_m2,k+1) = sum(DY(m11,j1).*DY(m12,j2),2);
        end

        if division
             DY(index_d,k+1) = B1*(DY(d1,k+1) - sum(DY(index_d,i1).*DY(d2,i2),2));
        end
        
        
        
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
        
        y(:,i) = sum(DY,2);

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


