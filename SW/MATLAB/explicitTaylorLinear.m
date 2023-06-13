function [t,y,ORD,DY_all]=explicitTaylorLinear(h,tspan,init,A,b,eps)
%  [t,y,ORD]=explicitTaylor(h,tspan,init,A,b,eps)
%  Numerical solution of ODEs with explicit Taylor series method
%  problem y' = A*y + b
%  where A is Jacobian matrix of constants and b is a vector of constants
%  ---------------------------------------------------------------
%  h ... integration step size
%  tmax ... interval of numerical solution [0,tmax]
%  ORD ... order of Taylor series - number of Taylor series Terms
%  DY ... cell of all derivatives in each time step DY [dy; ddy; d3y; ...]

ne=size(A,1); % number of equations
steps = round((tspan(2)-tspan(1))/h); % number of integration steps
t=tspan(1):h:steps*h;
%y=zeros([steps+1,ne]);
%DY=zeros(ORD,ne); % Taylor series terms - for linear system of n ODEs should be size [ORDxne]

DY_all=[];

% stoping rule: maximum number DY terms are smaller then eps
stopping=3;
% max ORD for double precision
maxORD=64;
% last step tolerance
ls_tol = 10^-10;


i=1; % timestep index
t(i)=tspan(1,1);
y(:,1)=init; % initial values

i=i+1;

while t(i-1)+ls_tol<tspan(1,2)
    DY=zeros(ne,maxORD);
    y(:,i)=y(:,i-1); % first term of Taylor series y_{i+1}=y_{i}+...

%     for n=1:ne
%        DY(1,n)=h*(sum(A(n,:).*y(i,:))+b(n));
%     end
    Ay=A*y(:,i);
    DY(:,1)=h*(Ay+b);
    
    y(:,i)=y(:,i)+DY(:,1); % first derivative
    
    k=1;
%     maxDY = ones(1,stopping)*10^10;
    maxDY = [ones(1,stopping)*10^10 zeros(1,maxORD+stopping)];

    k=2;
    while norm(maxDY(k-1:stopping+k-2))>eps
%         for n=1:ne
%             DY(k,n)=(h/k)*(sum(A(n,:).*DY(k-1,:)));
%         end
        ADy=A*DY(:,k-1);
%         ADy=DY(k-1,:)*AT;
        DY(:,k)=(h/k)*ADy;
        
        y(:,i)=y(:,i)+DY(:,k);
        
        maxDY(k+stopping-1)=abs(max(DY(:,k)));
        k=k+1;
        
        if k > maxORD % max ORD was reached
            break;
        end
        
    end
    
    if k > maxORD % max ORD was reached
        h=h/2; % we use half-size integration step
    else
        DY_all{i-1}=DY;
        ORD(i)=k;
        t(i)=t(i-1)+h;
        i=i+1;
    end
end




