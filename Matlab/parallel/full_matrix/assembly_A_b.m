A=sparse(nRLC*2+2,nRLC*2+2); % matrix A  
b=zeros(nRLC*2+2,1);      

% A = ( A11  A12 )
%     ( A21  A22 )
A11=sparse(nRLC,nRLC);
A12=sparse(nRLC,nRLC+2);
A21=sparse(nRLC+2,nRLC);
A22=sparse(nRLC+2,nRLC+2);

%% load coefficients
i=1; % row
for i=1:nRLC
    if (i~=1) && (i~=nRLC)
        A12(i,i:i+1)=[1/C,-1/C];
        A21(i,i-1:i)=[1/L,-1/L];
    elseif i==1 % first row
        A12(i,i:i+1)=[1/C,-1/C];
        A21(i,i)=-1/L;
        A22(i,i)=-R1/L;
    elseif i==nRLC % last row
        A11(i,i)=-1/(R2*C); 
        A12(i,i)=1/C;
        A21(i,i-1:i)=[1/L,-1/L];
    end
end

%% load the function for u0
A22(1,nRLC+1) = 1/L;

A22(nRLC+1,nRLC+2) = om;
A22(nRLC+2,nRLC+1) = -om;

A=[A11, A12;
  A21, A22];

%full(A)
%spy(A)
