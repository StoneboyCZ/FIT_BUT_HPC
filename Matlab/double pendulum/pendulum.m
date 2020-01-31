function dz = pendulum(~,z,data)
%% equations of the system
    % z1dot = phi1dot
    % z2dot = phi2dot
    % z3dot = 
    % (
    % -m2*L1*z3*z3*sin(z1-z2)*cos(z1-z2)+m2*g*sin(z2)*cos(z1-z2)-m2*L2*z4*z4*sin(z1-z2)
    % - (m1+m2)*g*sin(z1) 
    % ) / (L1* (m1+m2)-m2*L1*cos(z1-z2)*cos(z1-z2) 
    % z4dot = 
    % (
    % m2*L2*z4*z4*sin(z1-z2)*cos(z1-z2)+g*(m1+m2)sin(z1)cos(z1-z2)+L1*z3*z3*(m1+m2)*sin(z1-z2)
    % -g*sin(z2)*(m1+m2) / (L2*(m1+m2) - m2*L2*cos(z1-z2)*cos(z1-z2))
    %
    % z(1) -- phi1; z(2) -- phi2; z(3) -- phi1dot; z(4) -- phi2dot
    %
    
       
    m1 = data.m1;
    m2 = data.m2;
    L1 = data.L1;
    L2 = data.L2;
    g = data.g;
    
    dz = zeros(4,1);
    dz(1) = z(3);
    dz(2) = z(4);
    dz(3) = (-m2*L1*z(3)^2*sin(z(1)-z(2))*cos(z(1)-z(2))+m2*g*sin(z(2))*cos(z(1)-z(2)) - m2*L2*z(4)^2*sin(z(1)-z(2)) - (m1+m2)*g*sin(z(1)))/(L1*(m1+m2)-m2*L1*cos(z(1)-z(2))^2);
    dz(4) = (m2*L2*z(4)^2*sin(z(1)-z(2))*cos(z(1)-z(2))+(m1+m2)*g*sin(z(1))*cos(z(1)-z(2))+L1*z(3)^2*(m1+m2)*sin(z(1)-z(2))-g*(m1+m2)*sin(z(2)))/(L2*(m1+m2)-m2*L2*cos(z(1)-z(2))^2);