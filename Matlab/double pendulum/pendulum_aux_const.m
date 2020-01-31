function dz = pendulum_aux_const(~,z,data)
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
    % Auxillary system
    % 
    %
    %
    % z(1) -- phi1; z(2) -- phi2; z(3) -- phi1dot; z(4) -- phi2dot
    %

    A = data.A;
    B = data.B;
    C = data.C;
    D = data.D;
    E = data.E;
    F = data.F;
    
    dz = zeros(10,1);
    dz(1) = z(3);
    dz(2) = z(4);
    dz(3) = (-B*(z(3)^2)*z(5)*z(6) + F*z(9)*z(6) - D*(z(4)^2)*z(5) - (E+F)*z(7)) / ((A+B)-B*z(6)^2);
    dz(4) = ( D*(z(4)^2)*z(5)*z(6) + (E+F)*z(7)*z(6) + (A+B)*(z(3)^2)*z(5) - (E+F)*z(9)) / ((C+D) - D*z(6)^2);
    dz(5) = z(6)*(z(3)-z(4)); % sin(z1-z2)
    dz(6) = z(5)*(z(4)-z(3)); % cos(z1-z2)
    dz(7) = z(8)*z(3); % sin(z1)
    dz(8) = -z(7)*z(3); % cos(z1)
    dz(9) = z(10)*z(4); % sin(z2)
    dz(10) = -z(9)*z(4); %cos(z2)
   