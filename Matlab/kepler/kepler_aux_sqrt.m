function dy = kepler_aux_sqrt(~,y)
    dy = zeros(6,1);
    dy(1) = y(3);
    dy(2) = y(4);
    dy(3) = -y(1)/y(5); 
    dy(4) = -y(2)/y(5);
    dy(5) = 3*y(6)*(y(1)*y(3)+y(2)*y(4));
    dy(6) = (y(1)*y(3)+y(2)*y(4))/y(6);
end