function dy = kepler_aux_div_brackets(~,y)
    dy = zeros(12,1);
    
    dy(1) = y(3);
    dy(2) = y(4);
    dy(3) = -y(1)*y(7); 
    dy(4) = -y(2)*y(7); 
    dy(5) = 3*y(6)*(y(9)+y(10));
    dy(6) = y(8)*(y(9)+y(10));
    dy(7) = -3*y(6)*y(7)*y(7)*(y(9)+y(10));
    dy(8) = -y(8)*y(8)*y(8)*(y(9)+y(10));
    dy(9) = y(3)*y(3) - y(7)*y(11);
    dy(10) = y(4)*y(4) -y(7)*y(12);
    dy(11) = 2*y(9);
    dy(12) = 2*y(10);
end