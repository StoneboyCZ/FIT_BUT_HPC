function dy = kepler_aux_div_full(~,y)
    dy = zeros(16,1);
    dy(1) = y(3);
    dy(2) = y(4);
    dy(3) = -y(1)*y(7); 
    dy(4) = -y(2)*y(7); 
    dy(5) = 3*y(15);
    dy(6) = y(16);
    dy(7) = -3*y(15)*y(7)^2;
    dy(8) = -y(8)^2*y(16);
    dy(9) = y(3)^2 - y(7)*y(11);
    dy(10) = y(4)^2 - y(7)*y(12);
    dy(11) = 2*y(9);
    dy(12) = 2*y(10);
    dy(13) = y(3)^2 + y(4)^2 - y(7)*y(14);
    dy(14) = 2*(y(9) + y(10));
    dy(15) = y(8)*y(13)^2 + y(6)*(y(3)^2 + y(4)^2 - y(7)*y(14));
    dy(16) = -y(8)*(y(16)^2 - y(3)^2 - y(4)^2 + y(7)*y(14));
end