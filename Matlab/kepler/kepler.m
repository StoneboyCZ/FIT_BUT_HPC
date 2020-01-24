function dy = kepler(~,y)
    r = sqrt(y(1)*y(1) + y(2)*y(2));
    r3 = r*r*r;

    dy = zeros(4,1);
    dy(1) = y(3);
    dy(2) = y(4);
    dy(3) = -(y(1)/r3); 
    dy(4) = -(y(2)/r3);
end