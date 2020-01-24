function dy = pendulum_aux(~,y,data)
    dy = zeros(4,1);
    dy(1) = -data.k*y(1)-data.a*y(3); % velocity
    dy(2) = y(1); % position
    dy(3) = y(4)*y(1);
    dy(4) = -y(3)*y(1);
end