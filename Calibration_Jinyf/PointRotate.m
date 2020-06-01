function [output] = PointRotate(c0,delta,c1,theta)
output = zeros(3,1);
x = c0(1)-delta(1);
y = c0(2)-delta(2);
z = c0(3)-delta(3);
vx = c1(1);
vy = c1(2);
vz = c1(3);
c = cos(theta);
s = sin(theta);
output(1) = (vx*vx*(1-c)+c)*x + (vx*vy*(1-c)-vz*s)*y + (vx*vz*(1-c)+vy*s)*z;
output(2) = (vy*vx*(1-c)+vz*s)*x + (vy*vy*(1-c)+c)*y + (vy*vz*(1-c)-vx*s)*z;
output(3) = (vz*vx*(1-c)-vy*s)*x + (vz*vy*(1-c)+vx*s)*y + (vz*vz*(1-c)+c)*z;
end

