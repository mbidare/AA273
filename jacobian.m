function [A,B] = jacobian(x,u,dt)
    A = zeros(6,6); A(1,1) = 1; A(3,3) = 1; A(5,5) = 1; A(5,6) = dt; A(6,6)=1;
    B = zeros(6,3); B(1,1) = dt*cos(u(2)); B(1,2) = -dt*u(1)*sin(u(2));
    B(3,1) = dt*sin(u(2)); B(3,2) = dt*u(1)*cos(u(2)); B(6,3) = dt/mass;
end