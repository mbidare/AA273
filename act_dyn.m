% Actual Dynamics of state of system (non-linear)
% State X = [x;xdot,y,ydot,z,zdot]
% input U = [V (const vel);phi (rotation rate); z vel rate]
function x_next = act_dyn(x,u,t,dt)
    mass = 1;
    theta = u(2);
    x_next = zeros(6,1);
    x_next(1) = x(1) + dt*(u(1)*cos(theta));
   % x_next(2) = x(2) + dt*(u(1)/mass);
    x_next(3) = x(3) + dt*(u(1)*sin(theta));
   % x_next(4) = x(4) + dt*(u(2)/mass);
    x_next(5) = x(5) + dt*x(6);
    x_next(6) = x(6) + dt*u(3)/mass; 
end