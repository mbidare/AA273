% Linearized Dynamics of state of system

function x_next = lin_dyn(x,u,dt, A, B)
    mass = 1;
    x_next = x + dt*(A*x + B*u);
end