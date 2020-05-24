% Linearized Dynamics of state of system

function x_next = lin_dyn(x,u,dt, A, B)
    x_next = x + dt*(A*x + B*u);
end