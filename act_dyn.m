% Actual Dynamics of state of system (can modify to be non-linear)

function x_next = act_dyn(x,u,dt, A, B)
    x_next = x + dt*(A*x + B*u);
end