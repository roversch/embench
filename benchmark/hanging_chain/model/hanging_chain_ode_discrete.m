function x_next = hanging_chain_ode_discrete(x, u, Ts, num_masses)

    k1 = Ts * hanging_chain_ode(x, u, num_masses);
    k2 = Ts * hanging_chain_ode(x + k1./2, u, num_masses);
    k3 = Ts * hanging_chain_ode(x + k2./2, u, num_masses);
    k4 = Ts * hanging_chain_ode(x + k3, u, num_masses);
    
    x_next = x + 1/6 * (k1 + 2*k2 + 2*k3 + k4);
    
end
    