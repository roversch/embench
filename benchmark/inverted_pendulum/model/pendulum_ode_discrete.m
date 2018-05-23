function x_next = pendulum_ode_discrete(x, u, Ts)

k1 = Ts * pendulum_ode(x, u);
k2 = Ts * pendulum_ode(x + k1./2, u);
k3 = Ts * pendulum_ode(x + k2./2, u);
k4 = Ts * pendulum_ode(x + k3, u);

x_next = x + 1/6 * (k1 + 2*k2 + 2*k3 + k4);

end

