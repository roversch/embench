function dx_dt = pendulum_ode(x, u)
%PENDULUM_ODE Explicit ODE that models the movement of the inverted
%pendulum on a cart. The states are position p [m], velocity v [m/s] of the
%cart and angle theta [rad] of the rod with the vertical and angular
%velocity omega [rad]. The input is the horizontal force F [N] exerted on
%the cart.

m = 0.1; % mass of the pendulum
M = 1.0; % mass of the cart
l = 0.8; % length of the rod
g = 9.81; % gravitational acceleration

v = x(3);
theta = x(2);
omega = x(4);
F = u;

denom = M+m-m*cos(theta).^2;

dx_dt = [
    v; ...
    omega; ...
    (-m*l*sin(theta).*omega.^2 + m*g*cos(theta).*sin(theta) + F)./denom; ... 
    (-m*l*cos(theta).*sin(theta).*omega.^2 + F.*cos(theta) + (M+m)*g*sin(theta))./(l*denom)
  ];

end

