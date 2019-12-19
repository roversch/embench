function model = pendulum_ode_model()
    %PENDULUM_ODE Explicit ODE that models the movement of the inverted
    %pendulum on a cart. The states are position p [m], velocity v [m/s] of the
    %cart and angle theta [rad] of the rod with the vertical and angular
    %velocity omega [rad]. The input is the horizontal force F [N] exerted on
    %the cart.
    
    import casadi.*

    m = 0.1; % mass of the pendulum
    M = 1.0; % mass of the cart
    l = 0.8; % length of the rod
    g = 9.81; % gravitational acceleration

    %% named symbolic variables
    p = SX.sym('p');         % horizontal displacement of cart [m]
    theta = SX.sym('theta'); % angle of rod with the vertical [rad]
    v = SX.sym('v');         % horizontal velocity of cart [m/s]
    omega = SX.sym('omega'); % angular velocity of rod [rad/s]
    F = SX.sym('F');         % horizontal force acting on cart [N]

    %% (unnamed) symbolic variables
    sym_x = vertcat(p, theta, v, omega);
    model.nx = length(sym_x);
    sym_xdot = SX.sym('xdot', model.nx, 1);
    sym_u = F;

    %% ODE

    denom = M+m-m*cos(theta).^2;

    sin_theta = sin(theta);
    cos_theta = cos(theta);
    denom = M + m - m*cos(theta).^2;

    expr_f_expl = vertcat(v, ...
                          omega, ...
                          (- l*m*sin_theta*omega.^2 + F + g*m*cos_theta*sin_theta)/denom, ...
                          (- l*m*cos_theta*sin_theta*omega.^2 + F*cos_theta + g*m*sin_theta + M*g*sin_theta)/(l*denom));


    model.sym_x = sym_x;
    model.sym_xdot = sym_xdot;
    model.sym_u = sym_u;
    model.expr_f_expl = expr_f_expl;
end

