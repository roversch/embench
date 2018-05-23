
import casadi.*

nx = 4;
nu = 1;

x = SX.sym('x', nx);
u = SX.sym('u', nu);

f = pendulum_ode(x, u);

pendulum = Function('pendulum', {x, u}, {f});
pendulum.generate('pendulum', struct('casadi_int', 'int'));

vec = SX.sym('vec', nx);

adj_df_dx = Function('adj_df_dx', {x, u, vec}, {jtimes(f, x, vec, true)});
adj_df_dx.generate('adj_df_dx', struct('casadi_int', 'int'));

adj_df_du = Function('adj_df_du', {x, u, vec}, {jtimes(f, u, vec, true)});
adj_df_du.generate('adj_df_du', struct('casadi_int', 'int'));



