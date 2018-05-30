
import casadi.*
num_free_masses = 1;
nx = (2*num_free_masses + 1) * 3;
nu = 3;

x = SX.sym('x', nx);
u = SX.sym('u', nu);

f = hanging_chain_ode(x, u, num_free_masses);

pendulum = Function('hanging_chain', {x, u}, {f});
pendulum.generate('hanging_chain', struct('casadi_int', 'int'));

vec = SX.sym('vec', nx);

adj_df_dx = Function('adj_df_dx', {x, u, vec}, {jtimes(f, x, vec, true)});
adj_df_dx.generate('adj_df_dx', struct('casadi_int', 'int'));

adj_df_du = Function('adj_df_du', {x, u, vec}, {jtimes(f, u, vec, true)});
adj_df_du.generate('adj_df_du', struct('casadi_int', 'int'));
