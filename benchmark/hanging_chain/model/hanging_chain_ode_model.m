function model = hanging_chain_ode_model(num_free_masses)    
    import casadi.*

    %% named symbolic variables
    p = SX.sym('p', 3*num_free_masses);
    p_end = SX.sym('p_end', 3);
    v = SX.sym('v', 3*num_free_masses);
    u = SX.sym('u', 3);

    expr_f_expl = hanging_chain_ode([p; p_end; v], u, num_free_masses);

    nx = (num_free_masses * 2 + 1) * 3;
    nu = 3;

    sym_x = vertcat( p, p_end, v);
    sym_xdot = SX.sym('xdot', nx, 1);

    model.sym_x = sym_x;
    model.sym_xdot = sym_xdot;
    model.sym_u = u;
    model.expr_f_expl = expr_f_expl;
    model.nx = nx;
    model.nu = nu;
end
