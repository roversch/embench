function ode_rhs = hanging_chain_ode(x, u, num_free_masses)

p = x(1:3*num_free_masses);
v = x(3*num_free_masses+1:end-3);
p_end = x(end-2:end);

p0 = [0; 0; 0];
L  = 0.033;
D  = 1.0;
m  = 0.03;

g = [0; 0; -9.81];
f = v * 0;
f = f + repmat(g, num_free_masses, 1);
for i = 1:num_free_masses+1
    if i == 1
        dist = p((i-1)*3+1:i*3) - p0;
    elseif( i <= num_free_masses )
        dist = p((i-1)*3+1:i*3) - p((i-2)*3+1:(i-1)*3);
    else
        dist = p_end - p((num_free_masses-1)*3+1:end);
    end
    
    scale = D / m * (1-L./sqrt(dist(1).^2 + dist(2).^2 + dist(3).^2));
    F = scale .* dist;
        
    % mass on the right
    if i < num_free_masses+1
        f((i-1)*3+1:i*3) = f((i-1)*3+1:i*3) - F;
    end
    % mass on the left
    if i > 1
        f((i-2)*3+1:(i-1)*3) = f((i-2)*3+1:(i-1)*3) + F;
    end
end

ode_rhs = [ v; ...
            f; ...
            u];

end

