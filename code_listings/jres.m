function fid = jres_seq(spin_system, parameters, H, R, K)
    % Compose Liouvillian
    L = H + 1i * R + 1i * K; clear('H', 'R', 'K');
    % Coherent evolution timestep
    d1 = 1 / parameters.sweep(1);
    d2 = 1 / parameters.sweep(2);
    % Number of points
    pts1 = parameters.npoints(1);
    pts2 = parameters.npoints(2);
    % Targeted nucleus
    nuc = parameters.spins{1}
    % Initial state
    rho = state(spin_system, 'Lz', nuc);
    % Detection state
    coil = state(spin_system, 'L+', nuc);
    % Get the pulse operators
    Lp = operator(spin_system, 'L+', nuc);
    Lx = (Lp + Lp') / 2;
    % First pulse: 90 about x
    rho = step(spin_system, Lx, rho, pi / 2);
    % First half of t1 evolution
    rho = evolution(spin_system, L, [], rho, d1 / 2, pts1 - 1, 'trajectory');
    % Select "-1" coherence
    rho = coherence(spin_system, rho, {{nuc, -1}});
    % Second pulse: 180 about x
    rho = step(spin_system, Lx, rho, pi);
    % Select "+1" coherence
    rho = coherence(spin_system, rho, {{nuc, +1}});
    % Second half of t1 evolution
    rho = evolution(spin_system, L, [], rho, d1 / 2, pts1 - 1, 'refocus');
    % Run the F2 evolution
    fid = evolution(spin_system, L, coil, rho, d2, pts2 - 1, 'observable');
end
