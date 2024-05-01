function Amat = AmatDynamics(r_GCRF, v_GCRF, Cd, Params)

    r_sun_GCRF = SunPosition(Params);
    r_moon_GCRF = MoonPosition(Params);

    syms rx ry rz vx vy vz Cdrag
    r_sat = [rx; ry; rz];
    v_sat = [vx; vy; vz];
    r = sqrt(rx^2 + ry^2 + rz^2);
    
    % Two body
    a_2body = -Params.mu_earth / r^3 * r_sat;
    % J2
    a_J2 = -1.5 * Params.J2 * Params.mu_earth * Params.R_earth^2 / r^5 .* [rx * (1 - 5 * rz^2 / r^2); ry * (1 - 5 * rz^2 / r^2); rz * (3 - 5 * rz^2 / r^2)];
    % J3
    a_J3 = -2.5 * Params.J3 * Params.mu_earth * Params.R_earth^3 / r^7 .* [rx * (3 * rz - 7 * rz^3 / r^2); ry * (3 * rz - 7 * rz^3 / r^2); 6 * rz^2 - 7 * rz^4 / r^2 - 0.6 * r^2];
    % J4
    a_J4 = 1.875 * Params.J4 * Params.mu_earth * Params.R_earth^4 / r^7 .* [rx * (1 - 14 * rz^2 / r^2 + 21 * rz^4 / r^4); ry * (1 - 14 * rz^2 / r^2 + 21 * rz^4 / r^4); rz * (5 - 70 * rz^2 / 3 / r^2 + 21 * rz^4 / r^4)];
    
    % Third body Sun & Moon
    a_3bodySun = Params.mu_sun * ((r_sun_GCRF - r_sat) / norm(r_sun_GCRF - r_sat)^3 - r_sun_GCRF / norm(r_sun_GCRF)^3);
    a_3bodyMoon = Params.mu_moon * ((r_moon_GCRF - r_sat) / norm(r_moon_GCRF - r_sat)^3 - r_moon_GCRF / norm(r_moon_GCRF)^3);

    % Drag
    v_rel = v_sat - cross([0; 0; Params.w_earth], r_sat);
    
    A_drag = Params.FaceX.area;
    rho = Params.rho0 * exp(-(r - Params.r0) / Params.H);

    a_drag = -0.5 * rho * Cdrag * A_drag / Params.mass * norm(v_rel) * v_rel;

    % SRP
    F = EclipseShadow(r_GCRF, Params);
    p_srp = 1367 / Params.c;

    r_sat_sun = r_sun_GCRF - r_sat;
    uhat = r_sat_sun / norm(r_sat_sun);
    d = norm(r_sat_sun) / Params.AU;
    
    nsp = uhat;
    
    cos_theta_nsp_pos = dot(uhat, nsp);
    
    A_srp = Params.SolarPanel.area * ((2 / 3 * Params.SolarPanel.Cd * cos_theta_nsp_pos + 2 * Params.SolarPanel.Cs * cos_theta_nsp_pos^2) .* nsp ...
                                   + ((1 - Params.SolarPanel.Cs / 2) * cos_theta_nsp_pos^2) .* uhat);
    
    a_srp = -p_srp * F .* A_srp ./ Params.mass ./ d^2;
    
    % Total Acceleration
    a_tot = a_2body + a_J2 + a_J3 + a_J4 + a_3bodySun + a_3bodyMoon + a_drag + a_srp;
    %a_num = double(subs(a_tot, [rx; ry; rz; vx; vy; vz; Cdrag], [r_GCRF; v_GCRF; Cd]))

    % Partial Derivatives
    dadrx = diff(a_tot, rx);    dadrx = subs(dadrx, [rx; ry; rz; vx; vy; vz; Cdrag], [r_GCRF; v_GCRF; Cd]);
    dadry = diff(a_tot, ry);    dadry = subs(dadry, [rx; ry; rz; vx; vy; vz; Cdrag], [r_GCRF; v_GCRF; Cd]);
    dadrz = diff(a_tot, rz);    dadrz = subs(dadrz, [rx; ry; rz; vx; vy; vz; Cdrag], [r_GCRF; v_GCRF; Cd]);
    dadvx = diff(a_tot, vx);    dadvx = subs(dadvx, [rx; ry; rz; vx; vy; vz; Cdrag], [r_GCRF; v_GCRF; Cd]);
    dadvy = diff(a_tot, vy);    dadvy = subs(dadvy, [rx; ry; rz; vx; vy; vz; Cdrag], [r_GCRF; v_GCRF; Cd]);
    dadvz = diff(a_tot, vz);    dadvz = subs(dadvz, [rx; ry; rz; vx; vy; vz; Cdrag], [r_GCRF; v_GCRF; Cd]);
    dadCd = diff(a_tot, Cdrag); dadCd = subs(dadCd, [rx; ry; rz; vx; vy; vz; Cdrag], [r_GCRF; v_GCRF; Cd]) / 1e3;

    Amat = [zeros(3), eye(3), zeros(3, 1); ...
            dadrx, dadry, dadrz, dadvx, dadvy, dadvz, dadCd; ...
            zeros(1, 7)];

    Amat = double(Amat);

end