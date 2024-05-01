function Htilde = HtildeMeasurements(X, t, JD_UTC, stat_id, Params)

    X_lt_corr = LightTimeCorrection(X, t, JD_UTC, Params);
    Params = UpdateParams(JD_UTC, Params);

    W = PolarMatrix(Params);
    R = SiderealMatrix(Params);
    N = NutationMatrix(Params);
    P = PrecessionMatrix(Params);

    if stat_id == 1
        r_stat_ITRF = Params.r_Kwaj_ITRF;
    elseif stat_id == 2
        r_stat_ITRF = Params.r_DG_ITRF;
    elseif stat_id == 3
        r_stat_ITRF = Params.r_Arecibo_ITRF;
    end
    r_stat_GCRF = P * N * R * W * r_stat_ITRF;
    v_stat_GCRF = P * N * R * (W * zeros(3,1) + cross([0; 0; Params.w_earth], W * r_stat_ITRF));
    sx = r_stat_GCRF(1);
    sy = r_stat_GCRF(2);
    sz = r_stat_GCRF(3);
    svx = v_stat_GCRF(1);
    svy = v_stat_GCRF(2);
    svz = v_stat_GCRF(3);

    syms rx ry rz vx vy vz

    range = sqrt((rx - sx)^2 + (ry - sy)^2 + (rz - sz)^2);
    rangeRate = ((rx - sx) * (vx - svx) + (ry - sy) * (vy - svy) + (rz - sz) * (vz - svz)) / range;

    Y_obs = [range; rangeRate];

    dYdrx = diff(Y_obs, rx); dYdrx = subs(dYdrx, [rx; ry; rz; vx; vy; vz], X_lt_corr(1:6));
    dYdry = diff(Y_obs, ry); dYdry = subs(dYdry, [rx; ry; rz; vx; vy; vz], X_lt_corr(1:6));
    dYdrz = diff(Y_obs, rz); dYdrz = subs(dYdrz, [rx; ry; rz; vx; vy; vz], X_lt_corr(1:6));
    dYdvx = [0; diff(rangeRate, vx)]; dYdvx = subs(dYdvx, [rx; ry; rz; vx; vy; vz], X_lt_corr(1:6));
    dYdvy = [0; diff(rangeRate, vy)]; dYdvy = subs(dYdvy, [rx; ry; rz; vx; vy; vz], X_lt_corr(1:6));
    dYdvz = [0; diff(rangeRate, vz)]; dYdvz = subs(dYdvz, [rx; ry; rz; vx; vy; vz], X_lt_corr(1:6));
    dYdCd = zeros(2, 1);

    Htilde = [dYdrx, dYdry, dYdrz, dYdvx, dYdvy, dYdvz, dYdCd];
    Htilde = double(Htilde);

end