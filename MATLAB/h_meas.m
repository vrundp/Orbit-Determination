function Y = h_meas(X, t, stat_id, JD0_UTC, Sensor_data, Params)
% Make sure to correct for light-time!
% t = meas time since epoch (current meas time)

    JD_UTC = JD0_UTC + t / 86400; % current true time
    [lt, X_lt_corr] = LightTimeCorrection(X, t, JD_UTC, Sensor_data, Params);
    Params = UpdateParams(JD_UTC, Params);

    rx = X_lt_corr(1);
    ry = X_lt_corr(2);
    rz = X_lt_corr(3);
    vx = X_lt_corr(4);
    vy = X_lt_corr(5);
    vz = X_lt_corr(6);

    W = PolarMatrix(Params);
    R = SiderealMatrix(Params);
    N = NutationMatrix(Params);
    P = PrecessionMatrix(Params);

    b = 0;

    if stat_id == 1
        r_stat_ITRF = Params.Kwaj.r_ITRF;
    elseif stat_id == 2
        r_stat_ITRF = Params.DG.r_ITRF;
    elseif stat_id == 3
        r_stat_ITRF = Params.Arecibo.r_ITRF;
        b = X(8);
        %b = X(7); % consider filter
        %b = 20;
    end
    r_stat_GCRF = P * N * R * W * r_stat_ITRF;
    v_stat_GCRF = P * N * R * (W * zeros(3,1) + cross([0; 0; Params.w_earth], W * r_stat_ITRF));
    sx = r_stat_GCRF(1);
    sy = r_stat_GCRF(2);
    sz = r_stat_GCRF(3);
    svx = v_stat_GCRF(1);
    svy = v_stat_GCRF(2);
    svz = v_stat_GCRF(3);

    % range = sqrt((rx - sx)^2 + (ry - sy)^2 + (rz - sz)^2);
    % rangeRate = ((rx - sx) * (vx - svx) + (ry - sy) * (vy - svy) + (rz - sz) * (vz - svz)) / range;
    range = sqrt((rx - sx + lt * svx)^2 + (ry - sy + lt * svy)^2 + (rz - sz + lt * svz)^2);
    rangeRate = ((rx - sx + lt * svx) * (vx - svx) + (ry - sy + lt * svy) * (vy - svy) + (rz - sz + lt * svz) * (vz - svz)) / range;

    Y = [range + b; rangeRate];

end