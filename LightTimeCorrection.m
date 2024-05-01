function [lt, X_lt_corr] = LightTimeCorrection(X, t, JD_UTC, Sensor_data, Params)

    Params = UpdateParams(JD_UTC, Params);

    r_GCRF = X(1:3);

    range_sensor = Sensor_data.LEO_DATA_Apparent(Sensor_data.LEO_DATA_Apparent(:,2) == t, 3) * 1e3;
    lt = range_sensor / Params.c;

    t_lt = t - lt;
    JD_UTC_lt = JD_UTC - lt / 86400;

    [~, X_lt] = ode113(@orbPropHF, [t; t_lt], X, Params.ode_options, JD_UTC, Params);
    Params = UpdateParams(JD_UTC_lt, Params);

    W = PolarMatrix(Params);
    R = SiderealMatrix(Params);
    N = NutationMatrix(Params);
    P = PrecessionMatrix(Params);

    stat_id = Sensor_data.LEO_DATA_Apparent(Sensor_data.LEO_DATA_Apparent(:,2) == t, 1);
    if stat_id == 1
        r_stat_ITRF = Params.Kwaj.r_ITRF;
    elseif stat_id == 2
        r_stat_ITRF = Params.DG.r_ITRF;
    elseif stat_id == 3
        r_stat_ITRF = Params.Arecibo.r_ITRF;
    end

    r_lt_corr_ITRF = W' * R' * N' * P' * X_lt(end,1:3)';

    range = norm(r_lt_corr_ITRF - r_stat_ITRF);

    lt = range / Params.c;

    delta = norm(r_lt_corr_ITRF - W' * R' * N' * P' * r_GCRF);
    tol = 1e-3;

    while delta > tol

        t_lt = t - lt;
        JD_UTC_lt = JD_UTC - lt / 86400;
    
        [~, X_lt] = ode113(@orbPropHF, [t; t_lt], X, Params.ode_options, JD_UTC, Params);
        Params = UpdateParams(JD_UTC_lt, Params);
    
        W = PolarMatrix(Params);
        R = SiderealMatrix(Params);
        N = NutationMatrix(Params);
        P = PrecessionMatrix(Params);
    
        stat_id = Sensor_data.LEO_DATA_Apparent(Sensor_data.LEO_DATA_Apparent(:,2) == t, 1);
        if stat_id == 1
            r_stat_ITRF = Params.Kwaj.r_ITRF;
        elseif stat_id == 2
            r_stat_ITRF = Params.DG.r_ITRF;
        elseif stat_id == 3
            r_stat_ITRF = Params.Arecibo.r_ITRF;
        end
    
        r_new_lt_corr_ITRF = W' * R' * N' * P' * X_lt(end,1:3)';
    
        range = norm(r_new_lt_corr_ITRF - r_stat_ITRF);
    
        lt = range / Params.c;

        delta = norm(r_new_lt_corr_ITRF - r_lt_corr_ITRF);
        r_lt_corr_ITRF = r_new_lt_corr_ITRF;

    end
 
    X_lt_corr = X_lt(end,:)';

end