function [xhatHist, PHist, tMeasHist, statIdHist, residHist, residCovHist, xhat_deliv, P_deliv] = ExtendedKalmanFilter...
         (X0, P0, tMeasMax, tPropMax, JD0_UTC, Sensor_data, fitCase, Params)

    dataMaxIdx = find(Sensor_data.LEO_DATA_Apparent(:,2) == tMeasMax);

    nx = length(X0);
    ny = 2;

    xhatHist = zeros(dataMaxIdx, nx);
    PHist = zeros(dataMaxIdx, nx^2);
    tMeasHist = zeros(dataMaxIdx, 1);
    statIdHist = zeros(dataMaxIdx, 1);
    residHist = zeros(dataMaxIdx, ny);
    residCovHist = zeros(dataMaxIdx, ny^2);

    % First iteration
    t = Sensor_data.LEO_DATA_Apparent(1, 2); % meas time since epoch (current meas time)
    stat_id = Sensor_data.LEO_DATA_Apparent(1, 1);
    if stat_id == 1
        Rk = diag([Params.Kwaj.sigma_r, Params.Kwaj.sigma_rr]).^2;
        r_stat_ITRF = Params.Kwaj.r_ITRF;
    elseif stat_id == 2
        Rk = diag([Params.DG.sigma_r, Params.DG.sigma_rr]).^2;
        r_stat_ITRF = Params.DG.r_ITRF;
    elseif stat_id == 3
        Rk = diag([Params.Arecibo.sigma_r, Params.Arecibo.sigma_rr]).^2;
        r_stat_ITRF = Params.Arecibo.r_ITRF;
    end

    if fitCase == 'A'
        Rk(2,2) = 1e9; % Fit range only
    elseif fitCase == 'B'
        Rk(1,1) = 1e9; % Fit range-rate only
    end

    % First iteration
    xbark = X0;
    Pbark = P0;

    JDt_UTC = JD0_UTC + t / 86400;
    [lt, X_lt_corr] = LightTimeCorrection(xbark, t, JDt_UTC, Sensor_data, Params);
    Params = UpdateParams(JDt_UTC, Params);
    W = PolarMatrix(Params);
    R = SiderealMatrix(Params);
    N = NutationMatrix(Params);
    P = PrecessionMatrix(Params);
    r_stat_GCRF = P * N * R * W * r_stat_ITRF;
    v_stat_GCRF = P * N * R * (W * zeros(3,1) + cross([0; 0; Params.w_earth], W * r_stat_ITRF));
    sx = r_stat_GCRF(1);
    sy = r_stat_GCRF(2);
    sz = r_stat_GCRF(3);
    svx = v_stat_GCRF(1);
    svy = v_stat_GCRF(2);
    svz = v_stat_GCRF(3);
    H = HmatMeas(X_lt_corr(1), X_lt_corr(2), X_lt_corr(3), X_lt_corr(4), X_lt_corr(5), X_lt_corr(6), sx, sy, sz ,svx, svy, svz, X_lt_corr(8));

    yk_obs = Sensor_data.LEO_DATA_Apparent(1,3:4)' .* 1e3;
    yk_comp = h_meas(xbark, t, stat_id, JD0_UTC, Sensor_data, Params);
    yk = yk_obs - yk_comp;

    Sk = H * Pbark * H' + Rk;
    Kk = Pbark * H' * inv(Sk);

    % need to add 'G' case (short-arc)
    if fitCase == 'C' && stat_id ~= 1
        Kk = zeros(nx, ny);
    elseif fitCase == 'D' && stat_id ~= 2
        Kk = zeros(nx, ny);
    elseif fitCase == 'E' && stat_id ~= 3
        Kk = zeros(nx, ny);
    end
    %Kk(7,:) = zeros(1,ny); % consider Cd

    xhatk = xbark + Kk * yk;
    Pk = (eye(nx) - Kk * H) * Pbark * (eye(nx) - Kk * H)' + Kk * Rk * Kk';

    xhatHist(1,:) = xhatk;
    PHist(1,:) = Pk(:);
    tMeasHist(1) = t;
    statIdHist(1) = stat_id;
    residHist(1,:) = yk;
    residCovHist(1,:) = Sk(:);
    
    for kk = 2 : dataMaxIdx

        t = Sensor_data.LEO_DATA_Apparent(kk, 2); % meas time since epoch (current meas time)
        tprev = Sensor_data.LEO_DATA_Apparent(kk-1, 2);
        stat_id = Sensor_data.LEO_DATA_Apparent(kk, 1);
        if stat_id == 1
            Rk = diag([Params.Kwaj.sigma_r, Params.Kwaj.sigma_rr]).^2;
            r_stat_ITRF = Params.Kwaj.r_ITRF;
        elseif stat_id == 2
            Rk = diag([Params.DG.sigma_r, Params.DG.sigma_rr]).^2;
            r_stat_ITRF = Params.DG.r_ITRF;
        elseif stat_id == 3
            Rk = diag([Params.Arecibo.sigma_r, Params.Arecibo.sigma_rr]).^2;
            r_stat_ITRF = Params.Arecibo.r_ITRF;
        end

        if fitCase == 'A'
            Rk(2,2) = 1e9; % Fit range only
        elseif fitCase == 'B'
            Rk(1,1) = 1e9; % Fit range-rate only
        end
    
        % Prediction
        Phik0 = eye(nx);
        Xk0 = [xhatk; Phik0(:)];
        [~, XkHist] = ode113(@orbPropSTM, [tprev; t], Xk0, Params.ode_options, JD0_UTC, Params);

        xbark = XkHist(end, 1:nx)';
        R = xbark(1:3) / norm(xbark(1:3));
        C = cross(xbark(1:3), xbark(4:6)) / norm(cross(xbark(1:3), xbark(4:6)));
        I = cross(C, R);
        T_GCRF_RIC = [R, I, C];
        Qk = T_GCRF_RIC * Params.Qk_RIC * T_GCRF_RIC';
        dt = t - tprev;
        if length(X0) == 6
            Gammak = [dt^2/2 * eye(3); dt * eye(3)];
        elseif length(X0) > 6
            Gammak = [dt^2/2 * eye(3); dt * eye(3); zeros(length(X0)-6,3)];
        end
        Phik = reshape(XkHist(end, nx+1:end), [nx,nx]);
        Pbark = Phik * Pk * Phik' + Gammak * Qk * Gammak';
    
        % Measurement update
        JDt_UTC = JD0_UTC + t / 86400;
        [lt, X_lt_corr] = LightTimeCorrection(xbark, t, JDt_UTC, Sensor_data, Params);
        Params = UpdateParams(JDt_UTC, Params);
        W = PolarMatrix(Params);
        R = SiderealMatrix(Params);
        N = NutationMatrix(Params);
        P = PrecessionMatrix(Params);
        r_stat_GCRF = P * N * R * W * r_stat_ITRF;
        v_stat_GCRF = P * N * R * (W * zeros(3,1) + cross([0; 0; Params.w_earth], W * r_stat_ITRF));
        sx = r_stat_GCRF(1);
        sy = r_stat_GCRF(2);
        sz = r_stat_GCRF(3);
        svx = v_stat_GCRF(1);
        svy = v_stat_GCRF(2);
        svz = v_stat_GCRF(3);
        H = HmatMeas(X_lt_corr(1), X_lt_corr(2), X_lt_corr(3), X_lt_corr(4), X_lt_corr(5), X_lt_corr(6), sx, sy, sz ,svx, svy, svz, X_lt_corr(8));
    
        yk_obs = Sensor_data.LEO_DATA_Apparent(kk,3:4)' .* 1e3;
        yk_comp = h_meas(xbark, t, stat_id, JD0_UTC, Sensor_data, Params);
        yk = yk_obs - yk_comp;
    
        Sk = H * Pbark * H' + Rk;
        Kk = Pbark * H' * inv(Sk);
        
        % need to add 'G' case (short-arc)
        if fitCase == 'C' && stat_id ~= 1
            Kk = zeros(nx, ny);
        elseif fitCase == 'D' && stat_id ~= 2
            Kk = zeros(nx, ny);
        elseif fitCase == 'E' && stat_id ~= 3
            Kk = zeros(nx, ny);
        end
        %Kk(7,:) = zeros(1,ny); % consider Cd

        xhatk = xbark + Kk * yk;
        Pk = (eye(nx) - Kk * H) * Pbark * (eye(nx) - Kk * H)' + Kk * Rk * Kk';
    
        xhatHist(kk,:) = xhatk;
        PHist(kk,:) = Pk(:);
        tMeasHist(kk) = t;
        statIdHist(kk) = stat_id;
        residHist(kk,:) = yk;
        residCovHist(kk,:) = Sk(:);

    end

    % % Propagate to delivery
    xhat_deliv = zeros(nx, 1);
    P_deliv = zeros(nx);

    Phik0 = eye(nx);
    Xk0 = [xhatk; Phik0(:)];
    [~, XkHist] = ode113(@orbPropSTM, tMeasMax:60:tPropMax, Xk0, Params.ode_options, JD0_UTC, Params);

    xbark = XkHist(end, 1:nx)';
    R = xbark(1:3) / norm(xbark(1:3));
    C = cross(xbark(1:3), xbark(4:6)) / norm(cross(xbark(1:3), xbark(4:6)));
    I = cross(C, R);
    T_GCRF_RIC = [R, I, C];
    Qk = T_GCRF_RIC * Params.Qk_RIC * T_GCRF_RIC';
    dt = t - tprev;
    if length(X0) == 6
        Gammak = [dt^2/2 * eye(3); dt * eye(3)];
    elseif length(X0) > 6
        Gammak = [dt^2/2 * eye(3); dt * eye(3); zeros(length(X0)-6,3)];
    end
    Phik = reshape(XkHist(end, nx+1:end), [nx,nx]);
    Pbark = Phik * Pk * Phik' + Gammak * Qk * Gammak';

    xhatk = xbark;
    Pk = Pbark;

    xhat_deliv = xhatk;
    P_deliv = Pk;

end