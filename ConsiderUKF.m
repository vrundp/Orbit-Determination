function [zhatHist, PzzHist, tMeasHist, statIdHist, residHist, residCovHist, xhat_deliv, P_deliv] = ConsiderUKF...
         (X0, Pxx0, c0, Pcc0, tMeasMax, tPropMax, JD0_UTC, Sensor_data, fitCase, Params)

    dataMaxIdx = find(Sensor_data.LEO_DATA_Apparent(:,2) == tMeasMax);

    nx = length(X0);
    nc = length(c0);
    nz = nx + nc;
    ny = 2;
    w = 1 / (2 * (nz));

    if fitCase == 'G'
        startIdx = 1290;
    else
        startIdx = 1;
    end

    zhatHist = zeros(dataMaxIdx-startIdx+1, nz);
    PzzHist = zeros(dataMaxIdx-startIdx+1, (nz)^2);
    tMeasHist = zeros(dataMaxIdx-startIdx+1, 1);
    statIdHist = zeros(dataMaxIdx-startIdx+1, 1);
    residHist = zeros(dataMaxIdx-startIdx+1, ny);
    residCovHist = zeros(dataMaxIdx-startIdx+1, ny^2);

    % First iteration
    t = Sensor_data.LEO_DATA_Apparent(startIdx, 2); % meas time since epoch (current meas time)
    stat_id = Sensor_data.LEO_DATA_Apparent(startIdx, 1);
    if stat_id == 1
        Rk = diag([Params.Kwaj.sigma_r, Params.Kwaj.sigma_rr]).^2;
    elseif stat_id == 2
        Rk = diag([Params.DG.sigma_r, Params.DG.sigma_rr]).^2;
    elseif stat_id == 3
        Rk = diag([Params.Arecibo.sigma_r, Params.Arecibo.sigma_rr]).^2;
    end

    if fitCase == 'A'
        Rk(2,2) = 1e9; % Fit range only
    elseif fitCase == 'B'
        Rk(1,1) = 1e9; % Fit range-rate only
    end
    
    z0 = [X0; c0];
    Pzz0 = [Pxx0, zeros(nx, nc); zeros(nc, nx), Pcc0];
    Spkm1 = chol(Pzz0, 'lower');
    sigmaZpbarkHist = zeros(2 * (nz), nz);
    % Initialize Sigma Points and Propagate
    for ii = 1 : 2 * (nz)
        jj = ii;
        sign = 1;
        if ii > nz
            jj = ii - nz;
            sign = -1;
        end
        sigmaZk = z0 + sign * sqrt(nz) * Spkm1(:,jj);
        sigmaZpbarkHist(ii,:) = sigmaZk;
    end
    % Initial prediction
    zbark = 0;
    for ii = 1 : 2 * (nz)
        zbark = zbark + w * sigmaZpbarkHist(ii,:)';
    end
    R = X0(1:3) / norm(X0(1:3));
    C = cross(X0(1:3), X0(4:6)) / norm(cross(X0(1:3), X0(4:6)));
    I = cross(C, R);
    T_GCRF_RIC = [R, I, C];
    Qk = T_GCRF_RIC * Params.Qk_RIC * T_GCRF_RIC';
    dt = 0;
    if length(X0) == 6
        Gammak = [dt^2/2 * eye(3); dt * eye(3); zeros(nc, 3)];
    elseif length(X0) > 6
        Gammak = [dt^2/2 * eye(3); dt * eye(3); zeros(length(X0)-6,3); zeros(nc, 3)];
    end
    Pzzbark = 0;
    for ii = 1 : 2 * (nz)
        Pzzbark = Pzzbark + w * (sigmaZpbarkHist(ii,:)' - zbark) * (sigmaZpbarkHist(ii,:)' - zbark)';
    end
    Pzzbark = Pzzbark + Gammak * Qk * Gammak';
    Pxxbark = Pzzbark(1:nx,1:nx);
    Pxcbark = Pzzbark(1:nx,nx+1:end);
    Pcxbark = Pzzbark(nx+1:end,1:nx);
    Pccbark = Pzzbark(nx+1:end,nx+1:end);
    % Initial measurement update
    Suk = chol(Pzzbark, 'lower');
    sigmaZubarkHist = zeros(2 * (nz), nz);
    YkHist = zeros(2 * (nz), ny);
    for ii = 1 : 2 * (nz)
        jj = ii;
        sign = 1;
        if ii > nz
            jj = ii - nz;
            sign = -1;
        end
        sigmaZk = zbark + sign * sqrt(nz) * Suk(:,jj);
        sigmaZubarkHist(ii,:) = sigmaZk;
        YkHist(ii,:) = h_meas(sigmaZk, t, stat_id, JD0_UTC, Sensor_data, Params);
    end
    ybark = 0;
    for ii = 1 : 2 * (nz)
        ybark = ybark + w * YkHist(ii,:)';
    end
    Pyyk = 0;
    Pzyk = 0;
    for ii = 1 : 2 * (nz)
        Pyyk = Pyyk + w * (YkHist(ii,:)' - ybark) * (YkHist(ii,:)' - ybark)';
        Pzyk = Pzyk + w * (sigmaZubarkHist(ii,:)' - zbark) * (YkHist(ii,:)' - ybark)';
    end
    Pyyk = Pyyk + Rk;
    Kzk = Pzyk * inv(Pyyk);
    Kxk = Kzk(1:nx,:);
    Kck = Kzk(nx+1:end,:);

    if fitCase == 'C' && stat_id ~= 1
        Kxk = zeros(nx, ny);
    elseif fitCase == 'D' && stat_id ~= 2
        Kxk = zeros(nx, ny);
    elseif fitCase == 'E' && stat_id ~= 3
        Kxk = zeros(nx, ny);
    end

    % Initial state and covariance update
    yk = Sensor_data.LEO_DATA_Apparent(startIdx,3:4)' .* 1e3;
    zhatk = zbark + [Kxk; zeros(nc,ny)] * (yk - ybark);
    Pzzk = [Pxxbark, Pxcbark; Pcxbark, Pccbark] - [Kxk * Pyyk * Kxk', Kxk * Pyyk * Kck'; Kck * Pyyk * Kxk', zeros(nc, nc)];

    zhatHist(1,:) = zhatk;
    PzzHist(1,:) = Pzzk(:);
    tMeasHist(1) = t;
    statIdHist(1) = stat_id;
    residHist(1,:) = yk - ybark;
    residCovHist(1,:) = Pyyk(:);

    for kk = (startIdx + 1) : dataMaxIdx

        t = Sensor_data.LEO_DATA_Apparent(kk, 2);
        tprev = Sensor_data.LEO_DATA_Apparent(kk-1, 2);
        stat_id = Sensor_data.LEO_DATA_Apparent(kk, 1);
 
        if stat_id == 1
            Rk = diag([Params.Kwaj.sigma_r, Params.Kwaj.sigma_rr]).^2;
        elseif stat_id == 2
            Rk = diag([Params.DG.sigma_r, Params.DG.sigma_rr]).^2;
        elseif stat_id == 3
            Rk = diag([Params.Arecibo.sigma_r, Params.Arecibo.sigma_rr]).^2;
        end

        if fitCase == 'A'
            Rk(2,2) = 1e9; % Fit range only
        elseif fitCase == 'B'
            Rk(1,1) = 1e9; % Fit range-rate only
        end

        Spkm1 = chol(Pzzk, 'lower');
        sigmaZpbarkHist = zeros(2 * (nz), nz);
        % Initialize Sigma Points and Propagate
        parfor ii = 1 : 2 * (nz)
            jj = ii;
            sign = 1;
            if ii > nz
                jj = ii - nz;
                sign = -1;
            end
            sigmaZk = zhatk + sign * sqrt(nz) * Spkm1(:,jj);
            [~, ZkHist] = ode113(@orbPropHF, [tprev; t], sigmaZk, Params.ode_options, JD0_UTC, Params); % ACCOUNT FOR CD IN PROPAGATOR!!!
            sigmaZpbarkHist(ii,:) = ZkHist(end,:);
        end
        % Initial prediction
        zbark = 0;
        for ii = 1 : 2 * (nz)
            zbark = zbark + w * sigmaZpbarkHist(ii,:)';
        end
        R = zbark(1:3) / norm(zbark(1:3));
        C = cross(zbark(1:3), zbark(4:6)) / norm(cross(zbark(1:3), zbark(4:6)));
        I = cross(C, R);
        T_GCRF_RIC = [R, I, C];
        Qk = T_GCRF_RIC * Params.Qk_RIC * T_GCRF_RIC';
        dt = t - tprev;
        if length(X0) == 6
            Gammak = [dt^2/2 * eye(3); dt * eye(3); zeros(nc, 3)];
        elseif length(X0) > 6
            Gammak = [dt^2/2 * eye(3); dt * eye(3); zeros(length(X0)-6,3); zeros(nc, 3)];
        end
        Pzzbark = 0;
        for ii = 1 : 2 * (nz)
            Pzzbark = Pzzbark + w * (sigmaZpbarkHist(ii,:)' - zbark) * (sigmaZpbarkHist(ii,:)' - zbark)';
        end
        Pzzbark = Pzzbark + Gammak * Qk * Gammak';
        Pxxbark = Pzzbark(1:nx,1:nx);
        Pxcbark = Pzzbark(1:nx,nx+1:end);
        Pcxbark = Pzzbark(nx+1:end,1:nx);
        Pccbark = Pzzbark(nx+1:end,nx+1:end);
        % Initial measurement update
        Suk = chol(Pzzbark, 'lower');
        sigmaZubarkHist = zeros(2 * (nz), nz);
        YkHist = zeros(2 * (nz), ny);
        parfor ii = 1 : 2 * (nz)
            jj = ii;
            sign = 1;
            if ii > nz
                jj = ii - nz;
                sign = -1;
            end
            sigmaZk = zbark + sign * sqrt(nz) * Suk(:,jj);
            sigmaZubarkHist(ii,:) = sigmaZk;
            YkHist(ii,:) = h_meas(sigmaZk, t, stat_id, JD0_UTC, Sensor_data, Params);
        end
        ybark = 0;
        for ii = 1 : 2 * (nz)
            ybark = ybark + w * YkHist(ii,:)';
        end
        Pyyk = 0;
        Pzyk = 0;
        for ii = 1 : 2 * (nz)
            Pyyk = Pyyk + w * (YkHist(ii,:)' - ybark) * (YkHist(ii,:)' - ybark)';
            Pzyk = Pzyk + w * (sigmaZubarkHist(ii,:)' - zbark) * (YkHist(ii,:)' - ybark)';
        end
        Pyyk = Pyyk + Rk;
        Kzk = Pzyk * inv(Pyyk);
        Kxk = Kzk(1:nx,:);
        Kck = Kzk(nx+1:end,:);

        if fitCase == 'C' && stat_id ~= 1
            Kxk = zeros(nx, ny);
        elseif fitCase == 'D' && stat_id ~= 2
            Kxk = zeros(nx, ny);
        elseif fitCase == 'E' && stat_id ~= 3
            Kxk = zeros(nx, ny);
        end
    
        % Initial state and covariance update
        yk = Sensor_data.LEO_DATA_Apparent(kk,3:4)' .* 1e3;
        zhatk = zbark + [Kxk; zeros(nc,ny)] * (yk - ybark);
        Pzzk = [Pxxbark, Pxcbark; Pcxbark, Pccbark] - [Kxk * Pyyk * Kxk', Kxk * Pyyk * Kck'; Kck * Pyyk * Kxk', zeros(nc, nc)];
    
        zhatHist(kk-startIdx+1,:) = zhatk;
        PzzHist(kk-startIdx+1,:) = Pzzk(:);
        tMeasHist(kk-startIdx+1) = t;
        statIdHist(kk-startIdx+1) = stat_id;
        residHist(kk-startIdx+1,:) = yk - ybark;
        residCovHist(kk-startIdx+1,:) = Pyyk(:);

    end

    % Propagate to delivery
    xhat_deliv = zeros(1, nz);
    P_deliv = zeros(1, (nz)^2);

    Spkm1 = chol(Pzzk, 'lower');
    sigmaZpbarkHist = zeros(2 * (nz), nz);
    % Initialize Sigma Points and Propagate
    parfor ii = 1 : 2 * (nz)
        jj = ii;
        sign = 1;
        if ii > nz
            jj = ii - nz;
            sign = -1;
        end
        sigmaZk = zhatk + sign * sqrt(nz) * Spkm1(:,jj);
        [~, ZkHist] = ode113(@orbPropHF, tMeasMax:60:tPropMax, sigmaZk, Params.ode_options, JD0_UTC, Params); % ACCOUNT FOR CD IN PROPAGATOR!!!
        sigmaZpbarkHist(ii,:) = ZkHist(end,:);
    end
    % prediction
    zbark = 0;
    for ii = 1 : 2 * (nz)
        zbark = zbark + w * sigmaZpbarkHist(ii,:)';
    end
    R = zbark(1:3) / norm(zbark(1:3));
    C = cross(zbark(1:3), zbark(4:6)) / norm(cross(zbark(1:3), zbark(4:6)));
    I = cross(C, R);
    T_GCRF_RIC = [R, I, C];
    Qk = T_GCRF_RIC * Params.Qk_RIC * T_GCRF_RIC';
    dt = tPropMax - tMeasMax;
    if length(X0) == 6
        Gammak = [dt^2/2 * eye(3); dt * eye(3); zeros(nc, 3)];
    elseif length(X0) > 6
        Gammak = [dt^2/2 * eye(3); dt * eye(3); zeros(length(X0)-6,3); zeros(nc, 3)];
    end
    Pzzbark = 0;
    for ii = 1 : 2 * (nz)
        Pzzbark = Pzzbark + w * (sigmaZpbarkHist(ii,:)' - zbark) * (sigmaZpbarkHist(ii,:)' - zbark)';
    end
    Pzzbark = Pzzbark + Gammak * Qk * Gammak';
    Pxxbark = Pzzbark(1:nx,1:nx);

    xhatk = zbark(1:nx);
    Pk = Pxxbark;

    xhat_deliv = xhatk;
    P_deliv = Pk;

end