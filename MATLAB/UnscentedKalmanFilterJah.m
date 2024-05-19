function [xhatHist, PHist, xbarHist, PbarHist, sigmaXpbarHist, tMeasHist, statIdHist, residHist, residCovHist, xhat_deliv, P_deliv] = UnscentedKalmanFilterJah...
         (X0, P0, tMeasMax, tPropMax, JD0_UTC, Sensor_data, fitCase, Params)
    
    dataMaxIdx = find(Sensor_data.LEO_DATA_Apparent(:,2) == tMeasMax);

    nx = length(X0);
    ny = 2;
    w = 1 / (2 * nx);

    if fitCase == 'G'
        startIdx = 1290;
    else
        startIdx = 1;
    end
    
    xhatHist = zeros(dataMaxIdx-startIdx+1, nx);
    PHist = zeros(dataMaxIdx-startIdx+1, nx^2);
    xbarHist = zeros(dataMaxIdx-startIdx+1, nx);
    PbarHist = zeros(dataMaxIdx-startIdx+1, nx^2);
    sigmaXpbarHist = cell(dataMaxIdx-startIdx+1, 1);
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
    
    Spkm1 = chol(P0, 'lower');
    sigmaXpbarkHist = zeros(2 * nx, nx);
    % Initialize Sigma Points and Propagate
    for ii = 1 : 2 * nx 
        jj = ii;
        sign = 1;
        if ii > nx
            jj = ii - nx;
            sign = -1;
        end
        sigmaXk = X0 + sign * sqrt(nx) * Spkm1(:,jj);
        sigmaXpbarkHist(ii,:) = sigmaXk;
    end
    % Initial prediction
    xbark = 0;
    for ii = 1 : 2 * nx
        xbark = xbark + w * sigmaXpbarkHist(ii,:)';
    end
    R = X0(1:3) / norm(X0(1:3));
    C = cross(X0(1:3), X0(4:6)) / norm(cross(X0(1:3), X0(4:6)));
    I = cross(C, R);
    T_GCRF_RIC = [R, I, C];
    Qk = T_GCRF_RIC * Params.Qk_RIC * T_GCRF_RIC';
    dt = 0;
    if length(X0) == 6
        Gammak = [dt^2/2 * eye(3); dt * eye(3)];
    elseif length(X0) > 6
        Gammak = [dt^2/2 * eye(3); dt * eye(3); zeros(length(X0)-6,3)];
    end
    Pbark = 0;
    for ii = 1 : 2 * nx
        Pbark = Pbark + w * (sigmaXpbarkHist(ii,:)' - xbark) * (sigmaXpbarkHist(ii,:)' - xbark)';
    end
    Pbark = Pbark + Gammak * Qk * Gammak';
    % Initial measurement update
    Suk = chol(Pbark, 'lower');
    sigmaXubarkHist = zeros(2 * nx, nx);
    YkHist = zeros(2 * nx, ny);
    for ii = 1 : 2 * nx
        jj = ii;
        sign = 1;
        if ii > nx
            jj = ii - nx;
            sign = -1;
        end
        sigmaXk = xbark + sign * sqrt(nx) * Suk(:,jj);
        sigmaXubarkHist(ii,:) = sigmaXk;
        YkHist(ii,:) = h_meas(sigmaXk, t, stat_id, JD0_UTC, Sensor_data, Params);
    end
    ybark = 0;
    for ii = 1 : 2 * nx
        ybark = ybark + w * YkHist(ii,:)';
    end
    Pyyk = 0;
    Pxyk = 0;
    for ii = 1 : 2 * nx
        Pyyk = Pyyk + w * (YkHist(ii,:)' - ybark) * (YkHist(ii,:)' - ybark)';
        Pxyk = Pxyk + w * (sigmaXubarkHist(ii,:)' - xbark) * (YkHist(ii,:)' - ybark)';
    end
    Pyyk = Pyyk + Rk;
    Kk = Pxyk / Pyyk;
    if fitCase == 'C' && stat_id ~= 1
        Kk = zeros(nx, ny);
    elseif fitCase == 'D' && stat_id ~= 2
        Kk = zeros(nx, ny);
    elseif fitCase == 'E' && stat_id ~= 3
        Kk = zeros(nx, ny);
    end
    Kk(7,:) = zeros(1,ny); % consider Cd
    % Initial state and covariance update
    yk = Sensor_data.LEO_DATA_Apparent(startIdx,3:4)' .* 1e3;
    xhatk = xbark + Kk * (yk - ybark);
    %Pk = Pbark - Kk * Pyyk * Kk';
    Pk = Pbark - Pxyk * Kk' - Kk * Pxyk' + Kk * Pyyk * Kk';

    xhatHist(1,:) = xhatk;
    PHist(1,:) = Pk(:);
    xbarHist(1,:) = xbark;
    PbarHist(1,:) = Pbark(:);
    sigmaXpbarHist{1} = sigmaXpbarkHist;
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

        Spkm1 = chol(Pk, 'lower');
        sigmaXpbarkHist = zeros(2 * nx, nx);
        % Initialize Sigma Points and Propagate
        parfor ii = 1 : 2 * nx
            jj = ii;
            sign = 1;
            if ii > nx
                jj = ii - nx;
                sign = -1;
            end
            sigmaXk = xhatk + sign * sqrt(nx) * Spkm1(:,jj);
            [~, XkHist] = ode113(@orbPropHF, [tprev; t], sigmaXk, Params.ode_options, JD0_UTC, Params);
            sigmaXpbarkHist(ii,:) = XkHist(end,:);
        end
        % prediction
        xbark = 0;
        for ii = 1 : 2 * nx
            xbark = xbark + w * sigmaXpbarkHist(ii,:)';
        end
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
        Pbark = 0;
        for ii = 1 : 2 * nx
            Pbark = Pbark + w * (sigmaXpbarkHist(ii,:)' - xbark) * (sigmaXpbarkHist(ii,:)' - xbark)';
        end
        Pbark = Pbark + Gammak * Qk * Gammak';
        % measurement update
        Suk = chol(Pbark, "lower");
        sigmaXubarkHist = zeros(2 * nx, nx);
        YkHist = zeros(2 * nx, ny);
        parfor ii = 1 : 2 * nx
            jj = ii;
            sign = 1;
            if ii > nx
                jj = ii - nx;
                sign = -1;
            end
            sigmaXk = xbark + sign * sqrt(nx) * Suk(:,jj);
            sigmaXubarkHist(ii,:) = sigmaXk;
            YkHist(ii,:) = h_meas(sigmaXk, t, stat_id, JD0_UTC, Sensor_data, Params);
        end
        ybark = 0;
        for ii = 1 : 2 * nx
            ybark = ybark + w * YkHist(ii,:)';
        end
        Pyyk = 0;
        Pxyk = 0;
        for ii = 1 : 2 * nx
            Pyyk = Pyyk + w * (YkHist(ii,:)' - ybark) * (YkHist(ii,:)' - ybark)';
            Pxyk = Pxyk + w * (sigmaXubarkHist(ii,:)' - xbark) * (YkHist(ii,:)' - ybark)';
        end
        Pyyk = Pyyk + Rk;
        Kk = Pxyk / Pyyk;
        if fitCase == 'C' && stat_id ~= 1
            Kk = zeros(nx, ny);
        elseif fitCase == 'D' && stat_id ~= 2
            Kk = zeros(nx, ny);
        elseif fitCase == 'E' && stat_id ~= 3
            Kk = zeros(nx, ny);
        end
        Kk(7,:) = zeros(1,ny); % consider Cd
        % state and covariance update
        yk = Sensor_data.LEO_DATA_Apparent(kk,3:4)' .* 1e3;
        xhatk = xbark + Kk * (yk - ybark);
        %Pk = Pbark - Kk * Pyyk * Kk';
        Pk = Pbark - Pxyk * Kk' - Kk * Pxyk' + Kk * Pyyk * Kk';
    
        xhatHist(kk-startIdx+1,:) = xhatk;
        PHist(kk-startIdx+1,:) = Pk(:);
        xbarHist(kk-startIdx+1,:) = xbark;
        PbarHist(kk-startIdx+1,:) = Pbark(:);
        sigmaXpbarHist{kk-startIdx+1} = sigmaXpbarkHist;
        tMeasHist(kk-startIdx+1) = t;
        statIdHist(kk-startIdx+1) = stat_id;
        residHist(kk-startIdx+1,:) = yk - ybark;
        residCovHist(kk-startIdx+1,:) = Pyyk(:);

    end

    % Propagate to delivery
    xhat_deliv = zeros(nx, 1);
    P_deliv = zeros(nx); 

    Spkm1 = chol(Pk, 'lower');
    sigmaXpbarkHist = zeros(2 * nx, nx);
    parfor ii = 1 : 2 * nx
        jj = ii;
        sign = 1;
        if ii > nx
            jj = ii - nx;
            sign = -1;
        end
        sigmaXk = xhatk + sign * sqrt(nx) * Spkm1(:,jj);
        [~, XkHist] = ode113(@orbPropHF, tMeasMax:60:tPropMax, sigmaXk, Params.ode_options, JD0_UTC, Params);
        sigmaXpbarkHist(ii,:) = XkHist(end,:);
    end
    % prediction
    xbark = 0;
    for ii = 1 : 2 * nx
        xbark = xbark + w * sigmaXpbarkHist(ii,:)';
    end
    R = xbark(1:3) / norm(xbark(1:3));
    C = cross(xbark(1:3), xbark(4:6)) / norm(cross(xbark(1:3), xbark(4:6)));
    I = cross(C, R);
    T_GCRF_RIC = [R, I, C];
    Qk = T_GCRF_RIC * Params.Qk_RIC * T_GCRF_RIC';
    dt = tPropMax - tMeasMax;
    if length(X0) == 6
        Gammak = [dt^2/2 * eye(3); dt * eye(3)];
    elseif length(X0) > 6
        Gammak = [dt^2/2 * eye(3); dt * eye(3); zeros(length(X0)-6,3)];
    end
    Pbark = 0;
    for ii = 1 : 2 * nx
        Pbark = Pbark + w * (sigmaXpbarkHist(ii,:)' - xbark) * (sigmaXpbarkHist(ii,:)' - xbark)';
    end
    Pbark = Pbark + Gammak * Qk * Gammak';

    xhatk = xbark;
    Pk = Pbark;

    xhat_deliv = xhatk;
    P_deliv = Pk;

end