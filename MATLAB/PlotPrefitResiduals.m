function PlotPrefitResiduals(XHist, tHist, tmax, JD0_UTC, Params)

    dataMaxIdx = find(Params.Sensor_data.LEO_DATA_Apparent(:,2) == tmax);

    residuals = zeros(dataMaxIdx, 4);

    for ii = 1 : dataMaxIdx

        t = Params.Sensor_data.LEO_DATA_Apparent(ii, 2);
        stat_id = Params.Sensor_data.LEO_DATA_Apparent(ii, 1);

        if stat_id == 1
            r_stat_ITRF = Params.Kwaj.r_ITRF;
        elseif stat_id == 2
            r_stat_ITRF = Params.DG.r_ITRF;
        elseif stat_id == 3
            r_stat_ITRF = Params.Arecibo.r_ITRF;
        end

        X = XHist(tHist == t, :)';
        JD_UTC = JD0_UTC + t / 86400;
        X_lt_corr = LightTimeCorrection(X, t, JD_UTC, Params);

        Params = UpdateParams(JD_UTC, Params);
        W = PolarMatrix(Params);
        R = SiderealMatrix(Params);
        N = NutationMatrix(Params);
        P = PrecessionMatrix(Params);

        r_stat_GCRF = P * N * R * W * r_stat_ITRF;
        v_stat_GCRF = P * N * R * (W * zeros(3,1) + cross([0; 0; Params.w_earth], W * r_stat_ITRF));
        
        r_GCRF = X_lt_corr(1:3);
        v_GCRF = X_lt_corr(4:6);

        range = sqrt((r_GCRF(1) - r_stat_GCRF(1))^2 + (r_GCRF(2) - r_stat_GCRF(2))^2 + (r_GCRF(3) - r_stat_GCRF(3))^2);
        rangeRate = ((r_GCRF(1) - r_stat_GCRF(1)) * (v_GCRF(1) - v_stat_GCRF(1)) ...
                   + (r_GCRF(2) - r_stat_GCRF(2)) * (v_GCRF(2) - v_stat_GCRF(2)) ...
                   + (r_GCRF(3) - r_stat_GCRF(3)) * (v_GCRF(3) - v_stat_GCRF(3))) / range;

        residuals(ii,1) = t;
        residuals(ii,2) = stat_id;
        residuals(ii,3) = Params.Sensor_data.LEO_DATA_Apparent(ii,3) * 1e3 - range;
        residuals(ii,4) = Params.Sensor_data.LEO_DATA_Apparent(ii,4) * 1e3 - rangeRate;

    end

    RMS_range = sqrt(sum(residuals(:,3).^2) / length(residuals)) / 1e3;
    RMS_rangeRate = sqrt(sum(residuals(:,4).^2) / length(residuals)) / 1e3;
    
    figure
    sgtitle('Pre-fit Measurement Residuals')
    stat_id = unique(residuals(:,2));
    subplot(2,1,1)
    hold on
    subplot(2,1,2)
    hold on
    for ii = 1 : length(stat_id)
        if stat_id(ii) == 1
            station = 'Kwajalein';
        elseif stat_id(ii) == 2
            station = 'Diego Garcia';
        elseif stat_id(ii) == 3
            station = 'Arecibo';
        end
        idx = residuals(:,2) == stat_id(ii);
        subplot(2,1,1)
        plot(residuals(idx,1)./3600, residuals(idx,3)./1e3, '.', 'MarkerSize', 12, 'DisplayName', station)
        subplot(2,1,2)
        plot(residuals(idx,1)./3600, residuals(idx,4)./1e3, '.', 'MarkerSize', 12, 'DisplayName', station)
    end
    subplot(2,1,1)
    title(['RMS = ' num2str(RMS_range) ' km'])
    xlabel('Time since epoch [hr]')
    ylabel('Range [km]')
    grid on
    legend('Location', 'northwest')
    hold off
    subplot(2,1,2)
    title(['RMS = ' num2str(RMS_rangeRate) ' km/s'])
    xlabel('Time since epoch [hr]')
    ylabel('Range-Rate [km/s]')
    grid on
    legend('Location', 'northwest')
    hold off

end