function PlotPostFitResiduals(tMeasHist, residHist, residCovHist, statIdHist, fitCase)

    RMS_range = sqrt(sum(residHist(:,1).^2) / length(residHist)) / 1e3;
    RMS_rangeRate = sqrt(sum(residHist(:,2).^2) / length(residHist)) / 1e3;
    figure
    sgtitle(['Post-fit Measurement Residuals: Case ' fitCase])
    stat_id = unique(statIdHist);
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
        idx = statIdHist == stat_id(ii);
        subplot(2,1,1)
        plot(tMeasHist(idx)./3600, residHist(idx,1)./1e3, '.', 'MarkerSize', 12, 'DisplayName', station)
        subplot(2,1,2)
        plot(tMeasHist(idx)./3600, residHist(idx,2)./1e3, '.', 'MarkerSize', 12, 'DisplayName', station)
    end
    subplot(2,1,1)
    plot(tMeasHist./3600, sqrt(residCovHist(:,1)).*3./1e3, '--k', 'LineWidth', 1, 'DisplayName', '+3\sigma')
    plot(tMeasHist./3600, sqrt(residCovHist(:,1)).*-3./1e3, '--k', 'LineWidth', 1, 'DisplayName', '-3\sigma')
    hold off
    title(['RMS = ' num2str(RMS_range) ' km'])
    xlabel('Time since epoch [hr]')
    ylabel('Range [km]')
    grid on
    legend('Location','northeastoutside')
    subplot(2,1,2)
    plot(tMeasHist./3600, sqrt(residCovHist(:,4)).*3./1e3, '--k', 'LineWidth', 1, 'DisplayName', '+3\sigma')
    plot(tMeasHist./3600, sqrt(residCovHist(:,4)).*-3./1e3, '--k', 'LineWidth', 1, 'DisplayName', '-3\sigma')
    hold off
    title(['RMS = ' num2str(RMS_rangeRate) ' km/s'])
    xlabel('Time since epoch [hr]')
    ylabel('Range-Rate [km/s]')
    grid on
    legend('Location', 'northeastoutside')

end