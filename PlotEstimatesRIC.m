function PlotEstimatesRIC(pos_long_RIC, pos_short_RIC, pos_range_RIC, pos_rangerate_RIC, pos_kwaj_RIC, pos_dg_RIC, pos_arecibo_RIC,...
                          posCov_long_RIC, posCov_short_RIC, posCov_range_RIC, posCov_rangerate_RIC, posCov_kwaj_RIC, posCov_dg_RIC, posCov_arecibo_RIC)

    % Plot R-C
    [ellipR_long, ellipC_long] = CovarianceEllipse([posCov_long_RIC(1,1) posCov_long_RIC(1,3); posCov_long_RIC(3,1), posCov_long_RIC(3,3)]);
    [ellipR_short, ellipC_short] = CovarianceEllipse([posCov_short_RIC(1,1) posCov_short_RIC(1,3); posCov_short_RIC(3,1), posCov_short_RIC(3,3)]);
    [ellipR_range, ellipC_range] = CovarianceEllipse([posCov_range_RIC(1,1) posCov_range_RIC(1,3); posCov_range_RIC(3,1), posCov_range_RIC(3,3)]);
    [ellipR_rangerate, ellipC_rangerate] = CovarianceEllipse([posCov_rangerate_RIC(1,1) posCov_rangerate_RIC(1,3); posCov_rangerate_RIC(3,1), posCov_rangerate_RIC(3,3)]);
    [ellipR_kwaj, ellipC_kwaj] = CovarianceEllipse([posCov_kwaj_RIC(1,1) posCov_kwaj_RIC(1,3); posCov_kwaj_RIC(3,1), posCov_kwaj_RIC(3,3)]);
    [ellipR_dg, ellipC_dg] = CovarianceEllipse([posCov_dg_RIC(1,1) posCov_dg_RIC(1,3); posCov_dg_RIC(3,1), posCov_dg_RIC(3,3)]);
    [ellipR_arecibo, ellipC_arecibo] = CovarianceEllipse([posCov_arecibo_RIC(1,1) posCov_arecibo_RIC(1,3); posCov_arecibo_RIC(3,1), posCov_arecibo_RIC(3,3)]);
    figure
    hold on
    plot(pos_long_RIC(1)-pos_long_RIC(1), pos_long_RIC(3)-pos_long_RIC(3), '.', 'Color', '#0072BD', 'MarkerSize', 15, 'DisplayName', 'Case F')
    plot(pos_long_RIC(1)-pos_long_RIC(1)+ellipR_long, pos_long_RIC(3)-pos_long_RIC(3)+ellipC_long, '--', 'LineWidth', 1.5, 'Color', '#0072BD', 'DisplayName', 'Case F')

    plot(pos_short_RIC(1)-pos_long_RIC(1), pos_short_RIC(3)-pos_long_RIC(3), '.', 'Color', '#A2142F', 'MarkerSize', 15, 'DisplayName', 'Case G')
    plot(pos_short_RIC(1)-pos_long_RIC(1)+ellipR_short, pos_short_RIC(3)-pos_long_RIC(3)+ellipC_short, '--', 'LineWidth', 1.5, 'Color', '#A2142F', 'DisplayName', 'Case G')
    
    plot(pos_range_RIC(1)-pos_long_RIC(1), pos_range_RIC(3)-pos_range_RIC(3), '.', 'Color', '#D95319', 'MarkerSize', 15, 'DisplayName', 'Case A')
    plot(pos_range_RIC(1)-pos_long_RIC(1)+ellipR_range, pos_range_RIC(3)-pos_range_RIC(3)+ellipC_range, '--', 'LineWidth', 1.5, 'Color', '#D95319', 'DisplayName', 'Case A')
    
    plot(pos_rangerate_RIC(1)-pos_long_RIC(1), pos_rangerate_RIC(3)-pos_long_RIC(3), '.', 'Color', '#EDB120', 'MarkerSize', 15, 'DisplayName', 'Case B')
    plot(pos_rangerate_RIC(1)-pos_long_RIC(1)+ellipR_rangerate, pos_rangerate_RIC(3)-pos_long_RIC(3)+ellipC_rangerate, '--', 'LineWidth', 1.5, 'Color', '#EDB120', 'DisplayName', 'Case B')
    
    plot(pos_kwaj_RIC(1)-pos_long_RIC(1), pos_kwaj_RIC(3)-pos_long_RIC(3), '.', 'Color', '#7E2F8E', 'MarkerSize', 15, 'DisplayName', 'Case C')
    plot(pos_kwaj_RIC(1)-pos_long_RIC(1)+ellipR_kwaj, pos_kwaj_RIC(3)-pos_long_RIC(3)+ellipC_kwaj, '--', 'LineWidth', 1.5, 'Color', '#7E2F8E', 'DisplayName', 'Case C')
    
    plot(pos_dg_RIC(1)-pos_long_RIC(1), pos_dg_RIC(3)-pos_long_RIC(3), '.', 'Color', '#77AC30', 'MarkerSize', 15, 'DisplayName', 'Case D')
    plot(pos_dg_RIC(1)-pos_long_RIC(1)+ellipR_dg, pos_dg_RIC(3)-pos_long_RIC(3)+ellipC_dg, '--', 'LineWidth', 1.5, 'Color', '#77AC30', 'DisplayName', 'Case D')
    
    plot(pos_arecibo_RIC(1)-pos_long_RIC(1), pos_arecibo_RIC(3)-pos_long_RIC(3), '.', 'Color', '#4DBEEE', 'MarkerSize', 15, 'DisplayName', 'Case E')
    plot(pos_arecibo_RIC(1)-pos_long_RIC(1)+ellipR_arecibo, pos_arecibo_RIC(3)-pos_long_RIC(3)+ellipC_arecibo, '--', 'LineWidth', 1.5, 'Color', '#4DBEEE', 'DisplayName', 'Case E')
    hold off
    title('Position Estimate Error Radial-Cross-track Plane')
    xlabel('Radial Error [km]')
    ylabel('Cross-track Error [km]')
    grid on
    legend('Location', 'northeastoutside')
    
    % Plot R-I
    [ellipR_long, ellipI_long] = CovarianceEllipse([posCov_long_RIC(1,1) posCov_long_RIC(1,2); posCov_long_RIC(2,1), posCov_long_RIC(2,2)]);
    [ellipR_short, ellipI_short] = CovarianceEllipse([posCov_short_RIC(1,1) posCov_short_RIC(1,2); posCov_short_RIC(2,1), posCov_short_RIC(2,2)]);
    [ellipR_range, ellipI_range] = CovarianceEllipse([posCov_range_RIC(1,1) posCov_range_RIC(1,2); posCov_range_RIC(2,1), posCov_range_RIC(2,2)]);
    [ellipR_rangerate, ellipI_rangerate] = CovarianceEllipse([posCov_rangerate_RIC(1,1) posCov_rangerate_RIC(1,2); posCov_rangerate_RIC(2,1), posCov_rangerate_RIC(2,2)]);
    [ellipR_kwaj, ellipI_kwaj] = CovarianceEllipse([posCov_kwaj_RIC(1,1) posCov_kwaj_RIC(1,2); posCov_kwaj_RIC(2,1), posCov_kwaj_RIC(2,2)]);
    [ellipR_dg, ellipI_dg] = CovarianceEllipse([posCov_dg_RIC(1,1) posCov_dg_RIC(1,2); posCov_dg_RIC(2,1), posCov_dg_RIC(2,2)]);
    [ellipR_arecibo, ellipI_arecibo] = CovarianceEllipse([posCov_arecibo_RIC(1,1) posCov_arecibo_RIC(1,2); posCov_arecibo_RIC(2,1), posCov_arecibo_RIC(2,2)]);
    figure
    hold on
    plot(pos_long_RIC(1)-pos_long_RIC(1), pos_long_RIC(2)-pos_long_RIC(2), '.', 'Color', '#0072BD', 'MarkerSize', 15, 'DisplayName', 'Case F')
    plot(pos_long_RIC(1)-pos_long_RIC(1)+ellipR_long, pos_long_RIC(2)-pos_long_RIC(2)+ellipI_long, '--', 'LineWidth', 1.5, 'Color', '#0072BD', 'DisplayName', 'Case F')
    
    plot(pos_short_RIC(1)-pos_long_RIC(1), pos_short_RIC(2)-pos_long_RIC(2), '.', 'Color', '#A2142F', 'MarkerSize', 15, 'DisplayName', 'Case G')
    plot(pos_short_RIC(1)-pos_long_RIC(1)+ellipR_short, pos_short_RIC(2)-pos_long_RIC(2)+ellipI_short, '--', 'LineWidth', 1.5, 'Color', '#A2142F', 'DisplayName', 'Case G')

    plot(pos_range_RIC(1)-pos_long_RIC(1), pos_range_RIC(2)-pos_long_RIC(2), '.', 'Color', '#D95319', 'MarkerSize', 15, 'DisplayName', 'Case A')
    plot(pos_range_RIC(1)-pos_long_RIC(1)+ellipR_range, pos_range_RIC(2)-pos_long_RIC(2)+ellipI_range, '--', 'LineWidth', 1.5, 'Color', '#D95319', 'DisplayName', 'Case A')
    
    plot(pos_rangerate_RIC(1)-pos_long_RIC(1), pos_rangerate_RIC(2)-pos_long_RIC(2), '.', 'Color', '#EDB120', 'MarkerSize', 15, 'DisplayName', 'Case B')
    plot(pos_rangerate_RIC(1)-pos_long_RIC(1)+ellipR_rangerate, pos_rangerate_RIC(2)-pos_long_RIC(2)+ellipI_rangerate, '--', 'LineWidth', 1.5, 'Color', '#EDB120', 'DisplayName', 'Case B')
    
    plot(pos_kwaj_RIC(1)-pos_long_RIC(1), pos_kwaj_RIC(2)-pos_long_RIC(2), '.', 'Color', '#7E2F8E', 'MarkerSize', 15, 'DisplayName', 'Case C')
    plot(pos_kwaj_RIC(1)-pos_long_RIC(1)+ellipR_kwaj, pos_kwaj_RIC(2)-pos_long_RIC(2)+ellipI_kwaj, '--', 'LineWidth', 1.5, 'Color', '#7E2F8E', 'DisplayName', 'Case C')
    
    plot(pos_dg_RIC(1)-pos_long_RIC(1), pos_dg_RIC(2)-pos_long_RIC(2), '.', 'Color', '#77AC30', 'MarkerSize', 15, 'DisplayName', 'Case D')
    plot(pos_dg_RIC(1)-pos_long_RIC(1)+ellipR_dg, pos_dg_RIC(2)-pos_long_RIC(2)+ellipI_dg, '--', 'LineWidth', 1.5, 'Color', '#77AC30', 'DisplayName', 'Case D')
    
    plot(pos_arecibo_RIC(1)-pos_long_RIC(1), pos_arecibo_RIC(2)-pos_long_RIC(2), '.', 'Color', '#4DBEEE', 'MarkerSize', 15, 'DisplayName', 'Case E')
    plot(pos_arecibo_RIC(1)-pos_long_RIC(1)+ellipR_arecibo, pos_arecibo_RIC(2)-pos_long_RIC(2)+ellipI_arecibo, '--', 'LineWidth', 1.5, 'Color', '#4DBEEE', 'DisplayName', 'Case E')
    hold off
    title('Position Estimate Error Radial-In-track Plane')
    xlabel('Radial Error [km]')
    ylabel('In-track Error [km]')
    grid on
    legend('Location', 'northeastoutside')
    
    % Plot C-I
    [ellipI_long, ellipC_long] = CovarianceEllipse([posCov_long_RIC(2,2) posCov_long_RIC(2,3); posCov_long_RIC(3,2), posCov_long_RIC(3,3)]);
    [ellipI_short, ellipC_short] = CovarianceEllipse([posCov_short_RIC(2,2) posCov_short_RIC(2,3); posCov_short_RIC(3,2), posCov_short_RIC(3,3)]);
    [ellipI_range, ellipC_range] = CovarianceEllipse([posCov_range_RIC(2,2) posCov_range_RIC(2,3); posCov_range_RIC(3,2), posCov_range_RIC(3,3)]);
    [ellipI_rangerate, ellipC_rangerate] = CovarianceEllipse([posCov_rangerate_RIC(2,2) posCov_rangerate_RIC(2,3); posCov_rangerate_RIC(3,2), posCov_rangerate_RIC(3,3)]);
    [ellipI_kwaj, ellipC_kwaj] = CovarianceEllipse([posCov_kwaj_RIC(2,2) posCov_kwaj_RIC(2,3); posCov_kwaj_RIC(3,2), posCov_kwaj_RIC(3,3)]);
    [ellipI_dg, ellipC_dg] = CovarianceEllipse([posCov_dg_RIC(2,2) posCov_dg_RIC(2,3); posCov_dg_RIC(3,2), posCov_dg_RIC(3,3)]);
    [ellipI_arecibo, ellipC_arecibo] = CovarianceEllipse([posCov_arecibo_RIC(2,2) posCov_arecibo_RIC(2,3); posCov_arecibo_RIC(3,2), posCov_arecibo_RIC(3,3)]);
    figure
    hold on
    plot(pos_long_RIC(3)-pos_long_RIC(3), pos_long_RIC(2)-pos_long_RIC(2), '.', 'Color', '#0072BD', 'MarkerSize', 15, 'DisplayName', 'Case F')
    plot(pos_long_RIC(3)-pos_long_RIC(3)+ellipC_long, pos_long_RIC(2)-pos_long_RIC(2)+ellipI_long, '--', 'LineWidth', 1.5, 'Color', '#0072BD', 'DisplayName', 'Case F')
    
    plot(pos_short_RIC(3)-pos_long_RIC(3), pos_short_RIC(2)-pos_long_RIC(2), '.', 'Color', '#A2142F', 'MarkerSize', 15, 'DisplayName', 'Case G')
    plot(pos_short_RIC(3)-pos_long_RIC(3)+ellipC_short, pos_short_RIC(2)-pos_long_RIC(2)+ellipI_short, '--', 'LineWidth', 1.5, 'Color', '#A2142F', 'DisplayName', 'Case G')

    plot(pos_range_RIC(3)-pos_long_RIC(3), pos_range_RIC(2)-pos_long_RIC(2), '.', 'Color', '#D95319', 'MarkerSize', 15, 'DisplayName', 'Case A')
    plot(pos_range_RIC(3)-pos_long_RIC(3)+ellipC_range, pos_range_RIC(2)-pos_long_RIC(2)+ellipI_range, '--', 'LineWidth', 1.5, 'Color', '#D95319', 'DisplayName', 'Case A')
    
    plot(pos_rangerate_RIC(3)-pos_long_RIC(3), pos_rangerate_RIC(2)-pos_long_RIC(2), '.', 'Color', '#EDB120', 'MarkerSize', 15, 'DisplayName', 'Case B')
    plot(pos_rangerate_RIC(3)-pos_long_RIC(3)+ellipC_rangerate, pos_rangerate_RIC(2)-pos_long_RIC(2)+ellipI_rangerate, '--', 'LineWidth', 1.5, 'Color', '#EDB120', 'DisplayName', 'Case B')
    
    plot(pos_kwaj_RIC(3)-pos_long_RIC(3), pos_kwaj_RIC(2)-pos_long_RIC(2), '.', 'Color', '#7E2F8E', 'MarkerSize', 15, 'DisplayName', 'Case C')
    plot(pos_kwaj_RIC(3)-pos_long_RIC(3)+ellipC_kwaj, pos_kwaj_RIC(2)-pos_long_RIC(2)+ellipI_kwaj, '--', 'LineWidth', 1.5, 'Color', '#7E2F8E', 'DisplayName', 'Case C')
    
    plot(pos_dg_RIC(3)-pos_long_RIC(3), pos_dg_RIC(2)-pos_long_RIC(2), '.', 'Color', '#77AC30', 'MarkerSize', 15, 'DisplayName', 'Case D')
    plot(pos_dg_RIC(3)-pos_long_RIC(3)+ellipC_dg, pos_dg_RIC(2)-pos_long_RIC(2)+ellipI_dg, '--', 'LineWidth', 1.5, 'Color', '#77AC30', 'DisplayName', 'Case D')
    
    plot(pos_arecibo_RIC(3)-pos_long_RIC(3), pos_arecibo_RIC(2)-pos_long_RIC(2), '.', 'Color', '#4DBEEE', 'MarkerSize', 15, 'DisplayName', 'Case E')
    plot(pos_arecibo_RIC(3)-pos_long_RIC(3)+ellipC_arecibo, pos_arecibo_RIC(2)-pos_long_RIC(2)+ellipI_arecibo, '--', 'LineWidth', 1.5, 'Color', '#4DBEEE', 'DisplayName', 'Case E')
    
    hold off
    title('Position Estimate Error Cross-track-In-track Plane')
    xlabel('Cross-track Error [km]')
    ylabel('In-track Error [km]')
    grid on
    legend('Location', 'northeastoutside')

end