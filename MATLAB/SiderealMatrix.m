function R = SiderealMatrix(Params)
% PEF to TOD

    T_TT = Params.T_TT;
    T_UT1 = Params.T_UT1;

    eps_bar_1980 = Params.deg2rad * (23.439291 - 0.0130042 * T_TT - 1.64e-7 * T_TT^2 + 5.04e-7 * T_TT^3);
    r = 360; % deg
    M_moon = mod(Params.deg2rad * (134.96298139 + (1325 * r + 198.8673981) * T_TT + 0.0086972 * T_TT^2 + 1.78e-5 * T_TT^3), 2*pi);
    M_sun = mod(Params.deg2rad * (357.52772333 + (99 * r + 359.0503400) * T_TT - 0.0001603 * T_TT^2 - 3.3e-6 * T_TT^3), 2*pi);
    u_bar_moon = mod(Params.deg2rad * (93.27191028 + (1342 * r + 82.0175381) * T_TT - 0.0036825 * T_TT^2 + 3.1e-6 * T_TT^3), 2*pi);
    D_sun = mod(Params.deg2rad * (297.85036306 + (1236 * r + 307.1114800) * T_TT - 0.0019142 * T_TT^2 + 5.3e-6 * T_TT^3), 2*pi);
    lambda_bar_ecl_moon = mod(Params.deg2rad * (125.04452222 - (5 * r + 134.1362608) * T_TT + 0.0020708 * T_TT^2 + 2.2e-6 * T_TT^3), 2*pi);
    
    a_n1 = Params.nut_data(:,1);
    a_n2 = Params.nut_data(:,2);
    a_n3 = Params.nut_data(:,3);
    a_n4 = Params.nut_data(:,4);
    a_n5 = Params.nut_data(:,5);
    A = Params.nut_data(:,6) .* 0.0001; % arcsec / T_TT
    B = Params.nut_data(:,7) .* 0.0001; % arcsec / T_TT
    
    a_p = a_n1 .* M_moon + a_n2 .* M_sun + a_n3 .* u_bar_moon + a_n4 .* D_sun + a_n5 .* lambda_bar_ecl_moon;
    
    delta_psi_1980 = sum(Params.arcsec2rad * ((A + B .* T_TT) .* sin(a_p)));
    delta_psi_1980_corr = delta_psi_1980 + Params.d_delta_psi_1980;
    
    
    Eq_equinox_1982 = delta_psi_1980_corr * cos(eps_bar_1980) + Params.arcsec2rad * (0.00264 * sin(lambda_bar_ecl_moon) + 0.000063 * sin(2 * lambda_bar_ecl_moon));
    theta_GMST_1982 = mod(Params.arcsec2rad * ((67310.54841 + (876600 * 3600 + 8640184.812866) * T_UT1 + 0.093104 * T_UT1^2 - 6.2e-6 * T_UT1^3) * 15), 2*pi);
    theta_GAST_1982 = theta_GMST_1982 + Eq_equinox_1982;
    
    %R = Rot3(-theta_GAST_1982);

    cos_neg_theta_GAST_1982 = cos(-theta_GAST_1982);
    sin_neg_theta_GAST_1982 = sin(-theta_GAST_1982);

    % R = [ cos(-theta_GAST_1982), sin(-theta_GAST_1982), 0.0; ...
    %      -sin(-theta_GAST_1982), cos(-theta_GAST_1982), 0.0; ...
    %            0.0,      0.0, 1.0];

    R = [cos_neg_theta_GAST_1982, sin_neg_theta_GAST_1982, 0.0;
        -sin_neg_theta_GAST_1982, cos_neg_theta_GAST_1982, 0.0;
         0.0, 0.0, 1.0];

end