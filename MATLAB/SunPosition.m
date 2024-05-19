function r_sun_GCRF = SunPosition(Params)

    %Params = UpdateParams(JD_UTC, Params);

    T_TDB = Params.T_TDB;
    %T_TDB = Params.T_UT1;

    % Vallado Ed. 5
    % lambda_bar_ecl_sun = mod(Params.deg2rad * (280.460 + 36000.771285 * T_TDB), 2*pi);
    % 
    % M_sun = mod(Params.deg2rad * (357.528 + 35999.050957 * T_TDB), 2*pi);
    % 
    % lambda_ecl_sun = mod(lambda_bar_ecl_sun + Params.deg2rad * (1.915 * sin(M_sun) + 0.020 * sin(2 * M_sun)), 2*pi);
    % 
    % r_sun = 1.00014 - 0.01671 * cos(M_sun) - 0.00014 * cos(2 * M_sun);
    % epsilon = mod(Params.deg2rad * (23.439291 - 0.01461 * T_TDB), 2*pi);
    % 
    % r_sun_vec = [r_sun * cos(lambda_ecl_sun); r_sun * cos(epsilon) * sin(lambda_ecl_sun); r_sun * sin(epsilon) * sin(lambda_ecl_sun)];
    % 
    % r_sun_TOD = r_sun_vec .* Params.AU;

    % Vallado Code:
    lambda_bar_ecl_sun = mod(Params.deg2rad * (280.460 + 36000.77 * T_TDB), 2*pi);
    M_sun = mod(Params.deg2rad * (357.5277233 + 35999.05034 * T_TDB), 2*pi);
    lambda_ecl_sun = mod(lambda_bar_ecl_sun + Params.deg2rad * (1.914666471 * sin(M_sun) + 0.019994643 * sin(2 * M_sun)), 2*pi);
    epsilon = mod(Params.deg2rad * (23.439291 - 0.0130042 * T_TDB), 2*pi);
    r_sun = 1.000140612 - 0.016708617 * cos(M_sun) - 0.000139589 * cos(2 * M_sun);
    r_sun_vec = [r_sun * cos(lambda_ecl_sun); r_sun * cos(epsilon) * sin(lambda_ecl_sun); r_sun * sin(epsilon) * sin(lambda_ecl_sun)];
    r_sun_TOD = r_sun_vec .* Params.AU;

    P = PrecessionMatrix(Params);
    N = NutationMatrix(Params);

    r_sun_GCRF = P * N * r_sun_TOD;

end