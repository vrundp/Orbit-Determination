function Params = UpdateParams(JD_UTC, Params)
% Initialize Params at top level!

    % date = JD2cal(JD_UTC);
    % year = date(1);
    % month = date(2);
    % day = date(3);
    
    % dataRow_IERS = Params.EOP_IERS_data(Params.EOP_IERS_data(:,2) == year & Params.EOP_IERS_data(:,3) == month & Params.EOP_IERS_data(:,4) == day, :);
    % Params.delta_UT1 = dataRow_IERS(15); % 15 32
    % Params.LOD = dataRow_IERS(17) / 1e3;
    % Params.xp = Params.arcsec2rad * (dataRow_IERS(6)); % 6 29
    % Params.yp = Params.arcsec2rad * (dataRow_IERS(8)); % 8 30
    % Params.d_delta_psi_1980 = Params.arcsec2rad * (dataRow_IERS(20) / 1e3); % 20 34
    % Params.d_delta_eps_1980 = Params.arcsec2rad * (dataRow_IERS(22) / 1e3); % 22 35
    
    % Interpolating -------------------------------------------------------
    % dataRow0_IERS = Params.EOP_IERS_data(Params.EOP_IERS_data(:,2) == year & Params.EOP_IERS_data(:,3) == month & Params.EOP_IERS_data(:,4) == day, :);
    % dataRow1_IERS = Params.EOP_IERS_data(Params.EOP_IERS_data(:,2) == year & Params.EOP_IERS_data(:,3) == month & Params.EOP_IERS_data(:,4) == day + 1, :);

    MJD_UTC = JD_UTC - 2400000.5; % start from midnight for fraction of day
    x0 = fix(MJD_UTC); 
    x1 = x0 + 1;

    dataRow0_IERS = Params.EOP_IERS_data(Params.EOP_IERS_data(:,1) == x0, :);
    dataRow1_IERS = Params.EOP_IERS_data(Params.EOP_IERS_data(:,1) == x1, :);
    
    Params.delta_UT1 = Params.linearInterpolation(MJD_UTC, x0, x1, dataRow0_IERS(15), dataRow1_IERS(15));
    Params.LOD = Params.linearInterpolation(MJD_UTC, x0, x1, dataRow0_IERS(17), dataRow1_IERS(17)) / 1e3;
    Params.xp = Params.arcsec2rad * (Params.linearInterpolation(MJD_UTC, x0, x1, dataRow0_IERS(6), dataRow1_IERS(6)));
    Params.yp = Params.arcsec2rad * (Params.linearInterpolation(MJD_UTC, x0, x1, dataRow0_IERS(8), dataRow1_IERS(8)));
    Params.d_delta_psi_1980 = Params.arcsec2rad * (Params.linearInterpolation(MJD_UTC, x0, x1, dataRow0_IERS(20), dataRow1_IERS(20)) / 1e3);
    Params.d_delta_eps_1980 = Params.arcsec2rad * (Params.linearInterpolation(MJD_UTC, x0, x1, dataRow0_IERS(22), dataRow1_IERS(22)) / 1e3);
    % ---------------------------------------------------------------------

    %dataRow_celestrak = Params.EOP_celestrak_data(Params.EOP_celestrak_data(:,1) == year & Params.EOP_celestrak_data(:,2) == month & Params.EOP_celestrak_data(:,3) == day, :);
    dataRow_celestrak = Params.EOP_celestrak_data(Params.EOP_celestrak_data(:,4) == x0, :);
    % Params.delta_UT1 = dataRow_celestrak(7);
    % Params.LOD = dataRow_celestrak(8);
    Params.delta_AT = dataRow_celestrak(13);
    % Params.xp = Params.arcsec2rad * (dataRow_celestrak(5));
    % Params.yp = Params.arcsec2rad * (dataRow_celestrak(6));
    % Params.d_delta_psi_1980 = Params.arcsec2rad * (dataRow_celestrak(9)); 
    % Params.d_delta_eps_1980 = Params.arcsec2rad * (dataRow_celestrak(10));

    Params.w_earth = 7.292115146706979e-5 * (1 - Params.LOD / 86400);
    %Params.w_earth = 7.292115146706979e-5;

    Params.JD_TT = JD_UTC + (Params.delta_AT / 86400) + (32.184 / 86400);
    Params.JD_UT1 = JD_UTC + Params.delta_UT1 / 86400;
    Params.T_TT = (Params.JD_TT - 2451545.0) / 36525;
    Params.T_UT1 = (Params.JD_UT1 - 2451545.0) / 36525;

    Params.JD_TDB = Params.JD_TT + (0.001657 * sin(628.3076 * Params.T_TT + 6.2401) ...
                  + 0.000022 * sin(575.3385 * Params.T_TT + 4.2970) ...
                  + 0.000014 * sin(1256.6152 * Params.T_TT + 6.1969) ...
                  + 0.000005 * sin(606.9777 * Params.T_TT + 4.0212) ...
                  + 0.000005 * sin(52.9691 * Params.T_TT + 0.4444) ...
                  + 0.000002 * sin(21.3299 * Params.T_TT + 5.5431) ...
                  + 0.000010 * Params.T_TT * sin(628.3076 * Params.T_TT + 4.2490)) / 86400;

    Params.T_TDB = (Params.JD_TDB - 2451545.0) / 36525;

end