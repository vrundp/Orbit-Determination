% Utility Parameters for OD Project

Params.deg2rad = pi / 180;
Params.arcsec2rad = pi / 180 / 3600;

Params.linearInterpolation = @(x,x0,x1,y0,y1)y0+(x-x0)*(y1-y0)/(x1-x0);
%Params.wrapTo2Pi = @(x) mod(x, 2 * pi);

EOP_celestrak_data_all = readmatrix('EOP_All_data.txt');
EOP_IERS_data_all = readmatrix('finals.all.csv');

Params.EOP_celestrak_data = EOP_celestrak_data_all;%(20558:20567+5, :);
Params.EOP_IERS_data = EOP_IERS_data_all;%(16516:16525+5, :);
Params.nut_data = load('nut80.dat');

% Call UpdateParams to update
% Params.delta_UT1 = 0;
% Params.LOD = 0;
% Params.xp = 0;
% Params.yp = 0;
% Params.d_delta_psi_1980 = 0;
% Params.d_delta_eps_1980 = 0;
% Params.delta_AT = 0;
% Params.w_earth = 7.292115146706979e-5;
% Params.JD_TT = 0;
% Params.JD_UT1 = 0;
% Params.T_TT = 0;
% Params.T_UT1 = 0;
% Params.JD_TDB = 0;
% Params.T_TDB = 0;

egm96coeffs = load('egm96coeffs.txt');
Params.Nmax = 20;
Params.Nnm = NormalizeCoeff(Params.Nmax);
Params.Nnmp1 = NormalizeCoeffNnmp1(Params.Nmax);
Params.Cnm = zeros(Params.Nmax+1);
Params.Snm = zeros(Params.Nmax+1);
for ii = 1 : 228
    Params.Cnm(egm96coeffs(ii,1)+1, egm96coeffs(ii,2)+1) = egm96coeffs(ii,3);
    Params.Snm(egm96coeffs(ii,1)+1, egm96coeffs(ii,2)+1) = egm96coeffs(ii,4);
end

Params.ode_options = odeset('RelTol', 3e-14, 'AbsTol', 1e-16);
%Params.ode_options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

Params.Kwaj.r_ITRF = [-6143584; 1364250; 1033743];      % m
Params.DG.r_ITRF = [1907295; 6030810; -817119];         % m
Params.Arecibo.r_ITRF = [2390310; -5564341; 1994578];   % m

Params.Kwaj.sigma_r = 10;                               % m
Params.Kwaj.sigma_rr = 0.5 / 1e3;                       % m/s
Params.DG.sigma_r = 5;                                  % m
Params.DG.sigma_rr = 1 / 1e3;                           % m/s
Params.Arecibo.sigma_r = 10;                            % m
Params.Arecibo.sigma_rr = 0.5 / 1e3;                    % m/s

Params.Qk_RIC = diag([1e-9 1e-9 1e-9]).^2;            % m/s^2
%Params.Qk_RIC = diag([0, 0, 0]).^2;

Params.mu_earth = 398600.4415 * 1e3^3;                  % m^3/s^2
Params.R_earth = 6378.1363 * 1e3;                       % m
Params.J2 = -Params.Cnm(3,1) * sqrt(2 * 2 + 1);
Params.J3 = -Params.Cnm(4,1) * sqrt(2 * 3 + 1);
Params.J4 = -Params.Cnm(5,1) * sqrt(2 * 4 + 1);
Params.J5 = -Params.Cnm(6,1) * sqrt(2 * 5 + 1);
Params.J6 = -Params.Cnm(7,1) * sqrt(2 * 6 + 1);
Params.mu_sun = 132712440018 * 1e3^3;                   % m^3/s^2
Params.AU = 149597870.7 * 1e3;                          % m
Params.mu_moon = 4902.800066 * 1e3^3;                   % m^3/s^2
Params.e_earth = 0.081819221456;
Params.c = 299792458;                                   % m/s
Params.R_sun = 695700 * 1e3;                            % m

Params.r0 = 700000.0 + Params.R_earth;                  % m
Params.H = 88667.0;                                     % m
Params.rho0 = 3.614e-13;                                % kg/m^3

Params.mass = 2000;                                     % kg
Params.FaceX.area = 6;                                  % m^2
Params.FaceY.area = 8;                                  % m^2
Params.FaceZpos.area = 12;                              % m^2
Params.FaceZneg.area = 12;                              % m^2
Params.SolarPanel.area = 15;                            % m^2

Params.FaceX.Cd = 0.04;
Params.FaceX.Cs = 0.59;
Params.FaceY.Cd = 0.04;
Params.FaceY.Cs = 0.59;
Params.FaceZpos.Cd = 0.80;
Params.FaceZpos.Cs = 0.04;
Params.FaceZneg.Cd = 0.28;
Params.FaceZneg.Cs = 0.18;
Params.SolarPanel.Cd = 0.04;
Params.SolarPanel.Cs = 0.04;



