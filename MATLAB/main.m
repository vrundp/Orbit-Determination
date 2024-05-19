%% Final Project - Methods of Orbit Determination
% Vrund Patel
format long g

% Plot settings
% set(0,'DefaultTextFontSize',12);
% set(0,'DefaultAxesFontSize',12);
% set(0,'DefaultLegendFontSize',12);
% set(0,'DefaultLegendInterpreter','Latex');
% set(0,'DefaultAxesTickLabelInterpreter','Latex');
% set(0,'DefaultTextInterpreter','Latex');

r0_GCRF = [6984.45711518852; 1612.2547582643; 13.0925904314402] .* 1e3; % m
v0_GCRF = [-1.67667852227336; 7.26143715396544; 0.259889857225218] .* 1e3; % m/s

JD0_UTC = cal2JD(2018, 3, 23, 8, 55, 3.0);
JDdv1_UTC = cal2JD(2018, 3, 30, 8, 55, 3.0);
tPropMax = (JDdv1_UTC - JD0_UTC) * 86400; % time of DV1 since initial epoch

InitParams;
Params = UpdateParams(JD0_UTC, Params);
Sensor_data_3days = load("LEO_DATA_Apparent_3Days.mat");
Sensor_data_6days = load("LEO_DATA_Apparent_Days4-6.mat");
Sensor_data_6days.LEO_DATA_Apparent(:,2) = Sensor_data_6days.LEO_DATA_Apparent(:,2) + 3*86400;
Sensor_data.LEO_DATA_Apparent = [Sensor_data_3days.LEO_DATA_Apparent; Sensor_data_6days.LEO_DATA_Apparent];

% a_20x20 = AccelGrav20x20(r0_GCRF, Params)
% a_3b = AccelThirdBody(r0_GCRF, Params)
% a_drag = AccelDrag(r0_GCRF, v0_GCRF, 1.88, Params)
% a_srp  = AccelSRP(r0_GCRF, v0_GCRF, Params)
% a_tot = -Params.mu_earth * r0_GCRF / norm(r0_GCRF)^3 + AccelGrav20x20(r0_GCRF, Params) + AccelThirdBody(r0_GCRF, Params) + AccelDrag(r0_GCRF, v0_GCRF, 1.88, Params) + AccelSRP(r0_GCRF, v0_GCRF, Params)

%tic
%tMeasMax = 86220;
%tMeasMax = 258780;
tMeasMax = 518400;
% X0 = [r0_GCRF; v0_GCRF; 1.88; 20];
% P0 = diag([10e3, 10e3, 10e3, 1, 1, 1, 0.01, 5]).^2;

% Short arc initial condition
% X0 = [-2308694.19239785; 6780309.29965407; 221201.361692046; -7048.61472970489; -2412.13027110639; -134.864992679877; 2.23958790382862; 19.4427942522296];
% X0 = [-2308694.19239785; 6780309.29965407; 221201.361692046; -7048.61472970489; -2412.13027110639; -134.864992679877; 2.00; 20];
% P0 = diag([0.0150747483971509, 0.00433058020276977, 0.000671670146447589, 1.4647321254137e-09, 1.68666616932259e-08, 7.30011805666739e-10, 5.84676917082706e-05, 0.331133920804175]).*2;
% P0 = diag([10, 10, 10, 1, 1, 1, 0.01, 5]).^2;

% Good initial condition
X0 = [6978573.77997675; 1616518.33237809; 19632.2800280167; -1662.92967859376; 7260.85500411775; 270.486436585235; 2.00; 20];
%P0 = diag([0.00855749100446701, 0.0342724449001253, 0.030238039791584, 9.6708430796788e-07, 1.07922862424559e-07, 1.35654633662829e-05, 3.52350210002855e-06, 0.167507424020204]);
P0 = diag([5, 5, 5, 1, 1, 1, 0.06125, 0.5]).^2;
% 
[xhatHist, PHist, xbarHist, PbarHist, sigmaXpbarHist, tMeasHist, statIdHist, residHist, residCovHist, xhat_deliv, P_deliv] = UnscentedKalmanFilterJah...
         (X0, P0, tMeasMax, tPropMax, JD0_UTC, Sensor_data, 'F', Params);
% [xhatHist, PHist, tMeasHist, statIdHist, residHist, residCovHist, xhat_deliv, P_deliv] = ExtendedKalmanFilter...
%          (X0, P0, tMeasMax, tPropMax, JD0_UTC, Sensor_data, 'F', Params);

xhat_deliv
P_deliv(1:3,1:3)./1e6

%[xstarHist, PstarHist] = UnscentedKalmanSmoother(xhatHist, PHist, xbarHist, PbarHist, sigmaXpbarHist);
% PlotPostFitResiduals(tMeasHist, residHist, residCovHist, statIdHist, 'F')
%[~, XHist] = ode113(@orbPropHF, [tMeasMax; 0], xhatHist(end,1:6)', Params.ode_options, JD0_UTC, Params);


% Consider UKF ------------------------------------------------------------
% X0 = [6978573.77997675; 1616518.33237809; 19632.2800280167; -1662.92967859376; 7260.85500411775; 270.486436585235; 20];
% %X0 = [-2308694.19239785; 6780309.29965407; 221201.361692046; -7048.61472970489; -2412.13027110639; -134.864992679877; 20];
% Pxx0 = diag([10, 10, 10, 1, 1, 1, 2]).^2; % tune?
% c0 = 2.00;
% Pcc0 = 0.075^2; % tune
% [zhatHist, PzzHist, tMeasHist, statIdHist, residHist, residCovHist, xhat_deliv, P_deliv] = ConsiderUKF...
%          (X0, Pxx0, c0, Pcc0, tMeasMax, tPropMax, JD0_UTC, Sensor_data, 'F', Params);
% xhat_deliv
% P_deliv(1:3,1:3)./1e6
%PlotPostFitResiduals(tMeasHist, residHist, residCovHist, statIdHist, 'F')
% -------------------------------------------------------------------------
%toc
% best from NAG: Cd 0.06, Q 1e-9 1e-9 1e-9
% best so far?: Cd 0.075, Q 1e-9 5e-9 1e-9
%% Filter all cases
%tMeasMax = 258780;
tMeasMax = 518400;
X0 = [6978573.77997675; 1616518.33237809; 19632.2800280167; -1662.92967859376; 7260.85500411775; 270.486436585235; 2.00; 20];
P0 = diag([5, 5, 5, 1, 1, 1, 0.06125, 2]).^2;
tic
for ii = ['F','A','B','C','D','E','G']
    if ii == 'G'
        X0 = [-2308694.19239785; 6780309.29965407; 221201.361692046; -7048.61472970489; -2412.13027110639; -134.864992679877; 2.00; 20];
    end
    [xhatHist, PHist, xbarHist, PbarHist, sigmaXpbarHist, tMeasHist, statIdHist, residHist, residCovHist, xhat_deliv, P_deliv] = UnscentedKalmanFilterJah...
    (X0, P0, tMeasMax, tPropMax, JD0_UTC, Sensor_data, ii, Params);
    if ii == 'A'
        Range.xhat = xhat_deliv;
        Range.P = P_deliv;
        save('Phase2Range6daysCUKF.mat', 'Range');
        patel_pos_caseA = xhat_deliv(1:3) / 1e3;
        patel_poscov_caseA = P_deliv(1:3,1:3) ./ 1e6;
        PlotPostFitResiduals(tMeasHist, residHist, residCovHist, statIdHist, ii)
        saveas(gcf, 'CaseA6daysCUKF.fig')
    elseif ii == 'B'
        RangeRate.xhat = xhat_deliv;
        RangeRate.P = P_deliv;
        save('Phase2RangeRate6daysCUKF.mat', 'RangeRate');
        patel_pos_caseB = xhat_deliv(1:3) / 1e3;
        patel_poscov_caseB = P_deliv(1:3,1:3) ./ 1e6;
        PlotPostFitResiduals(tMeasHist, residHist, residCovHist, statIdHist, ii)
        saveas(gcf, 'CaseB6daysCUKF.fig')
    elseif ii == 'C'
        Kwaj.xhat = xhat_deliv;
        Kwaj.P = P_deliv;
        save('Phase2Kwaj6daysCUKF.mat', 'Kwaj');
        patel_pos_caseC = xhat_deliv(1:3) / 1e3;
        patel_poscov_caseC = P_deliv(1:3,1:3) ./ 1e6;
        PlotPostFitResiduals(tMeasHist, residHist, residCovHist, statIdHist, ii)
        saveas(gcf, 'CaseC6daysCUKF.fig')
    elseif ii == 'D'
        DG.xhat = xhat_deliv;
        DG.P = P_deliv;
        save('Phase2DG6daysCUKF.mat', 'DG');
        patel_pos_caseD = xhat_deliv(1:3) / 1e3;
        patel_poscov_caseD = P_deliv(1:3,1:3) ./ 1e6;
        PlotPostFitResiduals(tMeasHist, residHist, residCovHist, statIdHist, ii)
        saveas(gcf, 'CaseD6daysCUKF.fig')
    elseif ii == 'E'
        Arecibo.xhat = xhat_deliv;
        Arecibo.P = P_deliv;
        save('Phase2Arecibo6daysCUKF.mat', 'Arecibo');
        patel_pos_caseE = xhat_deliv(1:3) / 1e3;
        patel_poscov_caseE = P_deliv(1:3,1:3) ./ 1e6;
        PlotPostFitResiduals(tMeasHist, residHist, residCovHist, statIdHist, ii)
        saveas(gcf, 'CaseE6daysCUKF.fig')
    elseif ii == 'F'
        Long.xhat = xhat_deliv;
        Long.P = P_deliv;
        save('Phase2Long6daysCUKF.mat', 'Long');
        patel_pos_caseF = xhat_deliv(1:3) / 1e3;
        patel_poscov_caseF = P_deliv(1:3,1:3) ./ 1e6;
        PlotPostFitResiduals(tMeasHist, residHist, residCovHist, statIdHist, ii)
        saveas(gcf, 'CaseF6daysCUKF.fig')
    elseif ii == 'G'
        Short.xhat = xhat_deliv;
        Short.P = P_deliv;
        save('Phase2Short6daysCUKF.mat', 'Short');
        patel_pos_caseG = xhat_deliv(1:3) / 1e3;
        patel_poscov_caseG = P_deliv(1:3,1:3) ./ 1e6;
        PlotPostFitResiduals(tMeasHist, residHist, residCovHist, statIdHist, ii)
        saveas(gcf, 'CaseG6daysCUKF.fig')
    end
end
toc
save('patel.mat', "patel_poscov_caseG", "patel_pos_caseG", "patel_poscov_caseF", "patel_pos_caseF", ...
     "patel_poscov_caseE", "patel_pos_caseE", "patel_poscov_caseD", "patel_pos_caseD", "patel_poscov_caseC", "patel_pos_caseC",...
     "patel_poscov_caseB", "patel_pos_caseB", "patel_poscov_caseA", "patel_pos_caseA");

%% Delivery stats

Long = load('Phase2Long6daysCUKF.mat');
Short = load('Phase2Short6daysCUKF.mat');
Range = load('Phase2Range6daysCUKF.mat');
RangeRate = load('Phase2RangeRate6daysCUKF.mat');
Kwaj = load('Phase2Kwaj6daysCUKF.mat');
DG = load('Phase2DG6daysCUKF.mat');
Arecibo = load('Phase2Arecibo6daysCUKF.mat');

pos_long = Long.Long.xhat(1:3) ./1e3
vel_long = Long.Long.xhat(4:6) ./1e3;
posCov_long = Long.Long.P(1:3,1:3) ./1e6;

R = pos_long / norm(pos_long);
C = cross(pos_long / norm(pos_long), vel_long / norm(vel_long)) / norm(cross(pos_long / norm(pos_long), vel_long / norm(vel_long)));
I = cross(C, R);
T_RIC_GCRF = orthodcm([R, I, C]');

pos_short = Short.Short.xhat(1:3) ./1e3
posCov_short = Short.Short.P(1:3,1:3) ./1e6;

pos_range = Range.Range.xhat(1:3) ./1e3
posCov_range = Range.Range.P(1:3,1:3) ./1e6;

pos_rangerate = RangeRate.RangeRate.xhat(1:3) ./1e3
posCov_rangerate = RangeRate.RangeRate.P(1:3,1:3) ./1e6;

pos_kwaj = Kwaj.Kwaj.xhat(1:3) ./1e3
posCov_kwaj = Kwaj.Kwaj.P(1:3,1:3) ./1e6;

pos_dg = DG.DG.xhat(1:3) ./1e3
posCov_dg = DG.DG.P(1:3,1:3) ./1e6;

pos_arecibo = Arecibo.Arecibo.xhat(1:3) ./1e3
posCov_arecibo = Arecibo.Arecibo.P(1:3,1:3) ./1e6;

% position estimates in RIC
pos_long_RIC = T_RIC_GCRF * pos_long;
pos_short_RIC = T_RIC_GCRF * pos_short;
pos_range_RIC = T_RIC_GCRF * pos_range;
pos_rangerate_RIC = T_RIC_GCRF * pos_rangerate;
pos_kwaj_RIC = T_RIC_GCRF * pos_kwaj;
pos_dg_RIC = T_RIC_GCRF * pos_dg;
pos_arecibo_RIC = T_RIC_GCRF * pos_arecibo;
% covariences in RIC
posCov_long_RIC = T_RIC_GCRF * posCov_long * T_RIC_GCRF';
posCov_short_RIC = T_RIC_GCRF * posCov_short * T_RIC_GCRF';
posCov_range_RIC = T_RIC_GCRF * posCov_range * T_RIC_GCRF';
posCov_rangerate_RIC = T_RIC_GCRF * posCov_rangerate * T_RIC_GCRF';
posCov_kwaj_RIC = T_RIC_GCRF * posCov_kwaj * T_RIC_GCRF';
posCov_dg_RIC = T_RIC_GCRF * posCov_dg * T_RIC_GCRF';
posCov_arecibo_RIC = T_RIC_GCRF * posCov_arecibo * T_RIC_GCRF';

PlotEstimatesRIC(pos_long_RIC, pos_short_RIC, pos_range_RIC, pos_rangerate_RIC, pos_kwaj_RIC, pos_dg_RIC, pos_arecibo_RIC,...
                 posCov_long_RIC, posCov_short_RIC, posCov_range_RIC, posCov_rangerate_RIC, posCov_kwaj_RIC, posCov_dg_RIC, posCov_arecibo_RIC)
