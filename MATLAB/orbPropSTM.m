function Xdot = orbPropSTM(t, X, JD0_UTC, Params)

    JD_UTC = JD0_UTC + t / 86400; % current time
    Params = UpdateParams(JD_UTC, Params);

    r_GCRF = X(1:3);
    v_GCRF = X(4:6);
    rx = X(1);
    ry = X(2);
    rz = X(3);
    vx = X(4);
    vy = X(5);
    vz = X(6);
    Cd = X(7);
    b = X(8);

    r_sun_GCRF = SunPosition(Params);
    r_moon_GCRF = MoonPosition(Params);

    A = AmatJacobian(rx, ry, rz, vx, vy, vz, Cd, b, r_sun_GCRF(1), r_sun_GCRF(2), r_sun_GCRF(3), r_moon_GCRF(1), r_moon_GCRF(2), r_moon_GCRF(3));

    Phi = reshape(X(9:72), [8,8]);
    Phidot = A * Phi;

    a_GCRF = -Params.mu_earth * r_GCRF / norm(r_GCRF)^3 ...
          + AccelGrav20x20(r_GCRF, Params) ...
          + AccelThirdBody(r_GCRF, Params) ...
          + AccelDrag(r_GCRF, v_GCRF, Cd, Params) ...
          + AccelSRP(r_GCRF, v_GCRF, Params);

    Xdot = [v_GCRF; a_GCRF; 0; 0; Phidot(:)];

end