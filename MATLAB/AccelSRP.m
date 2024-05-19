function a_GCRF = AccelSRP(r_GCRF, v_GCRF, Params)

    %Params = UpdateParams(JD_UTC, Params);
    
    FracSun = EclipseShadow(r_GCRF, Params);

    p_srp = 1367 / Params.c;

    r_sun_GCRF = SunPosition(Params);
    r_sat_sun = r_sun_GCRF - r_GCRF;
    uhat = r_sat_sun / norm(r_sat_sun);
    d = norm(r_sat_sun) / Params.AU;

    nz = r_GCRF / norm(r_GCRF);
    ny = cross(r_GCRF, v_GCRF) / norm(cross(r_GCRF, v_GCRF));
    nx = cross(ny, nz);
    nsp = uhat;

    cos_theta_nx_pos = max(0, dot(uhat, nx));
    cos_theta_nx_neg = max(0, dot(uhat, -nx));
    cos_theta_ny_pos = max(0, dot(uhat, ny));
    cos_theta_ny_neg = max(0, dot(uhat, -ny));
    cos_theta_nz_pos = max(0, dot(uhat, nz));
    cos_theta_nz_neg = max(0, dot(uhat, -nz));
    cos_theta_nsp_pos = max(0, dot(uhat, nsp));

    A = Params.FaceX.area * ((2 / 3 * Params.FaceX.Cd * cos_theta_nx_pos + 2 * Params.FaceX.Cs * cos_theta_nx_pos^2) .* nx ...
                          + ((1 - Params.FaceX.Cs / 2) * cos_theta_nx_pos^2) .* uhat) ...
      + Params.FaceX.area * ((2 / 3 * Params.FaceX.Cd * cos_theta_nx_neg + 2 * Params.FaceX.Cs * cos_theta_nx_neg^2) .* -nx ...
                          + ((1 - Params.FaceX.Cs / 2) * cos_theta_nx_neg^2) .* uhat) ...
      + Params.FaceY.area * ((2 / 3 * Params.FaceY.Cd * cos_theta_ny_pos + 2 * Params.FaceY.Cs * cos_theta_ny_pos^2) .* ny ...
                          + ((1 - Params.FaceY.Cs / 2) * cos_theta_ny_pos^2) .* uhat) ...
      + Params.FaceY.area * ((2 / 3 * Params.FaceY.Cd * cos_theta_ny_neg + 2 * Params.FaceY.Cs * cos_theta_ny_neg^2) .* -ny ...
                          + ((1 - Params.FaceY.Cs / 2) * cos_theta_ny_neg^2) .* uhat) ...
      + Params.FaceZpos.area * ((2 / 3 * Params.FaceZpos.Cd * cos_theta_nz_pos + 2 * Params.FaceZpos.Cs * cos_theta_nz_pos^2) .* nz ...
                             + ((1 - Params.FaceZpos.Cs / 2) * cos_theta_nz_pos^2) .* uhat) ...
      + Params.FaceZneg.area * ((2 / 3 * Params.FaceZneg.Cd * cos_theta_nz_neg + 2 * Params.FaceZneg.Cs * cos_theta_nz_neg^2) .* -nz ...
                             + ((1 - Params.FaceZneg.Cs / 2) * cos_theta_nz_neg^2) .* uhat) ...
      + Params.SolarPanel.area * ((2 / 3 * Params.SolarPanel.Cd * cos_theta_nsp_pos + 2 * Params.SolarPanel.Cs * cos_theta_nsp_pos^2) .* nsp ...
                               + ((1 - Params.SolarPanel.Cs / 2) * cos_theta_nsp_pos^2) .* uhat);

    a_GCRF = -p_srp * FracSun .* A ./ Params.mass ./ d^2;

end