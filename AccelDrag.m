function a_GCRF = AccelDrag(r_GCRF, v_GCRF, Cd, Params)

    %Params = UpdateParams(JD_UTC, Params);

    v_rel = v_GCRF - cross([0; 0; Params.w_earth], r_GCRF); % ECEF transformation???

    nz = r_GCRF / norm(r_GCRF); 
    ny = cross(r_GCRF, v_GCRF) / norm(cross(r_GCRF, v_GCRF));
    nx = cross(ny, nz);

    r_sun_GCRF = SunPosition(Params);
    r_sat_sun = r_sun_GCRF - r_GCRF;
    nsp = r_sat_sun / norm(r_sat_sun);

    cos_theta_nx_pos = max(0, dot(v_rel / norm(v_rel), nx));
    cos_theta_nx_neg = max(0, dot(v_rel / norm(v_rel), -nx));
    cos_theta_ny_pos = max(0, dot(v_rel / norm(v_rel), ny));
    cos_theta_ny_neg = max(0, dot(v_rel / norm(v_rel), -ny));
    cos_theta_nz_pos = max(0, dot(v_rel / norm(v_rel), nz));
    cos_theta_nz_neg = max(0, dot(v_rel / norm(v_rel), -nz));
    cos_theta_nsp_pos = max(0, dot(v_rel / norm(v_rel), nsp));
    cos_theta_nsp_neg = max(0, dot(v_rel / norm(v_rel), -nsp));

    A = Params.FaceX.area * cos_theta_nx_pos ...
      + Params.FaceX.area * cos_theta_nx_neg ...
      + Params.FaceY.area * cos_theta_ny_pos ...
      + Params.FaceY.area * cos_theta_ny_neg ...
      + Params.FaceZpos.area * cos_theta_nz_pos ...
      + Params.FaceZneg.area * cos_theta_nz_neg ...
      + Params.SolarPanel.area * cos_theta_nsp_pos ...
      + Params.SolarPanel.area * cos_theta_nsp_neg;

    r = norm(r_GCRF);
    rho = Params.rho0 * exp(-(r - Params.r0) / Params.H);

    a_GCRF = -0.5 * rho * Cd * A / Params.mass * norm(v_rel) * v_rel;

end
