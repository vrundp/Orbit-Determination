function a_GCRF = AccelThirdBody(r_GCRF, Params)

    %Params = UpdateParams(JD_UTC, Params);

    r_sun_GCRF = SunPosition(Params);
    r_moon_GCRF = MoonPosition(Params);

    r_sat_sun = r_sun_GCRF - r_GCRF;
    r_sat_moon = r_moon_GCRF - r_GCRF;

    a_GCRF = Params.mu_sun * (r_sat_sun / norm(r_sat_sun)^3 - r_sun_GCRF / norm(r_sun_GCRF)^3) ...
           + Params.mu_moon * (r_sat_moon / norm(r_sat_moon)^3 - r_moon_GCRF / norm(r_moon_GCRF)^3);

end