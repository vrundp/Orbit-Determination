function FracSun = EclipseShadow(r_GCRF, Params)

    r_sun_GCRF = SunPosition(Params);

    r_sat_sun = r_sun_GCRF - r_GCRF;
    eta = acos(dot(-r_GCRF, r_sat_sun) / (norm(-r_GCRF) * norm(r_sat_sun)));
    rho_e = asin(Params.R_earth / norm(-r_GCRF));
    rho_s = asin(Params.R_sun / norm(r_sat_sun));

    if eta >= (rho_s + rho_e)

        FracSun = 1;

    elseif eta > (rho_e - rho_s) && eta < (rho_e + rho_s)

        L = 0.5 * (rho_e + rho_s + eta);
        q = 2 / eta * sqrt(L * (L - eta) * (L - rho_e) * (L - rho_s));

        if imag(q) ~= 0

            FracSun = 1; % antumbra case?

        else

            if eta^2 >= abs(rho_s^2 - rho_e^2)

                A = rho_s^2 * asin(q / rho_s) + rho_e^2 * asin(q / rho_e) - q * eta;

            elseif eta^2 < abs(rho_s^2 - rho_e^2)

                A = rho_e^2 * asin(q / rho_e) + (pi - asin(q / rho_s)) * rho_s^2 - q * eta;

            end

            FracSun = 1 - A / (pi * rho_s^2);

        end

    else

        FracSun = 0;
        
    end

end