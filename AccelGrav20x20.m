function a_GCRF = AccelGrav20x20(r_GCRF, Params)

    %Params = UpdateParams(JD_UTC, Params);

    W = PolarMatrix(Params);
    R = SiderealMatrix(Params);
    N = NutationMatrix(Params);
    P = PrecessionMatrix(Params);

    r_ITRF = W' * R' * N' * P' * r_GCRF;
    r = norm(r_ITRF);
    ri = r_ITRF(1);
    rj = r_ITRF(2);
    rk = r_ITRF(3);

    phi = asin(rk / r);
    lambda = atan2(rj, ri);

    Pnm = AssociatedLegendrePolys(sin(phi), Params.Nmax);

    sum_dVdr = 0;
    sum_dVdphi = 0;
    sum_dVdlambda = 0;
    for n = 2 : Params.Nmax

        for m = 0 : n

            Pbar = Pnm(n+1,m+1) / Params.Nnm(n+1,m+1);
            if m == n
                Pbar_nmp1 = 0;
            else
                Pbar_nmp1 = Pnm(n+1,m+2) / Params.Nnm(n+1,m+2) * Params.Nnmp1(n+1,m+1);
            end
            Cbar = Params.Cnm(n+1,m+1);
            Sbar = Params.Snm(n+1,m+1);

            sum_dVdr = sum_dVdr + (Params.R_earth / r)^n * (n + 1) * Pbar * (Cbar * cos(m * lambda) + Sbar * sin(m * lambda));
            sum_dVdphi = sum_dVdphi + (Params.R_earth / r)^n * (Pbar_nmp1 - m * tan(phi) * Pbar) * (Cbar * cos(m * lambda) + Sbar * sin(m * lambda));
            sum_dVdlambda = sum_dVdlambda + (Params.R_earth / r)^n * m * Pbar * (Sbar * cos(m * lambda) - Cbar * sin(m * lambda));

        end

    end

    dVdr = -(Params.mu_earth / r^2) * sum_dVdr;
    dVdphi = (Params.mu_earth / r) * sum_dVdphi;
    dVdlambda = (Params.mu_earth / r) * sum_dVdlambda;

    ai = (dVdr / r - dVdphi * rk / r^2 / sqrt(ri^2 + rj^2)) * ri - (dVdlambda / (ri^2 + rj^2)) * rj;
    aj = (dVdr / r - dVdphi * rk / r^2 / sqrt(ri^2 + rj^2)) * rj + (dVdlambda / (ri^2 + rj^2)) * ri;
    ak = (dVdr / r) * rk + dVdphi * sqrt(ri^2 + rj^2) / r^2;

    a_ITRF = [ai; aj; ak];

    a_GCRF = P * N * R * W * a_ITRF;

end

