function P = PrecessionMatrix(Params)
% MOD to GCRF

    T_TT = Params.T_TT;

    zeta = Params.arcsec2rad * (2306.2181 * T_TT + 0.30188 * T_TT^2 + 0.017998 * T_TT^3);
    theta = Params.arcsec2rad * (2004.3109 * T_TT - 0.42665 * T_TT^2 - 0.041833 * T_TT^3);
    z = Params.arcsec2rad * (2306.2181 * T_TT + 1.09468 * T_TT^2 + 0.018203 * T_TT^3);
    
    %P = Rot3(zeta) * Rot2(-theta) * Rot3(z);

    cos_neg_theta = cos(-theta);
    sin_neg_theta = sin(-theta);
    cos_zeta = cos(zeta);
    sin_zeta = sin(zeta);
    cos_z = cos(z);
    sin_z = sin(z);

    % P = [  cos(-theta)*cos(zeta)*cos(z) - sin(zeta)*sin(z), cos(z)*sin(zeta) + cos(-theta)*cos(zeta)*sin(z), -cos(zeta)*sin(-theta);
    %      - cos(zeta)*sin(z) - cos(-theta)*cos(z)*sin(zeta), cos(zeta)*cos(z) - cos(-theta)*sin(zeta)*sin(z),  sin(-theta)*sin(zeta);
    %                              cos(z)*sin(-theta),                                    sin(-theta)*sin(z),             cos(-theta)];

    P = [cos_neg_theta * cos_zeta * cos_z - sin_zeta * sin_z, cos_z * sin_zeta + cos_neg_theta * cos_zeta * sin_z, -cos_zeta * sin_neg_theta;
        -cos_zeta * sin_z - cos_neg_theta * cos_z * sin_zeta, cos_zeta * cos_z - cos_neg_theta * sin_zeta * sin_z, sin_neg_theta * sin_zeta;
         cos_z * sin_neg_theta, sin_neg_theta * sin_z, cos_neg_theta];

end