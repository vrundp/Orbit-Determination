function W = PolarMatrix(Params)
% ITRF to PEF
    
    %W = Rot1(Params.yp) * Rot2(Params.xp);
    %W = [1 0 -xp; 0 1 yp; xp -yp 1];

    cos_xp = cos(Params.xp);
    sin_xp = sin(Params.xp);
    cos_yp = cos(Params.yp);
    sin_yp = sin(Params.yp);

    % W = [           cos(Params.xp),           0,           -sin(Params.xp);
    %      sin(Params.xp)*sin(Params.yp),  cos(Params.yp), cos(Params.xp)*sin(Params.yp);
    %      cos(Params.yp)*sin(Params.xp), -sin(Params.yp), cos(Params.xp)*cos(Params.yp)];

    W = [cos_xp, 0, -sin_xp;
         sin_xp * sin_yp, cos_yp, cos_xp * sin_yp;
         cos_yp * sin_xp, -sin_yp, cos_xp * cos_yp];
    
end