function Xdot = orbPropHF(t, X, JD0_UTC, Params)

    JD_UTC = JD0_UTC + t / 86400; % current time
    Params = UpdateParams(JD_UTC, Params);

    r_GCRF = X(1:3);
    v_GCRF = X(4:6);
    
    if length(X) >= 7
        Cd = X(7);
    else
        Cd = 1.88;
        %Cd = 2.00;
    end

    %Cd = 1.88;
    %Cd = 2.00;

    %Cd = X(8); % consider filter
    
    a_GCRF = -Params.mu_earth * r_GCRF / norm(r_GCRF)^3 ...
             + AccelGrav20x20(r_GCRF, Params) ...
             + AccelThirdBody(r_GCRF, Params) ...
             + AccelDrag(r_GCRF, v_GCRF, Cd, Params) ...
             + AccelSRP(r_GCRF, v_GCRF, Params);

    Xdot = [v_GCRF; a_GCRF];
    if length(X) == 7
        Xdot = [v_GCRF; a_GCRF; 0];
    elseif length(X) > 7
        Xdot = [v_GCRF; a_GCRF; zeros(length(X)-6,1)];
    end

end