function [xhist, yhist] = CovarianceEllipse(P)

    [V, D] = eig(P);
    
    t = 0 : 0.01 : 2*pi;
    xhist = zeros(length(t), 1);
    yhist = zeros(length(t), 1);

    for ii = 1 : length(t)

        xy = V * sqrt(D) * [cos(t(ii)); sin(t(ii))];

        xhist(ii) = xy(1);
        yhist(ii) = xy(2);

    end

end