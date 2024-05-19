function [xstarHist, PstarHist] = UnscentedKalmanSmoother(xhatHist, PHist, xbarHist, PbarHist, sigmaXpbarHist)

    maxIdx = length(xhatHist);

    nx = length(xhatHist(1,:));
    w = 1 / (2 * nx);

    xstarHist = zeros(maxIdx, nx);
    PstarHist = zeros(maxIdx, nx^2);

    xstarkp1 = xhatHist(end,:)';
    Pstarkp1 = reshape(PHist(end,:), [nx,nx]);
    xstarHist(maxIdx,:) = xstarkp1;
    PstarHist(maxIdx,:) = Pstarkp1(:);

    for kk = maxIdx:-1:2
        % Index kk corresponds to k+1, kk-1 corresponds to k

        xhatk = xhatHist(kk-1,:)';
        Pk = reshape(PHist(kk-1,:), [nx,nx]);
        Shatk = chol(Pk, 'lower');
        sigmaXhatkHist = zeros(2 * nx, nx);
        for ii = 1 : 2 * nx
            jj = ii;
            sign = 1;
            if ii > nx
                jj = ii - nx;
                sign = -1;
            end
            sigmaXhatk = xhatk + sign * sqrt(nx) * Shatk(:,jj);
            sigmaXhatkHist(ii,:) = sigmaXhatk;
        end

        xbarkp1 = xbarHist(kk,:)';
        Pbarkp1 = reshape(PbarHist(kk,:), [nx,nx]);

        Ck = 0;
        for ii = 1 : 2 * nx
            Ck = Ck + w * (sigmaXhatkHist(ii,:)' - xhatk) * (sigmaXpbarHist{kk}(ii,:)' - xbarkp1)';
        end

        Ak = Ck / Pbarkp1;

        xstark = xhatk + Ak * (xstarkp1 - xbarkp1);
        Pstark = Pk + Ak * (Pstarkp1 - Pbarkp1) * Ak';

        xstarkp1 = xstark;
        Pstarkp1 = Pstark;

        xstarHist(kk-1,:) = xstark;
        PstarHist(kk-1,:) = Pstark(:);

    end
    
end