function Pnm = AssociatedLegendrePolys(x, n)

    Pnm = zeros(n + 1, n + 1); 
    Pnm(1, 1) = 1;
    
    if n >= 1

        Pnm(2, 1) = x;
        Pnm(2, 2) = sqrt(1 - x^2);

    end

    for ii = 2 : n

        for jj = 0 : ii

            if jj == 0

                Pnm(ii+1,jj+1) = ((2 * ii-1) * x * Pnm(ii,jj+1) - (ii - 1) * Pnm(ii-1,jj+1)) / ii;

            elseif jj == ii

                Pnm(ii+1,jj+1) = sqrt(1 - x^2) * (2*ii - 1) * Pnm(ii,jj);

            else

                Pnm(ii+1,jj+1) = ((2 * ii - 1) * x * Pnm(ii,jj+1) - (ii + jj - 1) * Pnm(ii-1,jj+1)) / (ii - jj);

            end

        end

    end

end