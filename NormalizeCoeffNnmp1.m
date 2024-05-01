function Nnmp1 = NormalizeCoeffNnmp1(n)
% Special Normalization coefficient for latitude partial of potential
% Index based on actual n and m value not m+1

    Nnmp1 = zeros(n + 1, n + 1);

    for ii = 0 : n

        for jj = 0 : ii

            if jj == 0
                d0m = 1;
            else
                d0m = 0;
            end

            Nnmp1(ii + 1, jj + 1) = sqrt((ii + jj + 1) * (ii - jj) * (2 - d0m) / 2);

        end

    end

end