function Nnm = NormalizeCoeff(n)

    % if m == 0
    %     dk = 1;
    % else
    %     dk = 2;
    % end
    % 
    % Nnm = sqrt(factorial(n - m) * dk * (2 * n + 1) / factorial(n + m));

    Nnm = zeros(n + 1, n + 1);

    for ii = 0 : n

        for jj = 0 : ii

            if jj == 0
                d0m = 1;
            else
                d0m = 0;
            end
        
            Nnm(ii + 1, jj + 1) = sqrt(factorial(ii + jj) / ((2 - d0m) * factorial(ii - jj) * (2 * ii + 1)));

        end

    end

end