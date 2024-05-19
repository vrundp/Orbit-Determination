function date = JD2cal(JD)

    % Algorithm valid 1900 - 2100
    
    LMonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    
    T_1900 = (JD - 2415019.5) / 365.25;
    year = 1900 + fix(T_1900);
    
    leap_years = fix((year - 1900 - 1) * 0.25);
    if mod(year, 4) == 0
        LMonth(2) = 29;
    end
    
    days = (JD - 2415019.5) - ((year - 1900) * 365.0 + leap_years);
    if days < 1.0
        year = year - 1;
        leap_years = fix((year - 1900 - 1) * 0.25);
        days = (JD - 2415019.5) - ((year - 1900) * 365.0 + leap_years);
    end
    
    DOY = fix(days);
    
    month = 0;
    day = 0;
    sum_days = 0;
    for i = 1 : length(LMonth)
        sum_days = sum_days + LMonth(i);
        if sum_days > DOY
            day = DOY - (sum_days - LMonth(i));
            month = i;
            break;
        end
    end
    
    tau = (days - DOY) * 24;
    hour = fix(tau);
    minute = fix((tau - hour) * 60);
    seconds = (tau - hour - (minute / 60)) * 3600;
    
    date = [year, month, day, hour, minute, seconds];

end
