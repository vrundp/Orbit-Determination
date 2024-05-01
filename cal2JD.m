function JD = cal2JD(year, month, day, hour, minute, seconds)

    JD = 367 * year ... 
         - fix((7 * (year + fix((month + 9) / 12))) / 4) + fix(275 * month / 9) ...
         + day + 1721013.5 + ((hour + (minute + (seconds / 60)) / 60) / 24);

end