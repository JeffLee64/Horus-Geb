function ztemp = VenusAtmosphereTemperature(altitude)
    % Return the temperature at a chosen altitude, in degrees Kelvin

    % Blumenthal, Kay, Palen, Smith (2012). Understanding Our Universe. 
    % New York: W.W. Norton & Company. p. 167. ISBN 9780393912104.

    if altitude >= 150000
        ztemp = 127;
    elseif altitude >= 110000
        ztemp = ((altitude-150000)/((150000-110000)/((127)-(-112))))+(127);
    elseif altitude >= 100000
        ztemp = -112;
    elseif altitude >= 90000
        ztemp = ((altitude-100000)/((100000-90000)/((-112)-(-104))))+(-112);
    elseif altitude >= 80000
        ztemp = ((altitude-90000)/((90000-80000)/((-104)-(-76))))+(-104);
    elseif altitude >= 70000
        ztemp = ((altitude-80000)/((80000-70000)/((-76)-(-43))))+(-76);
    elseif altitude >= 60000
        ztemp = ((altitude-70000)/((70000-60000)/((-43)-(-10))))+(-43);
    elseif altitude >= 55000
        ztemp = ((altitude-60000)/((60000-55000)/((-10)-(27))))+(-10);
    elseif altitude >= 50000
        ztemp = ((altitude-55000)/((55000-50000)/((27)-(75))))+(27);
    elseif altitude >= 45000
        ztemp = ((altitude-50000)/((50000-45000)/((75)-(110))))+(75);
    elseif altitude >= 40000
        ztemp = ((altitude-45000)/((45000-40000)/((110)-(143))))+(110);
    elseif altitude >= 30000
        ztemp = ((altitude-40000)/((40000-30000)/((143)-(222))))+(143);
    elseif altitude >= 20000
        ztemp = ((altitude-30000)/((30000-20000)/((222)-(306))))+(222);
    elseif altitude >= 10000
        ztemp = ((altitude-20000)/((20000-10000)/((306)-(385))))+(306);
    elseif altitude >= 0
        ztemp = ((altitude-10000)/((10000-0)/(385-(462))))+(385);
    end
    % ---------------------------
    % ... convert temperature from celcius to kelvin (K)
    ztemp = (ztemp + 273.15);

end