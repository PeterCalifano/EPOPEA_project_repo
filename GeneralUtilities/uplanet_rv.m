function [ rr, vv ] = uplanet_rv( TSPAN_mjd2000, planet )

% Calculate analytical ephemeris for planets
% Function is similar to uplanet.m, but output is already in cartesian coordinates
%
% PROTOTYPE:
%   [ rr, vv ] = uplanet_rv( TSPAN_mjd2000, planet )
%
% INPUT:
%   TSPAN_mjd2000[1x1 or 1xTSPAN]     MJD2000 time: Can be a scalar, or a vector of times
%                                                   TSPAN_mjd2000 = [T0, ..., TF]
%   planet                      String containing name of planet for ephemeris
%                                         1:   'Mercury'
%                                         2:   'Venus'
%                                         3:   'Earth'
%                                         4:   'Mars'
%                                         5:   'Jupiter'
%                                         6:   'Saturn'
%                                         7:   'Uranus'
%                                         8:   'Neptune'
%                                         9:   'Pluto'
%                                         10:  'Sun'
% OUTPUT:
%   rr[3x1 or TSPANx3]                Position vector of planet in HECI frame (heliocentric)
%   vv[3x1 or TSPANx3]                Velocity of planet in HECI frame (heliocentric)

len = length(TSPAN_mjd2000);

if len == 1

    [Kep, mu_sun] = uplanet ( TSPAN_mjd2000, planet );
    flags.generateMatrix = false;
    [ rr, vv, ~ ] = kp2rv( Kep, mu_sun, flags );

else

    Kep = zeros(len,6);

    for k = 1:len

        [Kep(k,:), mu_sun] = uplanet ( TSPAN_mjd2000(k), planet );

    end

    flags.generateMatrix = false;
    [ rr, vv, ~ ] = kp2rv( Kep, mu_sun, flags );

end

end