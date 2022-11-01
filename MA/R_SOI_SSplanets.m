function [Planet_R_SOI] = R_SOI_SSplanets(planetid)
%% PROTOTYPE
% [Planet_R_SOI] = SOI_radius(planetid)
% -------------------------------------------------------------------------
%% DESCRIPTION
% Function computing the approximate SOI radius of the SS planets identified by
% "planetid". CRTBP model with SUN as main attractor is assumed.
% -------------------------------------------------------------------------
%% INPUT
% planetid: [integer] from 1 to 8 (Mercury to Neptun)
% -------------------------------------------------------------------------
%% OUTPUT
% Planet_R_SOI
% -------------------------------------------------------------------------
%% CHANGELOG
% 01/11/2022 - Pietro Califano - Function recovered from OM codes
% -------------------------------------------------------------------------

% PLANETS LEGEND:
% Integer number identifying the celestial body (< 11)
%  1:   Mercury
%  2:   Venus
%  3:   Earth
%  4:   Mars
%  5:   Jupiter
%  6:   Saturn
%  7:   Uranus
%  8:   Neptun


switch planetid

    case 1
        rSOI_over_r_P = 46.1;
        Planet_R_SOI = rSOI_over_r_P * astroConstants(11);

    case 2
        rSOI_over_r_P = 101.7;
        Planet_R_SOI = rSOI_over_r_P * astroConstants(12);

    case 3
        rSOI_over_r_P = 145.3;
        Planet_R_SOI = rSOI_over_r_P * astroConstants(13);
       
    case 4
        rSOI_over_r_P = 170.0;
        Planet_R_SOI = rSOI_over_r_P * astroConstants(14);

    case 5
        rSOI_over_r_P = 675.1;
        Planet_R_SOI = rSOI_over_r_P * astroConstants(15);

    case 6
        rSOI_over_r_P = 906.9;
        Planet_R_SOI = rSOI_over_r_P * astroConstants(16);

    case 7
        rSOI_over_r_P = 2024.6;
        Planet_R_SOI = rSOI_over_r_P * astroConstants(17);

    case 8
        rSOI_over_r_P = 3494.8;
        Planet_R_SOI = rSOI_over_r_P * astroConstants(18);
end

end