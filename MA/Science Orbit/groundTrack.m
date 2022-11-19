function [alpha, delta, lat, lon] = groundTrack(T, Y, theta0_g, w_p)
% PROTOTYPE: [alpha, delta, lat, lon] = groundTrack(T, Y, theta0_g, w_p)
%
% DESCRIPTION: Determines the ground track of a S/C orbiting Earth
%
% INPUT:
%          T             [Nx1] Time vector
%          Y             [Nx6] Evolution of the state (cartesian coordinates)
%          theta0_g      Longitude of the Greenwich meridian at t0 [deg]
%          w_p           Angular velocity of the planet (default: Earth)
%
% OUTPUT:
%          alpha    Right Ascension [deg]
%          delta    Declination [deg]
%          lat      Latitude [deg]
%          lon      Longitude[deg]
%
% CONTRIBUTORS: Matteo Menessini
%               Roberta Pecchioli
%               Elena Pilo
%               Fabrizio Maccari
%
% VERSIONS: 2021-10-10: First version

% Initial Setup
if nargin == 3
    w_p = deg2rad(15.04)/3600;
end

t0 = T(1);
theta_g = @(t) theta0_g + w_p*(t-t0);

delta = [];
alpha = [];
lon = [];

for t = 1:length(T)
    x = Y(t,1);
    y = Y(t,2);
    z = Y(t,3);
    r_norm = norm(Y(t,1:3));
    
    % Declination
    D = asin(z/r_norm);
    delta = [delta rad2deg(D)];

    % Right Ascension
    A = acos(x/r_norm/cos(D));
    
    if y/r_norm<=0
        A = 2*pi - A;
    end
    alpha = [alpha rad2deg(A)];

    % Longitude
    L = A - theta_g(t*T(end)/length(T));
    
    lon = [lon rad2deg(wrapToPi(L))];
end

lat = delta;

end