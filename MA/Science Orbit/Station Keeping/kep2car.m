function [r,v] = kep2car(Kep,mu)
% PROTOTYPE: [r,v] = kep2car(Kep,mu)
%
% DESCRIPTION: Conversion from Keplerian elements to Cartesian coordinates. Angles in 
% degrees
%
% INPUT:
%          Kep    [1x6] orbital elements:
%                                          a      [1x1] Semi major axis [km]
%                                          emod   [1x1] Eccentricity [-]
%                                          i      [1x1] Inclination [deg]
%                                          Omega  [1x1] RAAN [deg]
%                                          w      [1x1] Pericentre anomaly [deg]
%                                          f  [1x1] True anomaly [deg]
%          mu     [1x1] Gravitational parameter [km^3/s^2]
%
% OUTPUT:
%          r [3x1] Position vector [km]
%          v [3x1] Velocity vector [km/s]
%
% CONTRIBUTORS: Matteo Menessini
%               Roberta Pecchioli
%               Elena Pilo
%               Fabrizio Maccari
%
% VERSIONS: 2021-11-12: Second version


a = Kep(1);
e = Kep(2);
i = Kep(3);
Omega = Kep(4);
w = Kep(5);
f = Kep(6);

R = Rmatrix(i,Omega,w);

t = deg2rad(f);

p = a*(1-e^2);
rmod = p/(1+e*cos(t));

%% Position and velocity vectors in PQW

r_pqw = rmod.*[cos(t);sin(t);0];
v_pqw = sqrt(mu/p).*[-sin(t);e+cos(t);0];

%% Position and velocity vectors in ECI

r = R' * r_pqw;
v = R' * v_pqw;

end


function [R] = Rmatrix(i,Omega,w)
% PROTOTYPE: [R] = KEP_CAR(i,Omega,w,deg)
%
% DESCRIPTION: Returns the rotational matrix from ECI [i,j,k] system to PQW [e,p,h].
% Angles in degrees
%
% CONTRIBUTORS: Matteo Menessini
%
% VERSIONS: 2021-10-09: First version


i = i*pi/180;
Omega = Omega*pi/180;
w = w*pi/180;

R1 = [cos(Omega),sin(Omega),0;-sin(Omega),cos(Omega),0;0,0,1];
R2 = [1,0,0;0,cos(i),sin(i);0,-sin(i),cos(i)];
R3 = [cos(w),sin(w),0;-sin(w),cos(w),0;0,0,1];

R = R3*R2*R1;

end
