function [r,v] = kep2car(a, e, i, OM, om, th, mu)

% kep2car.m - Conversion from Keplerian elements to Cartesian coordinates
%
% PROTOTYPE:
% [r, v] = kep2car(a, e, i, OM, om, th, mu)
%
% DESCRIPTION:
% Conversion from Keplerian elements to Cartesian coordinates. Angles in
% degrees.
%
% INPUT:
% a [1x1] Semi-major axis [km]
% e [1x1] Eccentricity [-]
% i [1x1] Inclination [rad]
% OM [1x1] RAAN [rad]
% om [1x1] Pericentre anomaly [rad]
% th [1x1] True anomaly [rad]
% mu [1x1] Gravitational parameter [km^3/s^2]
%
% OUTPUT:
% r [3x1] Position vector [km]
% v [3x1] Velocity vector [km/s]

% R1=rotazione di OM intorno al versore k
% R2=rotazione di i intorno al versore i'
% R3=rotazione di om intorno al versore k''
% R= prodotto matrici di rotazione (DA TRASPORRE)
%%
if nargin == 2
    mu = e;
    e = a(2);
    i = a(3);
    OM = a(4);
    om = a(5);
    th = a(6);
    a = a(1);
end
%% conversione in radianti
% i=deg2rad(i);
% om=deg2rad(om);
% OM=deg2rad(OM);
% th=deg2rad(th);

%% calcolo p ed |r|
p=a*(1-e^2);
r_nor=p/(1+e*cos(th));

%% calcolo r,v nel sistema perifocale (PF)
r_pf=r_nor.*[cos(th);sin(th);0];
v_pf=sqrt(mu/p)*[-sin(th);e+cos(th);0];

%% calcolo matrici di rotazione per arrivare al PF

R1=[cos(OM), sin(OM), 0; -sin(OM), cos(OM), 0; 0, 0, 1];
R2=[1, 0, 0; 0, cos(i), sin(i); 0, -sin(i), cos(i)];
R3=[cos(om), sin(om), 0; -sin(om), cos(om), 0; 0, 0, 1];
R=R3*R2*R1;

%% calcolo r e v nel sistema geocentrico equatoriale (ECI)

r=R'*r_pf;
v=R'*v_pf;

end