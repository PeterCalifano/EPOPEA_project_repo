function [a,e,i,OM,om,th] = car2kep(r,v,mu)
% car2kep.m - Conversion from Cartesian coordinates to Keplerian elements
%
% PROTOTYPE:
% [a, e, i, OM, om, th] = car2kep(r, v, mu)
%
% DESCRIPTION:
% Conversion from Cartesian coordinates to Keplerian elements. Angles in
% degrees.
%
% INPUT:
% r [3x1] Position vector [km]
% v [3x1] Velocity vector [km/s]
% mu [1x1] Gravitational parameter [km^3/s^2]
%
% OUTPUT:
% a [1x1] Semi-major axis [km]
% e [1x1] Eccentricity [-]
% i [1x1] Inclination [rad]
% OM [1x1] RAAN [rad]
% om [1x1] Pericentre anomaly [rad]
% th [1x1] True anomaly [rad]

%% Modulo vettore posizione r e vettore velocità v 

r_nor = norm(r);
v_nor = norm(v);

%% Momento angolare specifico h

h = cross(r,v);
h_nor = norm(h);

%% Inclinazione i
if h_nor == 0
    h = [0;0;1];
    h_nor = 1;
    i = 0;
else
    i = acos(h(3)/h_nor);
end

%% Vettore eccentricità ed eccentricità e

e_vec = 1/mu*((v_nor^2 - mu/r_nor)*r - dot(r,v)*v);
e = norm(e_vec);

%% Energia meccanica specifica E e semiasse maggiore a

E = 0.5*v_nor^2 - mu/r_nor;
a = -mu/(2*E);

%% Linea dei nodi N

K = [0;0;1];
N = cross(K,h);
N_nor = norm(N);

%% Ascensione retta del nodo ascendente OM
if N_nor == 0
    N = [1;0;0];
    N_nor = 1;
    OM = acos(N(1)/N_nor);
else
    if N(2) >= 0
        OM = acos(N(1)/N_nor);
    else
        OM = 2*pi - acos(N(1)/N_nor);
    end
    
end

%% Anomalia del pericentro om

if e_vec(3) >= 0
    om = acos(dot(N,e_vec)/(N_nor*e));
else
    om = 2*pi - acos(dot(N,e_vec)/(N_nor*e));
end


%% Velocità radiale vr

vr = dot(r,v)/r_nor;

% se vr > 0 ci stiamo allontanando dal pericentro
% se vr < 0 ci stiamo avvicinando dal pericentro

%% Anomalia vera th

if vr >= 0
    th = acos(dot(e_vec,r)/(e*r_nor));
else
    th = 2*pi - acos(dot(e_vec,r)/(e*r_nor));
end


end

