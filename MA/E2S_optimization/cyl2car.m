function coord_car = cyl2car(coord_cyl)
%
%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Conversion from cylindrical to cartesian coordinates
%
%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% coord_cyl = [r,phi,z] - [3,1] vector with the cylindrical coordinates
%
%%% OUTPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% coord_car = [x,y,z] - [3,1] vector with the cartesian coordinates
%
% -------------------------------------------------------------------------

r = coord_cyl(1);
phi = coord_cyl(2);
z = coord_cyl(3);

x = r*cos(phi);
y = r*sin(phi);

coord_car = [x; y; z];

end