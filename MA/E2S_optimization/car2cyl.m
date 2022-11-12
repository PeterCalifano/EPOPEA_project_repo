function coord_cyl = car2cyl(coord_car)
%
%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Conversion from cartesian to cylindrical coordinates
%
%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% coord_car = [x,y,z] - [3,1] vector with the cartesian coordinates
%
%%% OUTPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% coord_cyl = [r,phi,z] - [3,1] vector with the cylindrical coordinates
%
%-------------------------------------------------------------------------


x = coord_car(1);
y = coord_car(2);
z = coord_car(3);

r = sqrt(x^2+y^2);

if x == 0
    
    if y == 0
        error('Singularity in the conversion (x = y = 0)')
    else
        phi = pi/2 * sign(y);
    end

elseif x > 0

    if y >= 0 
        phi = atan(y/x);
    else
        phi = atan(y/x) + 2*pi;
    end

else

    phi = atan(y/x) + pi;

end

coord_cyl = [r;phi;z];

end