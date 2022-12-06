function  [F, JF] = land_objfunpp(var, state_i_rot, latlon, par, N)
%
%
%
%
%----------------------------------------------------------------------------

%Initial and final time
t1 = var(end-1);
tN = var(end);

% Point latitude and longitude
lat = deg2rad(latlon(1));
lon = deg2rad(latlon(2));
we = par(7);
Re = par(5);

% From latitudinal to cartesian
r_xz = Re*cos(lat);
Xrot = r_xz*sin(lon);
Yrot = Re*sin(lat);
Zrot = r_xz*cos(lon);
vec_rot = [Xrot; Yrot; Zrot];

% Enceladus rotation
th_e = we*(tN-t1);

% From rotating enceladus to IAU_Enceladus
A_rot2IAU = [cos(th_e)    0     sin(th_e)
                 0        1         0
            -sin(th_e)    0     cos(th_e)];
% desired final state
vec_pp = A_rot2IAU*vec_rot;

% final state
var_N = var(end-12:end-2);
r_N = var_N(1:3);

% minimize distance final point - target
F_pos = abs(r_N(1)-vec_pp(1))+abs(r_N(2)-vec_pp(2))+abs(r_N(3)-vec_pp(3));

% minimize propellant mass
F_mass = -var(end-6);

% weights ?????
w_pos = 0.1;
w_mass = 0.9;

% obj fcn
F = w_pos*F_pos+w_mass*F_mass;

%add derivative
if nargout>1
    JF = [];
end
