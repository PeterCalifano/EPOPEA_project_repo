function Vinfout_IRF = FindVinfOut(vinfinvec_i, beta_i, Rp_i, mu_pl, Vpl)
% Vinfout = v infinity of outbound leg 

% vinfinvec_i = v infinity of inbound leg (vectorial)
% beta_i = angle of hyperbolic plane
% Rp_i = pericentre radius of the hyperbola
% mu_pl = grav. parameter of planet

vinfin_i = norm(vinfinvec_i);

delta = 2*asin(1/(1 + (Rp_i*vinfin_i.^2)./mu_pl));

if delta > rad2deg(150)
    warning('Computed turning angle is very high!')
end

% Reference frame:
I = vinfinvec_i./vinfin_i;
Jvec = cross(I, Vpl);
J = Jvec./norm(Jvec);
K = cross(I, J);

Vinfout = vinfin_i.*[cos(delta); % I 
   cos(beta_i)*sin(delta); % J
   sin(beta_i)*sin(delta)]; % K

RotMatrixLocal2IRF = [I, J, K]; 

% Rotate the velocity from the Local "flyby" RF to the Inertial RF
Vinfout_IRF = RotMatrixLocal2IRF * Vinfout;


end