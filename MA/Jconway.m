function [ vettore_conway,theta_f ] = Jconway( RI,VI,RF,VF,TOF,N,theta_i,n_integrator )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%% algorithm---------------------------------------------------------------

% radii
ri = norm(RI);
rf = norm(RF);

% velocity modulus
vi = norm(VI);
vf = norm(VF);

% radii versors
RI_vers = RI/ri;
RF_vers = RF/rf;

% velocities versors
VI_vers = VI/vi;
VF_vers = VF/vf;

% initial velocity decomposition
VI_rad = RI_vers*dot(VI,RI_vers);
VI_tra = VI-VI_rad;

% final velocity decomposition
VF_rad = RF_vers*dot(VF,RF_vers);
VF_tra = VF-VF_rad;

% initial and final flight path angles
gam_i = atan(dot(VI,RI_vers)/norm(VI_tra));
gam_f = atan(dot(VF,RF_vers)/norm(VF_tra));

% initial and final angular velocities
rate_theta_i = norm(VI_tra)/ri;
rate_theta_f = norm(VF_tra)/rf;

% determine the angle between RI and RF 
RIvcRFv = cross(RI_vers,RF_vers);
rivDrfv =   dot(RI_vers,RF_vers);

if norm(RIvcRFv) < 1e-10
    if rivDrfv > 0 
        psy = 0;
    else 
        psy = pi;
    end 
else 
    if RIvcRFv(3) > 0 
        psy = acos(rivDrfv);
    else
        psy = 2*pi - acos(rivDrfv);
    end 
end

% compute final theta
theta_f = psy + theta_i + 2*N*pi ;

% compute a,b,c

a = 1/ri;
b = -(tan(gam_i)) / ri ;
c = (1/(2*ri)) * ((1/( (ri^3) * rate_theta_i^2 )) - 1) ;

% solving for d 
fun = @(d) find_d( d , a , b , c , theta_i , theta_f , rf , gam_f , rate_theta_f , TOF , n_integrator);
d = 0;
opt = optimset('Display','off', 'UseParallel', true);
d=fzero(fun,d,opt);
% hhh = 1e-5;
% er = 34567;
% 
% while er> 1e-4
%     
%     dfun = 1.0/hhh * (fun(d+hhh) - fun(d-hhh)) ;
%     do = d;
%     d = do - fun(d) / dfun ;
%     er = abs(d-do);
%     
% end


% calculation of e f g
AAA = [  30*theta_f^2 , -10*theta_f^3 ,    theta_f^4  ;
        -48*theta_f   ,  18*theta_f^2 , -2*theta_f^3  ;
         20           ,  -8*theta_f   ,    theta_f^2 ];

bbb = [ 1/rf - (a + b*theta_f + c*theta_f^2 + d*theta_f^3 );
       -(tan(gam_f))/rf - (b + 2*c*theta_f + 3*d*theta_f^2);
       1/((rf^4)*(rate_theta_f^2)) - (1/rf + 2*c + 6*d*theta_f)];
   
e = 0.5/theta_f^6*AAA(1,:)*bbb;
f = 0.5/theta_f^6*AAA(2,:)*bbb;
g = 0.5/theta_f^6*AAA(3,:)*bbb;

rf_c = 1./(a + b*theta_f + c*theta_f^2 + d*theta_f.^3 + e*theta_f^4 + f*theta_f^5 + g*theta_f^6) ;
gam_fc = atan(-rf_c.*(b + 2*c*theta_f + 3*d*theta_f^2 + 4*e*theta_f^3 + 5*f*theta_f.^4 + 6*g*theta_f^5));
gam_ic = atan(-1/a.*(b ));

if abs(gam_ic - gam_i ) > 0.01
    theta_f = -1;
end
if abs(gam_fc - gam_f ) > 0.01
    theta_f = -1;
end



vettore_conway = [a,b,c,d,e,f,g];
end

