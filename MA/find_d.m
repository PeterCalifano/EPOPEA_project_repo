function [ error_TOF ] = find_d( d , a , b , c , theta_i , theta_f , rf , gam_f , rate_theta_f , TOF , n_integrator)

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

% defining theta vector for the integration
t = linspace(theta_i,theta_f,n_integrator+1);

for i = 1:n_integrator
    tm(i) = 0.5*(t(i+1)+t(i));
end
% calculation r at each theta and theta_m
r = 1./(a + b*t + c*t.^2 + d*t.^3 + e*t.^4 + f*t.^5 + g*t.^6) ;
rm = 1./(a + b*tm + c*tm.^2 + d*tm.^3 + e*tm.^4 + f*tm.^5 + g*tm.^6) ;

% calculation of dt at each theta
dTOF  = ( abs( r.^4  .* (1./r  + 2*c + 6*d*t  + 12*e*t.^2  + 20*f*t.^3  + 30*g*t.^4))  ).^0.5 ;
dTOFm = ( abs( rm.^4 .* (1./rm + 2*c + 6*d*tm + 12*e*tm.^2 + 20*f*tm.^3 + 30*g*tm.^4)) ).^0.5 ;

% step size
h = t(2)-t(1);

% Cavalieri method 
I = 0;
for i = 2:n_integrator+1
    I = I + dTOF(i-1) + 4*dTOFm(i-1) + dTOF(i);
end

I = I*h/6;

% computation of TOF error (residual)

error_TOF =TOF-I;

   
end

