function [l, time, T, m, r, flag] = lambLTj_plan(RI, RF, VI, VF, TOF, N, M_end, n_sol, n_integrator, Is)
%% PROTOTYPE
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% Date, User, brief summary of the modification
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades


%% Function code
flag = 0;

DU = 149597870.7;      % length of 1 DU[km]
muuu = 132712440018;    % actractor parameter of the sun[km3/s2]
TU = (DU^3/muuu )^0.5;  % duration of 1 TU[s]

TOF = TOF/TU; % Normalized time;

RCRRv  = cross(RI,RF)/norm(cross(RI,RF)) ;


if RCRRv(3) < 0

    RCRRv = -RCRRv;

end


% decompongo velocita' nel piano e fuori dal piano

v_perp_i = dot(VI, RCRRv);
V_PLAN_I = VI-v_perp_i*RCRRv;
v_plan_i = norm(V_PLAN_I);

v_perp_f = dot(VF,RCRRv);
V_PLAN_F = VF-v_perp_f*RCRRv;
v_plan_f = norm(V_PLAN_F);


[vettore_conway, L] = Jconway( RI,V_PLAN_I,RF,V_PLAN_F,TOF,N,0 ,n_integrator );

a = vettore_conway(1);
b = vettore_conway(2);
c = vettore_conway(3);
d = vettore_conway(4);
e = vettore_conway(5);
f = vettore_conway(6);
g = vettore_conway(7);

if isnan(a*b*c*d*e*f*g) 
    flag = 1;
end

% generating theta vector for the solution
l = linspace(0, L, n_sol);
% radius at each theta
r = 1./(a + b*l + c*l.^2 + d*l.^3 + e*l.^4 + f*l.^5 + g*l.^6) ;
% angular velocity
rate_theta = ((1./r.^4)./(1./r + 2*c + 6*d*l + 12*e*l.^2 + 20*f*l.^3 + 30*g*l.^4)).^0.5;

% flight path angle at each theta
gam = atan(-r.*(b + 2*c*l + 3*d*l.^2 + 4*e*l.^3 + 5*f*l.^4 + 6*g*l.^5));

% adimensional and dimensional acceleration
acc = -1./(2*r.^3.*cos(gam)).*(6*d+24*e*l+60*f*l.^2+120*g*l.^3-tan(gam)./r)./(1./r + 2*c + 6*d*l + 12*e*l.^2 + 20*f*l.^3 + 30*g*l.^4).^2;

% time vector
d_time = 1./rate_theta;
dt = l(2)-l(1);
time = zeros(n_sol,1);

for i = 1:n_sol-1
    time(i+1) = time(i) + d_time(i)*dt;
end

% thrust and mass profile for TOF
m = zeros(n_sol,1);
T = zeros(n_sol,1);
rate_m_dim = zeros(n_sol,1);
m(end) = M_end;

for i = n_sol:-1:2

    T(i) = acc(i)*m(i) * 1000* DU/TU^2;
    rate_m_dim(i) = abs(T(i))/(Is*9.81);
    m(i-1) = m(i) + rate_m_dim(i)*(time(i)-time(i-1))*TU;

end

T(1) = acc(1)*m(1) * 1000* DU/TU^2;


end

