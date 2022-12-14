function dxdt = SCR3BP2_dyn(t,x,mu_tbp,mu_v,R_v,J2_v)
% Shape-based formulation of Circular Restricted Three Body Problem, that
% means taking into account also the j2 effects of both Saturn and
% Enceladus. Plus the acceleration of the Sun is considered.

% Extract constants of the two bodies
mu_S = mu_v(1);
mu_E = mu_v(2);
R_Saturn = R_v(1);
R_Enceladus = R_v(2);
J2_S = J2_v(1);
J2_E = J2_v(2);

% Extract radius from the state
r = x(1:3);

R_RF2IN = [cos(t),-sin(t), 0;
    sin(t), cos(t), 0;
    0 0 1];

r1 = [r(1) + mu_tbp;
    r(2);
    r(3)];
r2 = [r(1) + mu_tbp - 1;
    r(2);
    r(3)];

r1_in = R_RF2IN*r1;
r2_in = R_RF2IN*r2;

%r1_norm = ((x(1)+mu_tbp).^2+x(2).^2+x(3).^2).^0.5;
%r2_norm = ((x(1)+mu_tbp-1).^2+x(2).^2+x(3).^2).^0.5;
r1_norm = norm(r1);
r2_norm = norm(r2);

aJ2_S = 1.5 * J2_S * mu_S * R_Saturn ^ 2 / r1_norm^4 * ...
    [r1_in(1)/r1_norm * (5*r1_in(3)^2 / r1_norm^2 - 1);...
     r1_in(2)/r1_norm * (5*r1_in(3)^2 / r1_norm^2 - 1); ...
     r1_in(3)/r1_norm * (5*r1_in(3)^2 / r1_norm^2 - 3)];

aJ2_E = 1.5 * J2_E * mu_E * R_Enceladus ^ 2 / r2_norm^4 * ...
    [r2_in(1)/r2_norm * (5*r2_in(3)^2 / r2_norm^2 - 1);...
     r2_in(2)/r2_norm * (5*r2_in(3)^2 / r2_norm^2 - 1); ...
     r2_in(3)/r2_norm * (5*r2_in(3)^2 / r2_norm^2 - 3)];


dxdt= [x(4);
    x(5);
    x(6);
    2*x(5)+x(1)-(1-mu_tbp)*(x(1)+mu_tbp)/r1_norm.^3 - mu_tbp/(r2_norm).^3*(x(1)+mu_tbp-1);
    -2*x(4)+x(2)-(1-mu_tbp)*x(2)/r1_norm.^3-mu_tbp/(r2_norm).^3*x(2);
    -(1-mu_tbp)/(r1_norm).^3*x(3)-mu_tbp/(r2_norm.^3)*x(3)];  %CRTBP dynamics

dxdt(4:6) = dxdt(4:6) + R_RF2IN'*aJ2_E + R_RF2IN'*aJ2_S;
a = 1;
end