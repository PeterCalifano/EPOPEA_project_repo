function [v1_t,v2_t] = Lambert(planet1,planet2,t1,t2)
%
% Description:
% The function computes the interplanetary leg through a Lambert solver
% between two planets in the solar system on two specified Julian days.
%
% Prototype:
% [v1_t,v2_t] = Lambert(planet1,planet2,t1,t2)
%
% Inputs:
% planet1 - Initial planet
% planet2 - Final planet
% t1 [1] - Julian Day of departure
% t2 [1] - Julian Day of arrival
%
% Outputs:
% v1_t [3x1] - Velocity at the beginning of the transfer arc
% v2_t [3x1] - Velocity at the end of the transfer arc
% -------------------------------------------------------------------------
%% Computation of initial and final positions and velocities of the planets

    mu = cspice_bodvrd('Sun','GM',1);
    x1 = cspice_spkezr(planet1,t1,'ECLIPJ2000','NONE','SUN');
    r1 = x1(1:3);
    v1 = x1(4:6);
    x2 = cspice_spkezr(planet2,t2,'ECLIPJ2000','NONE','SUN');
    r2 = x2(1:3);
    v2 = x2(4:6);

%% Computation of the lambert arc

    tof = (t2 - t1) * 3600 * 24;
    [a_t,~,e_t,~,v1_t,v2_t,~,~] = lambertMR( r1, r2, tof, mu);

%% Find output quantities
% 
%     DV = norm(v1_t' - v1) + norm(v2_t' - v2);
%     
%     [~,~,i_t,OM_t,om_t,th1_t] = car2kep(r1,v1_t',mu);
%     [~,~,~,~,~,th2_t] = car2kep(r2,v2_t',mu);
%     kep_t_start = [a_t,e_t,i_t,OM_t,om_t,th1_t];
%     kep_t_end = [a_t,e_t,i_t,OM_t,om_t,th2_t];
end