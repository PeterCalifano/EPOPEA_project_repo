function [v1_t,v2_t] = Lambert(id_1,id_2,t1,t2)
%
% Description:
% The function computes the interplanetary leg through a Lambert solver
% between two planets in the solar system on two specified Julian days.
%
% Prototype:
% [v1_t,v2_t] = Lambert(planet1,planet2,t1,t2)
%
% Inputs:
% id_1 - Identifier of initial planet
% id_2 - Identifier of final planet
% t1 [1] - Time of departure [SECONDS]
% t2 [1] - Time of arrival [SECONDS]
%
% Outputs:
% v1_t [3x1] - Velocity at the beginning of the transfer arc
% v2_t [3x1] - Velocity at the end of the transfer arc
% -------------------------------------------------------------------------
%% Computation of initial and final positions and velocities of the planets
    jd_1 = t1/(3600*24);
    jd_2 = t2/(3600*24);
    mu = astroConstants(4);
    [kep1,~] = uplanet(jd_1, id_1);
    [kep2,~] = uplanet(jd_2, id_2);
    [r1,v1] = kep2car(kep1, mu);
    [r2,v2] = kep2car(kep2, mu);

%% Computation of the lambert arc

    tof = (t2 - t1);
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