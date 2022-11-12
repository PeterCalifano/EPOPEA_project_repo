function [DeltaV_LEGS, v_inf_minus, v_inf_plus, DeltaV_dep, DeltaV_arr] = evaluate_legs_cost(x, n_dep_planet, n_flyby_planet, n_arr_planet)
%% PROTOTYPE
% [DeltaV_LEGS, v_inf_minus, v_inf_plus, DeltaV_dep, DeltaV_arr] = evaluate_legs_cost(x, n_dep_planet, n_flyby_planet, n_arr_planet)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Evaluates the DV cost of each interplanetary leg from departure planet to
% flyby planet and from flyby planet to arrive planet, required for the
% Patched Conic Method. The excess velocities for the flyby phase inside
% the encounter planet SOI are calculated too.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% x: [1x3] array containing departure (1), flyby (2) and arrive (3) dates in MJD2000 format
% n_dep_planet: [integer] number of the departure planet
% n_flyby_planet: [integer] number of the flyby planet
% n_arr_planet: [integer] number of the arrive planet
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% DeltaV_LEGS: [scalar] sum of the costs of the two legs [Km/s]
% v_inf_minus: [3x1] excess velocity of ingoing hyperbola
% v_inf_plus: [3x1] excess velocity of outgoing hyperbola
% DeltaV_dep: [scalar] deltaV cost of the departure burn
% DeltaV_arr: [scalar] deltaV cost of the arrive burn
% -------------------------------------------------------------------------------------------------------------
%% CONTRIBUTORS
% Gennaro Rizzo, Pietro Califano
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% V1: documentation created for already existing code, 26/12/2021
% V2: added output of each DV component, 01/01/2022
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% astroConstants()
% from_j20002carthesian()
% lambertMR()
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% 1) Add tof_parabolic as output
% 2) Add a feature to reject all DVs that are not physical or feasible
%    (occurs for small TOFs)

%% Function code

mu_Sun = astroConstants(4);
departure = x(1);
flyby = x(2);
arrive = x(3);

tof1 = (flyby - departure)*3600*24;
tof2 = (arrive - flyby)*3600*24;

%% Departure
[RR_DEP,VV_DEP] = from_j20002carthesian(departure, n_dep_planet);
%% FLYBY
[RR_FLYBY,VV_FLYBY] = from_j20002carthesian(flyby, n_flyby_planet);
%% Arrive 
[RR_ARR,VV_ARR] = from_j20002carthesian(arrive, n_arr_planet);

%% Helio Legs and cost computation
[~,~,~, ERROR_1,VI_1,VF_1, t_parabolic, ~] = lambertMR(RR_DEP,RR_FLYBY,tof1,mu_Sun,0,0,0,0);

DeltaV_dep = norm(VI_1' - VV_DEP);
v_inf_minus = VF_1' - VV_FLYBY;

if ERROR_1 ~= 0
    error('1st Leg is not possible');
end

[~,~,~,ERROR_2,VI_2,VF_2,~,~] = lambertMR(RR_FLYBY,RR_ARR,tof2,mu_Sun,0,0,0,0);
v_inf_plus = VI_2' - VV_FLYBY;
DeltaV_arr = norm(VF_2' - VV_ARR);

if ERROR_2 ~= 0
    error('2nd Leg is not possible');
end
 
DeltaV_LEGS =  abs(DeltaV_arr) + abs(DeltaV_dep);
end