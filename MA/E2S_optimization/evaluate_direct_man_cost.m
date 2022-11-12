function DIRECT_DELTA_V = evaluate_direct_man_cost(x, n_dep_planet, n_arr_planet)
%% PROTOTYPE
% DIRECT_DELTA_V = evaluate_direct_man_cost(x, n_dep_planet, n_arr_planet)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Evaluates the DV cost of the a Direct Manoeuvre using a Lambert arc from
% the departure planet to the arrive one, specified by their number
% (counted from the Sun outward). Departure and arrive dates must be given
% as MJD2000 in the array x
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% x: [1x2] array containing departure (1) and arrive (2) dates in MJD2000 format
% n_dep_planet: [integer] number of the departure planet
% n_arr_planet: [integer] number of the arrive planet
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% DIRECT_DELTA_V: [scalar] total cost of the direct transfer [Km/s]
% -------------------------------------------------------------------------------------------------------------
%% CONTRIBUTORS
% Gennaro Rizzo, Pietro Califano
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% V1: documentation created for already existing code, 26/12/2021
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% from_j20002carthesian()
% astroConstants()
% lambertMR()
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades


%% Function code
mu_Sun = astroConstants(4);
departure_date = x(1);
arrive_date = x(2);
TOF = (arrive_date - departure_date)*24*3600; % In seconds

if TOF > 0
    % Direct orbit
    %% Departure
    [RR_DEP, VV_DEP] = from_j20002carthesian(departure_date, n_dep_planet);
    %% Arrive
    [RR_ARR, VV_ARR] = from_j20002carthesian(arrive_date, n_arr_planet);

    %% Direct leg Lambert
    [~,~,~,~,VI,VF,~,~] = lambertMR(RR_DEP, RR_ARR, TOF, mu_Sun, 0,0,0,0);

    %% Cost computation
    DeltaV_dep_dir = norm(VI' - VV_DEP);
    DeltaV_arr_dir = norm(VF' - VV_ARR);
    DIRECT_DELTA_V = DeltaV_arr_dir + DeltaV_dep_dir;

else
    warning('Departure date cannot be later than arrive one')
    DIRECT_DELTA_V = nan;
end

end