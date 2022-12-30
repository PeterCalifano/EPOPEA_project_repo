close all
clear
clc

%% Load SPICE
cspice_kclear();
try
    cspice_furnsh('spice_kernels/pck00010.tpc')
    cspice_furnsh('spice_kernels/naif0012.tls')
    cspice_furnsh('spice_kernels/gm_de431.tpc')
    cspice_furnsh('spice_kernels/de440s.bsp')
catch
    cspice_furnsh('..\..\spice_kernels/pck00010.tpc')
    cspice_furnsh('..\..\spice_kernels/naif0012.tls')
    cspice_furnsh('..\..\spice_kernels/gm_de431.tpc')
    cspice_furnsh('..\..\spice_kernels/de440s.bsp')
end

% Normalization coefficients
% TU = 24*3600; % 1 day
% DU = 1.495978707e+8; % 1 AU

%% Load solution workspace
% Baseline
load('E_EEJ_S_bestvalue1.016.mat')

clearvars -except min_at_iter NLPoptset_local planet_seq min_pos
% Main backup
% load('E_EEJ_S_bestvalue1.14backup2035.mat')

% Secondary backup
% load('E_EEJ_S_bestvalue1.02backup2043.mat')

%% Compute trajectory
R_Saturn = astroConstants(26); % [km]
Ra_target = 200*R_Saturn;
Rp_target = 2.55*R_Saturn;

% Extract best solution
[~, index_iter] = min(min_at_iter(min_at_iter > 0));
index_pos = min_pos(index_iter);
Xsol = NLPoptset_local(index_pos, :, index_iter);

% Evaluate 
[DVsum, data] = objfun_EStransfer_analysis(Xsol, planet_seq, Ra_target, Rp_target);

Xplanets = data.Xplanets;
EventsEpochsJD = data.EventsEpochsJD;
vinf_minunsplus = data.vinf_minusplus;
DV_breakdown = data.DVs;

%% Evaluate flybys
howmanyfb = numel(EventEpochsJD);
fbplanets = planet_seq(2:end-1);

R_SOI = zeros(1, howmanyfb);
GM_planets = zeros(1, howmanyfb);
MeanR_planets = zeros(1, howmanyfb);

GM_Sun = cspice_bodvrd('SUN', 'GM', 1);

for idfb = 1:howmanyfb
    % Determine Planet properties
    switch fbplanets
        case 1
            bodynm = 'MERCURY';
        case 2
            bodynm = 'VENUS';
        case 3
            bodynm = 'EARTH';
        case 4
            bodynm = 'MARS';
        case 5
            bodynm = 'JUPITER';
        case 6
            bodynm = 'SATURN';
    end

    GM_planets(idfb) = cspice_bodvrd(bodynm, 'GM', 1);
    MeanR_planets(idfb) = mean(cspice_bodvrd(bodynm, 'RADII', 3));
    R_SOI(idfb) = MeanR_planets(idfb) * (GM_planets(idfb)/GM_Sun)^(2/5);

%     rp = Xsol();
    % [] = EvalFlyBy(vinf_minunsplus(idfb, :), rp, GM_planets(idfb), R_SOI(idfb), Xplanets(idfb, :));

end







