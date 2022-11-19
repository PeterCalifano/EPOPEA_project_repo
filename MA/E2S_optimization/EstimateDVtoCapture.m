function [dV_capture] = EstimateDVtoCapture(Vinf_entry, mu_main, Ra_target, Rp_target)
%% PROTOTYPE
% [DV_capture] = EstimateDVtoCapture(Vinf_entry, mu_main, Ra_target, Rp_target)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function estimating the Impulsive DeltaV required to achieve capture
% orbit specified by Apoapsis and Periapsis radii from hyperbolic entry at
% Vinf in the SOI of the main body. Two body problem assumed: be sure it is
% valid!
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% Vinf_entry: [scalar] Excess speed with respect to main attractor
% mu_main: [scalar] Gravitational Parameter of the main attractor (TBP model)
% Ra_target: [scalar] Apoapsis radius of capture orbit
% Rp_target: [scalar] Periapsis radius of capture orbit
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dV_capture: [scalar] Impulsive DeltaV at hyperbola pericentre to achieve
%                      specified capture orbit.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 01/11/2022 - Pietro Califano - Coded
% 17/11/2022 - Pietro Califano - Time computation deleted
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades


%% Function code
% Hypothesis: Planar orbits, entry hyperbolic leg shaped at will

e_target = (Ra_target - Rp_target)/(Ra_target + Rp_target);

% Velocity of hyperbolic orbit at Rp target
Vp_hyp = sqrt(Vinf_entry^2 + 2*mu_main/(Rp_target)); 

% Velocity at pericentre of target orbit
V_cap = sqrt(mu_main*(1 + e_target)/Rp_target);

% DV required to achieve elliptical capture orbit
dV_capture = Vp_hyp - V_cap;


end