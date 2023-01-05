function [dV_capture] = HaloCapture(Vinf_entry, mu_main, Rp_target, V_pericenter)
%% PROTOTYPE
% [DV_capture] = EstimateDVtoCapture(Vinf_entry, mu_main, Ra_target, Rp_target)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function estimating the Impulsive DeltaV required to achieve capture
% orbit specified by the velocity at periapsis radii from hyperbolic entry at
% Vinf in the SOI of the main body. Two body problem assumed: be sure it is
% valid!
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% Vinf_entry: [scalar] Excess speed with respect to main attractor
% mu_main: [scalar] Gravitational Parameter of the main attractor (TBP model)
% Rp_target: [scalar] Target Pericenter radius of orbit
% V_pericenter: [scalar] Target velocity at capture 
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dV_capture: [scalar] Impulsive DeltaV at hyperbola pericentre to achieve
%                      specified capture orbit.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 04/01/2023 - Matteo Lusvarghi - Coded
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades


%% Function code
% Hypothesis: Planar orbits, entry hyperbolic leg shaped at will

% Velocity of hyperbolic orbit at Rp target
Vp_hyp = sqrt(Vinf_entry^2 + 2*mu_main/(Rp_target)); 

% DV required to achieve elliptical capture orbit
dV_capture = Vp_hyp - V_pericenter;


end