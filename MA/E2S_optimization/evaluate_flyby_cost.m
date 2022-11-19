function [DeltaV_power, rp, turning_angle] = evaluate_flyby_cost(v_inf_minus, v_inf_plus, mu_flyby_planet, rp_guess)
%% PROTOTYPE
% [DeltaV_power, rp, turning_angle] = evaluate_flyby_cost(v_inf_minus, v_inf_plus, mu_flyby_planet, rp_guess)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Evaluates the quantities to characterize a powered flyby, assuming that
% the burn is given at pericentre of both hyperbolic branches coming from 
% the interplanetary legs. Computes DV powered, radius of pericentre where 
% the two branches are patched together, total turning_angle for the flyby.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% v_inf_minus: [3x1] excess velocity of ingoing hyperbola
% v_inf_plus: [3x1] excess velocity of outgoing hyperbola
% mu_flyby_planet: [scalar] gravitational parameter of the main body
% rp_guess: [scalar] initial guess for hyperbola periapsis
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% DeltaV_power: [scalar] module of the DV powered required
% rp: [scalar] radius of the pericentre of the flyby hyperbola
% turning_angle: [scalar] angles through which the asymptotes are rotared,
%                considering the DV powered [rad]
% -------------------------------------------------------------------------------------------------------------
%% CONTRIBUTORS
% Gennaro Rizzo, Pietro Califano
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% V1: documentation created for already existing code, 26/12/2021
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% fzero_proj()
% astroConstants()
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% 1) Modify the code such that if rp is not found the output does not stop
%    the minimization process in which it is used

%% Function code

turning_angle = (acos(dot(v_inf_minus/norm(v_inf_minus), v_inf_plus/norm(v_inf_plus))));
f_delta = @(rp) turning_angle - asin( 1 ./ (1 + rp*(norm(v_inf_minus).^2)/mu_flyby_planet)) - asin( 1 ./ (1 + rp*(norm(v_inf_plus).^2)/mu_flyby_planet));

rp = fzero_proj(f_delta, rp_guess); %fzero_proj = fzero without frpintf

if isnan(rp) || ~isreal(rp) %check if r_p is a real number

    %error('Error: fzero could not find a real solution'); % goes to catch
    DeltaV_power = 1000;

elseif rp < astroConstants(23) + 300

    %error('Error: the flyby is performed under 300 km of Height'); % goes to catch
    DeltaV_power = 1000;

else

    v_p_1 = sqrt((norm(v_inf_minus)^2 + 2*mu_flyby_planet/rp));
    v_p_2 = sqrt((norm(v_inf_plus)^2 + 2*mu_flyby_planet/rp));
    
    DeltaV_power = v_p_2 - v_p_1;

end

end


