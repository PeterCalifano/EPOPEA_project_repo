function [sys2BP_J2] = odesys_2BP_J2(t, state_vec)
%% PROTOTYPE
% [sys2BP_J2] = odesys_2BP_J2(r, v)
% --------------------------------------------------------------------------
%% INPUT
% --------------------------------------------------------------------------

%% OUTPUT

% --------------------------------------------------------------------------

%% CONTRIBUTORS
% Pietro Califano
% --------------------------------------------------------------------------

%% CHANGELOG
% v1: builds the 2BP ODEs system with J2 term, 31/10/2021
% --------------------------------------------------------------------------


%% Function 
mu = astroConstants(13); % [Km^3/s^2] planetary costant
R_E = astroConstants(23); % [Km] radius of Eath
J2 = astroConstants(9); % for the perturbed two body problem

r_norm = norm(state_vec(1:3));
aJ2 = 1.5 * ((J2*mu*R_E^2)./(r_norm^4)) * [(state_vec(1)/r_norm) * (5*(state_vec(3)/r_norm)^2 - 1);
                                           (state_vec(2)/r_norm) * (5*(state_vec(3)/r_norm)^2 - 1);
                                           (state_vec(3)/r_norm) * (5*(state_vec(3)/r_norm)^2 - 3)];

sys2BP_J2 = [state_vec(4:6);...
             (-mu./r_norm^3)*state_vec(1:3) + aJ2];


end