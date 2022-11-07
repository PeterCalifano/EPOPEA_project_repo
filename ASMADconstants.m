function const = ASMADconstants(constName)

%% PROTOTYPE
%   const = ASMADconstants(constName)
%
%% DESCRIPTION
%   Function used to retrieve any constants used in codes. Constants are
%   hardcoded for efficiency (no calls to SPICE)
%
%% INPUT
% constName [string] Name of constant to retrieve from list (case-sensitive):
%       'AU'                         Astronomical Unit [km]
%       'SolarIrradiance_1AU'        Solar irradiance @ 1 AU [W/m^2]
%       
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% const [1x1] Constant that is retrieved
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% Date, User, brief summary of the modification
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% More constants to be added

%% Function code

switch constName

    case 'AU'
        const = 149597870.613689 ; % From SPICE

    case 'SolarIrradiance_1AU'
        const = 1360 ;

    otherwise

        error( 'Constant not found, check name or add constant to this file' ) ;

end


end