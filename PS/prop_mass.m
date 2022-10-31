function [mp,mp_v,m0_v] = prop_mass(dV_v,m0)
% Computes total propellant mass required for the entire dV budget
% and variation of s/c mass for each maneuvre
% INPUTS
%       dV[1xn] = array containing the dVi required for all maneuvres in
%       order of time at which impulsive maneuvre is applied [m/s] 
%       Margins in dV must be included (according to CDF margins
%       standard)
%       m0[1x1] = total wet mass at launch                    [kg]
%
% OUTPUTS
%       mp[1x1] = total propellant mass required for the entire sequence 
%       of maneuvres [kg]
%       mp_v[1xn] = propellant mass required for each impulsive manuevre dV
%       [kg]
%       m0_v[1xn] = s/c wet mass along the sequence of maneuvres [kg]
%       

% Constants 
g0 = 9.81;                                          % [m/s^2]
n = length(dV_v);
% Propellant characteristics
Is = 300;   % to change                             % [s]

% Apply Tsiolkowsky equation
m0_v = zeros(1,n);
mf = zeros(1,n);
mp_v = zeros(1,n);
m0_v(1,1) = m0;
for i = 1:n
    mf(1,i) = m0_v(i)/exp(dV_v(i)/(Is*g0));
    mp_v(1,i) = m0_v(1,i) - mf(1,i);
    m0_v(1,i+1) = mf(1,i);
end

mp = sum(mp_v);

end