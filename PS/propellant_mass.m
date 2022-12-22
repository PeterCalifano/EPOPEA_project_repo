function [ m_prop_prim, m_prop_sec, m_prop_lan ] = propellant_mass( Isp_prim, dV_prim, dV_disp, Isp_sec, dV_sec_c, dV_sec_orb, m_dry_orb, Isp_lan, Isp_att, dV_lan, dV_att_lan, m_dry_lan, m_hazard )

% INPUTS

% Isp_prim           [s]   - Specific impulse of primary engines (ORBITER)
% dV_prim            [m/s] - Delta v for primary maneuvers (ORBITER)
% dV_disp            [m/s] - Delta v for disposal (ORBITER)
% Isp_sec            [s]   - Specific impulse of secondary engines (ORBITER)
% dV_sec_c           [m/s] - Delta v for secondary maneuvers when clamped configuration is considered (ORBITER)
% dV_sec_orb         [m/s] - Delta v for secondary maneuvers when only orbiter is considere (ORBITER)
% m_dry_orb          [kg]  - ORBITER dry mass
% Isp_lan            [s]   - Specific impulse of primary engines (LANDER)
% Isp_att            [s]   - Specific impulse of attitude thrusters (LANDER)
% dV_lan             [m/s] - Delta v for landing (LANDER)
% dV_att             [m/s] - Delta v for attitude maneuvers (LANDER)
% m_dry_lan          [kg]  - LANDER dry mass
% m_prop_lan         [kg]  - LANDER propellant mass neede in case of harazd maneuver (MAR_hazard)

% OUTPUTS - masses already include 20% MARGIN

% m_prop_prim        [kg]  - Propellant mass required for primary maneuvers
% m_prop_sec         [kg]  - Propellant mass required for secondary maneuvers
% m_prop_lan         [kg]  - Propellant mass required by lander - already includes HAZARD MARGIN
% -----

g0 = 9.81 ; % [m/s] - gravitational acceleration

%% LANDER PROPELLANT MASS
mf = m_dry_lan + m_hazard ;
mi = mf * exp( dV_lan / ( Isp_lan * g0 ) ) ;
m_prop_des = mi - mf ; % [kg] - propellant mass needed for lander descent

mf = m_dry_lan + m_hazard ;
mi = mf * exp( dV_att_lan / ( Isp_att * g0 ) ) ;
m_prop_att =  mi - mf ; % [kg] - propellant mass needed for attitude maneuvers during landing

m_prop_lan = m_prop_att + m_prop_des + m_hazard ;

%% ORBITER PROPELLANT MASS

% Disposal - only orbiter
mf = m_dry_orb ;
mi = mf * exp( dV_disp / ( Isp_prim * g0 ) ) ;
m_prop_disp = mi - mf ; % [kg] - propellant mass needed for SK when only orbiter is considered

% SK - only orbiter
mf = m_dry_orb + m_prop_disp ;
mi = mf * exp( dV_sec_orb / ( Isp_sec * g0 ) ) ;
m_prop_sec_orb = mi - mf ; % [kg] - propellant mass needed for SK when only orbiter is considered

% SK - clamped configuration
mf = m_dry_orb + m_dry_lan + m_prop_lan + m_prop_sec_orb + m_prop_disp ;
mi = mf * exp( dV_sec_c / ( Isp_sec * g0 ) ) ;
m_prop_sec_c = mi - mf ; % [kg] - propellant mass needed for SK when only orbiter is considered


m_prop_sec = m_prop_sec_c + m_prop_sec_c ; % [kg] - propellant mass needed for secondary maneuvers

% Interplanetary transfer - clamped configuration
mf = m_dry_orb + m_dry_lan + m_prop_lan + m_prop_sec + m_prop_disp ;
mi = mf * exp( dV_prim / ( Isp_prim * g0 ) ) ;
m_prop_int = mi - mf ; % [kg] - propellant mass needed for primary maneuvers

m_prop_prim = m_prop_int + m_prop_disp ;

% INCLUDE MARGINS:
m_prop_prim = m_prop_prim * 1.2 ;
m_prop_sec = m_prop_sec * 1.2 ;
m_prop_lan = m_prop_lan * 1.2 ;

end

