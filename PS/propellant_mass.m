function [ m_prop_orb, m_prop_lan, m_prop_main, m_prop_sk, m_prop_att, m_prop_main_lan, m_prop_att_lan ] = propellant_mass( orbiter, lander )

% INPUTS
% orbiter  -  structure:        
        % dV_disp                [m/s] - Delta v for disposal of the orbiter 
        % dV_sk_o                [m/s] - Delta v for sk when only orbiter
        % dV_att_o               [m/s] - Delta v for attitude maneuvers when only orbiter
        % dV_sk_cl               [m/s] - Delta v for sk when clamped configuration
        % dV_EOI                 [m/s] - Delta v for moon tour and orbit insertion
        % dV_cap                 [m/s] - Delta v for capture at Saturn
        % dV_int                 [m/s] - Delta v for interplanetary leg
        % Isp_main               [s]   - Specific impulse of primary engines
        % Isp_sk                 [s]   - Specific impulse of sk
        % Isp_att                [s]   - Specific impulse of attitude thrusters 

        % m_att_moon_cl          [kg] - Propellant mass for attitude maneuvers during moon tour
        % m_prop_att_science_cl  [kg] - Propellant mass for attitude maneuvers when clamped configuration (science orbit)
        % m_prop_att_science_o   [kg] - Propellant mass for attitude maneuvers when only orbiter (science orbit)
        % m_prop_att_dsm_cl      [kg] - Propellant mass for attitude maneuvers when clamped configuration (deep space maneuver)
        % dry                    [kg] - Dry mass

% lander  -  structure:
        % dV_des             [m/s] - Delta v for descending phase 
        % dV_att             [m/s] - Delta v for attitude maneuvers
        % Isp_main           [s]   - Specific impulse of primary engines (LANDER)
        % Isp_att            [s]   - Specific impulse of attitude thrusters (LANDER)
        % dry                [kg]  - Dry mass
        % hazard             [kg]  - Propellant mass needed for hazard maneuvers (HAZARD MARGIN)

% OUTPUTS - masses already include 2.5% MARGIN

% m_prop_orb         [kg]  - Propellant mass required by orbiter
% m_prop_lan         [kg]  - Propellant mass required by lander - already includes HAZARD MARGIN

% -----

%% Recover data from structures
dV_lan = lander.dV_des ;
dV_att_lan = lander.dV_att ;
Isp_lan = lander.Isp_main ;
Isp_att_lan = lander.Isp_att ;

dV_disp = orbiter.dV_disp ;
dV_sk_orb = orbiter.dV_sk_o ;

dV_sk_cl = orbiter.dV_sk_cl ;
dV_EOI = orbiter.dV_EOI ;
dV_cap = orbiter.dV_cap ;
dV_int = orbiter.dV_int ;
Isp_main_orb = orbiter.Isp_main ;
Isp_sk_orb = orbiter.Isp_sk ;
% Isp_att_orb = orbiter.Isp_att ;

m_dry_orb = orbiter.dry ;
m_dry_lan = lander.dry ;
m_hazard = lander.hazard ;

m_prop_att_science_o = orbiter.m_att_science_o ;
m_prop_att_science_cl = orbiter.m_att_science_cl ;
m_prop_att_dsm_cl = orbiter.m_att_dsm_cl ;
m_prop_att_moon = orbiter.m_att_moon_cl ;

g0 = 9.81 ; % [m/s] - gravitational acceleration

%% LANDER

% Descending
mf = m_dry_lan + m_hazard ;
m_prop_des = mf * ( exp( dV_lan / ( Isp_lan * g0 ) ) - 1 ) ; % [kg] - propellant mass needed for lander descent

% Attitude
mf = m_dry_lan + m_hazard ;
m_prop_att_lan = mf * ( exp( dV_att_lan / ( Isp_att_lan * g0 ) ) - 1 ) ; % [kg] - propellant mass needed for attitude maneuvers during landing

m_prop_lan = m_prop_att_lan + m_prop_des + m_hazard ; % [kg] Total propellant mass on lander

%% ORBITER 

% Disposal
mf = m_dry_orb ;
m_prop_disp = mf * ( exp( dV_disp / ( Isp_main_orb * g0 ) ) - 1 ) ; % [kg] - propellant mass needed for disposal of the orbiter

% SK and Attitude - science orbit
% 
% mf = m_dry_orb + m_prop_disp + m_prop_att_science_o ;

mf = m_dry_orb + m_prop_disp ;
m_prop_sk_orb = mf * ( exp( dV_sk_orb / ( Isp_sk_orb * g0 ) ) - 1 ) ; % [kg] - propellant mass needed for sk of orbiter

% m_prop_temp = mf * ( exp( dV_sk_orb / ( Isp_sk_orb * g0 ) ) - 1 ) ; % [kg] - propellant mass needed for sk of orbiter
% m_prop_sk_orb = m_prop_temp - m_prop_att_science_o ; % [kg] - propellant mass needed for sk of orbiter

%% CLAMPED (ORBITER + LANDER)

% SK - science orbit

% mf = m_dry_orb + m_dry_lan + m_prop_lan + m_prop_sk_orb + m_prop_disp + m_prop_att_science_o + m_prop_att_science_cl ;
mf = m_dry_orb + m_dry_lan + m_prop_lan + m_prop_sk_orb + m_prop_disp + m_prop_att_science_o ;
m_prop_sk_cl = mf * ( exp( dV_sk_cl / ( Isp_sk_orb * g0 ) ) - 1 ) ; % [kg] - propellant mass needed for sk of clamped configuration

% m_prop_temp = mf * ( exp( dV_sk_cl / ( Isp_sk_orb * g0 ) ) - 1 ) ; % [kg] - propellant mass needed for sk of clamped configuration
% m_prop_sk_cl = m_prop_temp - m_prop_att_science_cl ; % [kg] - propellant mass needed for sk of clamped configuration


m_prop_att = m_prop_att_science_cl + m_prop_att_science_o ; % [kg] - total propellant mass needed for attitude control
m_prop_sk = m_prop_sk_orb + m_prop_sk_cl ; % [kg] - total propellant mass needed for sk

% Moon tour and enceladus orbit insertion
mf = m_dry_orb + m_dry_lan + m_prop_lan + m_prop_disp + m_prop_att + m_prop_sk ;
m_prop_EOI = mf * ( exp( dV_EOI / ( Isp_main_orb * g0 ) ) - 1 ) ; % [kg] - propellant mass needed for moon tour and enceladus orbit insertion

% Capture at Saturn
mf = m_dry_orb + m_dry_lan + m_prop_lan + m_prop_disp + m_prop_att + m_prop_sk + m_prop_EOI + m_prop_att_moon ;
m_prop_cap = mf * ( exp( dV_cap / ( Isp_main_orb * g0 ) ) - 1 ) ; % [kg] - propellant mass needed for capture at Saturn
 
% Interplanetary leg
% mf = m_dry_orb + m_dry_lan + m_prop_lan + m_prop_disp + m_prop_att + m_prop_sk + m_prop_EOI + m_prop_cap + m_prop_att_dsm_cl ;
mf = m_dry_orb + m_dry_lan + m_prop_lan + m_prop_disp + m_prop_att + m_prop_sk + m_prop_EOI + m_prop_cap + m_prop_att_moon ;
m_prop_int = mf * ( exp( dV_int / ( Isp_main_orb * g0 ) ) - 1 ) ; % [kg] - propellant mass needed for interplanetary transfer
 
m_prop_main = m_prop_int + m_prop_cap + m_prop_EOI + m_prop_disp ;

m_prop_orb = m_prop_att + m_prop_sk + m_prop_EOI + m_prop_cap + m_prop_int + m_prop_att_dsm_cl + m_prop_att_moon ; % [kg] - total propellant mass needed on the orbiter

%% INCLUDE MARGINS:
m_prop_orb = m_prop_orb * 1.025 ;
m_prop_lan = m_prop_lan * 1.025 ;

m_prop_main = m_prop_main ;
m_prop_sk = m_prop_sk ;
m_prop_att = m_prop_att + m_prop_att_dsm_cl + m_prop_att_moon ;
m_prop_main_lan = (m_prop_des + m_hazard ) ;
m_prop_att_lan = m_prop_att_lan ;

end

