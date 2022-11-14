% MANEUVERING TIME

%% S.O & S.L
clear all, close all, clc

m_dry_tot = 1279 ;
m_dry_land = 746 ;
m_dry_orb = 533 ;

mdot_bip = 0.135 ;

% Interplanetary transfer
Dv_int = 1000 ; % [m/s]
Isp_bip = 321 ;
m_prop_int = preliminary_prop_mass(Dv_int,m_dry_tot,Isp_bip) ;
dt_int = m_prop_int / mdot_bip ;

% Capture at Saturn
Dv_Sat = 900 ;
m_prop_Sat = preliminary_prop_mass(Dv_Sat,m_dry_tot,Isp_bip) ;
dt_Sat = m_prop_Sat / mdot_bip ;

% Moon tour and Enceladus orbit insertion
Dv_EOI = 2200 ;
m_prop_EOI = preliminary_prop_mass(Dv_EOI,m_dry_tot,Isp_bip) ;
dt_EOI = m_prop_EOI / mdot_bip ;

% Landing
Dv_land = 514 ;
mdot_land = 0.03925 ;
Isp_land = 228 ;
m_prop_land = preliminary_prop_mass(Dv_land,m_dry_land,Isp_land) ;
dt_land = m_prop_land / mdot_land ;

% Station keeping
Dv_SK = 370 ;
mdot_SK = 0.00905 ; % [kg/s]
Isp_SK = 235 ;
m_prop_SK_orb_land = preliminary_prop_mass(Dv_SK,m_dry_tot,Isp_SK) ;
dt_SK_orb_land = m_prop_SK_orb_land / mdot_SK ;

m_prop_SK_orb = preliminary_prop_mass(Dv_SK,m_dry_orb,Isp_SK) ;
dt_SK_orb = m_prop_SK_orb / mdot_SK ;

dt_SK_tot = dt_SK_orb + dt_SK_orb_land ;


%% NS.O & S.L
clear all, close all, clc

m_dry_tot = 928 ;
m_dry_land = 746 ;
m_dry_orb = 182 ;

mdot_bip = 0.135 ;

% Interplanetary transfer
Dv_int = 1000 ; % [m/s]
Isp_bip = 321 ;
m_prop_int = preliminary_prop_mass(Dv_int,m_dry_tot,Isp_bip) ;
dt_int = m_prop_int / mdot_bip ;

% Capture at Saturn
Dv_Sat = 900 ;
m_prop_Sat = preliminary_prop_mass(Dv_Sat,m_dry_tot,Isp_bip) ;
dt_Sat = m_prop_Sat / mdot_bip ;

% Moon tour and Enceladus orbit insertion
Dv_EOI = 2200 ;
m_prop_EOI = preliminary_prop_mass(Dv_EOI,m_dry_tot,Isp_bip) ;
dt_EOI = m_prop_EOI / mdot_bip ;

% Landing
Dv_land = 514 ;
mdot_land = 0.03925 ;
Isp_land = 228 ;
m_prop_land = preliminary_prop_mass(Dv_land,m_dry_land,Isp_land) ;
dt_land = m_prop_land / mdot_land ;

% Station keeping
Dv_SK = 550 ;
mdot_SK = 0.00905 ; % [kg/s]
Isp_SK = 229 ;
m_prop_SK_orb_land = preliminary_prop_mass(Dv_SK,m_dry_tot,Isp_SK) ;
dt_SK_orb_land = m_prop_SK_orb_land / mdot_SK ;

m_prop_SK_orb = preliminary_prop_mass(Dv_SK,m_dry_orb,Isp_SK) ;
dt_SK_orb = m_prop_SK_orb / mdot_SK ;

dt_SK_tot = dt_SK_orb + dt_SK_orb_land ;