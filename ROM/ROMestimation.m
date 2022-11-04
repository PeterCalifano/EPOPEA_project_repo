clearvars ; close all ; clc ;

% ----------------------------------------------------
% Script computes ROM estimation for each architecture
% ----------------------------------------------------



% ONLY CHANGE PARAMETERS IN THE BOX BELOW
% ---------------------------------------

Isp = 300 ; % [s  ] - Hydrazine specific impulse
DV = 2400 ; % [m/s] - Assuming Cassini's deltaV with NO margins added

% ---------------------------------------
% ONLY CHANGE PARAMETERS IN THE BOX ABOVE






% DO NOT TOUCH BELOW HERE

% Fixed constants
g0 = 9.81 ;

% ------------------------------------------------------------------------

% Derive empirical mass law for missions with ONE module. Similar to Lavagna's plots
m_orbitDry_singlemoduleRegression = [ 235, 255, 400, 550, 585, 655, 840, 900, 1050, 1025, 1135 ] ;
m_pl_singlemoduleRegression = [ 10, 15, 58, 68, 90, 75, 117, 85, 94, 130, 158 ] ;

% Derive empirical mass law for mission with TWO modules (orbiter+lander)   
%   --> If you assume that lander mass (with their payloads) is part of on-orbit dry mass, regression does not provide realistic results (masses are extremely low)
%   --> To overcome this, ONLY PAYLOAD MASSES (of both orbiter and lander) ARE CONSIDERED AS PAYLOADS, SO STRUCTURE OF LANDER IS NOT INCLUDED IN THIS MASS
%   --> Overall, enter plot with TOTAL PAYLOAD MASS, SUMMING ORBITER'S AND LANDER'S PAYLOADS   
%    
    % Rosetta:
      m_orbiterDry_Rosetta = 1230 ; m_lander_Rosetta = 100 ; m_pl_orbiter_Rosetta = 165 ; m_pl_lander_Philae = 21 ;
      m_pl_Rosetta = m_pl_lander_Philae + m_pl_orbiter_Rosetta ;
      m_orbitDry_Rosetta = m_orbiterDry_Rosetta + m_lander_Rosetta ;
    % Curiosity:
      m_spacecraftDry_Curiosity = 935 ; m_rover_Curiosity = 899 ; m_pl_rover_Curiosity = 80;
      m_pl_Curiosity = m_pl_rover_Curiosity ;
      m_orbitDry_Curiosity = m_spacecraftDry_Curiosity + m_rover_Curiosity ;
    % Cassini:
      m_orbiterDry_Cassini = 2523 ; m_lander_Huygens = 320 ; m_pl_Huygens = 6.3+17.3+6.3+8.1+1.9+3.9 ; m_pl_orbiter_Cassini = 43+56.5+15.5+37.1+43.3+14.4+29+23.8+16.8+8.8+10.3+37.7 ;
      m_pl_Cassini = m_pl_Huygens + m_pl_orbiter_Cassini ;
      m_orbitDry_Cassini = m_orbiterDry_Cassini + m_lander_Huygens ;
    % Mars Express
      m_orbiterDry_MarsExpress = 606 ; m_pl_orbiter_MarsExpress = 113 ; m_pl_lander_MarsExpress = 9 ; m_lander_MarsExpress = 60 ;
      m_pl_MarsExpress = m_pl_orbiter_MarsExpress + m_pl_lander_MarsExpress ;
      m_orbitDry_MarsExpress = m_orbiterDry_MarsExpress + m_lander_MarsExpress ;

m_orbitDry_multimoduleRegression = [ m_orbitDry_Rosetta, m_orbitDry_Curiosity, m_orbitDry_Cassini, m_orbitDry_MarsExpress ] ;
m_pl_multimoduleRegression = [ m_pl_Rosetta, m_pl_Curiosity, m_pl_Cassini, m_pl_MarsExpress ] ;

% Create regression for single module and multi-module missions
p_singlemodule = polyfit(m_orbitDry_singlemoduleRegression, m_pl_singlemoduleRegression, 1 ) ;
p_multimodule = polyfit( m_orbitDry_multimoduleRegression, m_pl_multimoduleRegression, 1 ) ;

% Plot regressions
x_singlemodule = linspace( 100, 1400, 2000 ) ;
y_singlemodule = polyval( p_singlemodule, x_singlemodule ) ;
x_multimodule = linspace( 100, 3500, 2000 ) ;
y_multimodule = polyval( p_multimodule, x_multimodule ) ;

figure(1) ;
hold on ; grid on ;
plot( x_singlemodule, y_singlemodule, 'linewidth', 1.5 ) ;
scatter( m_orbitDry_singlemoduleRegression, m_pl_singlemoduleRegression, 30, 'filled' ) ;
xlabel( 'On Orbit Dry Mass [Kg]' ) ;
ylabel( 'Payload Mass [Kg]' ) ;
title( 'Planetary S/C ROM Mass Estimation: Single-Module' ) ;

figure(2) ;
hold on ; grid on ;
plot( x_multimodule, y_multimodule, 'linewidth', 1.5 ) ;
scatter( m_orbitDry_multimoduleRegression, m_pl_multimoduleRegression, 30, 'filled' ) ;
xlabel( 'On Orbit Dry Mass [Kg]' ) ;
ylabel( 'Payload Mass [Kg]' ) ;
title( 'Planetary S/C ROM Mass Estimation: Multi-Module' ) ;

% ------------------------------------------------------------------------

% Orbilander + Service Module
pl_OLSM = 96.9 ;
m_dry_OLSM = (pl_OLSM-p_singlemodule(2))/p_singlemodule(1) ;
m_wet_OLSM = m_dry_OLSM*exp(DV/(Isp*g0)) ;
m_prop_OLSM = m_wet_OLSM-m_dry_OLSM ;

% NS Orbiter + S Lander
pl_SL = 75.4 ;
pl_NSO = 30.5 ;
pl_NSOSL = pl_SL+pl_NSO ;
m_dry_NSOSL = (pl_NSOSL-p_multimodule(2))/p_multimodule(1) ;
m_wet_NSOSL = m_dry_NSOSL*exp(DV/(Isp*g0)) ;
m_prop_NSOSL = m_wet_NSOSL-m_dry_NSOSL ;

% S Orbiter + S Lander
pl_SO = 75.4 ;
pl_SL = 85.2 ;
pl_SOSL = pl_SL+pl_SO ;
m_dry_SOSL = (pl_SOSL-p_multimodule(2))/p_multimodule(1) ;
m_wet_SOSL = m_dry_SOSL*exp(DV/(Isp*g0)) ;
m_prop_SOSL = m_wet_SOSL-m_dry_SOSL ;

% S Orbiter + n Lander
pl_SO = 63.7 ;
pl_nL = 108.8 ;
pl_SOnL = pl_nL+pl_SO ;
m_dry_SOnL = (pl_SOnL-p_multimodule(2))/p_multimodule(1) ;
m_wet_SOnL = m_dry_SOnL*exp(DV/(Isp*g0)) ;
m_prop_SOnL = m_wet_SOnL-m_dry_SOnL ;

figure(1) ;
hold on ; grid on ;
scatter(m_dry_OLSM,pl_OLSM,40,'filled')
legend( 'Regression', 'Previous Missions', 'Orbiter-Lander + SM', 'location', 'northwest' ) ;

figure(2) ;
hold on ; grid on ;
scatter(m_dry_NSOSL,pl_NSOSL,40,'filled') ;
scatter(m_dry_SOSL,pl_SOSL,40,'filled') ;
scatter(m_dry_SOnL,pl_SOnL,40,'filled') ;
legend( 'Regression', 'Previous Missions', 'NS-O + SL', 'SO + SL', 'SO + nL', 'location', 'northwest' ) ;

xtitle = { '', 'OL + SM', 'NS-O S-L', 'SO + SL' } ;
ytitle = { 'Payload mass', 'On-orbit dry', 'Wet mass', 'propellant' }' ;
Results = [ pl_OLSM, m_dry_OLSM, m_wet_OLSM, m_prop_OLSM ; pl_NSOSL, m_dry_NSOSL, m_wet_NSOSL, m_prop_NSOSL; pl_SOSL, m_dry_SOSL, m_wet_SOSL, m_prop_SOSL ] ;


%% Power %%%%%
% Empirical power law
P = @(P_pl) 332.93*log(P_pl) - 1046.6;

% Orbilander + Service Module
p_pl_1 = 326;
p1 = P(p_pl_1);

% Small Orbiter + Lander
p_pl_2a = 269;
p_pl_2b = 64;
p_singlemodule = P(p_pl_2a+p_pl_2b);

% Orbiter with sampling + Lander 
p_pl_3a = 269;
p_pl_3b = 303.4;
p3 = P(p_pl_3a+p_pl_3b);

x_singlemodule = linspace(150,700);
y_singlemodule = P(x_singlemodule);
figure
plot(x_singlemodule,y_singlemodule,'k','linewidth',1.5)
grid on
hold on
scatter(p_pl_1,p1,40,'filled')
scatter(p_pl_2a+p_pl_2b,p_singlemodule,40,'filled')
scatter(p_pl_3a+p_pl_3b,p3,40,'filled')
xlabel('Payload Power [W]')
ylabel('Overall Power [Kg]')
title('Planetary S/C ROM Power Estimation')
legend('','Orbiter-Lander + SM','NS Orbiter + S Lander',...
    'S Orbiter + S Lander','location','northwest')
