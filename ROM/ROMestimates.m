clearvars ; close all ; clc ;
set( 0, 'defaultlegendinterpreter', 'latex' ) ;
set( 0, 'defaulttextinterpreter', 'latex' ) ;
color = { [0 0.4470 0.7410],[0.8500 0.3250 0.0980], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330], [0.9290 0.6940 0.1250], [0.6350 0.0780 0.1840] } ; % Blue, Red, Purple, Green, Cyan, Yellow
fontsize = 15 ;


% ----------------------------------------------------
% Script computes ROM estimation for each architecture
% ----------------------------------------------------



% ONLY CHANGE PARAMETERS IN THE BOX BELOW
% ---------------------------------------

Isp = 300 ; % [s  ] - Hydrazine specific impulse
DV = 2400 ; % [m/s] - Assuming Cassini's deltaV with NO margins added

% ---------------------------------------
% ONLY CHANGE PARAMETERS IN THE BOX ABOVE






% !!! DO NOT TOUCH BELOW HERE !!!





% Fixed constants
g0 = 9.81 ;

% ------------------------------------------------------------------------


% Past mission data    
    % Rosetta:
      m_orbiterDry_Rosetta = 1230 ; m_lander_Philae = 100 ; m_pl_orbiter_Rosetta = 165 ; m_pl_lander_Philae = 21 ;
      m_pl_Rosetta = m_pl_lander_Philae + m_pl_orbiter_Rosetta ;
      m_orbitDry_Rosetta = m_orbiterDry_Rosetta + m_lander_Philae ;
    % Curiosity:
      m_spacecraftDry_Curiosity = 935 ; m_rover_Curiosity = 899 ; m_pl_rover_Curiosity = 80;
      m_pl_Curiosity = m_pl_rover_Curiosity ;
      m_orbitDry_Curiosity = m_spacecraftDry_Curiosity + m_rover_Curiosity ;
    % Cassini:
      m_orbiterDry_Cassini = 2068 ; m_lander_Huygens = 320 ; m_pl_Huygens = 6.3+17.3+6.3+8.1+1.9+3.9 ; m_pl_orbiter_Cassini = 43+56.5+15.5+37.1+43.3+14.4+29+23.8+16.8+8.8+10.3+37.7 ;
      m_pl_Cassini = m_pl_Huygens + m_pl_orbiter_Cassini ;
      m_orbitDry_Cassini = m_orbiterDry_Cassini + m_lander_Huygens ;
    % Mars Express
      m_orbiterDry_MarsExpress = 606 ; m_pl_orbiter_MarsExpress = 113 ; m_pl_lander_MarsExpress = 9 ; m_lander_MarsExpress = 60 ;
      m_pl_MarsExpress = m_pl_orbiter_MarsExpress + m_pl_lander_MarsExpress ;
      m_orbitDry_MarsExpress = m_orbiterDry_MarsExpress + m_lander_MarsExpress ;
    % Perseverance
      m_lander_Perseverance = 1025 ; m_pl_Perseverance = 61+45 ;
    % InSight
      m_lander_InSight = 358 ; m_pl_InSight = 50 ;
    % ExoMars
      m_rover_ExoMars = 310 ; m_pl_ExoMars = 2.13 + 1.74 + 1 + 1.7 + 2 + 2.4 + 11.5 + 5 ;
    % Opportunity
      % Not enough info ------ m_rover_Opportunity = 185 ; m_pl_Opportunity =  ;

% Single-module: Considers that on-orbit dry mass is the sum of orbiter + lander
m_orbitDry_singlemoduleRegression = [ 235, 255, 400, 550, 585, 655, 840, 900, 1050, 1025, 1135 ] ;
m_pl_singlemoduleRegression = [ 10, 15, 58, 68, 90, 75, 117, 85, 94, 130, 158 ] ;

% Multi-module: Considers that lander is a payload of the orbiter
m_orbitDry_multimoduleRegression = [ m_orbitDry_Rosetta, m_orbitDry_Curiosity, m_orbitDry_Cassini, m_orbitDry_MarsExpress ] ;
m_pl_multimoduleRegression = [ m_pl_Rosetta, m_pl_Curiosity, m_pl_Cassini, m_pl_MarsExpress ] ;

% Does not consider lander dry mass or lander payloads
m_dry_orbiter = [ m_orbiterDry_Rosetta, m_orbiterDry_Cassini, m_orbiterDry_MarsExpress ] ;
m_pl_orbiter = [ m_pl_orbiter_Rosetta, m_pl_orbiter_Cassini, m_pl_orbiter_MarsExpress ] ;

% Does not consider orbiter dry mass or orbiter payloads
m_lander = [ m_lander_Philae, m_rover_Curiosity, m_lander_Huygens, m_lander_MarsExpress, m_lander_Perseverance, m_lander_InSight, m_rover_ExoMars ] ;
m_pl_lander = [ m_pl_lander_Philae, m_pl_rover_Curiosity, m_pl_Huygens, m_pl_lander_MarsExpress, m_pl_Perseverance, m_pl_InSight, m_pl_ExoMars ] ;

% Create regression for single module and multi-module missions
p_singlemodule = polyfit(m_orbitDry_singlemoduleRegression, m_pl_singlemoduleRegression, 1 ) ;
p_multimodule = polyfit( m_orbitDry_multimoduleRegression, m_pl_multimoduleRegression, 1 ) ;
p_orbiter = polyfit( m_dry_orbiter, m_pl_orbiter, 1 ) ;
p_lander = polyfit( m_lander, m_pl_lander, 1 ) ;

% Plot regressions
x_singlemodule = linspace( 0, 1600, 2000 ) ;
y_singlemodule = polyval( p_singlemodule, x_singlemodule ) ;
x_multimodule = linspace( 0, 3500, 2000 ) ;
y_multimodule = polyval( p_multimodule, x_multimodule ) ;
x_orbiter = linspace( 0, 2500, 2000 ) ;
y_orbiter =  polyval( p_orbiter, x_orbiter ) ;
x_lander = linspace( 0, 1800, 2000 ) ;
y_lander = polyval( p_lander, x_lander ) ;

figure(1) ; % Single-module
hold on ; grid on ;
plot( x_singlemodule, y_singlemodule, 'linewidth', 1.5 ) ;
scatter( m_orbitDry_singlemoduleRegression, m_pl_singlemoduleRegression, 40, 'filled' ) ;
xlabel( '\textbf{On-orbit dry mass [kg]}', 'fontsize', fontsize ) ;
ylabel( '\textbf{Payload mass [kg]}', 'fontsize', fontsize ) ;
title( '\textbf{Mass ROM: single-module architectures}', 'fontsize', fontsize ) ;

% This has to be fixed
% figure(2) ; % Multi-module
% hold on ; grid on ;
% plot( x_multimodule, y_multimodule, 'linewidth', 1.5 ) ;
% scatter( m_orbitDry_multimoduleRegression, m_pl_multimoduleRegression, 30, 'filled' ) ;
% xlabel( 'On Orbit Dry Mass [Kg]' ) ;
% ylabel( 'Payload Mass [Kg]' ) ;
% title( 'Planetary S/C ROM Mass Estimation: Multi-Module' ) ;

figure(3) ; % Orbiters
hold on ; grid on ;
plot( x_orbiter, y_orbiter, 'linewidth', 1.5 ) ;
scatter( m_dry_orbiter, m_pl_orbiter, 40, 'filled' ) ;
xlabel( '\textbf{Orbiter dry mass [kg]}', 'fontsize', fontsize ) ;
ylabel( '\textbf{Orbiter payload mass [kg]}', 'fontsize', fontsize ) ;
title( '\textbf{Mass ROM: multi-module architectures (orbiter)}', 'fontsize', fontsize ) ;

figure(4) ; % Landers/Rovers
hold on ; grid on ;
plot( x_lander, y_lander, 'linewidth', 1.5 ) ;
scatter( m_lander, m_pl_lander, 40, 'filled' ) ;
xlabel( '\textbf{Lander mass [kg]}', 'fontsize', fontsize ) ;
ylabel( '\textbf{Lander payload mass [kg]}', 'fontsize', fontsize ) ;
title( '\textbf{Mass ROM: multi-module architectures (lander)}', 'fontsize', fontsize ) ;

% ------------------------------------------------------------------------

% Orbilander + Service Module
pl_OLSM = 96.9 ;
m_dry_OLSM = (pl_OLSM-p_singlemodule(2))/p_singlemodule(1) ;
m_wet_OLSM = m_dry_OLSM*exp(DV/(Isp*g0)) ;
m_prop_OLSM = m_wet_OLSM-m_dry_OLSM ;

% NS Orbiter + S Lander
pl_SL_NSOSL = 75.4 ;
pl_NSO_NSOSL = 30.5 ;
pl_NSOSL = pl_SL_NSOSL+pl_NSO_NSOSL ;
m_dry_NSOSL = (pl_NSOSL-p_singlemodule(2))/p_singlemodule(1) ;
m_wet_NSOSL = m_dry_NSOSL*exp(DV/(Isp*g0)) ;
m_prop_NSOSL = m_wet_NSOSL-m_dry_NSOSL ;
m_lander_NSOSL = (pl_SL_NSOSL-p_lander(2))/p_lander(1) ;
m_orbiter_NSOSL = (pl_NSO_NSOSL-p_orbiter(2))/p_orbiter(1) ;

% S Orbiter + S Lander
pl_SO_SOSL = 85.2 ;
pl_SL_SOSL = 75.4 ;
pl_SOSL = pl_SL_SOSL+pl_SO_SOSL ;
m_dry_SOSL = (pl_SOSL-p_singlemodule(2))/p_singlemodule(1) ;
m_wet_SOSL = m_dry_SOSL*exp(DV/(Isp*g0)) ;
m_prop_SOSL = m_wet_SOSL-m_dry_SOSL ;
m_lander_SOSL = (pl_SL_SOSL-p_lander(2))/p_lander(1) ;
m_orbiter_SOSL = (pl_SO_SOSL-p_orbiter(2))/p_orbiter(1) ;

% S Orbiter + 2 Lander
pl_SO_SO2L = 85.2 ;
pl_2L_SO2L = 150.8 ;
pl_SO2L = pl_2L_SO2L+pl_SO_SO2L ;
m_dry_SO2L = (pl_SO2L-p_singlemodule(2))/p_singlemodule(1) ;
m_wet_SOnL = m_dry_SO2L*exp(DV/(Isp*g0)) ;
m_prop_SOnL = m_wet_SOnL-m_dry_SO2L ;
m_lander_SO2L = (pl_2L_SO2L-p_lander(2))/p_lander(1) ;
m_orbiter_SO2L = (pl_SO_SO2L-p_orbiter(2))/p_orbiter(1) ;

figure(1) ;
hold on ; grid on ;
plot(m_dry_OLSM,pl_OLSM, 'p', 'markersize', 15, 'markerfacecolor', color{7}, 'markeredgecolor', color{7})
legend( '\textbf{Regression}', '\textbf{Previous Missions}', '\textbf{Orbiter-Lander + SM}', 'location', 'northwest', 'fontsize', fontsize ) ;

% This has to be fixed
% figure(2) ;
% hold on ; grid on ;
% scatter(m_dry_NSOSL,pl_NSOSL,40,'filled') ;
% scatter(m_dry_SOSL,pl_SOSL,40,'filled') ;
% scatter(m_dry_SOnL,pl_SOnL,40,'filled') ;
% legend( 'Regression', 'Previous Missions', 'NS-O + SL', 'SO + SL', 'SO + nL', 'location', 'northwest' ) ;

figure(3) ;
hold on ; grid on ;
plot(m_orbiter_NSOSL,pl_NSO_NSOSL, 'd', 'markersize', 10, 'markerfacecolor', color{6}, 'markeredgecolor', color{6} ) ;
plot(m_orbiter_SOSL,pl_SO_SOSL, '^', 'markersize', 15, 'markerfacecolor', color{3}, 'markeredgecolor', color{3} ) ;
plot(m_orbiter_SO2L,pl_SO_SO2L, 's', 'markersize', 10, 'markerfacecolor', color{4}, 'markeredgecolor', color{4} ) ;
legend( '\textbf{Regression}', '\textbf{Previous Missions}', '\textbf{NSO + SL}', '\textbf{SO + SL}', '\textbf{SO + nL}', 'location', 'northwest', 'fontsize', fontsize ) ;

figure(4) ;
hold on ; grid on ;
plot(m_lander_NSOSL,pl_SL_NSOSL, 'd', 'markersize', 15, 'markerfacecolor', color{6}, 'markeredgecolor', color{6} ) ;
plot(m_lander_SOSL,pl_SL_SOSL, '^', 'markersize', 10, 'markerfacecolor', color{3}, 'markeredgecolor', color{3} ) ;
plot(m_lander_SO2L,pl_2L_SO2L, 's', 'markersize', 10, 'markerfacecolor', color{4}, 'markeredgecolor', color{4} ) ;
legend( '\textbf{Regression}', '\textbf{Previous Missions}', '\textbf{NSO + SL}', '\textbf{SO + SL}', '\textbf{SO + nL}', 'location', 'northwest', 'fontsize', fontsize ) ;


%% Power %%%%%

% It does not make sense to create new regressions from past missions, as
% the payload powers in some cases are higher than the overall power
% budget. This is because the payloads are not used contemporariliy. As a
% first approximation, the planetary s/c power budget regression (Lavagna's
% slides) is used, separately for orbiter and lander.

% Empirical power law
P = @(P_pl) 332.93*log(P_pl) - 1046.6;

% Orbilander + Service Module
p_pl_1 = 326;
p1 = P(p_pl_1);

% Small Orbiter + Lander
p_pl_2a = 269;
p_pl_2b = 64;
P_orbiter_ = P(p_pl_2a+p_pl_2b);

% Orbiter with sampling + Lander 
p_pl_3a = 269;
p_pl_3b = 303.4;
p3 = P(p_pl_3a+p_pl_3b);

x_singlemodule = linspace(150,700);
y_singlemodule = P(x_singlemodule);
figure() ;
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
