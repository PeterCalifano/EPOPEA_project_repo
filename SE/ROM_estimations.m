%%%%%% Rough Order of Magnitude of the architectures %%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc
%% Mass %%%%%
g0 = 9.81;
Isp = 300; %hydrazine specific impulse
DV = 3000; % DV of Cassini with 50% safety margin
% Empirical mass law
m_orb = [235,255,400,550,585,655,840,900,1050,1025,1135];
m_pl = [10,15,58,68,90,75,117,85,94,130,158];
%p = polyfit(m_pl,m_orb,1);

p2 = polyfit(m_orb,m_pl,1);
x = linspace(100,1400);
y = polyval(p2,x);

% Orbilander + Service Module
pl_1 = 96.9;
m_dry_1 = (pl_1-p2(2))/p2(1);
m_wet_1 = m_dry_1*exp(DV/(Isp*g0));
m_prop_1 = m_wet_1-m_dry_1;

% Small Orbiter + Lander
pl_2a = 75.4;
pl_2b = 30.5;
pl_2 = pl_2a+pl_2b;
m_dry_2 = (pl_2-p2(2))/p2(1);
m_wet_2 = m_dry_2*exp(DV/(Isp*g0));
m_prop_2 = m_wet_2-m_dry_2;

% Orbiter with sampling + Lander
pl_3a = 75.4;
pl_3b = 85.2;
pl_3 = pl_3a+pl_3b;
m_dry_3 = (pl_3-p2(2))/p2(1);
m_wet_3 = m_dry_3*exp(DV/(Isp*g0));
m_prop_3 = m_wet_3-m_dry_3;

plot(x,y,'k','linewidth',1.5)
hold on
grid on
scatter(m_orb,m_pl,30,'filled')
scatter(m_dry_1,pl_1,40,'filled')
scatter(m_dry_2,pl_2a+pl_2b,40,'filled')
scatter(m_dry_3,pl_3a+pl_3b,40,'filled')
xlabel('On Orbit Dry Mass [Kg]')
ylabel('Payload Mass [Kg]')
title('Planetary S/C ROM Mass Estimation')
legend('','Previous Missions','Orbiter-Lander + SM','NS Orbiter + S Lander',...
    'S Orbiter + S Lander','location','northwest')
%% Power %%%%%
% Empirical power law
P = @(P_pl) 332.93*log(P_pl) - 1046.6;

% Orbilander + Service Module
p_pl_1 = 326;
p1 = P(p_pl_1);

% Small Orbiter + Lander
p_pl_2a = 269;
p_pl_2b = 64;
p2 = P(p_pl_2a+p_pl_2b);

% Orbiter with sampling + Lander 
p_pl_3a = 269;
p_pl_3b = 303.4;
p3 = P(p_pl_3a+p_pl_3b);

x = linspace(150,700);
y = P(x);
plot(x,y,'k','linewidth',1.5)
grid on
hold on
scatter(p_pl_1,p1,40,'filled')
scatter(p_pl_2a+p_pl_2b,p2,40,'filled')
scatter(p_pl_3a+p_pl_3b,p3,40,'filled')
xlabel('Payload Power [W]')
ylabel('Overall Power [Kg]')
title('Planetary S/C ROM Power Estimation')
legend('','Orbiter-Lander + SM','NS Orbiter + S Lander',...
    'S Orbiter + S Lander','location','northwest')
