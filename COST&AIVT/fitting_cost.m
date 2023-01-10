
%% Fitting - Launch mass / Mission cost
clear
clc 
close all

% Mission Cost
MC = [3135.6e+6, 3.9e+9, 4.25e+9, 1.49e+9, 1.16e+9, 420e+6, 199e+6];

% Launch mass
LM = [6610, 5863, 6065, 3000, 2110, 670, 385];

M_epo = 12407;
% Plot
c = polyfit(LM,MC,2); xfit = linspace(380,15000,1000);
yfit = polyval(c,xfit); cost_epo = polyval(c,M_epo);

C = {'Orbilander','Cassini','Europa Clipper','Rosetta','Osiris-REX',...
    'Phoenix','Stardust'};
col =  [0 0.4470 0.7410
       0.8500 0.3250 0.0980
       0.9290 0.6940 0.1250
       0.4940 0.1840 0.5560
       0.4660 0.6740 0.1880
       0.3010 0.7450 0.9330
       0.6350 0.0780 0.1840];

figure 
for i = 1:length(MC)
    scatter(LM(i),MC(i)/10^9,40,col(i,:),'o',...
        'LineWidth',1.5); hold on;
end
legend(C);
scatter(M_epo,cost_epo/10^9,40,'o','k','LineWidth',1.5,...
    'DisplayName','EPOPEA'); legend();
hold on;
plot(xfit,yfit/10^9,'LineWidth',1.5,'Color',"#D95319",'DisplayName',...
    'Fitting curve'); grid on;
labelpoints(M_epo,cost_epo/10^9,'EPOPEA','SE',0.2,1)
xlabel('Launch mass [kg]'); ylabel('Mission cost [B$]');
title('Fitting curve Mass-Cost for past missions')

cost_epo = cost_epo*0.94;

%% Fitting - PL cost / Mission cost

% Mission Cost
MC = [3135.6e+6, 3.9e+9, 4.25e+9, 1.49e+9, 1.16e+9, 420e+6];

% Payload / mission operation cost
PC = [459500, 87.5e+6, 110e+6, 227.5e+6, 283e+6, 12.6e+6];

% Plot
c = polyfit(PC,MC,2); xfit = linspace(4e+5,13e+6,1000);
yfit = polyval(c,xfit);

figure 
scatter(PC/1000,MC/1000,40,'o','filled'); hold on; grid on
xlabel('k$'); ylabel('k$');
%plot(xfit,yfit/1000,'LineWidth',1.5);

