%% COVERAGE
clearvars; clc; cspice_kclear();

% Define Kernel
%cd('..')
cspice_furnsh('spice_kernels/naif0012.tls')
cspice_furnsh('spice_kernels/de440s.bsp')
cspice_furnsh('spice_kernels/sat441.bsp')
cspice_furnsh('spice_kernels/pck00010_mod.tpc')
cspice_furnsh('spice_kernels/gm_de431.tpc')
cspice_furnsh('spice_kernels/dss.bsp')
cspice_furnsh('spice_kernels/dss.tf')
cspice_furnsh('spice_kernels/pck00010.tpc')
cspice_furnsh('spice_kernels/plu058.bsp')
cspice_furnsh ('spice_kernels\LATs.bsp');
cspice_furnsh ('spice_kernels\LATs.tf');
cspice_furnsh('spice_kernels\earthstns_fx_201023.bsp')
cspice_furnsh('spice_kernels\earthstns_itrf93_201023.bsp')
cspice_furnsh('spice_kernels\earth_fixed.tf')
cspice_furnsh('spice_kernels\earth_latest_high_prec.bpc')
cspice_furnsh('spice_kernels\earth_topo_201023.tf')
cspice_furnsh('spice_kernels\earth_200101_990628_predict.bpc')
cd('MA')
%% Preliminary operations
clear; clc;
% Define constants
mu_Sun=cspice_bodvrd('SUN','GM',1);

% Create vector with the times of interest

%%%%%%%%% SET AT NEED OF THE USER %%%%%%%%%%%%%%%
et0 = cspice_str2et('2052 JAN 01 00:00:00 UTC');
etf = cspice_str2et('2052 FEB 01 00:00:00 UTC');
delta_et = 3600; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

et_vec = et0:delta_et:etf;

% Names of the stations
stations = {'DSS-14','DSS-43','DSS-63'};

%% Compute position of Enceladus in ECI

Earth2Enc = zeros(6,length(et_vec));

for ii = 1:length(et_vec)
    Earth2Enc(:,ii) = cspice_spkezr('602', et_vec(ii), 'ECLIPJ2000', 'NONE', 'Earth');
end

%% Compute the elevations of Saturn System

elevations = zeros(length(stations),length(et_vec));

for id = 1:length(stations)
        
    for i = 1:length(et_vec)

        % Save the actual time
        t_i = et_vec(i);

        % Save the state of Enceladus in ECI
        rv_eci = Earth2Enc(:,i);
        
        % Measurement Model
        [R, Az, El] = measurement_model(stations{id},t_i,rv_eci);

        % Save Azimuth and Elevation
        elevations(id,i) = El*cspice_dpr();
    end
end

%% Plot the elevations for each station
for id = 1:3
    figure
    plot(et_vec/cspice_spd(), elevations(id,:),'b')
    grid minor
    hold on
    xlabel('Time since SOI [days]')
    ylabel('Elevation [deg]')
    title(stations{id})
end








