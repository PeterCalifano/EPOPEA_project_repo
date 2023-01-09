%% COVERAGE
clearvars; clc; cspice_kclear();

% Define Kernel
cd('..')
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
cd('MA')
%% Preliminary operations

% Define constants
mu_Sun=cspice_bodvrd('SUN','GM',1);

% Create vector with the times of interest
et0 = cspice_str2et('2047 AUG 27 21:19:55 UTC');
etf = cspice_str2et('2055 AUG 27 21:19:55 UTC');
et_vec = linspace(et0,etf,1000);

% Names of the stations
stations = {'DSS-14','DSS-43','DSS_63'};

%% Compute position of Enceladus in ECI

Earth2Enc = zeros(6,length(et_vec));

for ii = 1:length(et_vec)
    Earth2Enc(:,ii) = cspice_spkezr('Sun', et_vec(ii), 'ECLIPJ2000', 'NONE', 'Earth');
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
        elevations(id,i) = El;
    end
end
%%
ROT_ECI2TOPO = cspice_sxform('J2000', 'DSS-13_TOPO' , et0);