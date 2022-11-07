%% Load SPICE Kernels
cspice_kclear();
cspice_furnsh('spice_kernels/pck00010.tpc')
cspice_furnsh('spice_kernels/naif0012.tls')
cspice_furnsh('spice_kernels/gm_de431.tpc')
cspice_furnsh('spice_kernels/de440s.bsp')
cspice_furnsh('spice_kernels/sat441.bsp')
cspice_furnsh('spice_kernels/dss.bsp')
%%
clear;close all;clc;
date_1 = '2022-11-07 00:00:0.000 UTC';
t1 = cspice_str2et(date_1);
t2 = t1 + 17*24*3600;
tv = linspace(t1,t2,200);
ind = 1;
Rhea = zeros(3,length(tv));
Dione = Rhea;
Tethys = Rhea;
Enceladus = Rhea;
Titan = Rhea;
Mimas = Rhea;
for t = tv
    Rhea(:,ind) = cspice_spkpos('Rhea',t,'ECLIPJ2000','None','Saturn');
    Dione(:,ind) = cspice_spkpos('Dione',t,'ECLIPJ2000','None','Saturn');
    Tethys(:,ind) = cspice_spkpos('Tethys',t,'ECLIPJ2000','None','Saturn');
    Enceladus(:,ind) = cspice_spkpos('Enceladus',t,'ECLIPJ2000','None','Saturn');
    Titan(:,ind) = cspice_spkpos('Titan',t,'ECLIPJ2000','None','Saturn');
    Mimas(:,ind) = cspice_spkpos('Mimas',t,'ECLIPJ2000','None','Saturn');
    ind = ind+1;
end