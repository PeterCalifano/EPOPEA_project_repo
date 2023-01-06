%ACQUISITION DATA RATES
clc; clear all;
%% Coarse imaging mode
AoI = 120/360*2*pi;
AltO = 60; %[km]
tOrbit = 12*3600; %[s]
tAcq = 1*3600; %[s]
tTMTC = 4*3600; %[s]

%WAC
WACDpM = 0.0042; %[Gb]
WACFoV = 40/360*2*pi; %[rad]

[WACmeasures_per_orbit,WACdata_per_orbit,WACorbit_tot,WoIWAC] = MeasureAnalysis(WACDpM,WACFoV,AoI,AltO)

%TES
TESDpM = 0.00031; %[Gb]
TESFoV = 2/360*2*pi; %[rad]

[TESmeasures_per_orbit,TESdata_per_orbit,TESorbit_tot,WoITES] = MeasureAnalysis(TESDpM,TESFoV,AoI,AltO,WoIWAC)

%Altimeter
ALTR = 9.89/8/1e6; %[Gb/s]
ALTdata_per_orbit = ALTR * tAcq; %[Gb]

%Total for Coarse imaging mode
MODE_CI_data_per_orbit = TESdata_per_orbit + WACdata_per_orbit + ALTdata_per_orbit;
MODE_CI_R_Acq = MODE_CI_data_per_orbit/tAcq*1e6; %[kb/s]
MODE_CI_R_TMTC = MODE_CI_data_per_orbit/tTMTC*1e6; %[kb/s]

%% NAC mode
AoI = 3.6/360*2*pi;
%NAC
NACDpM = 0.0042; %[Gb]
NACFoV = 0.3/360*2*pi; %[rad]

dist1 = 70.8 + 96.2;
dist2 = 40.0 + 68.6;

%Angleofinterest = atan(max([dist1,dist2])/(510/2))
Angleofinterest = atan(10/(510/2));

[NACmeasures_per_orbit,NACdata_per_orbit,NACorbit_tot,WoINAC] = MeasureAnalysis(NACDpM,NACFoV,Angleofinterest,AltO)


%% Functions
function [measures_per_orbit,data_per_orbit,orbit_tot,footprint_theta_enc_tot] = MeasureAnalysis(DpM,FoV,AoI,AltO,WoI)
% DpM = Data per Measure of instrument - from Payload description document
% FoV = Field of View of instrument - from Payload description document
% AoI = Arc of sphere of interest - length in deg of data measurement for one orbit (should = 120deg for SPR)
% WoI = Width of interest - width in deg of data measurement for one  orbit - to determine swath width & number
% AltO = Altitude of the orbit [km]

% MpO = # of measurement per orbit
% DpO = Data per orbit
% Otot = # total of orbit for mapping of AoI

%DATA 
R_enc = 510/2; %[km]
overlap = 2;

footprint = 2*AltO*atan(FoV/2);
footprint_woverlap = footprint/overlap;
footprint_theta_enc = 2*atan((footprint_woverlap/2)/R_enc)

if exist("WoI","var")
    measures_per_swath = ceil(WoI/footprint_theta_enc);
else
    measures_per_swath = 1;
end
footprint_theta_enc_tot = footprint_theta_enc*measures_per_swath;
measures_per_orbit = ceil(AoI/footprint_theta_enc_tot);


data_per_orbit = measures_per_orbit * DpM;
orbit_tot = ceil(AoI/(footprint_theta_enc*measures_per_swath));
end