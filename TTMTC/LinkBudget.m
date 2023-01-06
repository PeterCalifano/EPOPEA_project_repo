clear;
close all;
clc;

%% Data
max_distance_Earth = 1658.6e6*1000; % considering Saturn [m]
max_distance_OL = 100*1000; % for now, average altitude of the orbiter, wrt the 
% lander (when they are one on top of the other)

% Since BER -> graph on TMTC presentation
min_EbNo_sci = 4.2; % [dB] science case -> BER = 1e-5
min_EbNo_tel = 5.1; % [dB] telemetry case -> BER = 1e-7

% Frequency selection
% f_ka = 32e9; % [Hz] between 27 and 40GHz, for science 
f_x = 8.45e9; % between 7 and 10.7 GHz, for science
f_s = 3e9; % between 2 and 4 GHz, for O-L comm
% f_opt = 1e12; %Optical, TBD what frequency is most commonly used

% Antenna data
D_hga = 3.5; % [m] took from Orbilander paper, page 21 -> WE WOULD NEED 0.165m of diameter to have 3dB margin in uplink
D_gs = 70; % [m] GS antenna, Cebreros considered (all the same, for ESTRACK DSA)
D_lga = 0.065; % [m] from datasheets of existing helix antennas (by RUAG space)
D_patch = 0.040; % [m] from datasheets of existing patch antenna (by ISIS – Innovative Solutions In Space)

% Power transmitted
P_on_board = 160; % [W] based on orbitander paper page 27
P_Lander = 15; % [W]
P_gs = 400; %[W] not 100% sure, taken from https://www.esa.int/Enabling_Support/Operations/ESA_Ground_Stations/Cebreros_-_DSA_2

% Antenna efficiency, depending on the shape
eta_gs = 0.55; % for parabolic reflectors
eta_hga = 0.55; % for parabolic reflectors
eta_lga = 0.7; % 0.7 if helix, 0.52 if horn

% Reference temperatures for noise computation
T_gs = 21; % [K] of the Ground Station, reference: DSN antenna temperature
% NB. also by setting this much higher (290K) the margins are still very
% large
T_hga = 250; % [K] of the HGA on board
T_lga = 250; % [K] of the LGA on board

%DATA RATES
R_sci = 41.3*1e3; % data rate [bit/s], as coming from the PL, Data rate for Deep space communication
R_LtoO = 272.5*1e3; % data rate [bit/s], as coming from the PL, Data rate for Lander - Orbiter communication
R_NN = 2.8*4 *1e3; % data rate [bit/s]as coming from the PL, Data rate for non-nominal (Housekeeping & Telemetry)
R_NN_OL = 2.8*4 *1e3; % data rate [bit/s]as coming from the PL, Data rate for non-nominal (Housekeeping & Telemetry)

R_up = 700; % [bps] from Orbilander paper, page 133

%% SL + NSO
R_SLNSO = R_sci; % bit/s science data downlink
[EbNo_down_SLNSO, ~] = link_budget('to_Earth', 'downlink', 'atmosphere', R_SLNSO, max_distance_Earth, D_hga, D_gs, f_x, P_on_board, eta_hga, eta_gs, T_gs);
[EbNo_up_SLNSO, ~] = link_budget('to_Earth', 'uplink', 'atmosphere', R_up, max_distance_Earth, D_gs, D_hga, f_x, P_gs, eta_gs, eta_hga, T_hga);

disp('SL + NSO')
marginEbNo_down_SLNSO = EbNo_down_SLNSO - min_EbNo_sci % must be > 3dB
marginEbNo_up_SLNSO = EbNo_up_SLNSO - min_EbNo_tel % must be > 3dB

%% SL + SO
R_SLSO = R_sci; % bit/s science data downlink
[EbNo_down_SLSO, ~] = link_budget('to_Earth', 'downlink', 'atmosphere', R_SLSO, max_distance_Earth, D_hga, D_gs, f_x, P_on_board, eta_hga, eta_gs, T_gs);
[EbNo_up_SLSO, ~] = link_budget('to_Earth', 'uplink', 'atmosphere', R_up, max_distance_Earth, D_gs, D_hga, f_x, P_gs, eta_gs, eta_hga, T_hga);

disp('SL + SO')
marginEbNo_down_SLSO = EbNo_down_SLSO - min_EbNo_sci % must be > 3dB
marginEbNo_up_SLSO = EbNo_up_SLSO - min_EbNo_tel % must be > 3dB

%% L to O
[EbNo_LtoO, ~] = link_budget('L_to_O', 'uplink', 'no_atmosphere', R_LtoO, max_distance_OL, D_lga, D_lga, f_s, P_Lander, eta_lga, eta_lga, T_lga);
marginEbNo_down_LtoO = EbNo_LtoO - min_EbNo_sci % must be > 3dB



%% FUNCTION
function [EbNo, SNR] = link_budget(condition, direction, atmosphere, R, max_distance, Dtx, Drx, f, Pinput, eta_tx, eta_rx, Ts)
% NB. Uplink and Downlink is managed outside the function! But I left the
% input here also, just in case

% Define the wavelength
c = 299792458; %[m/s] speed of light
lambda = c/f;

% Adjust the data rate based on the encoding and modulation techniques
alpha_mod = 2; % QSPK modulation
alpha_enc_convolution = 2; % convolutional encoding

R = R * alpha_enc_convolution / alpha_mod;

% Cable losses
Lc = -2; % [dB], between -1 and -3

if strcmp(condition, 'to_Earth')
    % Tx antenna gain computation
    Gtx = 10*log10(pi*Dtx^2*eta_tx/lambda^2);

    % Ptx computation
    eta_TWTA = 0.58; % based on the input power at 150 W
    Ptx = Pinput * eta_TWTA;

    % EIRP evaluation
    % EIRP = 20*log(Ptx + Gtx + Lc);
    EIRP = 10*log10(Ptx) + Gtx + Lc;
    % Free space losses
    Ls = 20*log10(lambda/(4*pi*max_distance)); % [dB]

    % Atmospheric losses
    if strcmp(atmosphere, 'atmosphere')
        La = -0.2; % [dB] based on 32GHz of frequency, dry air
    else
        La = 0; % [dB] if no atmospheric losses are considere (L-O comm)
    end

    % Reciever gain
    Grx = 74.6;%10*log10(pi*Drx^2*eta_rx/lambda^2);

    % Compute beamwidth (parabolic reflector case)
    theta_tx = 65.3*lambda/Dtx; % [deg]
    theta_rx = 65.3*lambda/Drx; % [deg]

    % Pointing losses <- comes from the pointing budget!
    err = 0.01; % typical value 0.1 deg as from the slides
    Lp = -12*(err/theta_tx)^2;

    % Power recieved
    Prx = EIRP + Ls + La + Grx + Lp;

    %Total losses
    Losses = Ls + Lc + La + Lp;

    % Reciever noise density (due to temperature)
    k = 1.38e-23; % Boltzmann's constant
    N0 = 10*log10(k*Ts);

    % Eb/No computation
    EbNo = Prx - N0 - 10*log10(R);

    % Compute carrier margin
    mu = 78*pi/180; % [rad] standard modulation index [45-90°]
    Pmod_loss = 20*log10(cos(mu)); % power loss due to modulation

    Pcarrier = Prx + Pmod_loss;

    % Signal to Noise Ratio (of the carrier)
    SNR = Pcarrier - N0 - 10*log10(R);

else
    L = 0.289; % [m] from datasheets of existing helix antennas
    % Tx antenna gain computation
    Ctx = pi*Dtx;
    Gtx = 10.3 + 10*log10(Ctx^2*L/(lambda^3)); % peak gain value (see slide 29 of TMTC)

    % Ptx computation
    eta_TWTA = 0.58; % based on the input power at 150 W
    Ptx = Pinput * eta_TWTA;

    % EIRP evaluation
    EIRP = 10*log10(Ptx) + Gtx + Lc;

    % Free space losses
    Ls = 20*log10(lambda/(4*pi*max_distance)); % [dB]

    % Atmospheric losses
    if strcmp(atmosphere, 'atmosphere')
        La = -0.2; % [dB] based on 32GHz of frequency, dry air
    else
        La = 0; % [dB] if no atmospheric losses are considere (L-O comm)
    end

    % Reciever gain (assume same helix antennas on board)
    Crx = pi*Drx;
    Grx = 10.3 + 10*log10(Crx^2*L/(lambda^3)); 
    
    % Compute beamwidth (parabolic reflector case)
    theta_tx = 52/sqrt(Ctx^2*L/lambda^2);
    theta_rx = 52/sqrt(Crx^2*L/lambda^2);

    % Pointing losses <- comes from the pointing budget!
    err = 0.1/100; % typical value
    Lp = -12*(err/theta_rx)^2;

    % Power recieved
    Prx = EIRP + Ls + La + Grx + Lp;

    % Reciever noise density (due to temperature)
    k = 1.38e-23; % Boltzmann's constant
    N0 = 10*log10(k*Ts);

    % Eb/No computation
    EbNo = Prx - N0 - 10*log10(R);

    % Compute carrier margin
    mu = 78*pi/180; % [rad] standard modulation index [45-90°]
    Pmod_loss = 20*log10(cos(mu)); % power loss due to modulation

    Pcarrier = Prx + Pmod_loss;

    % Signal to Noise Ratio (of the carrier)
    SNR = Pcarrier - N0 - 10*log10(R);
end

end