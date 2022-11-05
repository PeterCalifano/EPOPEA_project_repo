clear;
close all;
clc;

% TO DO:
% -try X band for science downlink -> see how antenna dimensions change

%% Data
max_distance_Earth = 1658.6e6; % considering Saturn [km]
% NOTE. I tried useing meters but then we would need like too much power to
% have positive margins, so I feel like it's correct in km
max_distance_OL = 100; % for now, average altitude of the orbiter, wrt the 
% lander (when they are one on top of the other) <- SHOULD COME FROM MA

% Since BER -> graph on TMTC presentation
min_EbNo_sci = 4.2; % [dB] science case -> BER = 1e-5
min_EbNo_tel = 5.1; % [dB] telemetry case -> BER = 1e-7

% Frequency selection
f_ka = 32e9; % [Hz] between 27 and 40GHz, for science 
f_x = 8e9; % between 7 and 10.7 GHz, for interplanetary
f_s = 3e9; % between 2 and 4 GHz, for O-L comm
f_opt = 1e12; %Optical, TBD what frequency is most commonly used

% Antenna data
D_hga = 2.1; % [m] took from Orbilander paper, page 21
D_gs = 35; % [m] GS antenna, Cebreros considered (all the same, for ESTRACK DSA)
D_lga = 0.065; % [m] from datasheets of existing helix antennas (by RUAG space)
D_patch = 0.040; % [m] from datasheets of existing patch antenna (by ISIS – Innovative Solutions In Space)

% Power transmitted
P_on_board = 170; % [W] based on orbitander paper page 27
P_gs = P_on_board; %??????????????????

% NB. Since after the first iteration we got margins for EbNo which were
% very high -> we can lower the input power! Or tune some other parameter.
% Note: we are overdesigning sth because even with a very low input power,
% the margins are very large

% Antenna efficiency, depending on the shape
eta_gs = 0.55; % for parabolic reflectors
eta_hga = 0.55; % for parabolic reflectors
eta_lga = 0.7; % 0.7 if helix, 0.52 if horn

% Reference temperatures for noise computation
T_gs = 21; % [K] of the Ground Station, reference: DSN antenna temperature
% NB. also by setting this much higher (290K) the margins are still very
% large
T_hga = 250; % [K] of the HGA on board (common value <- from SSEO slides)
T_lga = 250; % [K] of the LGA on board 

%DATA RATES
% NOTE. This is just a first draft. It influences A LOT the margins we
% obtain. Once we have an idea of the visibility windows -> determine the
% real data rate -> margins will be better (no overdesigning) 
R_sci = 258.14*1e3; % data rate [bit/s], as coming from the PL, Data rate for Deep space communication
R_sci_OL = 774.42*1e3; % data rate [bit/s], as coming from the PL, Data rate for Deep space communication if HGA is on lander
R_LtoO = 8.51*1e3; % data rate [bit/s], as coming from the PL, Data rate for Lander - Orbiter communication
R_NN = 2.8*4 *1e3; % data rate [bit/s]as coming from the PL, Data rate for non-nominal (Housekeeping & Telemetry)
R_NN_OL = 2.8*4 *1e3; % data rate [bit/s]as coming from the PL, Data rate for non-nominal (Housekeeping & Telemetry)


%% SCIENCE PHASE -> COMMUNICATION MODE, with EARTH

R_down = R_sci * 1.05; % to include also attitude
[EbNo_down, SNR_down] = link_budget('to_Earth', 'downlink', 'atmosphere', R_down, max_distance_Earth, D_hga, D_gs, f_ka, P_on_board, eta_hga, eta_gs, T_gs);

marginEbNo_down = EbNo_down - min_EbNo_sci; % must be > 3dB
% marginSNR_down = SNR_down - min_SNR; % THIS IS COMMENTED BECAUSE I HAVE
% NO IDEA HOW TO RETRIEVE THE MINIMUM SNR

% if UPLINK 
R_up = 700; % [bps] from Orbilander paper, page 133
[EbNo_up, SNR_up] = link_budget('to_Earth', 'uplink', 'atmosphere', R_up, max_distance_Earth, D_gs, D_hga, f_ka, P_gs, eta_gs, eta_hga, T_hga);

marginEbNo_up = EbNo_up - min_EbNo_tel; % must be > 3dB
% marginSNR_up = SNR_up - min_SNR; 

%% SCIENCE PHASE -> COMMUNICATION MODE, with EARTH FROM LANDER

R_down_OL = R_sci_OL * 1.05; % to include also attitude and HK
[EbNo_down_OL, SNR_down_OL] = link_budget('to_Earth', 'downlink', 'atmosphere', R_down_OL, max_distance_Earth, D_hga, D_gs, f_ka, P_on_board, eta_hga, eta_gs, T_gs);

marginEbNo_down_OL = EbNo_down_OL - min_EbNo_sci; % must be > 3dB
% marginSNR_down_OL = SNR_down_OL - min_SNR; % THIS IS COMMENTED BECAUSE I HAVE
% NO IDEA HOW TO RETRIEVE THE MINIMUM SNR

% if UPLINK 
R_up_OL = 700*3; % [bps] from Orbilander paper, page 133, TIMES 3 FOR VISIBILITY WINDOW OF LANDER
[EbNo_up_OL, SNR_up_OL] = link_budget('to_Earth', 'uplink', 'atmosphere', R_up_OL, max_distance_Earth, D_gs, D_hga, f_ka, P_gs, eta_gs, eta_hga, T_hga);

marginEbNo_up_OL = EbNo_up_OL - min_EbNo_tel; % must be > 3dB
% marginSNR_up_OL = SNR_up_OL - min_SNR; 

%% SCIENCE PHASE -> COMMUNICATION BETWEEN ORBITER AND LANDER
% FROM LANDER TO ORBITER
[EbNo_LtoO, SNR_LtoO] = link_budget('L_to_O', 'uplink', 'no_atmosphere', R_LtoO, max_distance_OL, D_lga, D_lga, f_s, P_on_board, eta_lga, eta_lga, T_lga);

marginEbNo_LtoO = EbNo_LtoO - min_EbNo_sci; % must be > 3dB
% marginSNR_down = SNR_down - min_SNR; 

%% SCIENCE PHASE PRINT OUT RESULTS:

fprintf('--- SCIENCE PHASE ---\n')
fprintf('ORBITER DEEP SPACE COMMUNICATION CASE\n')
fprintf('The Eb/No for Deep Space downlink is: %.2f \nIts corresponding margin is: %.2f\n',EbNo_down,marginEbNo_down)
fprintf('The Eb/No for Deep Space uplink is: %.2f \nIts corresponding margin is: %.2f\n',EbNo_up,marginEbNo_up)
fprintf('The Eb/No Orbiter-lander communication is: %.2f \nIts corresponding margin is: %.2f\n',EbNo_LtoO,marginEbNo_LtoO)
fprintf('LANDER DEEP SPACE COMMUNICATION CASE\n')
fprintf('The Eb/No for Deep Space downlink is: %.2f \nIts corresponding margin is: %.2f\n',EbNo_down_OL,marginEbNo_down_OL)
fprintf('The Eb/No for Deep Space uplink is: %.2f \nIts corresponding margin is: %.2f\n',EbNo_up_OL,marginEbNo_up_OL)


%% NON-NOMINAL PHASE -> DEEP SPACE COMMUNICATION W/ HGA AT 1/10 POWER
% ORBITER DEEP SPACE COMMUNICATION CASE
P_NN = 1/10 * P_on_board;
%Downlink
[EbNo_down_NN_Ploss, SNR_down_NN_Ploss] = link_budget('to_Earth', 'downlink', 'atmosphere', R_NN, max_distance_Earth, D_hga, D_gs, f_ka, P_NN, eta_hga, eta_gs, T_gs);

marginEbNo_down_NN_Ploss = EbNo_down_NN_Ploss - min_EbNo_sci; % must be > 3dB
% marginSNR_down_NN_Ploss = SNR_down_NN_Ploss - min_SNR;

%Uplink 
[EbNo_up_NN_Ploss, SNR_up_NN_Ploss] = link_budget('to_Earth', 'uplink', 'atmosphere', R_up, max_distance_Earth, D_gs, D_hga, f_ka, P_gs, eta_gs, eta_hga, T_hga);

marginEbNo_up_NN_Ploss = EbNo_up_NN_Ploss - min_EbNo_tel; % must be > 3dB
% marginSNR_up_NN_Ploss = SNR_up_NN_Ploss - min_SNR; 

% LANDER DEEP SPACE COMMUNICATION CASE
%Downlink
[EbNo_down_NN_Ploss_OL, SNR_down_NN_Ploss_OL] = link_budget('to_Earth', 'downlink', 'atmosphere', R_NN_OL, max_distance_Earth, D_hga, D_gs, f_ka, P_NN, eta_hga, eta_gs, T_gs);

marginEbNo_down_NN_Ploss_OL = EbNo_down_NN_Ploss - min_EbNo_sci; % must be > 3dB
% marginSNR_down_NN_Ploss_OL = SNR_down_NN_Ploss_OL - min_SNR;

%Uplink 
[EbNo_up_NN_Ploss_OL, SNR_up_NN_Ploss_OL] = link_budget('to_Earth', 'uplink', 'atmosphere', R_up_OL, max_distance_Earth, D_gs, D_hga, f_ka, P_gs, eta_gs, eta_hga, T_hga);

marginEbNo_up_NN_Ploss_OL = EbNo_up_NN_Ploss_OL - min_EbNo_tel; % must be > 3dB
% marginSNR_up_NN_Ploss = SNR_up_NN_Ploss_OL - min_SNR; 

fprintf('\n\n--- NON-NOMINAL PHASE ---\n')
fprintf('\n- CASE : 1/10 POWER HGA DEEP SPACE COMM.\n')
fprintf('ORBITER DEEP SPACE COMMUNICATION CASE\n')
fprintf('The Eb/No for Deep Space downlink is: %.2f \nIts corresponding margin is: %.2f\n',EbNo_down_NN_Ploss,marginEbNo_down_NN_Ploss)
fprintf('The Eb/No for Deep Space uplink is: %.2f \nIts corresponding margin is: %.2f\n',EbNo_up_NN_Ploss,marginEbNo_up_NN_Ploss)
fprintf('LANDER DEEP SPACE COMMUNICATION CASE\n')
fprintf('The Eb/No for Deep Space downlink is: %.2f \nIts corresponding margin is: %.2f\n',EbNo_down_NN_Ploss_OL,marginEbNo_down_NN_Ploss_OL)
fprintf('The Eb/No for Deep Space uplink is: %.2f \nIts corresponding margin is: %.2f\n',EbNo_up_NN_Ploss_OL,marginEbNo_up_NN_Ploss_OL)

%% NON-NOMINAL PHASE -> DEEP SPACE COMMUNICATION W/ LGA 
% ORBITER DEEP SPACE COMMUNICATION CASE
%Downlink
[EbNo_down_NN_LGA, SNR_down_NN_LGA] = link_budget('to_Earth', 'downlink', 'atmosphere', R_NN, max_distance_Earth, D_lga, D_gs, f_ka, P_on_board, eta_lga, eta_gs, T_gs);

marginEbNo_down_NN_LGA = EbNo_down_NN_LGA - min_EbNo_sci; % must be > 3dB
% marginSNR_down_NN_LGA = SNR_down_NN_LGA - min_SNR;

%Uplink 
[EbNo_up_NN_LGA, SNR_up_NN_LGA] = link_budget('to_Earth', 'uplink', 'atmosphere', R_up, max_distance_Earth, D_gs, D_lga, f_ka, P_gs, eta_gs, eta_lga, T_lga);

marginEbNo_up_NN_LGA = EbNo_up_NN_LGA - min_EbNo_tel; % must be > 3dB
% marginSNR_up_NN_LGA = SNR_up_NN_LGA - min_SNR; 

% LANDER DEEP SPACE COMMUNICATION CASE
%Downlink
[EbNo_down_NN_LGA_OL, SNR_down_NN_LGA_OL] = link_budget('to_Earth', 'downlink', 'atmosphere', R_NN_OL, max_distance_Earth, D_lga, D_gs, f_ka, P_on_board, eta_lga, eta_gs, T_gs);

marginEbNo_down_NN_LGA_OL = EbNo_down_NN_LGA - min_EbNo_sci; % must be > 3dB
% marginSNR_down_NN_LGA_OL = SNR_down_NN_LGA_OL - min_SNR;

%Uplink 
[EbNo_up_NN_LGA_OL, SNR_up_NN_LGA_OL] = link_budget('to_Earth', 'uplink', 'atmosphere', R_up_OL, max_distance_Earth, D_gs, D_lga, f_ka, P_gs, eta_gs, eta_lga, T_lga);

marginEbNo_up_NN_LGA_OL = EbNo_up_NN_LGA_OL - min_EbNo_tel; % must be > 3dB
% marginSNR_up_NN_Ploss = SNR_up_NN_LGA_OL - min_SNR; 

fprintf('\n- CASE : LGA DEEP SPACE COMM.\n')
fprintf('ORBITER DEEP SPACE COMMUNICATION CASE\n')
fprintf('The Eb/No for Deep Space downlink is: %.2f \nIts corresponding margin is: %.2f\n',EbNo_down_NN_LGA,marginEbNo_down_NN_LGA)
fprintf('The Eb/No for Deep Space uplink is: %.2f \nIts corresponding margin is: %.2f\n',EbNo_up_NN_LGA,marginEbNo_up_NN_LGA)
fprintf('LANDER DEEP SPACE COMMUNICATION CASE\n')
fprintf('The Eb/No for Deep Space downlink is: %.2f \nIts corresponding margin is: %.2f\n',EbNo_down_NN_LGA_OL,marginEbNo_down_NN_LGA_OL)
fprintf('The Eb/No for Deep Space uplink is: %.2f \nIts corresponding margin is: %.2f\n',EbNo_up_NN_LGA_OL,marginEbNo_up_NN_LGA_OL)

%% NON-NOMINAL PHASE -> ORBITER-LANDER COMMUNICATION W/ PATCH 
% FROM LANDER TO ORBITER
[EbNo_LtoO_NN, SNR_LtoO_NN] = link_budget('L_to_O', 'uplink', 'no_atmosphere', R_LtoO, max_distance_OL, D_patch, D_lga, f_s, P_on_board, eta_lga, eta_lga, T_lga);

marginEbNo_LtoO_NN = EbNo_LtoO_NN - min_EbNo_sci; % must be > 3dB
% marginSNR_LtoO_NN = SNR_LtoO_NN - min_SNR; 
fprintf('\n- CASE : PATCH COMM FOR O-L\n')
fprintf('The Eb/No Orbiter-lander communication is: %.2f \nIts corresponding margin is: %.2f\n',EbNo_LtoO_NN,marginEbNo_LtoO_NN)


%% FOR ARCHITECTURE SELECTION (ONLY NOMINAL CASE)
clc;

% GENERAL NOTES COMING FROM PLAYING A BIT WITH PARAMETERS:
% - changing Tgs (temperature of the GS) lowers the margin (Downlink case
% for Orbilander goes from 74 to 49 dB if Tgs goes from 21K (DSN reference)
% to 250 (-20°C)). NB. The reference temperature used is the one of the
% cryogenically cooled amplifiers! So, it is more or less the same as for
% DSN (21K). -> CONCLUSION: I WOULD KEEP THE Tgs=21K.
% - if downlink from the HGA is in X-Band instead of Ka-band, the margin
% lowers from 74 to 71 dB (NOT RELEVANT, also because the science data are
% sent, typically, in Ka band for Deep Space Communication and the chosen
% GS, apart from NNO, are able to work in that frequency range). 
% - Power on board does not seem to be relevant -> if 170ish W power is used, based
% on Orbilander Reference, the margins are very large. But even if the
% power is lowered to 1/10 of this reference, margins still are very large (63 dB 
% instead of 74 dB for downlink in Ka band). 
% - if diameters are changed: nothing major changes in terms of margins.
% Moreover, the LGA dimension is taken from existing helix antenna, so it
% is NOT arbitrary, even if we select anothe one, dimensions will NOT
% change drastically. HGA instead is taken from the orbilander, if its D is
% lowered even to 0.5m (physically not reasonable, just look at past
% missions, the further you go, the biggere should be the antenna), the
% only margin that changes is the uplink. So the antenna that enters into
% the margin is really the recieving one! -> we could use this to lower the
% dimensions and still get good margins. 

% CONCLUSION: TUNE THE DATA RATES, WHICH ARE THE MOST RELEVANT PARAMETER TO
% HAVE BETTER RESULTS, AND THEN PLAY WITH THE DIMENSIONS OF THE ON-BOARD
% ANTENNAS and THE POWER NEEDED.

%-------------ORBILANDER---------------------------------------------------
% Downlink
[EbNo_down, ~] = link_budget('to_Earth', 'downlink', 'atmosphere', R_down, max_distance_Earth, D_hga, D_gs, f_ka, P_on_board, eta_hga, eta_gs, T_gs);
% Uplink (most critical is X-band, Ka returns much better margins)
[EbNo_up, ~] = link_budget('to_Earth', 'uplink', 'atmosphere', R_up, max_distance_Earth, D_gs, D_hga, f_x, P_gs, eta_gs, eta_hga, T_hga);

% Margins:
disp('ORBILANDER')
marginEbNo_down = EbNo_down - min_EbNo_sci 
marginEbNo_up = EbNo_up - min_EbNo_tel 

% ORBILANDER USES HGA AS MAIN COMMUNICATION ELEMENT, BOTH ON GROUND AND ON
% ORBIT.

%-------------S.O + S.L----------------------------------------------------

% Downlink
[EbNo_down, ~] = link_budget('to_Earth', 'downlink', 'atmosphere', R_down, max_distance_Earth, D_hga, D_gs, f_ka, P_on_board, eta_hga, eta_gs, T_gs);
% Uplink (most critical is X-band, Ka returns much better margins)
[EbNo_up, ~] = link_budget('to_Earth', 'uplink', 'atmosphere', R_up, max_distance_Earth, D_gs, D_hga, f_x, P_gs, eta_gs, eta_hga, T_hga);
% From Lander to Orbiter
[EbNo_LtoO, ~] = link_budget('L_to_O', 'uplink', 'no_atmosphere', R_LtoO, max_distance_OL, D_lga, D_lga, f_s, P_on_board, eta_lga, eta_lga, T_lga);

% Margins:
disp('S.O. + S.L.')
marginEbNo_down = EbNo_down - min_EbNo_sci 
marginEbNo_up = EbNo_up - min_EbNo_tel 
marginEbNo_LtoO = EbNo_LtoO - min_EbNo_sci

% LANDER and ORBITER have LGA to communicate in S-BAND (due to low
% distance), while ORBITER has HGA to communicate with Earth in Ka/X-BAND in nominal
% conditions.

%-------------NS.O + S.L---------------------------------------------------

% Downlink
[EbNo_down, ~] = link_budget('to_Earth', 'downlink', 'atmosphere', R_down, max_distance_Earth, D_hga, D_gs, f_ka, P_on_board, eta_hga, eta_gs, T_gs);
% Uplink (most critical is X-band, Ka returns much better margins)
[EbNo_up, ~] = link_budget('to_Earth', 'uplink', 'atmosphere', R_up, max_distance_Earth, D_gs, D_hga, f_x, P_gs, eta_gs, eta_hga, T_hga);
% From Lander to Orbiter
[EbNo_LtoO, ~] = link_budget('L_to_O', 'uplink', 'no_atmosphere', R_LtoO, max_distance_OL, D_lga, D_lga, f_s, P_on_board, eta_lga, eta_lga, T_lga);

% Margins:
disp('NS.O. + S.L.')
marginEbNo_down = EbNo_down - min_EbNo_sci 
marginEbNo_up = EbNo_up - min_EbNo_tel 
marginEbNo_LtoO = EbNo_LtoO - min_EbNo_sci

% LANDER and ORBITER have LGA to communicate in S-BAND (due to low
% distance), while ORBITER has HGA to communicate with Earth in Ka/X-BAND in nominal
% conditions (when the two modules are attached, it's the orbiter that comm
% with the Earth, using the data coming from the sampling of the lander attached).

%--------------N LANDERS---------------------------------------------------

% Downlink
[EbNo_down, ~] = link_budget('to_Earth', 'downlink', 'atmosphere', R_down, max_distance_Earth, D_hga, D_gs, f_ka, P_on_board, eta_hga, eta_gs, T_gs);
% Uplink (most critical is X-band, Ka returns much better margins)
[EbNo_up, ~] = link_budget('to_Earth', 'uplink', 'atmosphere', R_up, max_distance_Earth, D_gs, D_hga, f_x, P_gs, eta_gs, eta_hga, T_hga);
% From Landers (just consider one, since all of them will be equal) to Orbiter
[EbNo_LtoO, ~] = link_budget('L_to_O', 'uplink', 'no_atmosphere', R_LtoO, max_distance_OL, D_lga, D_lga, f_s, P_on_board, eta_lga, eta_lga, T_lga);

% Margins:
disp('MULTIPLE LANDERS')
marginEbNo_down = EbNo_down - min_EbNo_sci 
marginEbNo_up = EbNo_up - min_EbNo_tel 
marginEbNo_LtoO = EbNo_LtoO - min_EbNo_sci

% ORBITER has HGA to comm with Earth. It recieves data from the N LANDERS
% with a LGA in S-Band. 
% HOW TO TUNE DATA RATE? ORBITER will speak with 1 lander at a time. Then
% gather data, store them, and once it speaks with all landers -> downlink
% everything to Earth. This way, we are multiplying the data coming from
% the landers -> more data to Earth in downlink, which can help choosing
% architectures.

%% FUNCTION
function [EbNo, SNR] = link_budget(condition, direction, atmosphere, R, max_distance, Dtx, Drx, f, Pinput, eta_tx, eta_rx, Ts)
% NB. Uplink and Downlink is managed outside the function! But I left the
% input here also, just in case

% Define the wavelength
c = 2.98e8; % speed of light
lambda = c/f;

% Adjust the data rate based on the encoding and modulation techniques
alpha_mod = 2; % QSPK modulation
alpha_enc_convolution = 2; % convolutional encoding

R = R * alpha_enc_convolution / alpha_mod;

% Cable losses
Lc = -2; % [dB], between -1 and -3

if strcmp(condition, 'to_Earth')
    % Tx antenna gain computation
    Gtx = 10*log(pi*Dtx^2*eta_tx/lambda^2);

    % Ptx computation
    eta_TWTA = 0.58; % based on the input power at 150 W
    Ptx = Pinput * eta_TWTA;

    % EIRP evaluation
    EIRP = 20*log(Ptx + Gtx + Lc);

    % Free space losses
    Ls = 20*log(lambda/(4*pi*max_distance)); % [dB]

    % Atmospheric losses
    if strcmp(atmosphere, 'atmosphere')
        La = -0.2; % [dB] based on 32GHz of frequency, dry air
    else
        La = 0; % [dB] if no atmospheric losses are considere (L-O comm)
    end

    % Reciever gain
    Grx = 10*log(pi*Drx^2*eta_rx/lambda^2);

    % Compute beamwidth (parabolic reflector case)
    theta_tx = 65.3*lambda/Dtx;
    theta_rx = 65.3*lambda/Drx;

    % Pointing losses <- comes from the pointing budget!
    err = 0.1/100; % typical value 0.1% as from the slides
    Lp = -12*(err/theta_rx)^2;

    % Power recieved
    Prx = EIRP + Ls + La + Grx + Lp;

    % Reciever noise density (due to temperature)
    k = 1.38e-23; % Boltzmann's constant
    N0 = 10*log(k*Ts);

    % Eb/No computation
    EbNo = Prx - N0 - 10*log(R);

    % Compute carrier margin
    mu = 78*pi/180; % [rad] standard modulation index [45-90°]
    Pmod_loss = 20*log(cos(mu)); % power loss due to modulation

    Pcarrier = Prx + Pmod_loss;

    % Signal to Noise Ratio (of the carrier)
    SNR = Pcarrier - N0 - 10*log(R);

else
    L = 0.289; % [m] from datasheets of existing helix antennas
    % Tx antenna gain computation
    Ctx = pi*Dtx;
    Gtx = 10.3 + 10*log(Ctx^2*L/(lambda^3)); % peak gain value (see slide 29 of TMTC)

    % Ptx computation
    eta_TWTA = 0.58; % based on the input power at 150 W
    Ptx = Pinput * eta_TWTA;

    % EIRP evaluation
    EIRP = 20*log(Ptx + Gtx + Lc);

    % Free space losses
    Ls = 20*log(lambda/(4*pi*max_distance)); % [dB]

    % Atmospheric losses
    if strcmp(atmosphere, 'atmosphere')
        La = -0.2; % [dB] based on 32GHz of frequency, dry air
    else
        La = 0; % [dB] if no atmospheric losses are considere (L-O comm)
    end

    % Reciever gain (assume same helix antennas on board)
    Crx = pi*Drx;
    Grx = 10.3 + 10*log(Crx^2*L/(lambda^3)); 
    
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
    N0 = 10*log(k*Ts);

    % Eb/No computation
    EbNo = Prx - N0 - 10*log(R);

    % Compute carrier margin
    mu = 78*pi/180; % [rad] standard modulation index [45-90°]
    Pmod_loss = 20*log(cos(mu)); % power loss due to modulation

    Pcarrier = Prx + Pmod_loss;

    % Signal to Noise Ratio (of the carrier)
    SNR = Pcarrier - N0 - 10*log(R);
end

end