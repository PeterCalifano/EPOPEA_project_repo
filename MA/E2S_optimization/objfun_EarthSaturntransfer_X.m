function DV = objfun_EarthSaturntransfer_X(nlpvar, planets_id, Ra_target, Rp_target)
%% PROTOTYPE
% DV = objfun_EarthSaturntransfer(var, N, planets_id, Ra_target, Rp_target, TU)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Ojective function for the Earth-Saturn interplanetary transfer
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% % var - variables vector
% planets_id - String of identifier of planets visited (dep and arrival included)
% Ra_target - Apoapsis Radius for capture orbit (SEE EstimateDVtoCapture.m)
% Rp_target - Periapsis Radius for capture orbit (SEE EstimateDVtoCapture.m)
% TU - Time Unit used to adimensionalize from seconds to days 
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% DV: [scalar] Sum of the DSMs and of the DVs given during the Powered Fly Bys
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 17/11/2022 - Matteo Lusvarghi - First (MATLAB) version
% 19/11/2022 - Pietro Califano - Test to use codegen
% -------------------------------------------------------------------------------------------------------------
%% EXPLANATION OF VAR VECTOR
%
% var = [
%       t_1 --> initial time [DAYS]
%       v_inf --> magnitude of the departure infinite velocity [Km/s]
%       u,v --> parameters related to RA and Dec at departure [ ]
%       ToF_1, ..., ToF_(N+1) --> ToF for each arc from 2 planets [DAYS]
%       alpha_1, ..., alpha_(N+1) --> adimensional number from 0 to 1 that
%           indicated the percentage of ToF_i at which the DSM is performed
%           (e.g. if alpha_2 = 0.5, the 2nd DSM is performed at ToF_2/2 of
%           the second arc)
%       rp_1, ... , rp_N --> Radii of pericenter for each FB, normalized by
%           the mean radius of the FB planet [Planetary Radii]
%       beta_1, ... , beta_N --> angle of hyperbolic plane of each FB [rad]
%       ];
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
%  astroConstants
%  uplanet
%  kep2car
%  EstimateDVtoCapture
%  lambertMR
%  PropagateorHelio_2BP
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% NEEDED: - propagator
%           - rotation of initial velocity
%           - FB function
% -------------------------------------------------------------------------------------------------------------

% FOR DEVELOPMENT ONLY
% Algorithm logic: 1) Input handling
%                  2) Evaluation of the planets positions at given time
%                  instants
%                  3) Computation of transfer arcs with DSMs in between
%                  4) Patching of the arcs via flyby evaluation at
%                  specified planets
% Currently evaluating the feasilibity of using parallelization (it may not
% be worth inside the obj function. Tests are required.) 

%% CONSTANTS DEFINITION
mu_S = astroConstants(4);

DU = astroConstants(2);
TU = 365.5;
TU2 = TU*3600*24;
VU = DU/TU2;

% Obtain the number of fly-bys from the string of visited planets
N = (length(nlpvar) - 6)/4;

%% EXTRACT VARIABLES FROM VAR
% Static allocation
t1 = nlpvar(1); % Extract launch date time
ephtimes = zeros(N+2,1); % ephemeris time in JULIAN DAYS
ephtimes(1) = t1; 

% Departure parameters (v_inf, u, v)
v_inf_dep = nlpvar(2);
u = nlpvar(3);
v = nlpvar(4);

% Time of flights of each arc and ephemerides times
tof = nlpvar(5 : 5 + N);
for i = 2 : N+2
    ephtimes(i) = ephtimes(i-1) + tof(i-1);
end

% Percentage of arc of DSM
alpha_DSM = nlpvar((N+6) : (N+6 + N));

% Pericenter radii and beta angles

rp_FB = nlpvar((2*N+7) : (2*N+6 + N));

beta_FB = nlpvar((3*N+7):(3*N+6 + N));

%% Retrieve positions of the planets
r_planets = zeros(3,N+2);
v_planets = zeros(3,N+2);
Rad_planets = zeros(1,N+2);
mu_planets = zeros(1,N+2);

for i = 1:(N+2)

    % Retrieve keplerian coordinates and convert them
    [kep1,~] = uplanet(ephtimes(i), planets_id(i));
    [r1,v1] = kep2car(kep1(1), kep1(2), kep1(3), kep1(4), kep1(5), kep1(6), mu_S);
    r_planets(:,i) = r1;
    v_planets(:,i) = v1;

    % Retrieve mean radius of each planet [Km]
    Rad_planets(i) = astroConstants(20 + planets_id(i));
    mu_planets(i) = astroConstants(10 + planets_id(i));
end

%% DEPARTURE 
% Position and velocity of the Earth
r_E = r_planets(:,1);
v_E = v_planets(:,1);

% Obtain angles to define the vector v_inf
th = 2*pi*u;
phi = acos(2*v-1)-pi/2;

% Departure v_inf in local frame
v_inf_dep_local = v_inf_dep*[cos(th)*cos(phi);sin(th)*cos(phi);sin(phi)]';

% Build rotation matrix
i1 = v_E / norm(v_E);
k1 = cross(r_E,v_E)/norm(cross(r_E,v_E));
j1 = cross(k1,i1);
RotMatECI2IRF = [i1,j1,k1];

% Rotate the vector in heliocentric inertial frame
v_inf_dep = RotMatECI2IRF*v_inf_dep_local';

v_dep = v_inf_dep + v_E;

%% ADIMENSIONALIZATION
% ephtimes = ephtimes/TU;
tof = tof/TU;
rp_FB = (rp_FB .* Rad_planets(2:end-1)) / DU;
v_dep = v_dep / VU;
r_planets = r_planets / DU;
v_planets = v_planets / VU;
mu_S = mu_S * TU2^2 / DU^3;
mu_planets = mu_planets * TU2^2 / DU^3;
Ra_target = Ra_target / DU;
Rp_target = Rp_target / DU;


%% Compute Transfer Arcs

% Set up the matrices to save the states at the DSMs 
v_DSM = zeros(3,N+1,2); 
r_DSM = zeros(3,N+1);

% Initialize heliocentric exit velocity from the FB
%v_FB = zeros(3,N+2,2); 
%v_FB(:,1,2) = v_dep;
v_FB_out = v_dep;

% Loop over the arcs
for i = 1:N+1

    % Initial state: position of the previous planet and outgoing velocity
    x_0 = [r_planets(:,i);v_FB_out];
    
    % Compute the tof until DSM and after DSM knowing its location alpha
    tof_1 = alpha_DSM(i)*tof(i);
    tof_2 = (1 - alpha_DSM(i)) * tof(i);

    % Propagation from initial state until the DSM
    %x_DSM = propagate(x_0,[ephtimes(i),ephtimes(i) + tof_1]);

    [r_propagated, v_propagated, ~] = PropagatorHelio_2BP(x_0, tof_1, mu_S);

    r_DSM(:,i) = r_propagated(:,end);
    v_DSM(:,i,1) = v_propagated(:,end);
    
    % Lambert arc from i-th DSM to (i+1)-th planet
    [~,~,~,~,v1_t,v2_t,~,~] = lambertMR(r_DSM(:,i), r_planets(:,i+1), tof_2, mu_S);
    v_DSM(:,i,2) = v1_t;
    v_FB_in = v2_t;
    
    % FLY_BY SEGMENT (performed only in the FB planets, not at Saturn)

    if i ~= (N+1)
        
        % Find gravitational constant of the planet
        mu_pl = mu_planets(i+1);

        % Obtain the incoming infinite velocity
        v_inf_in = v_FB_in' - v_planets(:,i+1);

        % Compute the fly-by
        %v_inf_out = perform_FB(v_inf_in,rp_FB(i),beta_FB(i));
        v_inf_out = FindVinfOut(v_inf_in, beta_FB(i), rp_FB(i), mu_pl, v_planets(:, i+1));

        % Obtain the outcoming infinite velocity which is the input for the
        % next iteration of the loop
        v_FB_out = v_planets(:,i+1) + v_inf_out;

    end


end

%% Compute CAPTURE DV: 

mu_end = mu_planets(end);

Vinf_entry = v_FB_in' - v_planets(:,end);

norm_Vinf_entry = norm(Vinf_entry);

[DV_capture] = EstimateDVtoCapture(norm_Vinf_entry, mu_end, Ra_target, Rp_target);

%% COMPUTE DV:

DV_DSM = zeros(1,N+1);

for i = 1:(N+1)

    DV_DSM(i) = norm(v_DSM(:,i,2) - v_DSM(:,i,1));

end

DV = (DV_capture + sum(DV_DSM)) * VU;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOCAL FUNCTIONS

function [dV_capture] = EstimateDVtoCapture(Vinf_entry, mu_main, Ra_target, Rp_target)
%% PROTOTYPE
% [DV_capture] = EstimateDVtoCapture(Vinf_entry, mu_main, Ra_target, Rp_target)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function estimating the Impulsive DeltaV required to achieve capture
% orbit specified by Apoapsis and Periapsis radii from hyperbolic entry at
% Vinf in the SOI of the main body. Two body problem assumed: be sure it is
% valid!
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% Vinf_entry: [scalar] Excess speed with respect to main attractor
% mu_main: [scalar] Gravitational Parameter of the main attractor (TBP model)
% Ra_target: [scalar] Apoapsis radius of capture orbit
% Rp_target: [scalar] Periapsis radius of capture orbit
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dV_capture: [scalar] Impulsive DeltaV at hyperbola pericentre to achieve
%                      specified capture orbit.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 01/11/2022 - Pietro Califano - Coded
% 17/11/2022 - Pietro Califano - Time computation deleted
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades


%% Function code
% Hypothesis: Planar orbits, entry hyperbolic leg shaped at will

e_target = (Ra_target - Rp_target)/(Ra_target + Rp_target);

% Velocity of hyperbolic orbit at Rp target
Vp_hyp = sqrt(Vinf_entry^2 + 2*mu_main/(Rp_target)); 

% Velocity at pericentre of target orbit
V_cap = sqrt(mu_main*(1 + e_target)/Rp_target);

% DV required to achieve elliptical capture orbit
dV_capture = Vp_hyp - V_cap;

end

function [r,v] = kep2car(a, e, i, OM, om, th, mu)

% kep2car.m - Conversion from Keplerian elements to Cartesian coordinates
%
% PROTOTYPE:
% [r, v] = kep2car(a, e, i, OM, om, th, mu)
%
% DESCRIPTION:
% Conversion from Keplerian elements to Cartesian coordinates. Angles in
% degrees.
%
% INPUT:
% a [1x1] Semi-major axis [km]
% e [1x1] Eccentricity [-]
% i [1x1] Inclination [rad]
% OM [1x1] RAAN [rad]
% om [1x1] Pericentre anomaly [rad]
% th [1x1] True anomaly [rad]
% mu [1x1] Gravitational parameter [km^3/s^2]
%
% OUTPUT:
% r [3x1] Position vector [km]
% v [3x1] Velocity vector [km/s]

% R1=rotazione di OM intorno al versore k
% R2=rotazione di i intorno al versore i'
% R3=rotazione di om intorno al versore k''
% R= prodotto matrici di rotazione (DA TRASPORRE)
%%
if nargin == 2
    mu = e;
    e = a(2);
    i = a(3);
    OM = a(4);
    om = a(5);
    th = a(6);
    a = a(1);
end
%% conversione in radianti
% i=deg2rad(i);
% om=deg2rad(om);
% OM=deg2rad(OM);
% th=deg2rad(th);

%% calcolo p ed |r|
p=a*(1-e^2);
r_nor=p/(1+e*cos(th));

%% calcolo r,v nel sistema perifocale (PF)
r_pf=r_nor.*[cos(th);sin(th);0];
v_pf=sqrt(mu/p)*[-sin(th);e+cos(th);0];

%% calcolo matrici di rotazione per arrivare al PF

R1=[cos(OM), sin(OM), 0; -sin(OM), cos(OM), 0; 0, 0, 1];
R2=[1, 0, 0; 0, cos(i), sin(i); 0, -sin(i), cos(i)];
R3=[cos(om), sin(om), 0; -sin(om), cos(om), 0; 0, 0, 1];
R=R3*R2*R1;

%% calcolo r e v nel sistema geocentrico equatoriale (ECI)

r=R'*r_pf;
v=R'*v_pf;

end

function Vinfout_IRF = FindVinfOut(vinfinvec_i, beta_i, Rp_i, mu_pl, Vpl)
% Vinfout = v infinity of outbound leg 

% vinfinvec_i = v infinity of inbound leg (vectorial)
% beta_i = angle of hyperbolic plane
% Rp_i = pericentre radius of the hyperbola
% mu_pl = grav. parameter of planet

vinfin_i = norm(vinfinvec_i);

delta = 2*asin(1/(1 + (Rp_i*vinfin_i.^2)./mu_pl));

% if delta > rad2deg(150)
%     warning('Computed turning angle is very high!')
% end

% Reference frame:
I = vinfinvec_i./vinfin_i;
Jvec = cross(I, Vpl);
J = Jvec./norm(Jvec);
K = cross(I, J);

Vinfout = vinfin_i.*[cos(delta); % I 
   cos(beta_i)*sin(delta); % J
   sin(beta_i)*sin(delta)]; % K

RotMatrixLocal2IRF = [I, J, K]; 

% Rotate the velocity from the Local "flyby" RF to the Inertial RF
Vinfout_IRF = RotMatrixLocal2IRF * Vinfout;


end

function  [kep, ksun] = uplanet(mjd2000, ibody)

% uplanet.m - Analytical ephemerides for planets
%
% PROTOTYPE:
%  [kep, ksun] = uplanet (mjd2000, ibody);
%
% DESCRIPTION:
%   Planetary orbital elements are restituited in a Sun-centred ecliptic 
%   system.
%   These ephemerides were succesfully compared with JPL/NAIF/SPICE
%   ephemerides using de405.bps.
%
%  INPUT :
%	mjd2000[1]  Time, modified Julian day since 01/01/2000, 12:00 noon
%               (MJD2000 = MJD-51544.5)
%	ibody[1]    Integer number identifying the celestial body (< 11)
%                   1:   Mercury
%                   2:   Venus
%                   3:   Earth
%                   4:   Mars
%                   5:   Jupiter
%                   6:   Saturn
%                   7:   Uranus
%                   8:   Neptune
%                   9:   Pluto
%                   10:  Sun
%
%  OUTPUT:
%	kep[6]    	Mean Keplerian elements of date
%                 kep = [a e i Om om theta] [km, rad]
%	ksun[1]     Gravity constant of the Sun [km^3/s^2]
%
%   Note: The ephemerides of the Moon are given by EphSS_kep, according to
%           to the algorithm in ephMoon.m
%
%  FUNCTIONS CALLED:
%   (none)
%
% AUTHOR:
%   P. Dysli, 1977
%
% PREVIOUS VERSION:
%   P. Dysli, 1977, MATLAB, uplanet.m
%       - Header and function name in accordance with guidlines.
%
%  CHANGELOG:
%   28/12/06, Camilla Colombo: tidied up
%   10/01/2007, REVISION, Matteo Ceriotti 
%   03/05/2008, Camilla Colombo: Case 11 deleted.
%   11/09/2008, Matteo Ceriotti, Camilla Colombo:
%       - All ephemerides shifted 0.5 days back in time. Now mjd2000 used
%           in this function is referred to 01/01/2000 12:00. In the old
%           version it was referred to 02/01/2000 00:00.
%       - Corrected ephemeris of Pluto.
%   04/10/2010, Camilla Colombo: Header and function name in accordance
%       with guidlines.
%
% -------------------------------------------------------------------------

DEG2RAD=pi/180;

G=6.67259e-20;
msun=1.988919445342813e+030;
ksun=msun*G;

KM  = 149597870.66;

% Static allocation
kep = zeros(1, 6);

%  T = JULIAN CENTURIES SINCE 31/12/1899 at 12:00
T   = (mjd2000 + 36525)/36525.00;
TT  = T*T;
TTT = T*TT;
kep(1)=T*0;

%
%  CLASSICAL PLANETARY ELEMENTS ESTIMATION IN MEAN ECLIPTIC OF DATE
%
switch round(ibody)
    
    %
    %  MERCURY
    %
    case 1
        kep(1) = 0.38709860;
        kep(2) = 0.205614210 + 0.000020460*T - 0.000000030*TT;
        kep(3) = 7.002880555555555560 + 1.86083333333333333e-3*T - 1.83333333333333333e-5*TT;
        kep(4) = 4.71459444444444444e+1 + 1.185208333333333330*T + 1.73888888888888889e-4*TT;
        kep(5) = 2.87537527777777778e+1 + 3.70280555555555556e-1*T +1.20833333333333333e-4*TT;
        XM   = 1.49472515288888889e+5 + 6.38888888888888889e-6*T;
        kep(6) = 1.02279380555555556e2 + XM*T;
    %
    %  VENUS
    %
    case 2
        kep(1) = 0.72333160;
        kep(2) = 0.006820690 - 0.000047740*T + 0.0000000910*TT;
        kep(3) = 3.393630555555555560 + 1.00583333333333333e-3*T - 9.72222222222222222e-7*TT;
        kep(4) = 7.57796472222222222e+1 + 8.9985e-1*T + 4.1e-4*TT;
        kep(5) = 5.43841861111111111e+1 + 5.08186111111111111e-1*T -1.38638888888888889e-3*TT;
        XM   = 5.8517803875e+4 + 1.28605555555555556e-3*T;
        kep(6) = 2.12603219444444444e2 + XM*T;
    %
    %  EARTH
    %
    case 3
        kep(1) = 1.000000230;
        kep(2) = 0.016751040 - 0.000041800*T - 0.0000001260*TT;
        kep(3) = 0.00;
    	kep(4) = 0.00;
      	kep(5) = 1.01220833333333333e+2 + 1.7191750*T + 4.52777777777777778e-4*TT + 3.33333333333333333e-6*TTT;
     	XM   = 3.599904975e+4 - 1.50277777777777778e-4*T - 3.33333333333333333e-6*TT;
     	kep(6) = 3.58475844444444444e2 + XM*T;
    %
    %  MARS
    %
    case 4
        kep(1) = 1.5236883990;
        kep(2) = 0.093312900 + 0.0000920640*T - 0.0000000770*TT;
        kep(3) = 1.850333333333333330 - 6.75e-4*T + 1.26111111111111111e-5*TT;
        kep(4) = 4.87864416666666667e+1 + 7.70991666666666667e-1*T - 1.38888888888888889e-6*TT - 5.33333333333333333e-6*TTT;
        kep(5) = 2.85431761111111111e+2 + 1.069766666666666670*T +  1.3125e-4*TT + 4.13888888888888889e-6*TTT;
        XM   = 1.91398585e+4 + 1.80805555555555556e-4*T + 1.19444444444444444e-6*TT;
        kep(6) = 3.19529425e2 + XM*T;
    %
    %  JUPITER
    %
    case 5
        kep(1) = 5.2025610;
        kep(2) = 0.048334750 + 0.000164180*T  - 0.00000046760*TT -0.00000000170*TTT;
        kep(3) = 1.308736111111111110 - 5.69611111111111111e-3*T +  3.88888888888888889e-6*TT;
        kep(4) = 9.94433861111111111e+1 + 1.010530*T + 3.52222222222222222e-4*TT - 8.51111111111111111e-6*TTT;
        kep(5) = 2.73277541666666667e+2 + 5.99431666666666667e-1*T + 7.0405e-4*TT + 5.07777777777777778e-6*TTT;
        XM   = 3.03469202388888889e+3 - 7.21588888888888889e-4*T + 1.78444444444444444e-6*TT;
        kep(6) = 2.25328327777777778e2 + XM*T;
    %
    %  SATURN
    %
	case 6
        kep(1) = 9.5547470;
        kep(2) = 0.055892320 - 0.00034550*T - 0.0000007280*TT + 0.000000000740*TTT;
        kep(3) = 2.492519444444444440 - 3.91888888888888889e-3*T - 1.54888888888888889e-5*TT + 4.44444444444444444e-8*TTT;
        kep(4) = 1.12790388888888889e+2 + 8.73195138888888889e-1*T -1.52180555555555556e-4*TT - 5.30555555555555556e-6*TTT;
        kep(5) = 3.38307772222222222e+2 + 1.085220694444444440*T + 9.78541666666666667e-4*TT + 9.91666666666666667e-6*TTT;
        XM   = 1.22155146777777778e+3 - 5.01819444444444444e-4*T - 5.19444444444444444e-6*TT;
        kep(6) = 1.75466216666666667e2 + XM*T;
    %
    %  URANUS
    %
    case 7
        kep(1) = 19.218140;
        kep(2) = 0.04634440 - 0.000026580*T + 0.0000000770*TT;
        kep(3) = 7.72463888888888889e-1 + 6.25277777777777778e-4*T + 3.95e-5*TT;
        kep(4) = 7.34770972222222222e+1 + 4.98667777777777778e-1*T + 1.31166666666666667e-3*TT;
        kep(5) = 9.80715527777777778e+1 + 9.85765e-1*T - 1.07447222222222222e-3*TT - 6.05555555555555556e-7*TTT;
        XM   = 4.28379113055555556e+2 + 7.88444444444444444e-5*T + 1.11111111111111111e-9*TT;
        kep(6) = 7.26488194444444444e1 + XM*T;
    %
    %  NEPTUNE
    %
    case 8
        kep(1) = 30.109570;
        kep(2) = 0.008997040 + 0.0000063300*T - 0.0000000020*TT;
        kep(3) = 1.779241666666666670 - 9.54361111111111111e-3*T - 9.11111111111111111e-6*TT;
        kep(4) = 1.30681358333333333e+2 + 1.0989350*T + 2.49866666666666667e-4*TT - 4.71777777777777778e-6*TTT;
        kep(5) = 2.76045966666666667e+2 + 3.25639444444444444e-1*T + 1.4095e-4*TT + 4.11333333333333333e-6*TTT;
        XM   = 2.18461339722222222e+2 - 7.03333333333333333e-5*T;
        kep(6) = 3.77306694444444444e1 + XM*T;
    %
    %  PLUTO
    %
    case 9
        kep(1) = 39.481686778174627;
        kep(2) = 2.4467e-001;
        kep(3) = 17.150918639446061;
        kep(4) = 110.27718682882954;
        kep(5) = 113.77222937912757;
        XM   = 4.5982945101558835e-008;
        kep(6) = 1.5021e+001 + XM*mjd2000*86400;
    %
    %  SUN
    %
    case 10
        kep = [0 0 0 0 0 0];
    %
    %  MOON (AROUND EARTH)
    %
%     case 11
%         msun=59.736e23;
%         ksun=msun*G;
%         kep(1) = 0.0025695549067660;
%         kep(2) = 0.0549004890;
%         kep(3) = 5.145396388888888890;
%         kep(4) = 2.59183275e+2  - 1.93414200833333333e+3*T + 2.07777777777777778e-3*TT + 2.22222222222222222e-6*TTT;
%         kep(4) = mod(kep(4) + 108e3, 360);
%         kep(5) = 7.51462805555555556e+1    + 6.00317604166666667e+3*T - 1.24027777777777778e-2*TT - 1.47222222222222222e-5*TTT;
%         kep(5) = mod(kep(5), 360);
%         XM   = 4.77198849108333333e+5    + 9.19166666666666667e-3*T  + 1.43888888888888889e-5*TT;
%         kep(6) = 2.96104608333333333e2     + XM                    *T;
    
    otherwise
	disp(ibody)
    if round(ibody)==11
        error('no planet in the list. For the Moon use EphSS_kep instead')
    else
        error('no planet in the list')
    end
end
  
%
%  CONVERSION OF AU INTO KM, DEG INTO RAD AND DEFINITION OF  XMU
%
 
kep(1)   = kep(1)*KM;       % a [km]
kep(3:6) = kep(3:6)*DEG2RAD;    % Transform from deg to rad
kep(6)   = mod(kep(6),2*pi);
% XMU  = (XM*DEG2RAD/(864000*365250))^2*kep(1)^3;
phi_uplanet    = kep(6);          % phi is the eccentric anomaly, uses kep(6)=M as a first guess
     
for j =1:5
    g       = kep(6)-(phi_uplanet-kep(2)*sin(phi_uplanet)); 
	g_primo = (-1+kep(2)*cos(phi_uplanet));
 	phi_uplanet     = phi_uplanet-g/g_primo;   % Computes the eccentric anomaly kep
end 

theta=2*atan(sqrt((1+kep(2))/(1-kep(2)))*tan(phi_uplanet/2));

kep(6)=theta;


end

function [r, v, timevec] = PropagatorHelio_2BP(initial_state, tf,mu_SUN)
%% PROTOTYPE
% [time, r, v] = twobp_solver(T, initial_state, mu, perturbation, dT)
% -------------------------------------------------------------------------
%% DESCRIPTION
% Solves the two body problem for a given planetary constant mu,
% initial position r and velocity v, using ode113
% -------------------------------------------------------------------------
%% INPUT
% initial_state: [6x1] Initial state of the orbiting body (rx ry rz vx vy vz) [L] [L/T^2]
% -------------------------------------------------------------------------
%% OUTPUT
% r: [Nx3] position of the body [L] 
% v: [Nx3] velocity of the body [L/T^2]
% time: [Nx1] vector of time instants
% Note: N is the number of points representing the orbit
% -------------------------------------------------------------------------
%% CONTRIBUTORS
% Pietro Califano
% -------------------------------------------------------------------------
%% CHANGELOG
% v1: 2BP solver for Heliocentric trajectories, 17/11/2022
% -------------------------------------------------------------------------
%% Next upgrades
%mu_SUN = 1.32712440017987e+11;

%% Function code

r_init = initial_state(1:3);
v_init = initial_state(4:6);


options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
sysF = @(t, state_vec) [state_vec(4:6);...
    (-mu_SUN/norm(state_vec(1:3)).^3)*state_vec(1:3)];

[timevec, output_state] = ode78(sysF, [0,tf], [r_init; v_init], options);

r = output_state(:, 1:3)';
v = output_state(:, 4:6)';

end

function [A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR(RI,RF,TOF,MU,orbitType,Nrev,Ncase,optionsLMR)

% lambertMR.m - Lambert's problem solver for all possible transfers
%   (multi-revolution transfer included).
%
% PROTOTYPE:
%   [A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR(RI,RF,TOF,MU,orbitType,Nrev,Ncase,optionsLMR)
%
% DESCRIPTION:
%   Lambert's problem solver for all possible transfers:
%       1- zero-revolution (for all possible types of orbits: circles, ellipses,
%       	parabolas and hyperbolas)
%       2- multirevolution case
%       3- inversion of the motion
%
%   1- ZERO-REVOLUTION LAMBERT'S PROBLEM
%
%   For the solution of Lambert's problem with number of revolution = 0 the
%   subroutine by Chris D'Souza is included here.
%   This subroutine is a Lambert algorithm which given two radius vectors
%   and the time to get from one to the other, it finds the orbit
%   connecting the two. It solves the problem using a new algorithm
%   developed by R. Battin. It solves the Lambert problem for all possible
%   types of orbits (circles, ellipses, parabolas and hyperbolas).
%   The only singularity is for the case of a transfer angle of 360 degrees,
%   which is a rather obscure case.
%   It computes the velocity vectors corresponding to the given radius
%   vectors except for the case when the transfer angle is 180 degrees
%   in which case the orbit plane is ambiguous (an infinite number of
%   transfer orbits exist).
% 
%   2- MULTIREVOLUTION LAMBERT'S PROBLEM
%
%   For the solution of Lambert's problem with Nrev>0 number of revolution,
%   Battin's formulation has been extended to accomodate N-revolution
%   transfer orbits, by following the paper: "Using Battin Mathod to obtain 
%   Multiple-revolution Lambert's Solutions" by Shen and Tsiotras.
%
%   When Nrev>0 the possible orbits are just ellipses.
%   If 0<=Nrev<=Nmax, there are two Nrev-revolution transfer orbits.
%   These two transfer orbits have different semi-major axis and they may 
%   be all combinations of large-e and small-e transfer orbits.
%   The Original Successive Substitution Method by Battin converges to one
%   of the two possible solution with a viable initial guest, however it
%   diverges from the other one. Then a Reversed Successive Substitution is
%   used to converge to the second solution.
%   A procedure is implemented in order to guarantee to provide initial
%   guesses in the convergence region. If Nrew exceeds the maximum number
%   of revolution an ERROR is given:
%   warning('off','lambertMR:SuccessiveSubstitutionDiverged') to take out
%   the warnings or use optionsLMR(1) = 0.
% 
%   3- INVERSION OF THE MOTION
% 
%   Direct or retrograde option can be selected for the transfer.
%   
%   The algorithm computes the semi-major axis, the parameter (semi-latus 
%   rectum), the eccentricity and the velocity vectors.
% 
%   NOTE: If ERROR occurs or the 360 or 180 degree transfer case is 
%   encountered. 
%
% INPUT:
%	RI[3]           Vector containing the initial position in Cartesian
%                   coordinates [L].
%	RF[3]           Vector containing the final position vector in
%                   Cartesian coordinates [L].
%	TOF[1]          Transfer time, time of flight [T].
%  	MU[1]           Planetary constant of the planet (mu = mass * G) [L^3/T^2]
%	orbitType[1]    Logical variable defining whether transfer is
%                       0: direct transfer from R1 to R2 (counterclockwise)
%                       1: retrograde transfer from R1 to R2 (clockwise)
%	Nrev[1]         Number of revolutions.
%                   if Nrev = 0 ZERO-REVOLUTION transfer is calculated
%                   if Nrev > 0 two transfers are possible. Ncase should be
%                          defined to select one of the two.
%	Ncase[1]        Logical variable defining the small-a or large-a option
%                   in case of Nrev>0:
%                       0: small-a option
%                       1: large-a option
%	optionsLMR[1]	lambertMR options:
%                    optionsLMR(1) = display options:
%                                    0: no display
%                                    1: warnings are displayed only when
%                                       the algorithm does not converge
%                                    2: full warnings displayed
%
% OUTPUT:
%	A[1]        Semi-major axis of the transfer orbit [L].
% 	P[1]        Semi-latus rectum of the transfer orbit [L].
%  	E[1]        Eccentricity of the transfer orbit.
%	ERROR[1]	Error flag
%                   0:	No error
%                   1:	Error, routine failed to converge
%                   -1:	180 degrees transfer
%                   2:  360 degrees transfer
%                   3:  the algorithm doesn't converge because the number 
%                       of revolutions is bigger than Nrevmax for that TOF
%                   4:  Routine failed to converge, maximum number of
%                       iterations exceeded.
%	VI[3]       Vector containing the initial velocity vector in Cartesian
%               coordinates [L/T].
%	VT[1]		Vector containing the final velocity vector in Cartesian
%               coordinates [L/T].
%	TPAR[1] 	Parabolic flight time between RI and RF [T].
%	THETA[1]	Transfer angle [radians].
%
% NOTE: The semi-major axis, positions, times, and gravitational parameter
%       must be in compatible units.
%
% CALLED FUNCTIONS:
%   qck, h_E (added at the bottom of this file)
%
% REFERENCES:
%   - Shen and Tsiotras, "Using Battin method to obtain Multiple-Revolution
%       Lambert's solutions".
%   - Battin R., "An Introduction to the Mathematics and Methods of
%       Astrodynamics, Revised Edition", 1999.
%
% FUTURE DEVELOPMENT:
%   - 180 degrees transfer indetermination
%   - 360 degrees transfer singularity
%   - Nmax number of max revolution for a given TOF:
%     work in progress - Camilla Colombo
%
% ORIGINAL VERSION:
%   Chris D'Souza, 20/01/1989, MATLAB, lambert.m
%       verified by Darrel Monroe, 10/25/90
%       - Labert.m solved only direct transfer, without multi-revolution
%         option
%
% AUTHOR:
%   Camilla Colombo, 10/11/2006, MATLAB, lambertMR.m
%
% CHANGELOG:
%   13/11/2006, Camilla Colombo: added ERROR = 3 if Nrev > NrevMAX
%	21/11/2006, Camilla Colombo: added another case of ERROR = 3 (index
%   	N3) corresponding to the limit case when small-a solution = large-a
%       solution. No solution is given in this case.
%	06/08/2007, Camilla Colombo: optionsLMR added as an input
%	28/11/2007, Camilla Colombo: minor changes
%   29/01/2009, Matteo Ceriotti:
%       - Introduced variable for maximum number of iterations nitermax.
%       - Corrected final check on maximum number of iterations exceeded, from
%           "==" to ">=" (if N1 >= nitermax || N >= nitermax).
%       - Increased maxumum number of iterations to 2000, not to lose some
%           solutions.
%       - In OSS loop, added check for maximum number of iterations exceeded,
%           which then sets checkNconvOSS = 0.
%       - Changed the way of coumputing X given Y1 in RSS. Now the
%           Newton-Raphson method with initial guess suggested by Shen,
%           Tsiotras is used. This should guarantee convergence without the
%           need of an external zero finder (fsolve).
%       - Changed absolute tolerance into relative tolerance in all loops X0-X.
%           Now the condition is: while "abs(X0-X) >= abs(X)*TOL+TOL".
%       - Added return immediately when any error is detected.
%       - Moved check on 4*TOF*LAMBDA==0 after computing LAMBDA.
%       - Moved check on THETA==0 || THETA==2*PI after computing THETA.
%       - Added error code 4 (number of iterations exceeded).
%       - Removed variable Nwhile, as not strictly needed.
%       - Removed variable PIE=pi.
%   29/01/2009, REVISION: Matteo Ceriotti
%   21/07/2009, Matteo Ceriotti, Camilla Colombo:
%       added condition to detect case 180 degrees transfer indetermination
%   30/01/2010, Camilla Colombo: Header and function name in accordance
%       with guidlines.
%
% Note: Please if you have got any changes that you would like to be done,
%   do not change the function, please contact the author.
%
% -------------------------------------------------------------------------

% Check inputs
if nargin < 8
    optionsLMR = 0;
    if nargin < 6
        Nrev = 0;
        if nargin < 5
            orbitType = 0;
            if nargin < 4
                error('Not enough input arguments. See lambertMR.');
            end
        end
    end
end

nitermax = 2000; % Maximum number of iterations for loops
TOL = 1e-14;

TWOPI=2*pi;

% Reset
A=0;P=0;E=0;VI=[0,0,0];VF=[0,0,0];

% ----------------------------------
% Compute the vector magnitudes and various cross and dot products

RIM2   = dot(RI,RI);
RIM    = sqrt(RIM2);
RFM2   = dot(RF,RF);
RFM    = sqrt(RFM2);
CTH    = dot(RI,RF)/(RIM*RFM);
CR     = cross(RI,RF);
STH    = norm(CR)/(RIM*RFM);

% Choose angle for up angular momentum
switch orbitType
    case 0 % direct transfer
        if CR(3) < 0 
            STH = -STH;
        end
    case 1 % retrograde transfer
        if CR(3) > 0 
            STH = -STH;
        end
    otherwise
		error('%d is not an allowed orbitType',orbitType);
end
        
THETA  = qck(atan2(STH,CTH));
% if abs(THETA - pi) >= 0.01
if THETA == TWOPI || THETA==0
    ERROR = 2;
    A=0; P=0; E=0; VI=[0,0,0]; VF=[0,0,0]; TPAR=0; THETA=0;
    return
end

B1     = sign(STH); if STH == 0; B1 = 1; end;

% ----------------------------------
% Compute the chord and the semi-perimeter

C= sqrt(RIM2 + RFM2 - 2*RIM*RFM*CTH);
S= (RIM + RFM + C)/2;
BETA   = 2*asin(sqrt((S-C)/S));
PMIN   = TWOPI*sqrt(S^3/(8*MU));
TMIN   = PMIN*(pi-BETA+sin(BETA))/(TWOPI);
LAMBDA = B1*sqrt((S-C)/S);

if 4*TOF*LAMBDA == 0 || abs((S-C)/S) < TOL
    ERROR = -1;
    A=0; P=0; E=0; VI=[0,0,0]; VF=[0,0,0]; TPAR=0; THETA=0;
    return
end

% ----------------------------------
% Compute L carefully for transfer angles less than 5 degrees

if THETA*180/pi <= 5
   W   = atan((RFM/RIM)^.25) - pi/4;
   R1  = (sin(THETA/4))^2;
   S1  = (tan(2*W))^2;
   L   = (R1+S1)/(R1+S1+cos(THETA/2));
else
   L   = ((1-LAMBDA)/(1+LAMBDA))^2;
end

M= 8*MU*TOF^2/(S^3*(1+LAMBDA)^6);
TPAR   = (sqrt(2/MU)/3)*(S^1.5-B1*(S-C)^1.5);
L1     = (1 - L)/2;

CHECKFEAS = 0;
N1 = 0;
N_lambMR = 0;

if Nrev == 0
    % ----------------------------------
    % Initialize values of y, n, and x

    Y= 1;
    N_lambMR= 0;
    N1=0;
    ERROR  = 0;
    % CHECKFEAS=0;

    if (TOF-TPAR) <= 1e-3
        X0  = 0;
    else
        X0  = L;
    end

    X= -1.e8;

    % ----------------------------------
    % Begin iteration
    
    % ---> CL: 26/01/2009, Matteo Ceriotti: 
    %       Changed absolute tolerance into relative tolerance here below.
    while (abs(X0-X) >= abs(X)*TOL+TOL) && (N_lambMR <= nitermax)
        N_lambMR   = N_lambMR+1;
        X   = X0;
        ETA = X/(sqrt(1+X) + 1)^2;
        CHECKFEAS=1;

        % ----------------------------------
        % Compute x by means of an algorithm devised by
        % Gauticci for evaluating continued fractions by the
        % 'Top Down' method
        
        DELTA = 1;
        U     = 1;
        SIGMA = 1;
        M1    = 0;

        while abs(U) > TOL && M1 <= nitermax
            M1    = M1+1;
            GAMMA = (M1 + 3)^2/(4*(M1+3)^2 - 1);
            DELTA = 1/(1 + GAMMA*ETA*DELTA);
            U     = U*(DELTA - 1);
            SIGMA = SIGMA + U;
        end

        C1 = 8*(sqrt(1+X)+1)/(3+1/(5 + ETA + (9*ETA/7)*SIGMA));

        % ----------------------------------
        % Compute H1 and H2
        
        if N_lambMR == 1
            DENOM = (1 + 2*X + L)*(3*C1 + X*C1 +4*X);
            H1 = (L+X)^2*(C1 + 1 + 3*X)/DENOM;
            H2 = M*(C1+X-L)/DENOM;
        else
            QR = sqrt(L1^2 + M/Y^2);
            XPLL = QR - L1;
            LP2XP1 = 2*QR;
            DENOM = LP2XP1*(3*C1 + X*C1+4*X);
            H1 = ((XPLL^2)*(C1 + 1 + 3*X))/DENOM;
            H2 = M*(C1+X-L)/DENOM;
        end
        
        B = 27*H2/(4*(1+H1)^3);
        U = -B/(2*(sqrt(B+1)+1));

        % ----------------------------------
        % Compute the continued fraction expansion K(u)
        % by means of the 'Top Down' method
        
        % Y can be computed finding the roots of the formula and selecting
        % the real one:
        % y^3 - (1+h1)*y^2 - h2 = 0     (7.113) Battin
        %
        % Ycami_ = roots([1 -1-H1 0 -H2])
        % kcami = find( abs(imag(Ycami_)) < eps );
        % Ycami = Ycami_(kcami)

        DELTA = 1;
        U0 = 1;
        SIGMA = 1;
        N1 = 0;

        while N1 < nitermax && abs(U0) >= TOL
            if N1 == 0
                GAMMA = 4/27;
                DELTA = 1/(1-GAMMA*U*DELTA);
                U0 = U0*(DELTA - 1);
                SIGMA = SIGMA + U0;
            else
                for I8 = 1:2
                    if I8 == 1
                        GAMMA = 2*(3*N1+1)*(6*N1-1)/(9*(4*N1 - 1)*(4*N1+1));
                    else
                        GAMMA = 2*(3*N1+2)*(6*N1+1)/(9*(4*N1 + 1)*(4*N1+3));
                    end
                    DELTA = 1/(1-GAMMA*U*DELTA);
                    U0 = U0*(DELTA-1);
                    SIGMA = SIGMA + U0;
                end
            end

            N1 = N1 + 1;
        end

        KU = (SIGMA/3)^2;
        Y = ((1+H1)/3)*(2+sqrt(B+1)/(1-2*U*KU));    % Y = Ycami
        
        X0 = sqrt(((1-L)/2)^2+M/Y^2)-(1+L)/2;
        % fprintf('n= %d, x0=%.14f\n',N,X0);
    end
    
% MULTIREVOLUTION
elseif (Nrev > 0) && (4*TOF*LAMBDA~=0) %(abs(THETA)-pi > 0.5*pi/180)

    checkNconvRSS = 1;
    checkNconvOSS = 1;
    N3 = 1;
    
    while N3 < 3
        
        if Ncase == 0 || checkNconvRSS == 0

            % - Original Successive Substitution -
            % always converges to xL - small a

            % ----------------------------------
            % Initialize values of y, n, and x
            
            Y= 1;
            N_lambMR= 0;
            N1=0;
            ERROR = 0;
            % CHECKFEAS = 0;
%             if (TOF-TPAR) <= 1e-3
%                 X0 = 0;
%             else
            if checkNconvOSS == 0
                X0 = 2*X0;
                checkNconvOSS = 1;
                % see p. 11 USING BATTIN METHOD TO OBTAIN 
                % MULTIPLE-REVOLUTION LAMBERT'S SOLUTIONS - Shen, Tsiotras
            elseif checkNconvRSS == 0;
                % X0 is taken from the RSS
            else
                X0 = L;
            end

            X = -1.e8;

            % ----------------------------------
            % Begin iteration
            			
            % ---> CL: 26/01/2009,Matteo Ceriotti 
            %   Changed absolute tolerance into relative tolerance here
            %   below.
            while (abs(X0-X) >= abs(X)*TOL+TOL) && (N_lambMR <= nitermax)
                N_lambMR   = N_lambMR+1;
                X   = X0;
                ETA = X/(sqrt(1+X) + 1)^2;
                CHECKFEAS = 1;

                % ----------------------------------
                % Compute x by means of an algorithm devised by
                % Gauticci for evaluating continued fractions by the
                % 'Top Down' method
                

                DELTA = 1;
                U     = 1;
                SIGMA = 1;
                M1    = 0;

                while abs(U) > TOL && M1 <= nitermax
                    M1    = M1+1;
                    GAMMA = (M1 + 3)^2/(4*(M1+3)^2 - 1);
                    DELTA = 1/(1 + GAMMA*ETA*DELTA);
                    U     = U*(DELTA - 1);
                    SIGMA = SIGMA + U;
                end

                C1 = 8*(sqrt(1+X)+1)/(3+1/(5 + ETA + (9*ETA/7)*SIGMA));

                % ----------------------------------
                % Compute H1 and H2
                
                if N_lambMR == 1
                    DENOM = (1 + 2*X + L)*(3*C1 + X*C1 +4*X);
                    H1 = (L+X)^2*(C1 + 1 + 3*X)/DENOM;
                    H2 = M*(C1+X-L)/DENOM;
                else
                    QR = sqrt(L1^2 + M/Y^2);
                    XPLL = QR - L1;
                    LP2XP1 = 2*QR;
                    DENOM = LP2XP1*(3*C1 + X*C1+4*X);
                    H1 = ((XPLL^2)*(C1 + 1 + 3*X))/DENOM;
                    H2 = M*(C1+X-L)/DENOM;
                end

                H3 = M*Nrev*pi/(4*X*sqrt(X));
                H2 = H3+H2;

                B = 27*H2/(4*(1+H1)^3);
                U = -B/(2*(sqrt(B+1)+1));

                % ----------------------------------
                % Compute the continued fraction expansion K(u)
                % by means of the 'Top Down' method
                
                % Y can be computed finding the roots of the formula and selecting
                % the real one:
                % y^3 - (1+h1)*y^2 - h2 = 0     (7.113) Battin
                %
                % Ycami_ = roots([1 -1-H1 0 -H2])
                % kcami = find( abs(imag(Ycami_)) < eps );
                % Ycami = Ycami_(kcami)

                DELTA = 1;
                U0 = 1;
                SIGMA = 1;
                N1 = 0;

                while N1 < nitermax && abs(U0) >= TOL
                    if N1 == 0
                        GAMMA = 4/27;
                        DELTA = 1/(1-GAMMA*U*DELTA);
                        U0 = U0*(DELTA - 1);
                        SIGMA = SIGMA + U0;
                    else
                        for I8 = 1:2
                            if I8 == 1
                                GAMMA = 2*(3*N1+1)*(6*N1-1)/(9*(4*N1 - 1)*(4*N1+1));
                            else
                                GAMMA = 2*(3*N1+2)*(6*N1+1)/(9*(4*N1 + 1)*(4*N1+3));
                            end
                            DELTA = 1/(1-GAMMA*U*DELTA);
                            U0 = U0*(DELTA-1);
                            SIGMA = SIGMA + U0;
                        end
                    end

                    N1 = N1 + 1;
                end

                KU = (SIGMA/3)^2;
                Y = ((1+H1)/3)*(2+sqrt(B+1)/(1-2*U*KU));	% Y = Ycami
                if Y > sqrt(M/L)
                    if optionsLMR(1) == 2
                        warning('lambertMR:SuccessiveSubstitutionDiverged',...
                                ['Original Successive Substitution is diverging\n'...
                                '-> Reverse Successive Substitution used to find the proper XO.\n']);
                    end
                    checkNconvOSS = 0;
                    break
                end
                
                X0 = sqrt(((1-L)/2)^2+M/Y^2)-(1+L)/2;
                % fprintf('N: %d X0: %.14f\n',N,X0);
            end
            
            % When 2 solutions exist (small and big a), the previous loop
            % must either converge or diverge because Y > sqrt(M/L) at some
            % point. Thus, the upper bound on the number of iterations
            % should not be necessary. Though, nothing can be said in the
            % case tof<tofmin and so no solution exist. In this case, an
            % upper bound on number of iterations could be needed.
            
            if N_lambMR >= nitermax % Checks if previous loop ended due to maximum number of iterations
                if optionsLMR(1) == 2
                    warning('lambertMR:SuccessiveSubstitutionExceedMaxIter',...
                            ['Original Successive Substitution exceeded max number of iteration\n'...
                            '-> Reverse Successive Substitution used to find the proper XO.\n']);
                end
                checkNconvOSS = 0;
            end
        end
        if (Ncase == 1 || checkNconvOSS == 0) && ~(checkNconvRSS == 0 && checkNconvOSS == 0)

            % - Reverse Successive Substitution -
            % always converges to xR - large a

            % ----------------------------------
            % Initialize values of y, n, and x
            
            N_lambMR = 0;
            N1 = 0;
            ERROR  = 0;
            % CHECKFEAS=0;
            if checkNconvRSS == 0;
                X0 = X0/2; % XL/2
                checkNconvRSS = 1;
                % see p. 11 USING BATTIN METHOD TO OBTAIN 
                % MULTIPLE-REVOLUTION LAMBERT'S SOLUTIONS - Shen, Tsiotras
            elseif checkNconvOSS == 0
                % X0 is taken from the OSS
            else
                X0 = L;
            end

            X = -1.e8;

            % ----------------------------------
            % Begin iteration
            
            % ---> CL: 26/01/2009, Matteo Ceriotti
            %   Changed absolute tolerance into relative tolerance here
            %   below.
            while (abs(X0-X) >= abs(X)*TOL+TOL) && (N_lambMR <= nitermax)
                N_lambMR = N_lambMR+1;
                X = X0;
                CHECKFEAS=1;

                Y = sqrt(M/((L+X)*(1+X))); % y1 in eq. (8a) in Shen, Tsiotras

                if Y < 1
                    if optionsLMR(1) == 2
                        warning('lambertMR:SuccessiveSubstitutionDiverged',...
                                ['Reverse Successive Substitution is diverging\n' ...
                                '-> Original Successive Substitution used to find the proper XO.\n']);
                    end
                    checkNconvRSS = 0;
                    break
                end
                
                % ---> CL: 27/01/2009, Matteo Ceriotti
                %   This is the Newton-Raphson method suggested by USING
                %   BATTIN METHOD TO OBTAIN MULTIPLE-REVOLUTION LAMBERT'S
                %   SOLUTIONS - Shen, Tsiotras
                
                % To assure the Newton-Raphson method to be convergent
                Erss = 2*atan(sqrt(X));
                while h_E(Erss,Y,M,Nrev) < 0
                    Erss = Erss/2;
                end
                
                Nnew = 1;
                Erss_old = -1.e8;
                
                % The following Newton-Raphson method should always
                % converge, given the previous first guess choice,
                % according to the paper. Therefore, the condition on
                % number of iterations should not be neccesary. It could be
                % necessary for the case tof < tofmin.
                while (abs(Erss-Erss_old) >= abs(Erss)*TOL+TOL) && Nnew < nitermax
                    Nnew = Nnew+1;
                    [h, dh] = h_E(Erss,Y,M,Nrev);
                    Erss_old = Erss;
                    Erss = Erss - h/dh;
                    % fprintf('Nnew: %d Erss: %.16f h_E: %.16f\n',Nnew,Erss,h);
                end
                if Nnew >= nitermax
                    if optionsLMR(1) ~= 0
                        warning('lambertMR:NewtonRaphsonIterExceeded', 'Newton-Raphson exceeded max iterations.\n');
                    end
                end
                X0 = tan(Erss/2)^2;
            end
        end
        if checkNconvOSS == 1 && checkNconvRSS == 1
            break
        end
        
        if checkNconvRSS == 0 && checkNconvOSS == 0
            if optionsLMR ~=0
                warning('lambertMR:SuccessiveSubstitutionDiverged',...
                        ['Both Original Successive Substitution and Reverse ' ...
                        'Successive Substitution diverge because Nrev > NrevMAX.\n' ...
                        'Work in progress to calculate NrevMAX.\n']);
            end
            ERROR = 3;
            A=0; P=0; E=0; VI=[0,0,0]; VF=[0,0,0]; TPAR=0; THETA=0;
            return
        end
        
        N3 = N3+1;
    end
    
    if N3 == 3
        if optionsLMR ~=0
            warning('lambertMR:SuccessiveSubstitutionDiverged',...
                    ['Either Original Successive Substitution or Reverse ' ...
                    'Successive Substitution is always diverging\n' ...
                    'because Nrev > NrevMAX or because large-a solution = small-a solution (limit case).\n' ...
                    'Work in progress to calculate NrevMAX.\n']);
        end
        ERROR = 3;
        A=0; P=0; E=0; VI=[0,0,0]; VF=[0,0,0]; TPAR=0; THETA=0;
        return
    end
end

% ----------------------------------
% Compute the velocity vectors

if CHECKFEAS == 0
    ERROR = 1;
    A=0; P=0; E=0; VI=[0,0,0]; VF=[0,0,0]; TPAR=0; THETA=0;
    return
end

if N1 >= nitermax || N_lambMR >= nitermax
    ERROR = 4;
    if optionsLMR ~=0
        disp('Lambert algorithm has not converged, maximum number of iterations exceeded.');
    end
    A=0; P=0; E=0; VI=[0,0,0]; VF=[0,0,0]; TPAR=0; THETA=0;
    return
end

CONST = M*S*(1+LAMBDA)^2;
A = CONST/(8*X0*Y^2);

R11 = (1 + LAMBDA)^2/(4*TOF*LAMBDA);
S11 = Y*(1 + X0);
T11 = (M*S*(1+LAMBDA)^2)/S11;

VI(1:3) = -R11*(S11*(RI(1:3)-RF(1:3))-T11*RI(1:3)/RIM);
VF(1:3) = -R11*(S11*(RI(1:3)-RF(1:3))+T11*RF(1:3)/RFM);

P = (2*RIM*RFM*Y^2*(1+X0)^2*sin(THETA/2)^2)/CONST;
E = sqrt(1 - P/A);

return

% -------------------------------------------------------------------------


function [angle] = qck(angle)

% qck.m - Reduce an angle between 0 and 2*pi
%
% PROTOTYPE:
%   [angle]=qck(angle)
%
% DESCRIPTION:
%   This function takes any angle and reduces it, if necessary,
% 	so that it lies in the range from 0 to 2 PI radians.
% 
% INPUTS:
%   ANGLE[1]    Angle to be reduced (in radians)
% 
% OUTPUTS:
%   QCK[1]      The angle reduced, if necessary, to the range
%               from 0 to 2 PI radians (in radians)
% 
% CALLED FUNCTIONS:
%   pi (from MATLAB)
%
% AUTHOR:
%   W.T. Fowler, July, 1978
%
% CHANGELOG:
%   8/20/90, REVISION: Darrel Monroe
%
% -------------------------------------------------------------------------

twopi = 2*pi;
 
diff = twopi * (fix(angle/twopi) + min([0,sign(angle)]));

angle = angle -diff;

end
% -------------------------------------------------------------------------


function [h, dh] = h_E(E, y, m, Nrev)

% h_E.m - Equation of multirevolution Lambert's problem h = h(E).
%
% PROTOTYPE:
%   [h, dh] = h_E(E, y, m, Nrev)
%
% DESCRIPTION:
%   Equation of multirevolution Lambert's problem:
%   h(E) = (Nrev*pi + E - sin(E)) / tan(E/2)^3 - 4/m * (y^3 - y^2)
%   See: "USING BATTIN METHOD TO OBTAIN MULTIPLE-REVOLUTION LAMBERT'S 
%      SOLUTIONS", Shen, Tsiotras, pag. 12
%
% INPUT
%   E, y, m, Nrev   See paper for detailed description.
%
% OUTPUT
%   h               Value of h(E).
%   dh              Value of dh(E)/dE.
%
% ORIGINAL VERSION:
%   Camilla Colombo, 20/02/2006, MATLAB, cubicN.m
%
% AUTHOR:
%   Matteo Ceriotti, 27/01/2009
%   - changed name of cubicN.m and added at the bottom of lambertMR.m file
%
% -------------------------------------------------------------------------

tanE2 = tan(E/2);
h = (Nrev*pi + E - sin(E)) / tanE2^3 - 4/m * (y^3 - y^2);

if nargout > 1  % two output arguments
    % h'(E)
    dh = (1-cos(E))/tanE2^3 - 3/2*(Nrev*pi+E-sin(E))*sec(E/2)^2 / tanE2^4;
end

end

end

end





