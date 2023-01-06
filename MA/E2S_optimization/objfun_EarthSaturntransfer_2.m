function DV = objfun_EarthSaturntransfer_2(nlpvar, planets_id, Ra_target, Rp_target)
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
% 17/11/2022 - Matteo Lusvarghi - First version
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
% 
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
t1 = nlpvar(1);
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
    [r1,v1] = kep2car(kep1, mu_S);
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
R = [i1,j1,k1];

% Rotate the vector in heliocentric inertial frame
v_inf_dep = R*v_inf_dep_local';

v_dep = v_inf_dep + v_E;

%% ADIMENSIONALIZATION
ephtimes = ephtimes/TU;
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


end





