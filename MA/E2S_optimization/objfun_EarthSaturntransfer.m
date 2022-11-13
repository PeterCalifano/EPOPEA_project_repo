function DV = objfun_EarthSaturntransfer(var, N, planets_id, Ra_target, Rp_target, TU)
%% PROTOTYPE
% DV = objfun_EarthSaturntransfer(var, N, planets_id, Ra_target, Rp_target, TU)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Ojective function for the Earth-Saturn interplanetary transfer
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% % var - variables vector
% N - Number of fly-bys
% planets_id - String of identifier of planets visited (dep and arrival included)
% Ra_target - Apoapsis Radius for capture orbit (SEE EstimateDVtoCapture.m)
% Rp_target - Periapsis Radius for capture orbit (SEE EstimateDVtoCapture.m)
% TU - Time Unit used to adimensionalize (It's considered to be 1 year)
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% DV: [scalar] Sum of the DSMs and of the DVs given during the Powered Fly Bys
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 12/11/2022 - Matteo Lusvarghi - First version
% -------------------------------------------------------------------------------------------------------------
%% EXPLANATION OF VAR VECTOR
%
% var = [t_1, --> time of departure
%       ToF_1,...,ToF_(N+1), --> time of flight of each arc
%       r_1, th_1, z_1, t_DSM_1, --> coordinates in cylindrical of each DSM
%       ...
%       r_(N+1), th_(N+1), z_(N+1), t_DSM_(N+1) ]
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% 
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
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

%% EXTRACT VARIABLES FROM VAR
% Static allocation
t1 = var(1);
tof = zeros(N+1,1);
ephtimes = zeros(N+2,1); % ephemeris time in [s]
ephtimes(1) = t1;
cyl_DSM = zeros(3,N+1);
times_DSM = zeros(N+1,1);


for i = 2:N+2
    % Time of flights
    tof(i-1) = var(i);
    ephtimes(i) = ephtimes(i-1) + tof(i-1);
end

for i = 1:(N+1) 
    % States at the DSMs in cylindrical coords
    cyl_DSM(:,i) = var((N+2 + 4*i-3):(N+2 + 4*i-1));
    % Corresponding times at DSMs
    times_DSM(i) = var(N+2 + 4*i);
end

%% Convert the times from years to seconds with TU
ephtimes = ephtimes * TU; 
% tof = tof * TU;
times_DSM = times_DSM * TU;

%% Convert the DSM coordinates from cylindrical to cartesian
car_DSM = zeros(size(cyl_DSM));

for i = 1:(N+1)
    car_DSM(:,i) = cyl2car(cyl_DSM(:,i));
end

%% Retrieve positions of the planets
r_planets = zeros(3,N+2);
v_planets = zeros(3,N+2);

for i = 1:(N+2)

    % Obtain the time in days for uplanet
    jd_i = ephtimes(i) / (24 * 3600);

    % Retrieve keplerian coordinates and convert them
    [kep1,~] = uplanet(jd_i, planets_id(i));
    [r1,v1] = kep2car(kep1, mu_S);
    r_planets(:,i) = r1;
    v_planets(:,i) = v1;

end

%% Compute Transfer Arcs

% Set up the matrices to save the velocities at the planets and at the DSMs 
v_FB = zeros(N+1,3,2); 
v_DSM = zeros(N+1,3,2); 

% Loop over the N+1 arcs 
for i = 1:(N+1) 
    
    % Compute ToFs in seconds to be given to Lambert
    tof_1 = times_DSM(i) - ephtimes(i); % can it be executed vectorially in one operation?
    tof_2 = ephtimes(i+1) - times_DSM(i);

    % Lambert arc from i-th planet to the position of the i-th DSM
    [~,~,~,~,v1_t,v2_t,~,~] = lambertMR( r_planets(:,i), car_DSM(:,i), tof_1, mu_S);
    v_FB(i,:,2) = v1_t;
    v_DSM(i,:,1) = v2_t;

    % Lambert arc from i-th DSM to (i+1)-th planet
    [~,~,~,~,v1_t,v2_t,~,~] = lambertMR( car_DSM(:,i), r_planets(:,i+1), tof_2, mu_S);
    v_DSM(i,:,2) = v1_t;
    v_FB(i+1,:,1) = v2_t;
end

% Compute DV_DSM
DV_DSM = 0;

% CAN BE MADE VECTORIAL WITH VECNORM AND SUM
for i = 1:(N+1)
    DV_DSM = DV_DSM + norm(v_DSM(i,:,2) - v_DSM(i,:,1));
end

%% Fly-bys

% Set up the vector to save the DVs of the Powered Fly-bys
DV_fb = zeros(N,1);

% Loop over the planets at which there is a FB (so from 2 to N+1)
for i = 2:(N+1)

    % Compute the infinite velocities before and after each FB
    v_inf_minus = v_FB(i,:,1)' - v_planets(:,i);
    v_inf_plus = v_FB(i,:,2)' - v_planets(:,i);

    % Retrieve gravitational constant of the FB planet
    mu_flyby_planet = astroConstants(20 + planets_id(i)); % THIS STEP CAN BE DONE OUT OF THE OBJ FUNCTION 

    % Use 110% of the planet radius as initial guess 
    rp_guess = 1.1 * astroConstants(10 + planets_id(i)); % THIS STEP CAN BE DONE OUT OF THE OBJ FUNCTION 

    % Compute the fly-by and the DV
    [DV_i, ~, ~] = evaluate_flyby_cost(v_inf_minus, v_inf_plus, mu_flyby_planet, rp_guess);
    DV_fb(i) = DV_i; % MORE AUTOMATIC CONTROL CHECKS ARE REQUIRED for consistency (with warning triggers)

end

DV_FB = sum(DV_fb);


%% Compute FINAL DV (NEED TO DECIDE WHETHER TO INCLUDE IT OR NOT IN OPTIMIZATION)

[kep_last, ~] = uplanet(ephtimes(end), planets_id(end));
[~,V_last] = kep2car(kep_last, mu_S);

mu_end = astroConstants(20 + planets_id(end)); % CAN BE MOVED OUT THE OBJ FUNCTION AND PASSED TO IT

Vinf_entry = v_planets(end,:,1) - V_last;

norm_Vinf_entry = norm(Vinf_entry);

[DV_capture, ~] = EstimateDVtoCapture(norm_Vinf_entry, mu_end, Ra_target, Rp_target);

%% COMPUTE OVERALL DV

DV = DV_DSM + DV_FB + DV_capture;


end



