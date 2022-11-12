function DV = objfun_EarthSaturntransfer(var,N,planets_id,Ra_target,Rp_target,TU)
%%%% INPUTS %%%%%%
% N - Number of fly-bys
% planets - String of identifier of planets visited (dep and arrival included)
% var - variables vector, with t1 first and then the ToF. These times are
%       adimensionalized with TU (at the current stage TU is 1 year, the
%       variables are passed to the different functions accordingly)
% Ra_target - Apoapsis Radius for capture orbit (SEE EstimateDVtoCapture.m)
% Rp_target - Periapsis Radius for capture orbit (SEE EstimateDVtoCapture.m)
% TU - Time Unit used to adimensionalize (It's considered to be 1 year)
%
%%% EXPLANATION OF VAR VECTOR %%%%%
% var = [t_1,
%       ToF_1,...,ToF_(N+1),
%       r_1, th_1, z_1, t_DSM_1,
%       ...
%       r_(N+1), th_(N+1), z_(N+1), t_DSM_(N+1) ]
%--------------------------------------------------------------------------
%% RETRIEVE CONSTANTS:
mu_S = astroConstants(4);

%% EXTRACT VARIABLES FROM VAR
t1 = var(1);
tof = zeros(N+1,1);
times = zeros(N+2,1);
times(1) = t1;
for i = 2:N+2
    tof(i-1) = var(i);
    times(i) = times(i-1) + tof(i-1);
end

cyl_DSM = zeros(3,N+1);
times_DSM = zeros(N+1,1);
for i = 1:(N+1) 
    cyl_DSM(:,i) = var((N+2 + 4*i-3):(N+2 + 4*i-1));
    times_DSM(i) = var(N+2 + 4*i);
end

%% Convert the times from years to seconds with TU
times = times * TU; 
tof = tof * TU;
times_DSM = times_DSM * TU;

%% Convert the DSM coordinates from cylindrical to cartesian
car_DSM = zeros(size(cyl_DSM));

for i = 1:(N+1)
    car_DSM(:,i) = cyl2car(cyl_DSM(:,i));
end

%% Retrieve positions of the planets
r_planets = zeros(3,N+2);
for i = 1:(N+2)

    % Obtain the time in days for uplanet
    jd_i = times(i) / (24 * 3600);

    % Retrieve keplerian coordinates and convert them
    [kep1,~] = uplanet(jd_i, planets_id(i));
    [r1,~] = kep2car(kep1, mu_S);
    r_planets(:,i) = r1;
end

%% Arcs:

% Set up the matrices to save the velocities at the planets and at the DSMs 
v_planets = zeros(N+1,3,2);
v_DSM = zeros(N+1,3,2);

% Loop over the N+1 arcs
for i = 1:(N+1) 
    
    % Compute ToFs in seconds to be given to Lambert
    tof_1 = times_DSM(i) - times(i);
    tof_2 = times(i+1) - times_DSM(i);

    % Lambert arc from i-th planet to the position of the i-th DSM
    [~,~,~,~,v1_t,v2_t,~,~] = lambertMR( r_planets(i), car_DSM(:,i), tof_1, mu_S);
    v_planets(i,:,2) = v1_t;
    v_DSM(i,:,1) = v2_t;

    % Lambert arc from i-th DSM to (i+1)-th planet
    [~,~,~,~,v1_t,v2_t,~,~] = lambertMR( car_DSM(:,i), r_planets(i+1), tof_2, mu_S);
    v_DSM(i,:,1) = v1_t;
    v_planets(i+1,:,1) = v2_t;
end

%%% FROM HERE ON IS TO BE DONE PROPERLY %%%%%%
%% Fly-bys
DV_fb = zeros(N,1);
for i = 1:N
    [DV_i, rp_i, turning_angle_i] = evaluate_flyby_cost(v_inf_minus, v_inf_plus, mu_flyby_planet, rp_guess);
    %[DV_i,rp_i] = Gravity_Assist(planets(i+1),v_interp(i,:,2),v_interp(i+1,:,1),times(i+1));
    DV_fb(i) = DV_i;
end


DV = sum(DV_fb);


%% Compute FINAL DV (NEED TO DECIDE WHETHER TO INCLUDE IT OR NOT IN OPTIMIZATION)
X_last = cspice_spkezr(planets(end),times(end),'ECLIPJ2000','NONE','SUN');
[kep,ksun] = uplanet(times(end), planets_id(end));
V_last = X_last(4:6);
mu_main = cspice_bodvrd(planets(end),'GM',1);
Vinf_entry = v_planets(end,:,2) - V_last;
[DV_capture, T_capture] = EstimateDVtoCapture(Vinf_entry, mu_main, Ra_target, Rp_target);




end



