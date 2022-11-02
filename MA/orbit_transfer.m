function [deltaV_arrival, deltaV_injection, deltaV_tot, transfer_arc] = orbit_transfer(mu, r_i1, r_f2, v_i1, v_f2, ToF, transfer_method, transfer_arc_switch)
%% PROTOTYPE
% [deltaV_arrival, deltaV_injection, deltaV_tot, transfer_arc] = orbit_transfer(mu, r_i1, r_f2, v_i1, v_f2, ToF, transfer_method)
% -------------------------------------------------------------------------
%% DESCRIPTION
% Computes the particular transfer arc between two completely specified
% orbits and the deltaV of the two manoeuvres 
% -------------------------------------------------------------------------
%% INPUT
% mu: [scalar] planetary constant of the main body [Km^3/s^2]
% r_i1: position vector at time t1 [Km]
% r_f2: position vector at time t2 = t1 + dt [Km]
% v_i1: [3x1] velocity vector of orbit 1 at time t1 [Km/s]
% v_f2: [3x1] velocity vector of orbit 2 at time t2 [Km/s]
% ToF: [scalar] time of flight (t2-t1) [s]
% transfer_method: [scalar] define the direction of motion on the orbit. Prograde = 0, retrograde = 1
% transfer_arc_switch : 0 or 1 to switch off or on the transfer arc propagation
% -------------------------------------------------------------------------
%% OUTPUT
% deltaV_injection: [3x1] deltaV required to begin the transfer; 1 --> t
% deltaV_arrival: [3x1] deltaV required to match final velocity on orbit two; t --> 2 
% deltaV_tot: [scalar] total deltaV required for the transfer manoeuvre
% transfer_arc: [struct] fields --> [a, e, r, v, tof_min]
% -------------------------------------------------------------------------
%% CONTRIBUTORS
% Pietro Califano
% -------------------------------------------------------------------------
%% CHANGELOG
% V1: transfer arc characterization, deltaV computation 16/11/2021
% -------------------------------------------------------------------------
%% Next upgrades


%% Function code
r_i1 = [r_i1(1); r_i1(2); r_i1(3)];
r_f2 = [r_f2(1); r_f2(2); r_f2(3)];

%% Lambert solver
rev_num = 0; % number of revolutions is 0

[a_t, p_t, e_t, error_code, v_t1, v_t2, ToF_min, delta_theta] = lambertMR(r_i1, r_f2, ToF, mu, transfer_method , rev_num, 0, 2);

%% Transfer arc 
% Convert to column vector
v_t1 = v_t1'; 
v_t2 = v_t2';

% Propagation to get [r,v] of the transfer arc
transfer_arc = struct();

if transfer_arc_switch == 1
    if ToF > ToF_min
        [r_t, v_t, time_vec] = twobp_solver(ToF, [r_i1; v_t1], mu, 0, 10);

        transfer_arc.r = r_t;
        transfer_arc.v = v_t;
        transfer_arc.transf_times = time_vec;

    else
        transfer_arc.r = nan(1,1);
        transfer_arc.v = nan(1,1);
    end
end

transfer_arc.a = a_t;
transfer_arc.e = e_t;
transfer_arc.tof_min = ToF_min;

%% DeltaV computation
% Compute the two deltaV as vectors
deltaV_injection = v_t1 - v_i1;
deltaV_arrival = v_f2 - v_t2;

deltaV_tot = norm(deltaV_arrival) + norm(deltaV_injection);

end