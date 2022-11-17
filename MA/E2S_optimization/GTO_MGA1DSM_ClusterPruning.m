close all;
clear;
clc;

% Load SPICE
cspice_kclear();
cspice_furnsh('..\..\spice_kernels/pck00010.tpc')
cspice_furnsh('..\..\spice_kernels/naif0012.tls')
cspice_furnsh('..\..\spice_kernels/gm_de431.tpc')
cspice_furnsh('..\..\spice_kernels/de440s.bsp')


%% Problem initialization
% Define algorithm parameters
N1 = 100; % n° Global optimization runs per iter (how many tentative solution are found by ga)
iter = 1; 
maxiter = 20; % Number of maximum allowed iteration of while loop
perc = 80; % Percentile of objective function for LB, UB update
cost_thr = 1; % DV cost in km/s ?
stoptime = 600; % Stop time for ga solver
maxtime =  13*3600; % Max allowable execution time

rng shuffle

% GA options
opts_ga = optimoptions('ga', 'FunctionTolerance', 1e-10, 'MaxTime', stoptime,...
    'UseParallel', true, 'PopulationSize', 400, 'Display', 'iter', 'MaxGenerations', 1e3,...
    'CrossoverFraction', 0.8, 'MaxStallGenerations', 200);
% Fminunc options
opts_fmincon = optimoptions('fmincon', 'Display', 'iter', 'FunctionTolerance', 1e-12,...
               'OptimalityTolerance', 1e-9, 'MaxFunctionEvaluations', 2e5, 'MaxIterations', 300, ...
               'StepTolerance', 1e-12);

% Number of flybys
N_fb = 3;
NLPvars = ones(4*N_fb + 6, 1);

% Define sequence and constant parameters
fb_sequence = [2, 3, 5];
planet_seq = [3, fb_sequence, 6];
Ra_target = 200*astroConstants(26);
Rp_target = 3*astroConstants(26);

% Static allocation
LB = zeros(maxiter, length(NLPvars));
UB = zeros(maxiter, length(NLPvars));
NLPoptset_ga = zeros(N1, length(NLPvars), maxiter);
feval_ga = zeros(N1, 1, maxiter);
exitflag_ga = zeros(N1, 1, maxiter);

NLPoptset_local = zeros(N1, length(NLPvars), maxiter);
feval_local = zeros(N1, 1, maxiter);
exitflag_local = zeros(N1, 1, maxiter);

% Define INITIAL [LB, UB] for all decision variables
% LB_benchmark = [];
% UB_benchmark = [];

LBt_launchdate = cspice_str2et('2030-01-01 00:00:00.000 UTC')./(3600*24); % [days]
UBt_launchdate = cspice_str2et('2050-01-01 00:00:00.000 UTC')./(3600*24); % [days]

LBvinfdep = 3; % [km/s]
UBvinfdep = 5; % [km/s]

LBu = 0;
UBu = 1;
LBv = 0;
UBv = 1;

LBtof = [30, 30, 400, 1000];
UBtof = [400, 400, 2000, 6000]; 
LBeta = [0.2, 0.2, 0.2, 0.2];
UBeta = [0.8, 0.8, 0.8, 0.8];

LBRp_seq = [1.05, 1.05, 80];
UBRp_seq = [10, 10, 160];
LBbeta = [-pi, -pi, -pi];
UBbeta = [pi, pi, pi];

LB(1, :) = [LBt_launchdate, LBvinfdep, LBu, LBv, LBtof, LBeta, LBRp_seq, LBbeta];
UB(1, :) = [UBt_launchdate, UBvinfdep, UBu, UBv, UBtof, UBeta, UBRp_seq, UBbeta];

% Define convergence criterion
converge_flag = 0;

time_limit = tic();

while converge_flag ~= 1 

    % DEVNOTE: CONVERT HEUR AND LOCAL IN MEX

    for idsol = 1:N1
        % Heuristic optimization: produce f_ga costs [N1x1] and corresponding
        % decision variables NLPoptset_ga(idsol, :, iter)
        [NLPoptset_ga(idsol, :, iter), feval_ga(idsol, 1, iter), exitflag_ga(idsol, 1, iter)] =...
            ga(@(NLPvars) objfun_EarthSaturntransfer_2(NLPvars, planet_seq, Ra_target, Rp_target),...
            length(NLPvars), [], [], [], [], LB(iter, :), UB(iter, :), [], opts_ga);
    end

    parfor idsol = 1:N1
        % Local optimization (refinement): produce f_local costs [N1x1] and 
        % corresponding decision variables NLPoptset_local(idsol, :, iter)
        [NLPoptset_local(idsol, :, iter), feval_local(idsol, 1, iter), exitflag_local(idsol, 1, iter)] = ...
            fmincon(@(NLPvars) objfun_EarthSaturntransfer_2(NLPvars, planet_seq, Ra_target, Rp_target),...
            NLPoptset_ga(idsol, :, iter), [], [], [], [], LB(iter, :), UB(iter, :), [], opts_fmincon);
    end

    for idsol = 1:N1
        % Evaluate convergence criterion
        if min(feval_local(idsol, 1, iter)) <=  cost_thr
            converge_flag = 1;
        end

        % Redefine LB for iter+1
        LB_min = min(NLPoptset_local(:, :, iter), [], 1); % Find minimum of decision variables

        LB_new = LB_min - (UB(iter, :) - LB(iter, :))*2*perc/100;

        % Redefine UB for iter+1
        UB_max = max(NLPoptset_local(:, :, iter), [], 1); % Find maximum of decision variables
        UB_new = UB_max + (UB(iter, :) - LB(iter, :))*2*perc/100;

        for idVar = 1:length(NLPvars)

            if LB_new(idVar) > LB(iter, idVar)
                % If new LB is within previous LB
                LB(iter+1, :) = LB_min - (UB(iter, :) - LB(iter, :))*2*perc/100;
            else
                % Assign previous LB
                LB(iter+1, idVar) = LB(iter, idVar);
            end


            if UB_new(idVar) < UB(iter, idVar)
                % If new UB is within previous UB
                UB(iter+1, idVar) = UB_max + (UB(iter, :) - LB(iter, :))*2*perc/100;
                else
                    % Assign previous UB
                    UB(iter+1, idVar) = UB(iter, idVar);
            end
        end
    end

    % Check if allowed maxiter reached
    elapsedtime = toc(time_limit);

    if iter == maxiter 
        break;
    end

    % Update iteration counter
    iter = iter + 1;
   
    % Check if maxtime reached
    if  elapsedtime >= maxtime
        break;
    end

end


% Best solution
for iter_id = 1:iter

    [min_at_iter(iter_id), min_pos(iter_id)] = min(feval_local(:, 1, iter_id), [], 1);
    fprintf('\nBest solution found at iter = %2g in position %4g', min_at_iter(iter_id), min_pos(iter_id));
end

if converge_flag == 1
    save('WorkspaceWithConvergedSolution.mat')
else
    save('Workspace13hrsNoConvergence.mat')
end
