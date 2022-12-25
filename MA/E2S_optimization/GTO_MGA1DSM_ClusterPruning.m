close all;
clear;
clc;

% Load SPICE
cspice_kclear();
try
    cspice_furnsh('spice_kernels/pck00010.tpc')
    cspice_furnsh('spice_kernels/naif0012.tls')
    cspice_furnsh('spice_kernels/gm_de431.tpc')
    cspice_furnsh('spice_kernels/de440s.bsp')
catch
    cspice_furnsh('C:\Users\pietr\OneDrive - Politecnico di Milano\PoliMi - LM\Space Engineering - Master\CODE FILES\SGN\default_metakernel.tm');
end

%% Problem initialization
% Define algorithm parameters
N1 = 50; % n° Global optimization runs per iter (how many tentative solution are found by ga)
iter = 1; 
maxiter = 4; % Number of maximum allowed iteration of while loop
perc = [70*ones(1, floor(maxiter/2)), 90*ones(1, ceil(maxiter/2))]; % Percentile of objective function for LB, UB update
cost_thr = 0.8; % DV cost in km/s ?
stoptime = 3*60; % Stop time for ga solver
maxtime = 10*3600; % Max allowable execution time

rng shuffle

% Sequences: 1) E-VEJ-S, 2) E-VEE-S, 3)E-VEVE-S, 4) E-VEEJ-S , % 5) E-J-S,
% 6) E-EEJ-S
Sequence_selector = 6;
% Solver selection: 0) ga, 1) SA
heursolver_selector = 0;

fminconopts_hyb = optimoptions('fmincon', 'FunctionTolerance', 1e-14, 'StepTolerance', 1e-12, 'Display', 'iter',...
    'MaxFunctionEvaluations', 1e5, 'MaxIterations', 350);

% GA options
% opts_ga = optimoptions('ga', 'FunctionTolerance', 1e-20, 'MaxTime', stoptime,...
%     'UseParallel', true, 'PopulationSize', 200, 'Display', 'iter', 'MaxGenerations', 1e3,...
%     'CrossoverFraction', 0.7, 'MaxStallGenerations', 3, 'MaxStallTime', 0.5*stoptime,...
%     'HybridFcn', {@fmincon, fminconopts_hyb});

opts_ga = optimoptions('ga', 'FunctionTolerance', 1e-20, 'MaxTime', stoptime,...
    'UseParallel', true, 'PopulationSize', 500, 'Display', 'iter', 'MaxGenerations', 1e3,...
    'CrossoverFraction', 0.7, 'MaxStallGenerations', 3, 'MaxStallTime', 0.5*stoptime);

% Fminunc options
opts_fmincon = optimoptions('fmincon', 'Display', 'iter', 'FunctionTolerance', 1e-14,...
    'OptimalityTolerance', 1e-12, 'MaxFunctionEvaluations', 5e4, 'MaxIterations', 300, ...
    'StepTolerance', 1e-11);



% Number of flybys
N_fb = 4;
NLPvars = ones(4*N_fb + 6, 1)';

% 2042-2050

% Define sequence and constant parameters
switch Sequence_selector
    case 1 % E-VEJ-S,

        fb_sequence = [2, 3, 5];
        LBt_launchdate = cspice_str2et('2042-01-01 00:00:00.000 UTC')./(3600*24); % [days]
        UBt_launchdate = cspice_str2et('2050-01-01 00:00:00.000 UTC')./(3600*24); % [days]

        LBvinfdep = 3; % [km/s]
        UBvinfdep = 5; % [km/s]

        LBu = 0;
        UBu = 1;
        LBv = 0.4;
        UBv = 0.8;

        LBtof = [100, 120, 400, 1000];
        UBtof = [500, 500, 2000, 3000];
        LBeta = [0.2, 0.2, 0.2, 0.2];
        UBeta = [0.8, 0.8, 0.8, 0.8];


        LBRp_seq = [1.035, 1.035, 20];
        UBRp_seq = [6, 6, 160];
        LBbeta = [-pi, -pi, -pi];
        UBbeta = [pi, pi, pi];

        N_fb = length(fb_sequence);
        NLPvars = ones(4*N_fb + 6, 1)';

    case 2 % E-VEE-S

        fb_sequence = [2, 3, 3];
        LBt_launchdate = cspice_str2et('2030-01-01 00:00:00.000 UTC')./(3600*24); % [days]
        UBt_launchdate = cspice_str2et('2040-01-01 00:00:00.000 UTC')./(3600*24); % [days]

        LBvinfdep = 3; % [km/s]
        UBvinfdep = 5; % [km/s]

        LBu = 0;
        UBu = 1;
        LBv = 0.4;
        UBv = 0.8;

        LBtof = [60, 60, 300, 2000];
        UBtof = [400, 400, 400, 5000];
        LBeta = [0.2, 0.2, 0.2, 0.2];
        UBeta = [0.8, 0.8, 0.8, 0.8];

        LBRp_seq = [1.04, 1.035, 1.035];
        UBRp_seq = [5, 5, 5];
        LBbeta = [-pi, -pi, -pi];
        UBbeta = [pi, pi, pi];

        N_fb = length(fb_sequence);
        NLPvars = ones(4*N_fb + 6, 1)';

case 3 % E-VEVE-S

        fb_sequence = [2, 3, 2, 3];
        LBt_launchdate = cspice_str2et('2036-01-01 00:00:00.000 UTC')./(3600*24); % [days]
        UBt_launchdate = cspice_str2et('2046-01-01 00:00:00.000 UTC')./(3600*24); % [days]

        LBvinfdep = 3; % [km/s]
        UBvinfdep = 6.5; % [km/s]

        LBu = 0;
        UBu = 1;
        LBv = 0.4;
        UBv = 0.8;

        LBtof = [60, 60, 260, 200, 300];
        UBtof = [400, 400, 400, 600, 700];
        LBeta = [0.2, 0.2, 0.2, 0.2, 0.2];
        UBeta = [0.8, 0.8, 0.8, 0.8, 0.8];

        LBRp_seq = [1.02, 1.035, 1.02, 1.035];
        UBRp_seq = [8, 8, 8, 8];
        LBbeta = [-pi, -pi, -pi, -pi];
        UBbeta = [pi, pi, pi, pi];

        N_fb = length(fb_sequence);
        NLPvars = ones(4*N_fb + 6, 1)';

    case 4 % E-VEEJ-S,

        fb_sequence = [2, 3, 3, 5];


        LBt_launchdate = cspice_str2et('2030-01-01 00:00:00.000 UTC')./(3600*24); % [days]
        UBt_launchdate = cspice_str2et('2045-01-01 00:00:00.000 UTC')./(3600*24); % [days]

        LBvinfdep = 3; % [km/s]
        UBvinfdep = 6.05; % [km/s]

        LBu = 0;
        UBu = 1;
        LBv = 0.4;
        UBv = 0.8;

        LBtof = [60, 60, 300, 800, 1000];
        UBtof = [500, 500, 500, 1600, 2300];
        LBeta = [0.2, 0.2, 0.2, 0.2, 0.2];
        UBeta = [0.8, 0.8, 0.8, 0.8, 0.8];


        LBRp_seq = [1.035, 1.035, 1.035, 25];
        UBRp_seq = [5, 5, 8, 50];
        LBbeta = [-pi, -pi, -pi, -p];
        UBbeta = [pi, pi, pi, pi];

        N_fb = length(fb_sequence);
        NLPvars = ones(4*N_fb + 6, 1)';
        
    case 5 % E-J-S
        fb_sequence = [5];

        LBt_launchdate = cspice_str2et('2030-01-01 00:00:00.000 UTC')./(3600*24); % [days]
        UBt_launchdate = cspice_str2et('2044-01-01 00:00:00.000 UTC')./(3600*24); % [days]

        LBvinfdep = 3; % [km/s]
        UBvinfdep = 6.5; % [km/s]

        LBu = 0;
        UBu = 1;
        LBv = 0;
        UBv = 0.4;

        LBtof = [365.5*2, 365.5*3];
        UBtof = [365.5*5, 365.5*10];

        LBeta = [0.1, 0.1];
        UBeta = [0.9, 0.9];

        LBRp_seq = [5];
        UBRp_seq = [50];
        LBbeta = [-pi];
        UBbeta = [pi];

        N_fb = length(fb_sequence);
        NLPvars = ones(4*N_fb + 6, 1)';

    case 6 % E-EEJ-S

        fb_sequence = [3, 3, 5];

        LBt_launchdate = cspice_str2et('2037-01-01 00:00:00.000 UTC')./(3600*24); % [days]
        UBt_launchdate = cspice_str2et('2044-01-01 00:00:00.000 UTC')./(3600*24); % [days]

        LBvinfdep = 5.5; % [km/s]
        UBvinfdep = 7.5; % [km/s]

        LBu = 0;
        UBu = 0.6;
        LBv = 0;
        UBv = 1;

        LBtof = [100, 300, 600, 1000];
        UBtof = [500, 800, 1500, 3000];

        LBeta = [0.05, 0.05, 0.05, 0.05];
        UBeta = [0.8, 0.8, 0.9, 0.9];

%         LBRp_seq = [1.8, 1.1, 30];
%         UBRp_seq = [2.05, 1.2, 32];
%         LBbeta = [-1.25, -1.3, -1.92];
%         UBbeta = [-1, -1.1, -1.82];
        LBRp_seq = [1.1, 1.1, 15];
        UBRp_seq = [9, 9, 70];
%         LBbeta = [-1.25, -1.3, -1.92];
%         UBbeta = [-1, -1.1, -1.82];

        LBbeta = [-pi, -pi,-pi];
        UBbeta = [pi,  pi, pi];

        N_fb = length(fb_sequence);
        NLPvars = ones(4*N_fb + 6, 1)';

    case 7 % E-EEEJ-S

        fb_sequence = [3, 3, 3, 5];

        LBt_launchdate = cspice_str2et('2035-01-01 00:00:00.000 UTC')./(3600*24); % [days]
        UBt_launchdate = cspice_str2et('2043-01-01 00:00:00.000 UTC')./(3600*24); % [days]

        LBvinfdep = 5.5; % [km/s]
        UBvinfdep = 7.2; % [km/s]

        LBu = 0;
        UBu = 0.5;
        LBv = 0;
        UBv = 1;

        LBtof = [100, 100, 100, 300, 1000];
        UBtof = [500, 500, 600, 1300, 3000];

        LBeta = [0.05, 0.05, 0.05, 0.05];
        UBeta = [0.9, 0.9, 0.9, 0.9];

        %         LBRp_seq = [1.8, 1.1, 30];
        %         UBRp_seq = [2.05, 1.2, 32];
        %         LBbeta = [-1.25, -1.3, -1.92];
        %         UBbeta = [-1, -1.1, -1.82];
        LBRp_seq = [1.1, 1.1, 1.1, 10];
        UBRp_seq = [10, 10, 10,  70];
        %         LBbeta = [-1.25, -1.3, -1.92];
        %         UBbeta = [-1, -1.1, -1.82];

        LBbeta = [-pi, -pi, -pi, -pi];
        UBbeta = [pi, pi, pi, pi];

        N_fb = length(fb_sequence);
        NLPvars = ones(4*N_fb + 6, 1)';

end

planet_seq = [3, fb_sequence, 6];
Ra_target = 200*astroConstants(26);
Rp_target = 2.55*astroConstants(26);

% Static allocation
LB = zeros(maxiter, length(NLPvars));
UB = zeros(maxiter, length(NLPvars));

NLPoptset_ga = zeros(N1, length(NLPvars), maxiter);
feval_ga = zeros(N1, 1, maxiter);
exitflag_ga = zeros(N1, 1, maxiter);

NLPoptset_simulanneal = zeros(N1, length(NLPvars), maxiter);
feval_simulanneal = zeros(N1, 1, maxiter);
exitflag_simulanneal = zeros(N1, 1, maxiter);

NLPoptset_local = zeros(N1, length(NLPvars), maxiter);
feval_local = zeros(N1, 1, maxiter);
exitflag_local = zeros(N1, 1, maxiter);


% Assign INITIAL [LB, UB] for all decision variable

LB(1, :) = [LBt_launchdate, LBvinfdep, LBu, LBv, LBtof, LBeta, LBRp_seq, LBbeta];
UB(1, :) = [UBt_launchdate, UBvinfdep, UBu, UBv, UBtof, UBeta, UBRp_seq, UBbeta];

if Sequence_selector == 1
% LB(1, :) = [18247.3952550500	3.26130504680723	0.0324887957649668	0.189517384867105	193.843108455067	195.149083580162...
%     1076.04934255018	2000	0.273248141971223	0.200000000000000	0.268231439504458	0.200000000000000	1.03500000000000	1.03500000000000	20	-2.96849811064606	-2.29537257370371	-3.14159265358979];
% 
% LB(1, :) = [19115.8626158843	3.26130504680723	0.0324887957649668	0.271301686637625	362.679289372678	321.136406208845	1284.15975992840	2336.71867169394	0.273248141971223	0.200000000000000	0.365356692727062	0.200000000000000	1.97193075293984	1.03500000000000	20	-2.31365697597635	-1.97940059660343	-2.07801819011155];
% 
% % LB(1, :) = [19309.4503293261	4.73919575702108	0.166301572030575	0.372404898833974	420.445385260946	326.703286978812	1719.97665009843	2336.71867169394	0.293620249145098	0.200000000000000	0.524295065983319	0.200000000000000	1.97193075293984	1.03500000000000	21.7499035228576	-1.84493640620792	-1.95474812899966	-1.52469084880116];
% LB(1, :) = [19124.9753673395	3.33079382834192	0.205169962256986	0.602118178873050	398.423415490033	321.136406208845	1604.26111724497	2336.71867169394	0.273248141971223	0.200000000000000	0.484017745471751	0.200000000000000	1.97193075293984	1.03500000000000	20.0249308754403	-2.31365697597635	-1.97940059660343	-1.67528011987852];
% 
% % LB(1, :) = [19154.6173331302	4.23609248223818	0.376397451709376	0.676346247775956	402.208524129756	360.872404278232	1944.13613013672	3136.21541598248	0.329579565153566	0.485867421741454	0.575855946619520	0.239162173203062	2.70372086376352	1.14631436071520	32.0653459929266	-1.64515278142009	-1.84809939262272	-1.34869281519231];
% 
% 
% UB(1, :) = [20088.5008007396	5	0.799013207324930	1	500	500	2000	3983.95312149453	0.660292846622087	0.595740136725434	0.745463688009700	0.479757256374570	5.74061410103485	3.45401057453000...
%     87.0058491653975	0.499998723341889	-0.537386178820938	-0.294587772674394];
% 
% UB(1, :) = [20088.5008007396	5	0.776731281897155	0.967113179061435	486.053509121667	500	2000	3940.28528958966	0.660292846622087	0.595740136725434	0.646857580574019	0.325493377968005	5.74061410103485	3.04574697897508	79.8681087476366	0.422314250478805	-0.978503910504675	-1.20258857081786];
% 
% % UB(1, :) = [19481.9586081478	5	0.733665862975332	0.824548471909893	469.111589320308	500	2000	3917.17612884743	0.658438147296497	0.595740136725434	0.621382648236797	0.325493377968005	3.82595662874976	1.45936367922071	50.7642710398457	-0.695990821093950	-1.44600683903938	-1.27626914095766];
% UB(1, :) = [19656.6504788324	5	0.478383371149862	0.676346247775956	480.747508435950	500	2000	3616.43886976651	0.660292846622087	0.523125543124660	0.646857580574019	0.300378164600519	4.51144112066214	1.86297029822060	53.0481322691609	-0.553725961591942	-1.23480886956356	-1.21451729564210];

end




% Simulated Annealing options
InitialTempVec = UB(1, :) - LB(1, :);

opts_simulanneal = optimoptions('simulannealbnd', 'Display', 'iter',...
    'FunctionTolerance', 1e-14, 'MaxTime', stoptime,...
    'MaxFunctionEvaluations', 1e6, 'MaxStallIterations', 500, ...
    InitialTemperature=InitialTempVec, ReannealInterval=300, TemperatureFcn='temperatureexp');

% Define convergence criterion
converge_flag = 0;

time_limit = tic();

% codegen -config:mex -lang::c++ -report ...
%        -args {N1, iter, NLPvars, planet_seq, Ra_target, Rp_target, LB, UB, params, NLPoptset_ga, feval_ga, exitflag_ga} ...
%        -nargout {1} HeuristicOptimization

% codegen -config:mex objfun_EarthSaturntransfer_X.m -args {NLPvars, planet_seq, Ra_target, Rp_target} -nargout 1 -lang::c -report

% cfg = coder.config('mex');
% cfg.IntegrityChecks = false;
% cfg.SaturateOnIntegerOverflow = false;
% cfg.DynamicMemoryAllocation = 'Off';

% codegen -config:mex LocalOptimization.m -args {NLPoptset_ga(:, :, 1), planet_seq, Ra_target, Rp_target, LB(iter, :), UB(iter, :), N1}...
%     -nargout 3 -O enable:openmp -lang::c -report


while converge_flag ~= 1

    %% Heuristic solver step
    if heursolver_selector == 0

        for idsol = 1:N1
            % Heuristic optimization: produce f_ga costs [N1x1] and corresponding
            % decision variables NLPoptset_ga(idsol, :, iter)
            [NLPoptset_ga(idsol, :, iter), feval_ga(idsol, 1, iter), exitflag_ga(idsol, 1, iter)] =...
                ga(@(NLPvars) objfun_EarthSaturntransfer_2(NLPvars, planet_seq, Ra_target, Rp_target),...
                length(NLPvars), [], [], [], [], LB(iter, :), UB(iter, :), [], opts_ga);

            disp("GA solver iter: " + num2str(idsol) + " of iter cycle " + num2str(iter));
        end

    elseif heursolver_selector == 1

        parfor idsol = 1:N1
            % Heuristic optimization: produce f_ga costs [N1x1] and corresponding
            % decision variables NLPoptset_ga(idsol, :, iter)
            [NLPoptset_simulanneal(idsol, :, iter), feval_simulanneal(idsol, 1, iter), exitflag_simulanneal(idsol, 1, iter)] =...
                simulannealbnd(@(NLPvars) objfun_EarthSaturntransfer_2(NLPvars, planet_seq, Ra_target, Rp_target),...
                NLPvars, LB(iter, :), UB(iter, :), opts_simulanneal);
            
        end
    end

    %% Local solver step
%     [NLPoptset_local(:, :, iter), feval_local(:, 1, iter), exitflag_local(:, 1, iter)] =...
%         LocalOptimization_mex(NLPoptset_ga(:, :, iter), ...
%         planet_seq, Ra_target, Rp_target, LB(iter, :), UB(iter, :), N1);

        parfor idsol = 1:N1
            % Local optimization (refinement): produce f_local costs [N1x1] and
            % corresponding decision variables NLPoptset_local(idsol, :, iter)
            [NLPoptset_local(idsol, :, iter), feval_local(idsol, 1, iter), exitflag_local(idsol, 1, iter)] = ...
                fmincon(@(NLPvars) objfun_EarthSaturntransfer_2(NLPvars, planet_seq, Ra_target, Rp_target),...
                NLPoptset_ga(idsol, :, iter), [], [], [], [], LB(iter, :), UB(iter, :), [], opts_fmincon);
        end


    % Evaluate convergence criterion
    if min(feval_local(:, 1, iter)) <= cost_thr
        converge_flag = 1;
    end


    %% First Pruning tentative

    % Find percentile of cost
    perc_value = prctile(feval_local(:, :, iter), perc(iter));
    feval_temp = feval_local(:, :, iter);
    % Create mask to extract only solutions below desired percentile
    maskid_perc = feval_temp < perc_value;
    % Extract decision vectors corresponding to these solutions
    NLPopts_perc = NLPoptset_local(maskid_perc == 1, :, iter);

    [min_feval, minpos] = min(feval_local(:, :, iter));
    [max_feval, maxpos] = max(feval_local(:, :, iter));

    % Find minimum of decision variables perc percentile
    LB_min = min(NLPopts_perc, [], 1);
    % Find maximum of decision variables perc percentile
    UB_max = max(NLPopts_perc, [], 1);

    % Redefine LB for iter+1
    %        LB_min = min(NLPoptset_local(:, :, iter), [], 1); % Find minimum of decision variables
    LB_new = LB_min - (UB(iter, :) - LB(iter, :))*perc(iter)/100/2;

    % Redefine UB for iter+1
    %        UB_max = max(NLPoptset_local(:, :, iter), [], 1); % Find maximum of decision variables
    UB_new = UB_max + (UB(iter, :) - LB(iter, :))*perc(iter)/100/2;

    for idVar = 1:length(NLPvars)

        if LB_new(idVar) > LB(iter, idVar)
            % If new LB is within previous LB
            LB(iter+1, idVar) = LB_new(idVar);
        else
            % Assign previous LB
            LB(iter+1, idVar) = LB(iter, idVar);
        end

        if UB_new(idVar) < UB(iter, idVar)
            % If new UB is within previous UB
            UB(iter+1, idVar) = UB_new(idVar);
        else
            % Assign previous UB
            UB(iter+1, idVar) = UB(iter, idVar);
        end
    end

    %% "Emergency" pruning
    % If 1 in mask --> no change both in LB and in UB for the ith variable
    if sum(or(LB(iter+1, :) == LB(iter, :), UB(iter+1, :) == UB(iter, :))) > floor((4*N_fb+6)./4)
        warning('Emergency pruning triggered')
        % Redefine LB for iter+1
        %        LB_min = min(NLPoptset_local(:, :, iter), [], 1); % Find minimum of decision variables
        LB_min = NLPoptset_local(minpos, :, iter);
        LB_new = LB_min - (UB(iter, :) - LB(iter, :))*perc(iter)/100/4;

        % Redefine UB for iter+1
        %        UB_max = max(NLPoptset_local(:, :, iter), [], 1); % Find maximum of decision variables
        UB_max = NLPoptset_local(minpos, :, iter);
        UB_new = UB_max + (UB(iter, :) - LB(iter, :))*perc(iter)/100/4;

        for idVar = 1:length(NLPvars)

            if LB_new(idVar) > LB(iter, idVar)
                % If new LB is within previous LB
                LB(iter+1, idVar) = LB_new(idVar);
            else
                % Assign previous LB
                LB(iter+1, idVar) = LB(iter, idVar);
            end

            if UB_new(idVar) < UB(iter, idVar)
                % If new UB is within previous UB
                UB(iter+1, idVar) = UB_new(idVar);
            else
                % Assign previous UB
                UB(iter+1, idVar) = UB(iter, idVar);
            end
        end
    end


   
    % Check if allowed maxiter reached
    elapsedtime = toc(time_limit);

    [min_at_iter(iter), min_pos(iter)] = min(feval_local(:, 1, iter), [], 1);
    fprintf('\nBest solution found at iter = %2g in position %4g \n', min_at_iter(iter), min_pos(iter));

    if iter == maxiter
        break;
    end

    % Update iteration counter
    iter = iter + 1;

    % Check if maxtime reached
    if  elapsedtime >= maxtime
        break;
    end
    pause(0.01);
end


% Best solution
for iter_id = 1:iter
    [min_at_iter(iter_id), min_pos(iter_id)] = min(feval_local(:, 1, iter_id), [], 1);
    fprintf('\nBest solution found at iter = %2g in position %4g \n', min_at_iter(iter_id), min_pos(iter_id));
end

if converge_flag == 1
    save('WorkspaceWithConvergedSolution.mat')
else
    save('Workspace13hrsNoConvergence.mat')
end
