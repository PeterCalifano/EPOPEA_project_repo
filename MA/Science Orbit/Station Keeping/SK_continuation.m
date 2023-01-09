%% SK continuation
clear;clc;close all


load('DV0.27_6days_4orb.mat')

clear sum;
SK_points_old = SK_points;

options = optimoptions('fmincon', 'Algorithm', 'active-set', 'Display', 'iter',...
    'OptimalityTolerance', 1e-8, 'StepTolerance', 1e-8, 'ConstraintTolerance', 1e-8,...
    'SpecifyObjectiveGradient', false, 'SpecifyConstraintGradient', false, ...
    'MaxFunctionEvaluations',5000000,'MaxIterations',500000,'FunctionTolerance',1e-8); 

options_ode_stop = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-13,'Events',@FirstZeroCrossing);
options_ode_event = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-13,'Events',@ApseLineCrossing);
options_ode = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-13);

% Define the number of orbits per day 
N_orbits = 4;
n_var = 7*N_orbits*2;

% Define the number of days of propagation
N_days_prev = ii;
N_days_new = 12; 

check_lb = initial_guess_opt - lb;
check_ub =  - initial_guess_opt + ub;
list = 1:length(initial_guess_opt);

list_lb = list(check_lb<0)
list_ub = list(check_ub<0)

%%

pericenter_0 = new_pericenter_0;


%%% Initialize storage variables
DV_days = [DV_days, zeros(1,N_days_new)];
SK_points = [SK_points,zeros(7,N_days_new)];


for ii = N_days_prev + 1 : N_days_prev + N_days_new

   % Optimization
    [XX_ii, DV_ii] = fmincon(@(var) objfcn_multiple_SK(var,mu_tbp,mu_v,R_v,J2_v,pericenter_0,t_0,N_orbits),initial_guess_opt,A,B,...
        Aeq,Beq,lb,ub,@(var) nlcon_multiple_SK_old(var,mu_tbp,mu_v,R_v,J2_v,pericenter_0,t_0,N_orbits,lb_peri,ub_peri,lb_apo,ub_apo,...
        r_max,r_min,v_max,v_min), options);
    
    % Save outputs
    DV_days(ii) = DV_ii;
    for jj = 1 : 2*N_orbits
        SK_points(1:6,2*N_orbits*(ii-1) + jj) = XX_ii(7*jj - 6 : 7*jj-1);
        SK_points(7,2*N_orbits*(ii-1) + jj) = XX_ii(7*jj);
    end

    % Update initial pericenter position and initial guess
    
    t1 = SK_points(7,2*N_orbits*ii);
    t2 = t1 + 1;
    [~,prop_arc_0,t_0_new,new_pericenter_0,i_e] = ode113(@SCR3BP_dyn,[t1 t2],SK_points(1:6,2*N_orbits*ii),options_ode_event,mu_tbp,mu_v,R_v,J2_v);
    prop_arc_0=prop_arc_0';
    pericenter_0 = new_pericenter_0';

    initial_guess_opt = XX_ii;

    for kk = 1:2*N_orbits*7
        if rem(kk,7) == 0
            initial_guess_opt(kk) = initial_guess_opt(kk) - t_0 + t_0_new;
        end
    end

    % Update the bounds
    lb = lb - t_0 + t_0_new;
    ub = ub - t_0 + t_0_new;
    
    t_0 = t_0_new;

    a = 1;
    
end
   


%% post processing - optimized
close all
% Initialize DV of each SK fire
DV_array = zeros(1,2*N_orbits*(N_days_prev+N_days_new));
times_peri = zeros(1,2*N_orbits*(N_days_prev+N_days_new));
states_peri = zeros(6,2*N_orbits*(N_days_prev+N_days_new));
states_peri(:,1) = new_peri0';

% First propagation
x_sk1=SK_points(1:6,1);
t_sk1=SK_points(7,1);
[~,prop_state] = ode113(@SCR3BP_dyn,[0 t_sk1],new_peri0',options_ode,mu_tbp,mu_v,R_v,J2_v);
prop_state=prop_state';
DV_array(1) = norm(prop_state(4:6,end) - SK_points(4:6,1));

% Plot
Enceladus_3D(R_Enceladus,[(1-mu_tbp)*DU,0,0]);
for k = 1 : 2*N_orbits*(N_days_prev+N_days_new)-1
%for k = 1 : 16

    % Propagation
    t1 = SK_points(7,k);
    t2 = SK_points(7,k+1);

    plot3(SK_points(1,k)*DU,SK_points(2,k)*DU,SK_points(3,k)*DU,'o','markersize',5,'linewidth',2);
    
    [~,prop_arc,t_peri_ii,peri_ii,~] = ode113(@SCR3BP_dyn,[t1 t2],SK_points(1:6,k),options_ode_event,mu_tbp,mu_v,R_v,J2_v);
    prop_state=[prop_state,prop_arc'];
    times_peri(k+1) = t_peri_ii;
    states_peri(:,k+1) = peri_ii';
    DV_array(k+1)=norm(prop_state(4:6,end)-SK_points(4:6,k+1));
    
end

t1 = SK_points(7,end);
t2 = t1 + 1;
P1=plot3(SK_points(1,end)*DU,SK_points(2,end)*DU,SK_points(3,end)*DU,'ob','markersize',5,'linewidth',2);  
[~,prop_arc_fin,t_e,x_e,i_e] = ode113(@SCR3BP_dyn,[t1 t2],SK_points(1:6,end),options_ode_stop,mu_tbp,mu_v,R_v,J2_v);

prop_state=[prop_state,prop_arc_fin'];

DV_array_dim = DV_array*DU*1000/TU;
DV_days_dim = DV_days*DU*1000/TU

P2=plot3(prop_state(1,:)*DU,prop_state(2,:)*DU,prop_state(3,:)*DU,...
    'k--','linewidth',0.5);
P3=plot3(state_fullHalo(:,1),state_fullHalo(:,2),state_fullHalo(:,3),...
    'r--','linewidth',0.5,'DisplayName','Free dynamics Halo');
% P4=plot3(state_vec_Halo(1,:),state_vec_Halo(2,:),state_vec_Halo(3,:),...
%     'k','linewidth',1.25,'DisplayName','Southern Halo Orbit');
grid minor
%legend()
legend([P1,P2,P3],'SK points','trajectory','free dynamics crash','CR3BP')

%% Visualize the points
index = SK_points(7,:) * TU / (3600*24);
figure
scatter(index,DV_array_dim,20,'filled')
xlabel('Time Elapsed [days]')
ylabel('$\Delta v\;[m/s]$')
grid minor
hold on

%% Find apocenter and pericenter

norm_r = sqrt((prop_state(1,:) - (1-mu_tbp)).^2 + (prop_state(2,:)).^2 +...
    (prop_state(3,:)).^2);
apo_max = max(norm_r)*DU - R_Enceladus
peri_min = min(norm_r)*DU - R_Enceladus

peri0_old = sqrt((old_pericenter_0(1) - (1-mu_tbp)).^2 + (old_pericenter_0(2)).^2 +...
    (old_pericenter_0(3)).^2)*DU - R_Enceladus

