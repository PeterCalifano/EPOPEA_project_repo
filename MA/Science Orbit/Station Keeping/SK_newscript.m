%% SK optimization
clear;clc;close all

% Load workspace for nominal time
load('FirstTrial_SK_1day.mat', 'XX_ii','x0_Sk2')

old_initial_guess = XX_ii;
old_pericenter_0 = x0_Sk2;

load('Workspace_30days_1.1avg.mat')

%load('Workspace_4orbit_optim.mat')

clear sum;
SK_points_old = SK_points;

% Backpropagate to find the initial pericenter
t1=SK_points_old(7,1);
t2=t1-1;
[~,~,t_0_init,new_peri0,~] = ode113(@SCR3BP_dyn,[t1 t2],SK_points_old(1:6,1),options_ode_event,mu_tbp,mu_v,R_v,J2_v);

options = optimoptions('fmincon', 'Algorithm', 'active-set', 'Display', 'iter',...
    'OptimalityTolerance', 1e-8, 'StepTolerance', 1e-8, 'ConstraintTolerance', 1e-8,...
    'SpecifyObjectiveGradient', false, 'SpecifyConstraintGradient', false, ...
    'MaxFunctionEvaluations',5000000,'MaxIterations',500000,'FunctionTolerance',1e-8); 


% Define the number of orbits per day 
N_orbits = 4;
n_var = 7*N_orbits*2;

% Define the number of days of propagation
N_days = 3; 

%%% Create bounds for SK position and velocity
    norm_r = sqrt((SK_points_old(1,:) - (1-mu_tbp)).^2 + (SK_points_old(2,:)).^2 +...
        (SK_points_old(3,:)).^2);
    norm_v = sqrt((SK_points_old(4,:)).^2 + (SK_points_old(5,:)).^2 +...
        (SK_points_old(6,:)).^2);
    
    r_max = max(norm_r)*1.2;
    r_min = min(norm_r)*0.8;
    v_max = max(norm_v)*1.2;
    v_min = min(norm_v)*0.8;

    lb_peri = (R_Enceladus+19)/DU;
    ub_peri = (R_Enceladus+60)/DU;
    lb_apo = (R_Enceladus+800)/DU;
    ub_apo = (R_Enceladus+1500)/DU;

%%% Create linear constraints and boundaries
    A = [];
    B = [];
    Aeq = [];
    Beq = [];
    ub = 1e+10*ones(n_var,1);
    lb = -1e+10*ones(n_var,1);
    t_ref = 0;
    for ii = 1:2*N_orbits
        
        if rem(ii,2) ~= 0
    
            % SK1 ( peri --> apo )
            ub(7*ii) = t_ref + ((ii-1)/2 * t_orb + tf_SK)*1.2;
            lb(7*ii) = t_ref + ((ii-1)/2 * t_orb + tf_CI)*0.8;
    
        else
    
            % SK2 ( apo --> peri )
            ub(7*ii) = t_ref + (ii/2 * t_orb - tf_CI)*1.2;
            lb(7*ii) = t_ref + (ii/2 * t_orb - tf_SK)*0.8;
            
        end
    
    end

%%% Create initial guess as:
    % 1) the one coming from the preliminary optimization
            % initial_guess_opt = old_initial_guess;
    
    % 2) the one coming from the 30 days run with decent results
    initial_guess_opt = zeros(n_var,1);
    for i = 1:N_orbits
    initial_guess_opt(14*i-13:14*i) = [SK_points_old(:,2*i-1);
                         SK_points_old(:,2*i)];
    end

% Check on boundaries

lb_check = initial_guess_opt - lb;
list = 1:length(lb_check);
ind_lb = list(lb_check<0)
ub_check = - initial_guess_opt + ub;
ind_ub = list(ub_check<0)



%%

%%% Define the initial pericenter state as:

    % 1) the one coming from the preliminary optimization
            %pericenter_0 = old_pericenter_0';

    % 2) the one corresponding to the initial guess given
            pericenter_0 = new_peri0;

%%% Initialize the times
    t_0=0;
    t_update = 0;

%%% Initialize storage variables
    DV_days = zeros(1,N_days);
    SK_points = zeros(7,N_days*2*N_orbits);

for ii = 1:N_days

   % Optimization
    [XX_ii, DV_ii] = fmincon(@(var) objfcn_multiple_SK(var,mu_tbp,mu_v,R_v,J2_v,pericenter_0,t_0,N_orbits),initial_guess_opt,A,B,...
        Aeq,Beq,lb,ub,@(var) nlcon_multiple_SK_old(var,mu_tbp,mu_v,R_v,J2_v,pericenter_0,t_0,N_orbits,lb_peri,ub_peri,lb_apo,ub_apo,...
        r_max,r_min,v_max,v_min), options);
    
    % Save outputs
    DV_days(ii) = DV_ii;
    for jj = 1 : 2*N_orbits
        SK_points(1:6,4*(ii-1) + jj) = XX_ii(7*jj - 6 : 7*jj-1);
        SK_points(7,4*(ii-1) + jj) = t_update + XX_ii(7*jj);
    end

    % Update initial pericenter position and initial guess
    
    t1 = SK_points(7,2*N_orbits*ii);
    t2 = t1 + 1;
    [~,prop_arc_0,t_0_init_new,new_pericenter_0,i_e] = ode113(@SCR3BP_dyn,[t1 t2],SK_points(1:6,2*N_orbits*ii),options_ode_event,mu_tbp,mu_v,R_v,J2_v);
    prop_arc_0=prop_arc_0';
    pericenter_0 = new_pericenter_0';
    initial_guess_opt = XX_ii;
    t_update = t_0_init_new;
    t_0 = t_0_init_new;
end
    
%% post processing - optimized
close all
% Initialize DV of each SK fire
DV_array = zeros(1,2*N_orbits*N_days);

% First propagation
x_sk1=SK_points(1:6,1);
t_sk1=SK_points(7,1);
[~,prop_state] = ode113(@SCR3BP_dyn,[0 t_sk1],x0_Sk2',options_ode,mu_tbp,mu_v,R_v,J2_v);
prop_state=prop_state';
DV_array(1) = norm(prop_state(4:6,end) - SK_points(4:6,1));

% Plot
Enceladus_3D(R_Enceladus,[(1-mu_tbp)*DU,0,0]);
for k = 1 : 2*N_orbits*N_days-1
%for k = 1 : 4*4

    % Propagation
    t1 = SK_points(7,k);
    t2 = SK_points(7,k+1);

    plot3(SK_points(1,k)*DU,SK_points(2,k)*DU,SK_points(3,k)*DU,'o','markersize',5,'linewidth',2);%,'DisplayName',['SK n ', num2str(k)]);
    
    [~,prop_arc] = ode113(@SCR3BP_dyn,[t1 t2],SK_points(1:6,k),options_ode,mu_tbp,mu_v,R_v,J2_v);
    prop_state=[prop_state,prop_arc'];
    DV_array(k+1)=norm(prop_state(4:6,end)-SK_points(4:6,k+1));

end

t1 = SK_points(7,end);
t2 = t1 + 1;
P1=plot3(SK_points(1,end)*DU,SK_points(2,end)*DU,SK_points(3,end)*DU,'ob','markersize',5,'linewidth',2);%,'DisplayName',['SK n ', num2str(4*N_days)]);   
[~,prop_arc_fin,t_e,x_e,i_e] = ode113(@SCR3BP_dyn,[t1 t2],SK_points(1:6,end),options_ode_event,mu_tbp,mu_v,R_v,J2_v);

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

%% Find apocenter and pericenter

norm_r = sqrt((prop_state(1,:) - (1-mu_tbp)).^2 + (prop_state(2,:)).^2 +...
    (prop_state(3,:)).^2);
apo_max = max(norm_r)*DU - R_Enceladus
peri_min = min(norm_r)*DU - R_Enceladus

peri0_old = sqrt((old_pericenter_0(1) - (1-mu_tbp)).^2 + (old_pericenter_0(2)).^2 +...
    (old_pericenter_0(3)).^2)*DU - R_Enceladus

