function [c,ceq] = nlcon_multiple_SK(var,mu_tbp,mu_v,R_v,J2_v,x_0,N_orbits,lb_peri,ub_peri,lb_apo,ub_apo,...
    r_max,r_min,v_max,v_min)
% INPUTS
% var - [14*N_orbitsx1] design variables
% lb_peri - [1x1] lower bound for the pericentre
% ub_peri - [1x1] upper bound for the pericentre
% lb_apo - [1x1] lower bound for the apocentre
% ub_apo - [1x1] upper bound for the apocentre
% 
% OUTPUTS
% c - non linear inequality constraints
% ceq non linear equality constraints

options_ode=odeset('RelTol',1e-13,'AbsTol',1e-13);
options_ode_cycle=odeset( 'RelTol', 1e-13, 'AbsTol', 1e-13,'Events',@ApseLineCrossing);
options_ode_stop=odeset( 'RelTol', 1e-13, 'AbsTol', 1e-13,'Events',@LastCrossing);



c=zeros(4*N_orbits,1);
ceq=zeros(6*N_orbits,1);

%algorithm
%1st propagation
x_sk1=var(1:6);
t_sk1=var(7);
[~,prop_state] = ode113(@SCR3BP_dyn,[0 t_sk1],x_0,options_ode,mu_tbp,mu_v,R_v,J2_v);

flow1=prop_state(end,:)';
ceq(1:3)=flow1(1:3)-x_sk1(1:3);


for k=1:2*N_orbits-1
    %inequality constraints
    %propagation
    t1=var(7*k);
    t2=var(7*(k+1));
    [~,prop_state] = ode113(@SCR3BP_dyn,[t1 t2],var(7*k-6:7*k-1),options_ode,mu_tbp,mu_v,R_v,J2_v);
%     if exist("state_aps") == 0
%         plot3(prop_state(:,1),prop_state(:,2),prop_state(:,3),'linewidth',2)
%         axis equal
%         hold on
%         xlabel('X')
%         ylabel('Y')
%         zlabel('Z')
%         grid minor
%     end
    if rem(k,2)~=0
        max_dist = 0;
        for i = 1:length(prop_state(1,:))
            dist = norm(prop_state(i,1:3));
            if dist > max_dist
                max_dist = dist;
                index_apse = i;
            end
        end
    else
        min_dist = 1e+10;
        for i = 1:length(prop_state(1,:))
            dist = norm(prop_state(i,1:3));
            if dist < min_dist
                min_dist = dist;
                index_apse = i;
            end
        end

    end
    state_aps = prop_state(index_apse,:);
    state_aps(1)=state_aps(1)-(1-mu_tbp);
    r_aps=norm(state_aps(1:3));
    
    %position matching (equality constraints)
    flow=prop_state(end,:)';
    ceq(3*(k+1)-2:3*(k+1)) = flow(1:3)-var(7*(k+1)-6:7*(k+1)-4);

    %inequality constraints - apse line bounds
    if rem(k,2)~=0
        c(6*k-5)=-r_aps+lb_apo;
        c(6*k-4)=r_aps-ub_apo;
    else
        c(6*k-5)=-r_aps+lb_peri;
        c(6*k-4)=r_aps-ub_peri;
    end

    % Limits on position and velocity of SK points
    r_norm = norm(var(7*k-6:7*k-4) - [1-mu_tbp;0;0]);
    v_norm = norm(var(7*k-3:7*k-1));

    c(6*k - 3) = r_norm - r_max;
    c(6*k - 2) = - r_norm + r_min;
    c(6*k - 1) = v_norm - v_max;
    c(6*k) = - v_norm + v_min;
end
%final propagation, pericentre check 
t1=var(end);
t2=var(end)+1;
[~,prop_state] = ode113(@SCR3BP_dyn,[t1 t2],var(end-6:end-1),options_ode,mu_tbp,mu_v,R_v,J2_v);
min_dist = 1e+10;
for i = 1:length(prop_state(1,:))
    dist = norm(prop_state(i,1:3));
    if dist < min_dist
        min_dist = dist;
        index_apse = i;
    end
end
state_aps = prop_state(index_apse,:);
state_aps(1)=state_aps(1)-(1-mu_tbp);
r_aps=norm(state_aps(1:3));

%cosntraint
c(6*k+1)=-r_aps+lb_peri;
c(6*k+2)=r_aps-ub_peri;

r_norm = norm(var(7*(k+1)-6:7*(k+1)-4) - [1-mu_tbp;0;0]);
v_norm = norm(var(7*(k+1)-3:7*(k+1)-1));
c(6*k+3) = r_norm - r_max;
c(6*k+4) = - r_norm + r_min;
c(6*k+5) = v_norm - v_max;
c(6*k+6) = - v_norm + v_min;


end