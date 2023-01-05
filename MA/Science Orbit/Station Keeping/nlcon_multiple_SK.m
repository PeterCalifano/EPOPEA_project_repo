function [c,c_eq] = nlcon_multiple_SK(var,mu_tbp,mu_v,R_v,J2_v,x_0,N_orbits,lb_peri,ub_peri,lb_apo,ub_apo)
% INPUTS
% var - [7N_orbitsx1] design variables
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
    [~,prop_state,~,state_aps,~] = ode113(@SCR3BP_dyn,[t1 t2],var(7*k-6:7*k-1),options_ode_cycle,mu_tbp,mu_v,R_v,J2_v);
    state_aps(1)=state_aps(1)-(1-mu_tbp);
    r_aps=norm(state_aps(1:3));
    
    %position matching (equality constraints)
    flow=prop_state(end,:)';
    ceq(3*(k+1)-2:3*(k+1))=flow(1:3)-var(7*(k+1)-6:7*(k+1)-4);
    %inequality constraints - apse line bounds
    if rem(k,2)~=0
        c(2*k-1)=-r_aps+lb_apo;
        c(2*k)=r_aps-ub_apo;
    else
        c(2*k-1)=-r_aps+lb_peri;
        c(2*k)=r_aps-ub_peri;
    end
end
%final propagation, pericentre check 
t1=var(end);
t2=var(end)+1;
[~,prop_state,~,state_aps,~] = ode113(@SCR3BP_dyn,[t1 t2],var(end-6:end-1),options_ode_stop,mu_tbp,mu_v,R_v,J2_v);
state_aps(1)=state_aps(1)-(1-mu_tbp);
r_aps=norm(state_aps(1:3));
%cosntraint
c(2*k+1)=-r_aps+lb_peri;
c(2*k+2)=r_aps-ub_peri;


end