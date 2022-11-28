function [ C, Ceq, JC, JCeq ] = land_nonlincon(var, state_i, par, N)

%Initial REFERENCE state on Halo orbit
r_i = state_i(1:3)';
v_i = state_i(4:6)';

m_i = state_i(7);

%Initial and final time
t1 = var(end-1);
tN = var(end);

cd('..');
cd('Science Orbit')
% Orbit Propagator
mass_ratio = par(6);
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
t_range = [0 t1];
[~, output] = ode113(@CRTBP_dyn, t_range, [r_i;v_i], options, mass_ratio);
state_rot_fin = output(end,1:6);     %CHECK WHETHER ITS COL OR ROW
cd('..')
cd('Landing')

% POST PROCESSING
% From rotating Saturn-Enceladus to IAU_Enceladus
state_in = rot2iau_enc(tN, state_rot_fin, mass_ratio);
r_in = state_in(1:3);
v_in = state_in(4:6);

% DCM Matrix: rigid rotation around X
A_rotx = [1  0  0
          0  0  1
          0 -1  0];
r1 = A_rotx*r_in;
v1 = A_rotx*v_in;
m1 = m_i;
R_orb = norm(r1);

% step time grid
h = (tN-t1)/(N-1);

% Retrieve par
Re = par(5);

step_st = length(state_i);           % 7: rr,vv,m
step_var = step_st+4;                % 11: rr,vv,m,u,ax,ay,az

% Non-linear equality constraints
% Initialization
vec_def = zeros(step_st*(N-1),1);
vec_alpha = zeros(N,1);
% C = zeros(2*N,1);
C = zeros(2*N-1,1);

% Cycle to fill constraints
for k = 1:(N-1)
    %defects
    t_k = t1 + h*k;
    t_next = t_k + h;
    t_c = (t_k + t_next)/2;

    var_k = var(step_var*(k-1)+1:step_var*k);
    var_next = var(step_var*k+1:step_var*(k+1));

    x_k = var_k(1:7);
    x_next = var_next(1:7);

    u_k = var_k(step_st+1:step_st+4);
    u_next = var_next(step_st+1:step_st+4);
    u_c = (u_k + u_next)/2;

    f_k = landing_dyn(t_k, x_k, u_k, par);
    f_next = landing_dyn(t_next, x_next, u_next, par);

    % Hermite-Simpson
    x_c = 0.5 * (x_k + x_next) + (h/8) * (f_k - f_next);
    f_c = landing_dyn(t_c, x_c, u_c, par);
    zk = x_k - x_next + (h/6) * (f_k + 4*f_c + f_next);

    %thrust versor
    q_k = norm(var_k(step_st+2:step_st+4)) - 1;

    %columns for Ceq
    vec_def(step_st*(k-1)+1:step_st*(k-1)+step_st) = zk;
    vec_alpha(k) = q_k;

    %Inequality constraints: CHECK (REASON) ON SECOND CONSTRAINT
    C((k-1)*2+1:2*k) = [Re - norm(x_k(1:3)); norm(x_k(1:3)) - R_orb];

end
vec_alpha(end) = norm(var(end-4:end-2)) - 1;

%initial and final conditions
var_N = var(end-12:end-2);
r_N = var_N(1:3);
v_N = var_N(4:6);
psi_i = var(1:step_st) - [r1; v1; m1]; 
psi_f = [norm(r_N); norm(v_N)] - [Re; 0];

%fill Ceq
Ceq = [vec_def; vec_alpha; psi_i; psi_f];              %check column vecs

%add Final constraint to C (inequality)
% lim defines semi-angle of possible landing cone 
lim = deg2rad(20);          
% C(end-1:end) = [abs(r_N(1)); r_N(2)] - [Re*sin(lim); -Re*cos(lim)];
C(end) = r_N(2) + Re*cos(lim);

% derivative
if nargout > 2
JCeq = [];
JC = [];
end



end