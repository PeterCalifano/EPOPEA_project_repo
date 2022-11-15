function [ C, Ceq, JC, JCeq ] = land_nonlincon(var, state_i, par, N)

%Initial state: probably apocenter
r_i = state_i(1,2);
v_i = state_i(3,4);
m_i = state_i(5);

%Initial and final time
t1 = var(end-1);
tN = var(end);

% Orbit Propagator
% t_range = [0 t1];
% [~, output] = ode113(@dyn, t_range, [r_i;v_i], options, altrevariabili);
% r1 = output(end,1:2);
% v1 = output(end,3:4);
% m1 = m_i;

% step time grid
h = (tN-t1)/(N-1);

% Retrieve par
Re = par(5);

step_st = length(state_i);           % 5: rr,vv,m
step_var = 5+3;                      % 8: rr,vv,m,u,ax,ay

% Non-linear equality constraints
% Initialization
vec_def = zeros(5*(N-1),1);
vec_alpha = zeros(N,1);
C = zeros((N-1),1);

% Cycle to fill constraints
for k = 1:(N-1)
    %defects
    t_k = t1 + h*k;
    var_k = var(step_var*(k-1)+1:step_var*k);
    x_k = var_k(1:5);
    x_next = var(step_var*k+1:step_var*k+step_st);
    uk = var_k(6:8);

    f = landing_dyn(~, x_k, uk,par);
    zk = x_next - x_k - h.*f;

    %thrust versor
    q_k = norm(var_k(7:8)) - 1;

    %columns for Ceq
    vec_def(step_st*(k-1)+1:step_st*(k-1)+5) = zk;
    vec_alpha(k) = q_k;

    %Inequality constraints
    C(k) = Re - norm(x_k(1:2));

end
vec_alpha(end) = norm(var(end-3:end-2)) - 1;

%initial and final conditions
var_N = var(end-9:end-2);
r_N = var_N(1:2);
v_N = var_N(3:4);
psi_i = var(1:5) - [r1; v1; m1];
psi_f = [r_N(1); r_N(2); v_N(1); v_N(2)] - [0; -Re; 0; 0]; 

%fill Ceq
Ceq = [vec_def; vec_alpha; psi_i; psi_f];              %check column vecs

% derivative
if nargout > 2
JCeq = [];
JC = [];
end



end