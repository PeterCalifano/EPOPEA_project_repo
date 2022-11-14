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
% [~, output] = ode113(@dyn, t_range, [r0;v0], options, altrevariabili);
% r1 = output(end,1:2);
% v1 = output(end,3:4);
% m1 = m_i;

% step time grid
h = (tN-t1)/(N-1);

% Retrieve par
Re = par.Re;

step_state = length(state_i);
step_var = 5+3;

% Non-linear equality constraints
% Initialization

% Cycle to fill constraints
for k = 1:(N-1)
    t_k = t1 + h*k;
    var =  var((k-1)*8+1:8*k);
    x_k = var(1:5);
    x_next = var(k*8+1:k*8+5);
    

end
%initial and final conditions

%Constraints on alpha

% derivative
if nargout > 2


end

% NB : many could be linear inequalities

C = zeros(8*N+1,1);
C(end) = -tN;
C(end-1) = t1-tN;

var = var(1:8);
x = var(1); 
y = var(2);
u = var(6);
alpha_x = var(7);
alpha_y = var(8);
C(1) = Re^2-x^2-y^2;
C(2) = -u;
C(3) = -u-1;
C(4) = -alpha_x;
C(5) = -alpha_x-1;
C(6) = -alpha_y;
C(7) = -alpha_y-1;
for k = 2:N
    var = var((k-1)*8+1:8*k);
    x = var(1); 
    y = var(2);
    m = var(5);
    u = var(6);
    alpha_x = var(7);
    alpha_y = var(8);
    C(8*(k-2)+8) = Re^2-x^2-y^2;
    C(8*(k-2)+9) = -u;
    C(8*(k-2)+10) = -u-1;
    C(8*(k-2)+11) = -alpha_x;
    C(8*(k-2)+12) = -alpha_x-1;
    C(8*(k-2)+13) = -alpha_y;
    C(8*(k-2)+14) = -alpha_y-1;
    C(8*(k-2)+15) = m_dry -m;
end
JC = [];
JCeq = [];
end