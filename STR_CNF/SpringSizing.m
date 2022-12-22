clear ; close all ; clc ;

% This function is used to size the springs at the orbiter/lander
% interface, used during separation

% For the sizing, a frame fixed on the initial clamped configuration is
% considered
% Once the two s/c are separated, one spacecraft will move in the positive
% x-direction of this frame, and the other will move in the opposite
% direction
% Hence, the relative velocity is considered with respect to the initial
% frame fixed in the clamped configuration

% Lander data
m_lan = 788.92 ; % [kg] - Mass at separation

% Requested relative velocity
v_lan_f = 13e-2 ; % [12 cm/s] - Relative velocity of lander with respect to frame initially fixed with clamped configuration

% Mean acceleration at release
ReleaseTime = 1 ; % [seconds]
a = v_lan_f / ReleaseTime ;
fprintf(['Acceleration a = ', num2str(a), ' m/s^2\n']  ) ;

% Compute required stiffness
DX = 2e-2 ; % [5cm] - Length of spring travel before release
keq = m_lan * a / DX ;
fprintf(['Equivalent stiffness k = ', num2str(keq*1e-3), ' N/mm\n']  ) ;

% Select number of springs
Spring_number = 4 ; 
k_singleSpring = Spring_number * keq ; % The springs are in parallel
fprintf(['Stiffness of single spring k = ', num2str(k_singleSpring*1e-3), ' N/mm\n']  ) ;


%% From website

% Compute preload of each spring
F0 = k_singleSpring * DX ;



function meter = inch2meter(inch)

meter = inch / 39.37 ;

end

function inch = meter2inch(meter)

inch = meter * 39.37 ;

end