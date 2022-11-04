clear;
close all
clc;

% TISSERAND PLANE FOR SATURN SYSTEM
% R_S = 58232;
% r_a_vec = R_S*linspace(21,200);
% r_p_vec = R_S*linspace(2,20);
% [RA,RP] = meshgrid(r_a_vec,r_p_vec);
% 
% % Compute Tisserand parameter
% C_tiss = 2./(RA+RP) + 2 * sqrt((RA+RP)/2 .* (1 - ((RA-RP)./(RA+RP)).^2));
% 
% % Conversion to Km/s
% C_tiss = C_tiss/1000;
% 
% % Compute v_inf 
% v_inf = sqrt(3-C_tiss);
% 
% contour(RA/R_S,RP/R_S,v_inf,10)
% grid on
% colorbar 
% xlabel('Apoapsis [R_S]')
% ylabel('Periapsis [R_S]')

addpath("..\GeneralUtilities\")
cspice_kclear();
cspice_furnsh('..\..\EPOPEA_project_repo\EPOPEA_metakernel.tm')

mu_Saturn = cspice_bodvrd('SATURN', 'GM', 1); % [km/s^2]

% Reference Units (Enceladus orbit)
DU = 2.302500200000000e+05; % [km] % Radius of the circular orbit
TU = 2*pi*sqrt(DU^3./mu_Saturn); % [s] % Period of revolution 
VU = DU/TU; % Circular velocity of the secondary body (Moons)

C_Tiss = @(v_inf) 3 - v_inf^2;

% NOTE: Turning angle = Alfa of the Tisserand Graph
% Moving along a contour line at v_inf approximately constant -->
% non-impulsive flyby. v_inf: non dimensional hyperbolic excess speed wrt
% main body
i_SC = 0;

% Post flyby orbit parameters
SMA_SC = @(v_inf, alpha) 1./( 1 - v_inf^2 - 2*v_inf*cos(alpha) );

% For fixed v_inf: alpha goes from 0 to pi
e_SC = @(i_SC, v_inf, alpha) sqrt( 1 - (1/SMA_SC(v_inf, alpha))*( (3 - 1/SMA_SC(v_inf, alpha) - v_inf^2)/(2*cos(i_SC)) )^2);

% 
% TT = 2 * np.pi * np.sqrt((A_SC * R_body) ** 3 / body.parent.k)
% EE = -body.parent.k / (2 * A_SC * R_body)

% SMA = SMA_SC(v_inf);
% e = e_SC(i_SC, v_inf);
% Ra_SC = SMA*(1+e);
% Rp_SC = SMA*(1-e);

v_inf_range = 1:0.5:10;
v_inf_range = v_inf_range./VU;
alpha_range = linspace(0, pi, 300);

Ra_SC = nan(length(v_inf_range), length(alpha_range));
Rp_SC = nan(length(v_inf_range), length(alpha_range));


idv = 1;
for v_inf = v_inf_range
    ida = 1;
    for alpha = alpha_range
        SMA = SMA_SC(v_inf, alpha);
        e = e_SC(i_SC, v_inf, alpha);

        if e < 1 || SMA > 0
            Ra_SC(idv, ida) = SMA*(1+e);
            Rp_SC(idv, ida) = SMA*(1-e);
        end

        ida = ida + 1;
    end
    idv = idv + 1;
end


figure;
hold on;
cmap = jet(length(v_inf_range));
for idv = 1:length(v_inf_range)
plot(Ra_SC(idv, :), Rp_SC(idv, :), '.-', 'Color', cmap(idv, :), 'LineWidth', 1.05);

end
grid minor
axis auto
xlabel('Apoapsis Radius [DU]');
ylabel('Periapsis Radius [DU]');
