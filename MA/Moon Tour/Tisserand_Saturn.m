%%%%%%%%%%%%%%%%%%%%%%%%%% MOON TOUR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Latex font
set(groot,'defaulttextinterpreter','latex')
set(groot,'defaultLegendInterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex')
%% Load SPICE Kernels
cspice_kclear();
try
    cspice_furnsh('spice_kernels/pck00010.tpc')
    cspice_furnsh('spice_kernels/naif0012.tls')
    cspice_furnsh('spice_kernels/gm_de431.tpc')
    cspice_furnsh('spice_kernels/de440s.bsp')
    cspice_furnsh('spice_kernels/sat441.bsp')
catch
    cspice_furnsh('..\..\spice_kernels/pck00010.tpc')
    cspice_furnsh('..\..\spice_kernels/naif0012.tls')
    cspice_furnsh('..\..\spice_kernels/gm_de431.tpc')
    cspice_furnsh('..\..\spice_kernels/de440s.bsp')
end

%% MAIN SCRIPT
clear; close all; clc;

%%% DEFINE INITIAL CONSTANTS AND FORMULAS %%%
R_Saturn = 58232;
mu_Saturn = cspice_bodvrd('SATURN', 'GM', 1); % [km^3/s^2]
mu_Enceladus = cspice_bodvrd('602', 'GM', 1);

C_Tiss = @(v_inf) 3 - v_inf^2;

% NOTE: Turning angle = Alfa of the Tisserand Graph
% Moving along a contour line at v_inf approximately constant -->
% non-impulsive flyby. v_inf: non dimensional hyperbolic excess speed wrt
% main body
i_SC = 0;

% Post flyby orbit parameters
SMA_SC = @(v_inf, alpha) 1./(1 - v_inf^2 - 2*v_inf*cosd(alpha));

% For fixed v_inf: alpha goes from 0 to pi
e_SC = @(i_SC, v_inf, alpha) sqrt(1 - 1/SMA_SC(v_inf, alpha) * ((3 - 1/SMA_SC(v_inf, alpha) - v_inf^2)/(2*cosd(i_SC)))^2);

% Orbital Radii of the Moons
r_Titan = 1221.83e+3;
r_Rhea = 527.04e+3;
r_Dione = 377.40e+3;
r_Tethys = 294.66e+3;
r_Enceladus = 238.02e+3;
RR = [r_Titan,r_Rhea,r_Dione,r_Tethys,r_Enceladus];
TT = sqrt(RR.^3/mu_Saturn); % Vector of orbital periods

%%% DEFINE PARAMETERS OF THE FB's AS TAKEN BY THE REFERENCE %%%
% man stands for maneuvers and it's a 1x5 struct, where each sub-struct is
% one of the moons. Each substruct has 7 fields:
% - N:M represent the resonance with respect to the moon of the orbits
%       between each FB
% - vinf is the infinite velocity reached at the orbit with correpondence
%       resonance
% - perc is the percentile length of the lines in the plot (NOT IMPORTANT)
% - typeFB is a cell array to indicate the type of VILT. If it's empty, the
%        VILT is external, if it's 'int' the VILT is internal
% - DV is a vector with the DV required by each FB (some are unpowered)
% - T is a vector with the times between each FB (empty for the moment)
%
% All the data are taken from the paper
moons = {'Titan','Rhea','Dione','Tethys','Enceladus'};

% Titan (3)
man(1).N = [2,1,0];
man(1).M = [1,1,0];
man(1).vinf = [1.27,1.27,1.27];
man(1).perc = [1.2,1.5];
man(1).typeFB = cell(size(man(1).N));

% Rhea (12)
man(2).N = [2,11,7,8,3,7,4,6,1,6,4,0];
man(2).M = [1,6,4,5,2,5,3,5,1,7,5,0];
man(2).vinf = [1.75,1.76,1.76,1.77,1.21,1.02,0.88,0.9,0.99,0.75,0.75,0.75];
man(2).perc = [1.05,1.06,1.07,1.08,1.09,1.13,1.17,1.21,1.25,1.28,1.3,0];
man(2).typeFB = cell(size(man(2).N));
man(2).typeFB{end-2} = 'int';

% Dione (7)
man(3).N = [5,6,9,1,9,6,0];
man(3).M = [4,5,8,1,10,7,0];
man(3).vinf = [0.82,0.7,0.7,0.77,0.7,0.65,0.65];
man(3).perc = [1.05,1.06,1.1,1.14,1.15,1.2,0];
man(3).typeFB = cell(size(man(3).N));
man(3).typeFB{end-2} = 'int';
man(3).typeFB{end-1} = 'int';

% Tethys (9)
man(4).N = [6,7,9,14,1,13,9,7,0];
man(4).M = [5,6,8,13,1,14,10,8,0];
man(4).vinf = [0.77,0.7,0.7,0.7,0.7,0.67,0.67,0.63,0.63];
man(4).perc = [1.06,1.06,1.06,1.06,1.06,1.06,1.06,1.06,0];
man(4).typeFB = cell(size(man(4).N));
man(4).typeFB{end-3} = 'int';
man(4).typeFB{end-1} = 'int';

% Enceladus (9)
man(5).N = [7,15,8,17,9,10,11,13,14];
man(5).M = [6,13,7,15,8,9,10,12,13];
man(5).vinf = [0.8,0.82,0.75,0.6,0.6,0.5,0.52,0.37,0.30];
man(5).perc = [1.06,1.06,1.06,1.06,1.06,1.06,1.06,1.06];
man(5).typeFB = cell(size(man(5).N));

for ii = 1:length(man)
    man(ii).DV = zeros(size(man(ii).N));
    man(ii).T = zeros(size(man(ii).N));
end

%%% PERICENTER RAISE MANEUVER %%%
% Computes the PRM to raise the pericenter from 2.55 SR to Titan's orbit,
% to begin the sequence of FBs

ra00 = 200*R_Saturn;
ra0 = ra00;
rp00 = 2.55*R_Saturn;
rp0 = RR(1);
a00 = (rp00+ra00)/2;
a0 = (ra0+rp0)/2;
DV_PRM = sqrt(mu_Saturn)*abs(sqrt(2/ra00 - 1/a00) - sqrt(2/ra0 - 1/a0));


% Set the v_infinite at the beginning of each moon to the value taken by
% the paper (this represents a check to verify that the FBs performed are
% coherent)
v_inf0_moon = [1.9245,1.75,0.9,0.77,0.8];

% Loop over the moons
for moon = 1:length(RR)

    figure;

    %%%% Define Reference Units (different for each moon)
    DU = RR(moon);
    VU = sqrt(mu_Saturn/DU);
    TU = DU/VU;
    
    %%%% Build Tisserand plane of each moon 
    
    % Define range of v_infinite
    if moon == 5 % Enceladus' Tisserand starts from a lower v_inf
         v_inf_range_dim = linspace(0.3,v_inf0_moon(moon),6);
    else
        v_inf_range_dim = linspace(0.5,v_inf0_moon(moon),11);
    end
    v_inf_range = v_inf_range_dim./VU;

    % Define range of pump angles
    alpha_range = 0:1:180;
    
    % Initialize storing variables for the curves
    Ra_SC = zeros(length(v_inf_range),length(alpha_range));
    Rp_SC = Ra_SC;
    
    % Given alpha and v_inf, compute the outer orbit
    idv = 1;
    for v_inf = v_inf_range
        ida = 1;
        for alpha = alpha_range
            SMA = SMA_SC(v_inf, alpha);
            e = e_SC(i_SC, v_inf, alpha);
            Ra_SC(idv, ida) = SMA*(1+e)*DU;
            Rp_SC(idv, ida) = SMA*(1-e)*DU;
            ida = ida + 1;
        end
        idv = idv + 1;
    end
    
    % Draw the Tisserand
    hold on;
    cmap = jet(length(v_inf_range));
    for idv = 1:length(v_inf_range)
        h = plot(Ra_SC(idv, :), Rp_SC(idv, :), '-', 'Color', cmap(idv, :), 'LineWidth', 1.05,...
            'DisplayName',['$v_{\infty} = $',num2str(v_inf_range_dim(idv))]);
    end
    plot(RR(moon)*ones(100,1),linspace(1.05*min(Rp_SC(end,:)),RR(moon)),'k','linewidth',1.5)
    plot(linspace(RR(moon),1.05*max(Ra_SC(end,:))),RR(moon)*ones(100,1),'k','linewidth',1.5)
    t = text(0.98*RR(moon),1.02*RR(moon),moons{moon});
    t.FontSize = 9 - moon;
    t.Interpreter = 'latex';
    colormap jet
    c = colorbar;
    caxis([v_inf_range_dim(1),v_inf_range_dim(end)])
    c.Label.String = '$v_{\infty}$ [Km/s]';
    c.Label.FontSize = 11;
    c.Label.Interpreter = 'latex';
    grid minor
    axis auto
    xlabel('Apoapsis Radius [Km]');
    ylabel('Periapsis Radius [Km]');
    scatter3(ra0,rp0,20,'r','filled')
    
    %%% Build resonant orbits' lines %%%

    % Given a certain resonance N:M with respect to a certain moon, the
    % orbital period is defined, and the semi-major axis with it through
    % Kepler's third law. This function retrieves the radius of pericenter
    % from the apocenter (or viceversa) to achieve that SMA.
    rp_iso = @(ra_iso,N,M) - ra_iso + 2*(mu_Saturn*(N/M * TT(moon))^2)^(1/3);
    
    % Extract the resonant orbits considered for the moon
    N_vec = man(moon).N(1:end-1);
    M_vec = man(moon).M(1:end-1);
    
    % Plot the lines
    for id_iso = 1:length(N_vec)
        ra_v = plot_isoline(rp_iso,N_vec(id_iso),M_vec(id_iso),RR(moon),man(moon).perc(id_iso));
        plot(ra_v,rp_iso(ra_v,N_vec(id_iso),M_vec(id_iso)),'linewidth',1.5)
    end
    
    %%% COMPUTE VILTS AND FBs %%%

    % Initialize the v_inf at the beginning of the FBs with a certain moon
    v_inf0 = v_inf0_moon(moon);

    % Loop over the number of GA for each moon
    for id_man = 1:length(man(moon).N)

        % Save the parameters of the maneuver
        N = man(moon).N(id_man);
        M = man(moon).M(id_man);
        vinf_tar = man(moon).vinf(id_man);

        % Build the iso-lines for the v_inf of interest
        if v_inf0 == vinf_tar
            v_inf_range_dim = vinf_tar;
        else
            v_inf_range_dim = linspace(vinf_tar,v_inf0,2);
        end
        v_inf_range = v_inf_range_dim./VU;
        alpha_range = 0:1:180;
        Ra_SC = zeros(length(v_inf_range),length(alpha_range));
        Rp_SC = Ra_SC;
        idv = 1;
        for v_inf = v_inf_range
            ida = 1;
            for alpha = alpha_range
                SMA = SMA_SC(v_inf, alpha);
                e = e_SC(i_SC, v_inf, alpha);
                Ra_SC(idv, ida) = SMA*(1+e)*DU;
                Rp_SC(idv, ida) = SMA*(1-e)*DU;
                ida = ida + 1;
            end
            idv = idv + 1;
        end

        % Compute the VILT
        [DV1,r_plus,r_minus] = VILT(rp_iso,N,M,[Ra_SC(1,:);Rp_SC(1,:)],[Ra_SC(end,:);Rp_SC(end,:)],...
            mu_Saturn,man(moon).typeFB{id_man});
        man(moon).DV(id_man) = DV1*1000;

        % Plot the point after the FB (and after the VILT) if needed
        scatter(r_plus(1),r_plus(2),40,'b','filled')
        if (r_minus(1) ~= r_plus(1)) || (r_minus(2) ~= r_plus(2)) 
            scatter(r_minus(1),r_minus(2),40,'g','filled')
        end

        % Update v_inf0
        v_inf0 = vinf_tar;
    end

    % Update the initial orbit for next moon
    ra0 = r_plus(1);
    rp0 = r_plus(2);
    
    % Sum DV
    man(moon).DV_tot = sum(man(moon).DV);
end

DV_tot = 0;
for ind = 1:length(RR)
    DV_tot = DV_tot + man(ind).DV_tot;
end

% Compare DV_tot computed and DV_target taken by the paper
DV_tot
DV_target = 27+251+90+28+96

%%% Capture in the Halo %%%

% Adimensionalization Constants for Halo
DUH = 238411468.296/1000;
TUH = 118760.57/(2*pi); 
VUH = DUH/TUH;

% Initial state for Halo
x0_Halo=1.000062853735440;
y0_Halo=0;
z0_Halo=-0.00117884381145460;
vx0_Halo=0;
vy0_Halo=0.0168877463349484;
vz0_Halo=0;
r_Halo_Enc = [(x0_Halo-1);y0_Halo;z0_Halo];
v_Halo = [vx0_Halo;vy0_Halo;vz0_Halo];

% Compute velocity and distance of the initial state
V_pericenter = norm(v_Halo) * VUH;
Rp_target = norm(r_Halo_Enc) * DUH;

% Compute DV capture
[DV_capture] = HaloCapture(v_inf0, mu_Enceladus, Rp_target, V_pericenter);


