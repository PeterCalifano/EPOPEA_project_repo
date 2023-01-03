% clear;
% close all
% clc;
% 
% %TISSERAND PLANE FOR SATURN SYSTEM
% R_S = 58232;
% r_a_vec = R_S*linspace(20,200);
% r_p_vec = R_S*linspace(2,4);
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
%
clear; close all;
clc;
%addpath("..\GeneralUtilities\")
%cspice_kclear();
%cspice_furnsh('..\..\EPOPEA_project_repo\EPOPEA_metakernel.tm')
%% Set Latex font
set(groot,'defaulttextinterpreter','latex')
set(groot,'defaultLegendInterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex')
%% ORBITAL FORMULAS
mu_Saturn = cspice_bodvrd('SATURN', 'GM', 1); % [km^3/s^2]

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

% SELECT PLANET

% Periods 
T_Enceladus = 1.370218*3600*24;
T_Tethys = 1.887802*3600*24; 
T_Dione = 2.736915*3600*24;
T_Rhea = 4.517500*3600*24;
T_Titan = 15.945421*3600*24;
% TT = [T_Titan,T_Rhea,T_Dione,T_Tethys,T_Enceladus];

% Orbital Radii
R_Saturn = 58232;
R_Titan = 1221.83e+3;
R_Rhea = 527.04e+3;
R_Dione = 377.40e+3;
R_Tethys = 294.66e+3;
R_Enceladus = 230.25e+3;
RR = [R_Titan,R_Rhea,R_Dione,R_Tethys,R_Enceladus];
TT = sqrt(RR.^3/mu_Saturn);
man(1).N = [2,1];
man(1).M = [1,1];
man(1).vinf = [1.27,1.27];
man(1).perc = [1.3,1.5];
moons = {'Titan','Rhea','Dione','Tethys','Enceladus'};

v_inf0 = 1.9245;

figure;

% Loop over the moons
for moon = 1%:2%length(RR)
    
    % Define Reference Units
    DU = RR(moon);
    VU = sqrt(mu_Saturn/DU);
    TU = DU/VU;
    
    % Build Tisserand plane once and for all
    v_inf_range_dim = linspace(0.5,v_inf0,11);
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
    
    hold on;
    cmap = jet(length(v_inf_range));
    for idv = 1:length(v_inf_range)
        h = plot(Ra_SC(idv, :), Rp_SC(idv, :), '-', 'Color', cmap(idv, :), 'LineWidth', 1.05,...
            'DisplayName',['$v_{\infty} = $',num2str(v_inf_range_dim(idv))]);
    end
    plot(RR(moon)*ones(100,1),linspace(1.05*min(Rp_SC(end,:)),RR(moon)),'k','linewidth',1.5)
    plot(linspace(RR(moon),1.05*max(Ra_SC(end,:))),RR(moon)*ones(100,1),'k','linewidth',1.5)
    t = text(RR(moon)-40000,RR(moon)+30000,moons{moon});
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
    scatter3(200*R_Saturn,RR(1),20,'r','filled')
    
    % Build isolines
    rp_iso = @(ra_iso,N,M) - ra_iso + 2*(mu_Saturn*(N/M * TT(moon))^2)^(1/3);
    
    N_vec = man(moon).N;
    M_vec = man(moon).M;

    for id_iso = 1:length(N_vec)
        ra_v = plot_isoline(rp_iso,N_vec(id_iso),M_vec(id_iso),RR(moon),man(moon).perc(id_iso));
        plot(ra_v,rp_iso(ra_v,N_vec(id_iso),M_vec(id_iso)),'linewidth',1.5)
    end

    % Loop over the number of GA for each moon
    for id_man = 1:length(man(moon).N)

        % Save the parameters of the maneuver
        N = man(moon).N(id_man);
        M = man(moon).M(id_man);
        vinf_tar = man(moon).vinf(id_man);

        % BUILD TISSERAND
        if v_inf0 == vinf_tar
            v_inf_range_dim = vinf_tar;
        else
            v_inf_range_dim = linspace(vinf_tar,v_inf0,10);
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
        
        [DV1,r_plus,r_minus] = VILT(rp_iso,N,M,[Ra_SC(1,:);Rp_SC(1,:)],[Ra_SC(end,:);Rp_SC(end,:)],mu_Saturn);
        scatter(r_plus(1),r_plus(2),40,'b','filled')
        if r_minus(2) ~= r_plus(2)
            scatter(r_minus(1),r_minus(2),40,'g','filled')
        end

        % Update v_inf0
        v_inf0 = vinf_tar;
    end

end



