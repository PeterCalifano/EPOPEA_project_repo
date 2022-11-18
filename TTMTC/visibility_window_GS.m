function [rho_GS_body, azi_GS_body, ele_GS_body]=visibility_window_GS(body,GS,t_in_GS,t_end_GS,n_points,angle_Sun)

%% PROTOTYPE
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Computes the plots of elevation of the chosen body wrt the selected
% gorund station in the time span. Sun interference taken in account
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% body [1x1] string - Spice compatible ID of the body to inspect
% GS [1x1] string - Kernel compatible ID of the selected GS
% t_in_GS [1x1] UTC string - starting time of observation, UTC format string
% t_end_GS [1x1] UTC string - final time of observation, UTC format string
%                          n_points [1x1] double
% angle_Sun [1x1] double [°] - relative angle between sun and body below
% which the observation is too disturbed
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% rho_GS_body [1xn_points] double [km] radius vector GS-body
% azi_GS_body [1xn_points] double [°] array of azimuth in time
% ele_GS_body [1xn_points] double [°] array of elevation in time
% plot of the elevation in time
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 03/11/22 - Fabri's first implementation 
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% spice
% -------------------------------------------------------------------------------------------------------------

% Code
%time conversion
t0_GS=cspice_str2et(t_in_GS);
tf_GS=cspice_str2et(t_end_GS);
t_span_GS=linspace(t0_GS,tf_GS,n_points);

%shielding angle
ecl_ang_Sun=deg2rad(angle_Sun);

%initialization
%body
rho_GS_body=zeros(3,n_points);
azi_GS_body=zeros(1,n_points);
ele_GS_body=zeros(1,n_points);

rho_GS_Sun=zeros(1,n_points);
azi_GS_Sun=zeros(1,n_points);
ele_GS_Sun=zeros(1,n_points);
flag_r=0; %flag for legend
flag_k=0; %flag for legend
RF=strcat(GS,'_TOPO');
figure
hold on
for k=1:n_points
    %r vectors recovery
    r_GS_body=cspice_spkpos(body,t_span_GS(k),RF, 'NONE', GS);
    r_GS_Sun=cspice_spkpos('SUN',t_span_GS(k),RF, 'NONE', GS);
    %rho, azimuth and elevation computation
    [rho_body, azi_body, ele_body] = cspice_reclat(r_GS_body);
    rho_GS_body(k)=rho_body;
    azi_GS_body(k)=azi_body;
    ele_GS_body(k)=ele_body;
    
    [rho_Sun, azi_Sun, ele_Sun] = cspice_reclat(r_GS_Sun);
    rho_GS_Sun(k)=rho_Sun;
    azi_GS_Sun(k)=azi_Sun;
    ele_GS_Sun(k)=ele_Sun;
    % plot
    if abs(azi_GS_body(k)-azi_GS_Sun(k))<ecl_ang_Sun && abs(ele_GS_body(k)-ele_GS_Sun(k))<ecl_ang_Sun
        P2=scatter((t_span_GS(k)-t0_GS)/3600,rad2deg(ele_GS_body(k)),'k','filled','linewidth',0.5);
        flag_k=1;
    else
        P1=scatter((t_span_GS(k)-t0_GS)/3600,rad2deg(ele_GS_body(k)),'r','filled','linewidth',0.5);
   end
end
grid on
grid minor
xlabel('time (hours)')
ylabel('elevation')
title(strcat('Elevation plot of ',body,'wrt ',GS))
if flag_k==0
    legend([P1],'Elevation in time of the chosen body')
else
    legend([P1 P2],'Elevation in time of the chosen body','Sun interference')
end