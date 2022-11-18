clear,clc;close all
cspice_kclear();




%% kernel creation - don't need to run it
%system('pinpoint -def estrack.def -pck pck00010.tpc -spk estrack.bsp -fk estrack.tf');

%%
%Kernel loading - keep the metakernel out of the kernels folder, if needed look into
%it to see what kernels it loads
cspice_furnsh('..\spice_kernels\naif0012.tls');
cspice_furnsh('..\spice_kernels\sat441.bsp');
cspice_furnsh('..\spice_kernels\de440s.bsp');
cspice_furnsh('..\spice_kernels\pck00010.tpc');
cspice_furnsh('..\spice_kernels\estrack.bsp');
cspice_furnsh('..\spice_kernels\estrack.tf');

%%
%Time span - play with the dates here
t_beginning_GS_UTC='2040-01-01-00:00:00.000 UTC';
t_end_GS_UTC='2040-01-02-00:00:00.000 UTC';
t0_GS=cspice_str2et(t_beginning_GS_UTC);
tf_GS=cspice_str2et(t_end_GS_UTC);
t_span_GS=linspace(t0_GS,tf_GS,1000);

%Sun shielding angle - choose here the angle between Sun and Saturn system
%in the sky from which it can be said that Saturn is not in view
ecl_ang_Sun= 5; % [°]

ecl_ang_Sun=deg2rad(ecl_ang_Sun);
%% 
%Visibility windows for the 3 GS's - in the call to the function  type 602 for Enceladus, 699 for Saturn
% The plots present the altitude in the sky in time for each ground station

% plots
t_span_GS=(t_span_GS-t0_GS)/3600;
[rho_cer_Enc, azi_cer_Enc, ele_cer_Enc] = cspice_reclat(r_cer_Enc);


% CEBREROS
[rho_GS_body, azi_GS_body, ele_GS_body]=visibility_window_GS('602','ESA-CEBREROS',t_beginning_GS_UTC,t_end_GS_UTC,1000,0);
%NNO
[rho_GS_body, azi_GS_body, ele_GS_body]=visibility_window_GS('602','ESA-NNO',t_beginning_GS_UTC,t_end_GS_UTC,1000,0);
%MALARGUE
[rho_GS_body, azi_GS_body, ele_GS_body]=visibility_window_GS('602','ESA-MALARGUE',t_beginning_GS_UTC,t_end_GS_UTC,1000,0);








% %% Eclipse - work in progresssss
% 
% %time span - play with the dates here
% t_beginning_UTC='2040-01-01-00:00:00.000 UTC';
% t_end_UTC='2041-01-01-00:00:00.000 UTC';
% t0=cspice_str2et(t_beginning_UTC);
% tf=cspice_str2et(t_end_UTC);
% t_span_whole=linspace(t0,tf,100000);
% 
% rE_Sun=cspice_spkpos('SUN',t_span_whole,'IAU_EARTH','NONE','EARTH'); %Sun wrt Earth
% rE_Saturn=cspice_spkpos('SATURN',t_span_whole,'IAU_EARTH','NONE','EARTH'); %saturn wrt Earth
% rE_Enceladus=cspice_spkpos('602',t_span_whole,'IAU_EARTH','NONE','EARTH'); %Enceladus wrt Earth
% 
% 
% 
% % rSun_Saturn=cspice_spkpos('699',t_span_whole,'ECLIPJ2000','NONE','SUN');  %Saturn wrt Sun
% % rSun_Earth=cspice_spkpos('EARTH',t_span_whole,'ECLIPJ2000','NONE','SUN');  %Earth wrt Sun
% % rSun_Enceladus=cspice_spkpos('602',t_span_whole,'ECLIPJ2000','NONE','SUN');  %Enceladus wrt Sun
% 
% 
% Radii_Sun=cspice_bodvrd('SUN','RADII',3); %Sun mean radius
% Radii_Saturn=cspice_bodvrd('SATURN','RADII',3); %Saturn mean radius
% R_Sun=Radii_Sun(1);
% R_Saturn=Radii_Saturn(1);
% 
% Vis_Enc=zeros(size(t_span_whole)); %visibility of Enceladus wrt Earth
% Vis_Satsys=ones(size(t_span_whole)); %visibility of the saturn system (sun eclipse check)
% 
% 
% for k=1:length(t_span_whole)
% 
%     dir_Saturn=rE_Saturn/norm(rE_Saturn);
%     dir_Enceladus=rE_Enceladus/norm(rE_Enceladus);
%     alpha=acos(dot(dir_Saturn,dir_Enceladus));
%     h(k)=abs(norm(rE_Saturn)*tan(alpha));
%     if h(k)>R_Saturn
%        Vis_Enc(k)=1;
%     elseif norm(rE_Enceladus)<norm(rE_Saturn)
%        Vis_Enc(k)=1;
%     else
%         Vis_Enc(k)=0;
%     end
% 
% end
% 
% % for k=1:length(t_span_whole)
% %     rE_Saturn=rSun_Saturn(:,k)-rSun_Earth(:,k);
% %     rE_Sun=-rSun_Earth(:,k);
% %     
% %     dir_Saturn=rE_Saturn/norm(rE_Saturn);
% %     dir_Sun=rE_Sun/norm(rE_Sun);
% %     alpha=acos(dot(dir_Saturn,dir_Sun));
% %     h(k)=abs(norm(rE_Sun)*tan(alpha));
% %     if h(k)>R_Sun
% %        Vis_Satsys(k)=1; 
% %     elseif norm(rE_Saturn)<norm(rE_Sun)
% %        Vis_Satsys(k)=1;
% %     else
% %         Vis_Satsys(k)=0;
% %     end
% %         
% % end
% % 
% 
%   
% figure
% grid on
% grid minor
% plot(t_span_whole,Vis_Satsys,'b','linewidth',1.5);
% xlabel('time')
% ylabel('visibility (boolean)')
% legend('visibility of Enceladus from the Earth')
% 
% 
% figure
% plot3(rSun_Saturn(1,:)-rSun_Earth(1,:),rSun_Saturn(2,:)-rSun_Earth(2,:),rSun_Saturn(3,:)-rSun_Earth(3,:))
% axis equal
% 
% %% prova
% rS_Enc=cspice_spkpos('602',t_span_whole,'J2000','NONE','6');
% 
% figure
% plot3(rS_Enc(1,:),rS_Enc(2,:),rS_Enc(3,:))
% axis equal



