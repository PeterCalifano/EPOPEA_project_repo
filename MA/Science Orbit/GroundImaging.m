function [track_length,Nptot,Data_tot,rev,Nporb,Data_Orb,left_swath,right_swath] = GroundImaging(t_vec,pos_Enc,R_Enc,theta0_g,w_p,FoV,data_pic,Overlap,Max_Pics)
%% PROTOTYPE
% [N_pics,Data_Orb] = GroundImaging(t_vec,pos_Enc,R_Enc,theta0_g, w_p,Cam,R_Enc)

%% DESCRIPTION
% Computes the number of measurements and the consequent data volume during
% a remote sensing arc

%% INPUT
% t_vec [1xN] - [s] time history
% pos_Enc  [3xN] - [km] array of positions wrt Enc
% R_Enc [1x1] - [km] Enceladus radius
% theta0_g [1x1] - [deg] initial longitude
% w_p [1x1] - [rad/s] rotational velocity of Enc (0)
% FoV [1x1] - [deg] field of view of the instrument
% data_pic [1x1] - [Gb] data per picture
% Overlap [1x1] - in 1,..,100, required % overlap between images 
% Max_pics [1x1] - max number of pictures per orbit (if TMTC and OBDH have the braccino corto)
%% OUTPUT
% mo 
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% Date, User, brief summary of the modification
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
%

%% Function code


% latitude and longitude computation
[~,~, lat_deg, lon_deg] = groundTrack(t_vec, pos_Enc',theta0_g,w_p);
lat=deg2rad(lat_deg);
lon=deg2rad(lon_deg);

% track length and steps initialization
track_length=0;
steps_array=zeros(1,length(t_vec)-1); %array of the distance on the surface covered between time instants 
track_array=zeros(3,length(t_vec));


%first footprint on the surface on Enceladus
track_old=R_Enc*[cos(lat(1))*cos(lon(1));
                   cos(lat(1))*sin(lon(1));
                   sin(lat(1))];
track_array(:,1)=track_old;

%imaging - initialization
FoV2=deg2rad(FoV/2);
Nptot=1; %total number of pictures, first picture immediatly taken
h=norm(pos_Enc(:,1))-R_Enc; %first altitude
l=2*h*tan(FoV2); %size of the first image

dist=-l/2; %distance covered without taking a picture, at the beginning of the first step
% a full picture is taken 

%arrays for plots of the swath
left_pos=zeros(3,length(t_vec)-1);
right_pos=zeros(3,length(t_vec)-1);
for k=1:length(t_vec)-1
    %recovery of the angle between two consecutive tracklets on the surface
    %Position of the k+1th footprint
    track_k1=R_Enc*[cos(lat(k+1))*cos(lon(k+1));
                   cos(lat(k+1))*sin(lon(k+1));
                   sin(lat(k+1))];
               
    path_angle=abs(acos(dot(track_old/norm(track_old),track_k1/norm(track_k1)))); %spherical arc covered in an integration step
    step=2*R_Enc*path_angle; %[km] - distance covered between two consecutive time instants 
    %data collection
    steps_array(k)=step;
    track_length=track_length+step;%total distance
    track_array(:,k+1)=track_k1; 

    %imaging study
    h=norm(pos_Enc(:,k))-R_Enc; %altitude
    l=2*h*tan(FoV2); %size of the image
    Ol=l*(1-Overlap/100); %size of the non-overlapped part
    
    dist=step+dist; %distance without taking a picture 
    Nptot=Nptot+floor(dist/Ol); %number of pictures needed to cover the distance
    dist=dist-floor(dist/Ol)*Ol; %residual unimaged area

    %chicchina plot - cit non funge
    step_dir=(track_k1-track_old)/norm(track_k1-track_old);
    dir_swath=cross(step_dir,pos_Enc(:,k)); %perpendicular direction to the step    
    dir_swath=dir_swath/norm(dir_swath);
    
    left_pos(:,k)=track_old+dir_swath*l/2;
    right_pos(:,k)=track_old-dir_swath*l/2;
    
    
    %update
    track_old=track_k1;
end
if dist>0
    Nptot=Nptot+1; %last picture to cover any residual distance
end

%check with the maximum data volume per orbit
Data_tot=data_pic*Nptot; %total amount of data
if Nptot>Max_Pics
   rev=ceil(Nptot/Max_Pics);
   Nporb=Nptot/rev;
   Data_Orb=Nporb*data_pic;
else
    Data_Orb=Data_tot;
    Nporb=Nptot;
    rev=1;
end

[~,~, left_lat, left_lon] = groundTrack(t_vec(1:end-1),left_pos',theta0_g,w_p);

left_swath=[left_lon;left_lat];

[~,~, right_lat, right_lon] = groundTrack(t_vec(1:end-1),right_pos',theta0_g,w_p);

right_swath=[right_lon;right_lat];


end