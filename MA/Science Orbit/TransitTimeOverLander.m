clear;close all;clc;
options_ode=odeset('RelTol',1e-13,'AbsTol',1e-13);
%units
mu=1.90095713928102*1e-7;
DU=238411468.296/1000; %km
TU=118760.57/(2*pi); 

R_enc0=252.1; %km, Enceladus mean radius
R_enc=R_enc0/DU;


% Times
t0=0;
FlightDays=10; %days of prapagation
tf=7.5*60/TU; %final time of propagation
t_span = t0:(1/TU):tf;

%%% Pinpoint parameters %%%
we = 0.191;
lonlat = [-80; 90]; th_e0 =deg2rad(0);
lat = deg2rad(lonlat(1));
lon = deg2rad(lonlat(2));
Re = 1;

% From latitudinal to cartesian
r_xy = Re*cos(lat);
Yrot = r_xy*sin(lon);
Zrot = Re*sin(lat);
Xrot = r_xy*cos(lon);
vec_rot = [Xrot; Yrot; Zrot];

% Enceladus rotation
th_e = @(tf) we*(tf-t0)*TU/3600 + th_e0;

% From rotating enceladus to IAU_Enceladus
A_rot2IAU = @(th_e) [cos(th_e)     sin(th_e) 0;
            -sin(th_e)    cos(th_e) 0     ;
            0        0         1];
% desired final state
vec_pp_0 = A_rot2IAU(th_e(t0))*vec_rot*R_enc0;
vec_pp_f = A_rot2IAU(th_e(tf))*vec_rot*R_enc0;


%%% Halo Initial State %%%
x0_Halo=1.000062853735440;
y0_Halo=0;
z0_Halo=-0.00117884381145460;
vx0_Halo=0;
vy0_Halo=0.0168877463349484;
vz0_Halo=0;

state0_Halo=[x0_Halo,y0_Halo,z0_Halo,vx0_Halo,vy0_Halo,vz0_Halo]';


%propagation - Halo
[t_vec_Halo,state_vec_Halo]=ode113(@(t,x) CR3BP_dyn(t,x,mu),t_span,state0_Halo,options_ode);
state_vec_Halo=state_vec_Halo';

%Position wrt Enceladus
pos_Halo_Enc=state_vec_Halo(1:3,:);
pos_Halo_Enc(1,:)=pos_Halo_Enc(1,:)+mu-1;

%in the inertial frame
%conversion of the trajectory to the inertial Enceladus Centred frame
r_vec_Halo_in=zeros(3,length(t_vec_Halo));
for k=1:length(t_vec_Halo)
    r_vec_Halo_in(1,k)=(state_vec_Halo(1,k)+mu-1)*cos(t_vec_Halo(k))-state_vec_Halo(2,k)*sin(t_vec_Halo(k));
    r_vec_Halo_in(2,k)=(state_vec_Halo(1,k)+mu-1)*sin(t_vec_Halo(k))+state_vec_Halo(2,k)*cos(t_vec_Halo(k));
    r_vec_Halo_in(3,k)=state_vec_Halo(3,k);
end

%dimension recovery
state_vec_Halo(1:3,:)=state_vec_Halo(1:3,:)*DU;
state_vec_Halo(4:6,:)=state_vec_Halo(4:6,:)*DU/TU;
r_vec_Halo_in=r_vec_Halo_in*DU;
pos_Halo_Enc=pos_Halo_Enc*DU;


% Compute angle at beginning and end
diff_pp_0 = - vec_pp_0 + r_vec_Halo_in(:,1);
diff_pp_f = - vec_pp_f + r_vec_Halo_in(:,end);
th_initial =  acosd(dot(r_vec_Halo_in(:,1),diff_pp_0)/(norm(r_vec_Halo_in(:,1))*norm(diff_pp_0)))
th_final =  acosd(dot(r_vec_Halo_in(:,end),diff_pp_f)/(norm(r_vec_Halo_in(:,end))*norm(diff_pp_f)))

% Plot
Enceladus_3D(R_enc*DU,[0,0,0])
plot3(r_vec_Halo_in(1,:),r_vec_Halo_in(2,:),r_vec_Halo_in(3,:),'linewidth',2)
hold on
scatter3(vec_pp_0(1),vec_pp_0(2),vec_pp_0(3),30,'filled')
scatter3(vec_pp_f(1),vec_pp_f(2),vec_pp_f(3),70,'filled')
grid minor 
xlabel('X')
ylabel('Y')
zlabel('Z')

