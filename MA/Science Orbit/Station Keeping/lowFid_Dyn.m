function dxdt=lowFid_Dyn(t,x,muSat,muSun)




%%OPEN POINTS:
%-ASK IF THE APPROACH MAKES SENSE
%-SPHERICAL HARMONICS TO BE IMPLEMENTED
%-ENCELADUS AN EPHEMERIDES TO BE IMPLEMENTED
%-VALIDATION

%initialization - Saturn Centred orbital frame
dxdt=zeros(6,1);

r=x(1:3);
v=x(4:6);


dxdt(1:3)=[v(1);
           v(2);
           v(3)];
 
%Sun to Saturn position

%Saturn's orbital elements in the ecliptic ref frame
[KepSat,~] = uplanet(t,6);
i_Sat=KepSat(3);
OM_Sat=KepSat(4);
om_Sat=KepSat(5);
theta_Sat=KepSat(6);

%z-Rotation, om+theta 
R_an_Sat=[cos(om_Sat+theta_Sat) sin(om_Sat+theta_Sat) 0;
          -sin(om_Sat+theta_Sat) cos(om_Sat+theta_Sat) 0;
          0 0 1];
%x-Rotation, i 
R_i_Sat=[1 0 0;
         0 cos(i_Sat) sin(i_Sat);
         0 -sin(i_Sat) cos(i_Sat)];
%z-Rotation OM
R_OM_Sat=[cos(OM_Sat) sin(OM_Sat) 0;
          -sin(OM_Sat) cos(OM_Sat) 0;
          0 0 1];
      
%Saturn's position in the ecpliptic Sun Centred RF
[rSat,~]=kep2car(KepSat,muSun);

%Sun's position in the ecpliptic Saturn Centred RF
rSun_Sat_Ecl=-rSat;

%Sun's state in Saturn's orbital frame
rSun_Sat_F=R_an_Sat*R_i_Sat*R_OM_Sat*r_Sun_Sat_Ecl;


%Enceladus' orbital elements and state in the saturn ref frame
[KepEnc, rEnc, vEnc] = Enc_eph(t);

i_Enc=KepEnc(3);
OM_Enc=KepEnc(4);
om_Enc=KepEnc(5);
theta_Enc=KepEnc(6);

%z-Rotation, om+theta 
R_an_Enc=[cos(om_Enc+theta_Enc) sin(om_Enc+theta_Enc) 0;
          -sin(om_Enc+theta_Enc) cos(om_Enc+theta_Enc) 0;
          0 0 1];
%x-Rotation, i 
R_i_Enc=[1 0 0;
         0 cos(i_Enc) sin(i_Enc);
         0 -sin(i_Enc) cos(i_Enc)];
%z-Rotation OM
R_OM_Enc=[cos(OM_Enc) sin(OM_Enc) 0;
          -sin(OM_Enc) cos(OM_Enc) 0;
          0 0 1];

%transformation in the Saturn centred Saturn-Enceladus RF
rSun=R_an_Enc*R_i_Enc*R_OM_Enc*rSun_Sat_F;

%accelerations
%Sun: 
rSun_SC=rSun-r;
a_Sun=muSun*rSun_SC/norm(rSun_SC^3);


%Saturn:
a_Saturn=-muSat*r/norm(r^3);

%Enceladus - spherical harmonics
a_Enc=0;


%acc assembly
dxdt(4:6)=a_Sun+a_Saturn+a_Enc;

end