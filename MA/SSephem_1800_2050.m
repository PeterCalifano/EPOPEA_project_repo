% -------------------------------------------------------------------------
% [ X_ecl , V_ecl , KP ] = SSephem_1800_2050( planet , JD )
% -------------------------------------------------------------------------
%
% © Jacopo Prinetto - Aerospace Science and Technology Dept. - PoliMi - 
% Version:      1.0
% Status:       VALIDATED
% Last Update:  11/01/2019
%
% -------------------------------DESCRIPTION-------------------------------
%
% In this function Ephemeris for major objects of Solar System are
% implemented following E.M.Standish document available on NASA-JPL Solar
% System Dynamics group.
% The ephemeris are acceptable within the error declared in the document
% between 1800 AD and 2050 AD
% The reference is the mean ecliptic and equinox JD2000.
%
% ---------------------------------INPUTS----------------------------------  
%
%   planet = ID of the planet (1 for Mercury, 2 for Venus...)
%   JD     = Julian Data
%
% ---------------------------------OUTPUTS---------------------------------
%
%   X_ecl = (3x1) x_y_z position in the ecliptic RF [AU]                   
%   V_ecl = (3x1) x_y_z velocity in the ecliptic RF [AU/TU] with TU 
%           such that sun parameter is equal to 1
%   KP    = (1x6) Keplerian parameters [AU, - , rad , rad , rad , rad]
%
% ----------------------------------NOTES----------------------------------
%
%   see: https://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf (Sept. 2017)
%

function [ X_ecl , V_ecl , KP ] = SSephem_1800_2050( planet , JD)
%% Select planetary data
switch planet
    
    case 1
        KE = [0.38709927      0.20563593      7.00497902      252.25032350     77.45779628     48.33076593];
        KR = [0.00000037      0.00001906     -0.00594749   149472.67411175      0.16047689     -0.12534081];
    case 2
        KE = [0.72333566      0.00677672      3.39467605      181.97909950    131.60246718     76.67984255];
        KR = [0.00000390     -0.00004107     -0.00078890    58517.81538729      0.00268329     -0.27769418];
    case 3
        KE = [1.00000261      0.01671123     -0.00001531      100.46457166    102.93768193      0.0];
        KR = [0.00000562     -0.00004392     -0.01294668    35999.37244981      0.32327364      0.0];
    case 4
        KE = [1.52371034      0.09339410      1.84969142       -4.55343205    -23.94362959     49.55953891];
        KR = [0.00001847      0.00007882     -0.00813131    19140.30268499      0.44441088     -0.29257343];
    case 5
        KE = [5.20288700      0.04838624      1.30439695       34.39644051     14.72847983    100.47390909];
        KR = [-0.00011607     -0.00013253     -0.00183714     3034.74612775      0.21252668      0.20469106];
    case 6
        KE = [ 9.53667594      0.05386179      2.48599187       49.95424423     92.59887831    113.66242448];
        KR = [-0.00125060     -0.00050991      0.00193609     1222.49362201     -0.41897216     -0.28867794];
    case 7
        KE = [19.18916464      0.04725744      0.77263783      313.23810451    170.95427630     74.01692503];
        KR = [-0.00196176     -0.00004397     -0.00242939      428.48202785      0.40805281      0.04240589];
    case 8
        KE = [30.06992276      0.00859048      1.77004347      -55.12002969     44.96476227    131.78422574];
        KR = [0.00026291      0.00005105      0.00035372      218.45945325     -0.32241464     -0.00508664];
    case 9
        KE = [39.48211675      0.24882730     17.14001206      238.92903833    224.06891629    110.30393684];
        KR = [-0.00031596      0.00005170      0.00004818      145.20780515     -0.04062942     -0.01183482];
        
end
%% Computing Parameters

% initializing keplerian parameters & positions
KP = zeros(1,6);
X_perifocal = zeros(3,1);

%X_eq = X_perifocal;
%X_ecl= X_perifocal;

%V_perifocal = X_perifocal;
%V_rt = X_perifocal;
%V_ecl = X_perifocal;

%converting from deg to rad
for i = 1:4
    KE(i+2) = KE(i+2)*pi/180;
    KR(i+2) = KR(i+2)*pi/180;
end

%computing century
T = ( JD - 2451545 )/36525 ;

% compute elements at current epoch
for i = 1:6
    
    KE(i) = KE(i) + KR(i)*T;
    
end

% computing keplerian param.
KP(1) = KE(1);
KP(2) = KE(2);
KP(3) = KE(3);
KP(4) = KE(6);

% computing argomentuum of perihelion
KP(5) = KE(5) - KE(6) ;

% computing Mean Anomaly
M = KE(4) - KE(5) ;

% solving kepler's eq for Eccentric anomaly
e = KP(2);

E = M + e*sin(M);

toll = 1e-6;
DE = toll + 1;

while DE>toll 
    
    DM = M - (E-e*sin(E));
    DE = DM/(1-e*cos(E));
    E = E + DE ;
end

% computing true anomaly
KP(6) = 2*atan(( (1+e)/(1-e) )^0.5 * tan(E/2));

% computing perifocal coordinates
X_perifocal(1) = KP(1)*(cos(E)-e) ;
X_perifocal(2) = KP(1)*(1-e^2)^0.5*sin(E);
X_perifocal(3) = 0;

%% computing ecliptic coordinates (J2000 frame)
R = [ cos(KP(5))*cos(KP(4))-sin(KP(5))*sin(KP(4))*cos(KP(3)) , -sin(KP(5))*cos(KP(4))-cos(KP(5))*sin(KP(4))*cos(KP(3)) ,  13 ;
      cos(KP(5))*sin(KP(4))+sin(KP(5))*cos(KP(4))*cos(KP(3)) , -sin(KP(5))*sin(KP(4))+cos(KP(5))*cos(KP(4))*cos(KP(3)) ,  23 ;
      sin(KP(5))*sin(KP(3))                                  ,  cos(KP(5))*sin(KP(3))                                  ,  33];


X_ecl = R*X_perifocal;


V_rt = (1/(KP(1)*(1-KP(2)^2)))^0.5 * [KP(2)*sin(KP(6)) ; 1+KP(2)*cos(KP(6)) ; 0];
V_perifocal = [cos(KP(6)),sin(KP(6)),0 ;-sin(KP(6)) , cos(KP(6)) , 0 ; 0 , 0 ,1]' * V_rt;

V_ecl = R*V_perifocal;




end

