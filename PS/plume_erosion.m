clear all, close all, clc

set( 0, 'defaulttextinterpreter', 'latex' )
set( 0, 'defaultlegendinterpreter', 'latex' )
set( groot,'defaultaxesfontSize', 17 )

% ESTIMATE THRUSTER PLUME EROSION AND CONTAMINATION
% MR107T 100 N thruster (hydrazine)

%% Mass rate deposition
m_dry = 746 ; % [ kg ] - Lander dry mass
m_prop = 65 * 2 ; % [ kg ] - Propellant mass coming from MA + 100% margin
m_lander = m_dry + m_prop ; % [ kg ] - Lander total mass (from architecture_sizing.m)
g = 0.1 ; % [ m/s^2 ] - Enceladus gravitational acceleration
T = m_lander * g ; % [ N ] - Thrust needed during landing (assumed to be equal to lander weight)

Isp = 228 ; % [ s ] - Specific impulse of MR107T thruster
g0 = 9.81 ; % [ m/s ] - Gravitational acceleration on Earth at sea level
mdot = T / ( Isp * g0 ) ; % [ kg/s ] - Mass flow raterequired to generate thrust T

alpha = 45 ; % [ deg ] - Half-cone angle representing exhaust plume expansion -> height and ardius of the cone are equal

% v = 0.75 ; % [ m / s ] - Velocity at landing (assumed from literature)

%h = 20 ;

v_profile = [1.94741212617684
1.88197616378036
1.83340990904269
1.78176620263761
1.74333092604841
1.71974067982790
1.71120929767070
1.71076489098799
1.71096377860360
1.71066356625206
1.71024446759916
1.71053391742766
1.71103747365681
1.71156392637115
1.71211326028898
1.71268547984342
1.71328058964721
1.71389859449450
1.71453949937146
1.71520328992923
1.71588997149468
1.71659956910169
1.71733208848240
1.71808753552807
1.71886591629146
1.71966723700921
1.72049150409718
1.72133872413876
1.72220890390257
1.72310205032639
1.72401817051687
1.72495727175604
1.72591936150605
1.72690444738501
1.72791253718951
1.72894363887958
1.72999776058462
1.73107491059831
1.73217509737463
1.73329832954121
1.73444461587116
1.73561396530814
1.73680638695202
1.73802189005517
1.73926048402370
1.74052217841988
1.74180698295296
1.74311490747844
1.74444596207264
1.74580015719451
1.74717753435726
1.74857810361344
1.75000184420005
1.75144876679726
1.75291888222001
1.75441220140775
1.75592873542547
1.75746849545555
1.75903149279591
1.76061773885712
1.76222724515907
1.76386002332660
1.76551608508372
1.76719544224983
1.76889810673716
1.77062409054353
1.77237340574903
1.77414606451037
1.77594207905466
1.77776146167577
1.77960422472712
1.78147038061548
1.78335994179634
1.78527292076663
1.78720933005659
1.78916918222579
1.79115248985399
1.79315926553367
1.79518952185938
1.79724327142487
1.79932052682824
1.80142130061956
1.80354560532903
1.80569345347890
1.80786485749715
1.81005982978873
1.81227838266961
1.81452052840053
1.81678627912750
1.81907564691945
1.82138864372530
1.82372528134717
1.82608557148713
1.82846952568179
1.83087715529760
1.83330847153495
1.83576348538392
1.83824220765715
1.84074464893469
1.84327081954126
1.84582072958001
1.84839438886521
1.85099180691515
1.85361299297448
1.85625795594177
1.85892670435727
1.86161924645647
1.86433559002998
1.86707574250137
1.86983971085735
1.87262750165943
1.87543912096476
1.87827457436758
1.88113386694905
1.88401700322351
1.88692398715833
1.88985482212123
1.89280951086149
1.89578805549784
1.89879045744412
1.90181671741445
1.90486683542022
1.90794081065884
1.91103864155966
1.91416032571068
1.91730585982552
1.92047523973640
1.92366846030924
1.92688551548545
1.93012639813181
1.93339110011342
1.93667961217980
1.93999192396508
1.94332804271262
1.94668794985100
1.95007161825970
1.95347903830919
1.95691018757483
1.96036504794715
1.96384359987221
1.96734582229918
1.97087169261030
1.97442118657491
1.97799427831465
1.98159094022812
1.98521114290468
1.98885485513723
1.99252202792322
1.99621264211361
1.99992667644274
2.00366407580799
2.00742479884435
2.01120877338510
2.01501595203050
2.01884631516716
2.02269981379713
2.02657639472042
2.03047600212989
2.03439857754403
2.03834405973958
2.04231238464220
2.04630348525144
2.05031729156789
2.05435371207050
2.05841268891645
2.06249416085753
2.06659802646660
2.07072419924518
2.07487258910197
2.07904309230940
2.08323562156645
2.08745006705947
2.09168631322556
2.09594426908928
2.10022382088829
2.10452485054329
2.10873086405365
2.11295827503972
2.11725761865199
2.12094358614758
2.11524015382024
2.09197555611576
2.05118480344737
1.99242215826399
1.91474302070104
1.81678102488119
1.69864176386125
1.57083452454181
1.44256160135077
1.31353088646239
1.18392918027191
1.05387355238637] * 0.0697 * 1e3 ;

% 0.923272156056618
% 0.792041452489496
% 0.660903497177164
% 0.530269428301020
% 0.398741931445934
% 0.265783359486333
% 0.132475554893756
% 9.23679803873558e-11

h_vect = linspace( 2, 30, length(v_profile) ) ; % [ m ] - Height from ground at which we turn off the thrusters

t_landing = linspace( 0, 2000, length(v_profile) ) ; 

m_fluence = zeros( length( v_profile ), 1 ) ; % Initialize vector

%% Plot if we consider mean landing velocity and different cut-off altituted
for k = 1 : length( h_vect )
    h = h_vect(k);
    %v = v_profile( k ) ;
    v=mean(v_profile) ;
    m_fluence(k) = T / ( Isp * g0 * pi * h * v ) ; % [ kg / m^2 ] - Derivative of mass flux within the cone over time

end

m_fluence = m_fluence * 1e6 ; % [ mg/ m ^ 2 ]

scale = 0.1 ; % Scaling factor (account for percentage of ammonia in thruster plume) - from literature (Phoenix lander)

% Plot the mass fluence as function of cutoff height
figure()
plot( h_vect, m_fluence, 'Linewidth', 2 )
hold, grid on
plot( h_vect, m_fluence * scale, 'Linewidth', 2 )
xlabel( 'Cutoff height [m]' ), ylabel( 'Mass fluence [$mg/m^2$]' )
title( 'Exhaust plume deposition' )
legend( 'Exhaust plume deposition', '$NH_3$ deposition', 'location', 'best' )


%% Plot if we consider the actual variation of the velocity during landing phase and fixed cut-off altitude
h = 20 ;
for k = 1 : length( v_profile )
    v = v_profile( k ) ;
    h1 = 50 ;
    h2 = 30;
    h3 = 15;
    h4 = 2;
    m_fluence1(k) = T / ( Isp * g0 * pi * h1 * v ) ; % [ kg / m^2 ] - Derivative of mass flux within the cone over time
    m_fluence2(k) = T / ( Isp * g0 * pi * h2 * v ) ; % [ kg / m^2 ] - Derivative of mass flux within the cone over time
    m_fluence3(k) = T / ( Isp * g0 * pi * h3 * v ) ; % [ kg / m^2 ] - Derivative of mass flux within the cone over time
    m_fluence4(k) = T / ( Isp * g0 * pi * h4 * v ) ; % [ kg / m^2 ] - Derivative of mass flux within the cone over time

end

m_fluence1 = m_fluence1 * 1e6 ; % [ mg/ m ^ 2 ]
m_fluence2 = m_fluence2 * 1e6 ; % [ mg/ m ^ 2 ]
m_fluence3 = m_fluence3 * 1e6 ; % [ mg/ m ^ 2 ]
m_fluence4 = m_fluence4 * 1e6 ; % [ mg/ m ^ 2 ]

scale = 0.1 ; % Scaling factor (account for percentage of ammonia in thruster plume) - from literature (Phoenix lander)

% Plot the mass fluence as function of cutoff height (plot brutto)
figure()
plot( v_profile, m_fluence1 * scale, 'o', 'Color',  [0 0.4470 0.7410],  'Linewidth', 1 )
hold, grid on
plot( v_profile, m_fluence2 * scale, 'square', 'Color',[0.4660 0.6740 0.1880],'Linewidth', 1 )
plot( v_profile, m_fluence3 * scale, '*', 'Color', [0.9290 0.6940 0.1250],'Linewidth', 1 )
plot( v_profile, m_fluence4 * scale, '^', 'Color', [0.6350 0.0780 0.1840], 'Linewidth', 1 )
xlabel( 'Descent velocity [m/s]' ), ylabel( 'Mass fluence [$mg/m^2$]' )
title( '$NH_3$ deposition' )
legend( '$h=50 \ m$', '$h=30 \ m$', '$h=15 \ m$', '$h=2 \ m$', 'location', 'best' )

figure()
scatter(v_profile, m_fluence1 * scale, [], t_landing, 'Linewidth', 1.2 )
hold on, grid on
s = scatter(v_profile, m_fluence2 * scale, [], t_landing, 'Linewidth', 1.2 )
s.Marker = 'square' ;
s = scatter(v_profile, m_fluence3 * scale, [], t_landing, 'Linewidth', 1.2 )
s.Marker = '^' ;
s = scatter(v_profile, m_fluence4 * scale, [], t_landing, 'Linewidth', 1.2 )
s.Marker = 'x' ;
colorbar
colormap jet
xlabel( 'Descent velocity [m/s]' ), ylabel( 'Mass fluence [$mg/m^2$]' )
title( '$NH_3$ deposition' )
legend( '$h=50 \ m$', '$h=30 \ m$', '$h=15 \ m$', '$h=2 \ m$', 'location', 'best' )