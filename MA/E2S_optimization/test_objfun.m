
t1 = 475848/TU;
var = [t1, t1+30*3600*24/TU, t1 + 70*3600*24/TU, DU, pi/8, 1000, t1 + 10*3600*24/TU,...
     DU, pi/4, 900, t1 + 50*3600*24/TU];
DV = objfun_EarthSaturntransfer(var,1,[3,2,3],20*R_Saturn,200*R_Saturn,TU);
