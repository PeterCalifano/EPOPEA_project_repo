clear
close all
clc



nlpvar = [468483,3,0.5,0.5,67777,49939,0.3,0.5,7000,pi/2];
planets_id = [3,3,4];

DV = objfun_EarthSaturntransfer_2(nlpvar, planets_id, 20*6378, 30*6378);