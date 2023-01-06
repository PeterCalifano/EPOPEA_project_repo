function T_FB = time_of_FB(v_inf,rp,id_planet,VU,DU)
%% DESCRIPTION
% Computes the time spent to perform an Unpowered Gravity Assist inside the
% SOI of the planet.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% v_inf - Incoming infinite velocity [adimensional]
% rp - Pericenter radius [adimensional]
% id_planet - Identifier of the FB planet 
% VU - Adimensionalization constant for velocity
% DU - Adimensionalization constant for position
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% T_FB - Time to perform the FB [s]
% -------------------------------------------------------------------------------------------------------------


rp = rp*DU;
v_inf_m = v_inf*VU;
v_inf_p = v_inf_m;
G = astroConstants(1);
AU = astroConstants(2);
mu_p = astroConstants(10+id_planet);
Mass_p = mu_p/G;
Mass_S = astroConstants(4)/G;

    r_SOI = AU * (Mass_p/Mass_S)^(2/5);
    
    a_m = - mu_p/norm(v_inf_m)^2;
    e_m = 1 + rp*norm(v_inf_m)^2 / mu_p;
    a_p = - mu_p/norm(v_inf_p)^2;
    e_p = 1 + rp*norm(v_inf_p)^2 / mu_p;
    delta_m = 2*asin(1/e_m);
    delta_p = 2*asin(1/e_p);
    aim_p = -a_p*e_p*cos(delta_p/2);
    aim_m = -a_m*e_m*cos(delta_m/2);
    p_p = a_p * (1-e_p^2);
    p_m = a_m * (1-e_m^2);

    th_inf_m = acos(1/e_m * (p_m/r_SOI - 1));
    th_inf_p = acos(1/e_p * (p_p/r_SOI - 1));
    

    E_inf_m = 2*atanh(sqrt((e_m-1)/(e_m+1)) * tan(th_inf_m/2));
    E_inf_p = 2*atanh(sqrt((e_p-1)/(e_p+1)) * tan(th_inf_p/2));

    M_p = e_p*sinh(E_inf_p) - E_inf_p ;
    M_m = e_m*sinh(E_inf_m) - E_inf_m;
    
    t_p = M_p * sqrt(-a_p^3 / mu_p);
    t_m = M_m * sqrt(-a_m^3 / mu_p);
    
    T_FB = t_p + t_m;

end