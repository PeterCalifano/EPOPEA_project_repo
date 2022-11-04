function [DV,rp] = Gravity_Assist(planet,Vm,Vp,time)
%
% Description:
%   The function computes the powered gravity-assist manoeuvre around planet to link the 
% incoming and outcoming heliocentric velocities Vm and Vp.
%   If the fly-by is not feasible, DV == NaN.
%
% Prototype:
% [DV,rp] = Gravity_Assist(i_FB,kep_t1,kep_t2,mjd2000,flag)
%
% Inputs:
% planet [1] - Planet around which the FB is performed
% Vm [3x1] - Heliocentric velocity before the FB
% Vp [3x1] - Heliocentric Velocity after the FB
% time [1] - Julian Day of the Fly-by
%
% Outputs:
% DV [1] - Delta velocity to give at the pericenter to perform the powered FB [km/s] 
% rp [1] - radius of pericenter of the hyperbola [Km]
%
% -------------------------------------------------------------------------

%% Find velocities in heliocentric reference frame

    mu_S = cspice_bodvrd('Sun','GM',1);
    mu_p = cspice_bodvrd(planet,'GM',1);
    R_p = cspice_bodvrd(planet,'RADII',3);
    R_p = R_p(1);
    G = astroConstants(1);
    Mass_p = mu_p/G;
    Mass_S = mu_S/G;
    AU = astroConstants(2);
    
    x = cspice_spkezr(planet,time,'ECLIPJ2000','NONE','SUN');
    r_GA = x(1:3);
    v_GA = x(4:6);
%     [~,Vm] = kep2car(kep_t1, mu_S);
%     [~,Vp] = kep2car(kep_t2, mu_S);
    
%% Solve the hyperbola
    
    v_inf_m = Vm - v_GA;
    v_inf_p = Vp - v_GA;
    DV_nat = norm(v_inf_p - v_inf_m);

    delta = acos(dot(v_inf_m,v_inf_p)/(norm(v_inf_p)*norm(v_inf_m)));
    
    equation = @(r) - asin(1/(1 + (r * norm(v_inf_m)^2 /mu_p)))...
        - asin(1/(1 + (r * norm(v_inf_p) ^ 2 /mu_p))) + delta;
    
    options = optimset('TolFun', 1e-13,'Display','off');
    
    [rp,~,Exitflag] = fsolve(equation,1.5*(R_p+500),options);
    
    
    vp_m = sqrt(norm(v_inf_m)^2 + 2*mu_p/rp);
    vp_p = sqrt(norm(v_inf_p)^2 + 2*mu_p/rp);
    
    DV = norm(vp_p - vp_m);
    

%% If something went wrong in the solution of the non-linear function, set the DV to NaN

    if rp < R_p 
        disp('Error: fly-by not possible because of impact')
        DV = NaN;
    end

    if Exitflag < 0
        disp('Error in the solution of non-linear function')
        DV = NaN;
    end
end