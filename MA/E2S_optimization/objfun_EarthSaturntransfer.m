function DV = objfun_EarthSaturntransfer(var,N,planets,Ra_target,Rp_target)
%%%% INPUTS %%%%%%
% N - Number of fly-bys
% planets - String of planets visited (dep and arrival included)
% var - variables vector, with t1 first and then the ToF
% Ra_target - SEE EstimateDVtoCapture.m
% Rp_target - SEE EstimateDVtoCapture.m
t1 = var(1);
tof = zeros(N+1,1);
times = zeros(N+2,1);
times(1) = t1;
for i = 2:N+2
    tof(i-1) = var(i);
    times(i) = times(i-1) + tof(i-1);
end

%% Interplanetary Legs
v_interp = zeros(N+1,2);
for i = 1:N+1 
    [v1_t,v2_t] = Lambert(planets(i),planets(i+1),times(i),times(i+1));
    v_interp(i,1) = v1_t;
    v_interp(i,2) = v2_t;
end
%% Fly-bys
DV_fb = zeros(N,1);
for i = 1:N
    [DV_i,rp_i] = Gravity_Assist(planet,v_interp(i,2),v_interp(i+1,1),times(i+1));
    DV_fb(i) = DV_i;
end
%% Compute DV
X_last = cspice_spkezr(planets(end),times(end,'ECLIPJ2000','NONE','SUN');
V_last = X_last(4:6);
mu_main = cspice_bodvrd(planets(end),'GM',1);
Vinf_entry = v_interp(end,2) - V_last;
[DV_capture, T_capture] = EstimateDVtoCapture(Vinf_entry, mu_main, Ra_target, Rp_target);

DV = DV_capture + sum(DV_fb);



end



