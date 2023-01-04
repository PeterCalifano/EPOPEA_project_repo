function zero = SolveForVinf(par,Ra,Rp)

vinf = par(1);
alpha = par(2);

i_SC = 0;

SMA_SC = @(v_inf, alpha) 1./(1 - v_inf^2 - 2*v_inf*cosd(alpha));
e_SC = @(v_inf, alpha) sqrt(1 - 1/SMA_SC(v_inf, alpha) * ((3 - 1/SMA_SC(v_inf, alpha) - v_inf^2)/(2*cosd(i_SC)))^2);

e = (Ra-Rp)/(Ra+Rp);
a = (Ra+Rp)/2;

zero = norm([SMA_SC(vinf,alpha) - a;
        e_SC(vinf,alpha) - e]);

end
