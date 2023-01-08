function DV = biell_cost(N,N_targ,mu_v,R_v,VU)

mu_S = mu_v(1);
R_Saturn = R_v(1);

r1 = 4*R_Saturn;
r2 = N*R_Saturn;
r3 = N_targ*R_Saturn;

v_E = sqrt(mu_S/r1);



a1 = (r1 + r2)/2;
v1 = sqrt(mu_S*(2/r1 - 1/a1));
v2a = sqrt(mu_S*(2/r2 - 1/a1));


a2 = (r2 + r3)/2;
v2b = sqrt(mu_S*(2/r2 - 1/a2));

DV1 = abs(v1 - v_E)*VU;
DV2 = abs(v2a-v2b)*VU;

DV = DV1+DV2;

end

