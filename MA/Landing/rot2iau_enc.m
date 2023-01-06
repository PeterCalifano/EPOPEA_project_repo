function state_in = rot2iau_enc (t, state, mu)

xr = state(1);
yr = state(2);
zr = state(3);
vxr = state(4);
vyr = state(5);
vzr = state(6);

X = (xr+mu-1)*cos(t) - yr*sin(t);
Y = (xr+mu-1)*sin(t) + yr*cos(t);
Z = zr;
VX = (vxr-yr)*cos(t) - (vyr+xr+mu-1)*sin(t);
VY = (vxr-yr)*sin(t) + (vyr+xr+mu-1)*cos(t);
VZ = vzr;

state_in = [X; Y; Z; VX; VY; VZ];
end