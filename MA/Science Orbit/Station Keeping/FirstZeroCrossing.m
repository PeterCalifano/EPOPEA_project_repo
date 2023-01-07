function [value,isterminal,direction] = FirstZeroCrossing(t,x,mu_tbp,mu_v,R_v,J2_v) %Event function 
%to stop the integration at the intersection point with the x-z plane 

value=x(2);
isterminal=1;
direction=0;

end