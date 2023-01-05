function [value,isterminal,direction] = LastCrossing(t,x) %Event function 
%to stop the integration at the intersection point with the x-z plane 

value=x(2);
isterminal=1;
direction=0;

end