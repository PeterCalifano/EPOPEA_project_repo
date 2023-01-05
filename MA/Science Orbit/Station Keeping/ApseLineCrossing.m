function [value,isterminal,direction] = ApseLineCrossing(t,x) %Event function 
%to stop the integration at the intersection point with the x-z plane 

value=x(2);
isterminal=0;
direction=0;

end