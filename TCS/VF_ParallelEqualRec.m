function F12 = VF_ParallelEqualRec(H,W1,W2)
% VF_ParallelEqualRec computes the view factor between parallel equal
% rectangular plates of size W1, W2, sepatated a distance H

x = W1/H;
y = W2/H;

x1 = sqrt(1+x^2);
y1 = sqrt(1+y^2);

F12 = 1/(pi*x*y)*(log((x1^2*y1^2)/(x1^2+y1^2-1))+2*x*(y1*atan(x/y1)-atan(x))...
    + 2*y*(x1*atan(y/x1)-atan(y)));


end