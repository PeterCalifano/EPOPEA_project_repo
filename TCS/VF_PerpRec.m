function F12 = VF_PerpRec(H,W,L)
% VF_PerpRec computes the view factor between perpendicular
% rectangular plates: from horizontal rectangle of WxL to adjacent vertical
% rectangle of HxL

h = H/L;
w = W/L;
a = (1+h^2)*(1+w^2)/(1+h^2+w^2);
b = w^2*(1+h^2+w^2)/((1+w^2)*(h^2+w^2));
c = h^2*(1+h^2+w^2)/((1+h^2)*(h^2+w^2));


F12 = 1/(pi*w)*(h*atan(1/h)+w*atan(1/w)-sqrt(h^2+w^2)*atan(1/sqrt(h^2+w^2))+...
    1/4*log(a*b^(w^2)*c^(h^2)));

end