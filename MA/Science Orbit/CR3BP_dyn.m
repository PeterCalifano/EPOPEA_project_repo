function dxdt = CR3BP_dyn(t,x,mu)
%INPUT - è IL CR3BP LA SAPETE
%OUTPUT dxdt

r1=@(t,x) ((x(1)+mu).^2+x(2).^2+x(3).^2).^0.5;
r2=@(t,x) ((x(1)+mu-1).^2+x(2).^2+x(3).^2).^0.5;

dxdt= [x(4);
    x(5);
    x(6);
    2*x(5)+x(1)-(1-mu)*(x(1)+mu)/(r1(t,x)).^3 - mu/(r2(t,x)).^3*(x(1)+mu-1);
    -2*x(4)+x(2)-(1-mu)*x(2)/r1(t,x).^3-mu/(r2(t,x)).^3*x(2);
    -(1-mu)/(r1(t,x)).^3*x(3)-mu/(r2(t,x).^3)*x(3)];  %CRTBP dynamics



end

