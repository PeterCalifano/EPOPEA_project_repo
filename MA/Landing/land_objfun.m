function  [F, JF] = ObjFun(var,~)

t_1 = 0; % ???????
t_N = var(end);
% Unpack the control

N = floor(length(var)/8);
F = 0;
h = (t_N-t_1)/(N-1);
for k = 1:(N-1)
   u_k = var(6+8*(k-1));
   u_next = var(6+8*(k));
   F = F+(u_k+u_next)*h/2;                                                 % trapezoidal. change in Gauss!
end
JF = [];
% add derivative
%     if nargout>1
%        
% 
% 
%     end

end