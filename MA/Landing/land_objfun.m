function  [F, JF] = land_objfun(var,state_i, par, N)
%
%
%
%
%
%----------------------------------------------------------------------------

%Initial and final time
t_1 = var(end-1);
t_N = var(end);

% Time grid
h = (t_N-t_1)/(N-1);

F = 0;
for k = 1:(N-1)
   u_k = var(6+8*(k-1));
   u_next = var(6+8*(k));
   F = F+(u_k+u_next)*h/2;                 % trapezoidal. change in Gauss!
end


% add derivative
if nargout>1
   JF = [];
end

end