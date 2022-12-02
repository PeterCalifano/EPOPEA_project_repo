function  [F, JF] = land_objfun(var, state_i, ~, N)
%
%
%
%
%
%----------------------------------------------------------------------------

%Initial and final time
t_1 = var(end-1);
t_N = var(end);

F = -var(end-6);

%add derivative
if nargout>1
    JF = [];
end

end