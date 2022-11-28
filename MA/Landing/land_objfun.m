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

% Time grid
h = (t_N-t_1)/(N-1);

step_st = length(state_i);           % 7: rr,vv,m
step_var = step_st+4;                % 11: rr,vv,m,u,ax,ay

F = -var(end-6);

%add derivative
if nargout>1
    JF = [];
end

end