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

step_st = length(state_i);           % 5: rr,vv,m
step_var = 5+3;                      % 8: rr,vv,m,u,ax,ay

F = 0;
for k = 1:(N-1)
   u_k = var((step_st+1)+step_var*(k-1));
   u_next = var((step_st+1)+step_var*(k));
   F = F+(u_k+u_next)*h/2;            % trapezoidal. change in Gauss!
end


% add derivative
JF = zeros(step_var*N+2,1);
if nargout>1
    for k = 6:step_var:(step_var*N+2)-4
        if k == 6 || k == (step_var*N+2)-4
            JF(k) = h/2;
        else
            JF(k) = h;
        end
    end
% u1 = var(6);
% u2 = var(14);
% uN_1 = var(end-12);
% uN = var(end-4);
% JF(end-1) = -(u1+u2)/2; %%%%%%%%%%%%%%%%%
% JF(end) = (uN_1+uN)/2;  %%%%%%%%%%%%%%%%%
end

end