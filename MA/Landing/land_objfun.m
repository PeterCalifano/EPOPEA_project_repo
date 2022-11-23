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
u_k = zeros(N,1);
% for k = 1:(N-1)
% %    u_k(k) = var((step_st+1)+step_var*(k-1));
% %    u_next = var((step_st+1)+step_var*(k));
% %    F = F+(u_k(k)+u_next)*h/2;            % trapezoidal. change in Gauss!
% end

F = -var(end-5);

%add derivative
if nargout>1
JF = zeros(step_var*N+2,1);
    for k = 6:step_var:(step_var*N+2)-4
        if k == 6 || k == (step_var*N+2)-4
            JF(k) = h/2;
        else
            JF(k) = h;
        end
    end

    for k =  1:N
        if k == 1 || k == N
            JF(end-1) = JF(end-1) - u_k(k);
            JF(end) = JF(end) + u_k(k);
        else
            JF(end-1) = JF(end-1) - 2*u_k(k);
            JF(end) = JF(end) + 2*u_k(k);
        end
    end
JF(end-1) = JF(end-1)/(2*(N-1));
JF(end) = JF(end)/(2*(N-1));
end

end