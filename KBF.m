function f = KBF(x)
% KBF.m
% Evaluate objective function: Keane's Bump Function
% Accepts matrix where each column is a solution vector x
% Returns row vector containing f(x) in each column

y = [1:size(x,1)];
f = zeros(1,size(x,2));

for i = 1:size(x,2)
    xt = x(:,i);
    if feasible(xt)
        f(1,i) = (sum(cos(xt).^4)-2*prod(cos(xt).^2))/sqrt(y*(xt.^2));
    else
        f(1,i) = -1; % If infeasible x gets passed to KBF(x)
    end
end