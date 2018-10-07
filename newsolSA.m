function y = newsolSA(x,g,D)
% newsolSA.m
% Generate new SA solution 
% Current solution x, Generation index g

% g=1: x_{i+1} = x_i + Cu (Crude method)
% g=2: x_{i+1} - x_i + Du (Parks)

u = 2*(rand(size(x))-0.5);

if g == 1
    C = 2*eye(length(x));
    y = x + C*u;
elseif g == 2
    y = x + D*u;
end

% If new solution is infeasible, try again
if ~ feasible(y)
    y = newsolSA(x,g,D);
end