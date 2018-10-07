function x = SD(L,n,LTM,regionmap)
% SD.m
% Search Diversification: Find a new x in a previously unexplored region

unexplored = [0:L^n-1]';
unexplored = unexplored(LTM); % List of unexplored regions
regionindex = unexplored(1+floor(rand(1)*length(unexplored))); % Choose random from list
regcount = 1; % To avoid getting stuck within an infeasible region
x = (10/L)*(floor(mod(regionindex./regionmap,L))+rand(1,n))'; % Convert to region

while (~ feasible(x)) && regcount < 20
    regcount = regcount + 1;
    x = (10/L)*(floor(mod(regionindex./regionmap,L))+rand(1,n))';
end

if ~feasible(x)
%     All 20 attempts at choosing an x in that region were infeasible
    x = SD(L,n,LTM,regionmap);
end