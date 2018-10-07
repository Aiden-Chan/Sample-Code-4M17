function c = feasible(x)
% feasible.m
% Checks to see if vectors in matrix is a feasible solution.
% Returns true/false in a row vector
c = (max(x)<=10).*(min(x)>=0).*(prod(x)>0.75).*(sum(x)<(7.5*size(x,1)));


