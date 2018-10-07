function [dev,acc,maxf,avgf] = diverse(bx,bf)
% diverse.m
% bx is a matrix whose columns are the best solutions of each run
% dev is a measure of deviation between solutions xi and xj
% acc is a measure of accuracy; proportion of solutions close to 'true
% optimum'

[~,r] = size(bx);
M = zeros(r,r);
t = 0.5; % Deviation tolerance for solutions being 'the same'

y = [1:r]';
brun = y(bf == max(bf)); % Find best run

for u = 1:r-1
    for v = (u+1):r % Upper triangular, with leading diagonal = 0
        M(u,v) = sqrt(sum((bx(:,u)-bx(:,v)).^2));
    end
end

A = [M(:,brun); M(brun,:)']; % Deviation between best run and all others
A = A(A~=0); % Remove zeroes
A = A(A<t); % Find solutions that are suitably close

dev = mean(M(M~=0)); % Average standard deviation
acc = (length(A)+1)/r; % Proportion of 'true optimum' solutions
maxf = max(bf);
avgf = mean(bf);