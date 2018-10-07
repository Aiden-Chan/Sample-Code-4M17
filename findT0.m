function T0 = findT0(x,startf)
% findT0.m
% Conducts initial search around starting point x to find average df

chi = 0.8; % Kirkpatrick probability of initial acceptance
s = 10; % Number of initial search points

y = x*ones(1,s)+(2*rand(length(x),s)-1); % Random number between -1 and 1
df = startf - KBF(y);
df = df(df>0); % Only valid for startf > KBF(y)

if ~ isempty(df)
    avgdf = mean(df);
else
    avgdf = 0.05;
end
    
T0 = -avgdf/log(chi);