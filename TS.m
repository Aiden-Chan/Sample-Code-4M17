% TS.m
% Tabu Search

n = 2; % Dimension of problem
iter = 10000; % Maximum number of iterations per run; should not be limiting
r = 100; % number of runs
t = 10; % number of tests

d0 = 1.5; % Initial step size
alpha = 0.8; % Step size reduction factor
dmin = d0*alpha^9; % Step size conversion threshold
evalmax = 5000; % Maximum number of objective function evaluations
N = 7; % Entries in short term memory
M = 4; % Entries in medium term memory
L = 3; % Divisor for LTM regions
cI = 10; % Threshold to intensify
cD = 15; % Threshold to diversify
cR = 25; % Threshold to reduce

regionmap = L.^[0:n-1];
% feasmap finds regions which are definitely feasible by checking midpoint
feasmap = (feasible((10/L)*(floor(mod([0:L^n-1]'./regionmap,L))+0.5)'))';

bx = zeros(n,r); % Store optimum x vectors
bf = zeros(r,1); % Store optimum f
rt = zeros(r,1); % Store runtimes

results = zeros(6,t); % param, dev, acc, maxf, avgf, runtime

for test = 1:t
% %     Change test parameters here
    L = test;
    regionmap = L.^[0:n-1];
    feasmap = (feasible((10/L)*(floor(mod([0:L^n-1]'./regionmap,L))+0.5)'))';

%     d0 = 0.25 + (test-1)*0.25;
%     dmin = d0*alpha^9;
   
%     alpha = 0.9 - (test-1)*0.1;
%     dmin = d0*alpha^9;
    
    param = L; % e.g. param = L if testing

    for run = 1:r
        STM = zeros(n,N); % Last successfully visited locations
        MTMx = zeros(n,M); % Best solutions located thus far
        MTMf = zeros(1,M); % Corresponding best maximums located thus far
        LTM = false(L^n,1); % Record explored regions
        d = d0; % Reset initial step size

        trackx = zeros(n,iter);

        x = 10*rand(n,1); % Generate first solution and ensure feasibility
        while ~ feasible(x)
            x = 10*rand(n,1);
        end

        tic
        currentf = KBF(x);
        STM(:,1) = x;
        MTMx(:,1) = x;
        MTMf(1) = currentf;
        eval = 1; % Reset number of evaluations of KBF

        x0 = x;
        f0 = currentf;
        c = 0; % Counter

        for i = 1:iter
            trackx(:,i) = x;
            %     Update STM
            STM(:,mod(i,N)+1) = x;

            % Test convergence
            if d <= dmin || eval >= evalmax
                break
            else
                step = d*[eye(n) -eye(n)];
                trialx = x*ones(1,2*n) + step; % Take steps

            %     Check STM; zero out columns that have already been visited
                STMfilter = true(1,2*n);
                for m = 1:N
                    STMfilter = max(trialx~=STM(:,m)).*STMfilter; % 0 if visited before, 1 otherwise
                end
                trialx = trialx.*STMfilter; % Remove locations already visited

                trialf = KBF(trialx);
                eval = eval + sum(trialf ~= -1);
                if max(trialf) == -1
            %         All steps are tabu; intensify
                    x = sum(MTMx,2)/M;
                    currentf = KBF(x);
                else
                    filter = (trialf == max(trialf));
                    trialx = trialx(logical(true(n,1).*filter)); % Isolate best move
                    trialf = trialf(filter); % Isolate new f

                    if trialf > currentf
                        patternx = x + 2*step(logical(true(n,1).*filter));
                        patternf = KBF(patternx);
                        if patternf > trialf
                            x = patternx;
                            currentf = patternf;
                        else
                            x = trialx;
                            currentf = trialf;
                        end
                    else
                        x = trialx;
                        currentf = trialf;
                    end

                %     Update LTM
                    regionx = regionmap*floor(x./(10/L));
                    LTM(regionx+1) = true; % Region has been visited

                    if currentf > max(MTMf) % New optimum found
                        c = 0;
                    else
                        c = c+1;
                        if c == cI
                    %       Intensify
                            x = sum(MTMx,2)/M;
                            currentf = KBF(x);

                        elseif c == cD
                    %       Diversify
                    %       Check whether all feasible regions have already been searched:
                            if max(feasmap-LTM) > 0
                                x = SD(L,n,LTM,regionmap); % Diversify                            
                            else
        %                       All feasible regions have been visited at least once; choose random  
                                x = 10*rand(n,1);
                                while ~ feasible(x)
                                    x = 10*rand(n,1);
                                end
                            end
                            currentf = KBF(x);

                        elseif c == cR
                    %         Reduce step size
                            d = alpha*d;
                            c = 0;
                        end
                    end

            %     Update MTM
                    if currentf > min(MTMf) && min(max(abs(MTMx-x))) > 0.0001
            %           Ensure no repeats within MTM
                %       Find MTM column to be overwritten; avoids multiple minima
                        for col = 1:M
                            if MTMf(col) == min(MTMf)
                                break
                            end
                        end
                        MTMf(col) = currentf;
                        MTMx(:,col) = x;
                    end
                end
            end
        end

        bfilter = (MTMf == max(MTMf));
        bx(:,run) = MTMx(logical(true(n,1).*bfilter)); % Isolate optimum
        bf(run) = MTMf(bfilter);
        rt(run) = toc;
    end

    if r > 1
        [dev,acc,maxf,avgf] = diverse(bx,bf);
        results(1,test) = param;
        results(2,test) = dev;
        results(3,test) = acc;
        results(4,test) = maxf;
        results(5,test) = avgf;
        results(6,test) = mean(rt);
    end
end

results

if n == 2
    PlotKBF2D
    hold on
    if r == 1
        scatter3(trackx(1,:),trackx(2,:),KBF(trackx),20,'b','filled')
        scatter3(x0(1),x0(2),f0,20,'w','filled')
        scatter3(MTMx(1,:),MTMx(2,:),MTMf,20,'r','filled')
        legend('KBF','Track','Start','Best')
    else
        scatter3(bx(1,:),bx(2,:),bf,20,'k','filled')
        legend('KBF','Best')
    end
end