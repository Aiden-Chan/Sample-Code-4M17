% SA.m
% Simulated Annealing

n = 10; % Dimension of problem
g = 2; % Generation method
r = 1000; % number of runs
t = 1; % number of tests

a = 0.95; % Temperature decrement factor
Lk = 50; % minimum number of trials at each temperature
eta = ceil(0.6*Lk); % minimum number of accepted trials at each temperature
cres = 4000; % Number of iterations without new optimum for restart
iter = 5000; % number of iterations per run

alpha = 0.1; % Parks damping constant
omega = 2.1; % Parks weighting of R into D

bx = zeros(n,r); % Store optimum x vectors
bf = zeros(r,1); % Store optimum f
rt = zeros(r,1); % Store runtimes
trackx = zeros(n,iter); % Store path

results = zeros(6,t); % param, dev, acc, maxf, avgf, runtime

for test = 1:t
% %     Change testing parameters here
%     a = 1 - (test-1)*0.1;
%     Lk = 100 - (test-1)*10;
%     eta = ceil(0.6*Lk);
%     cres = 5000 - (test-1)*500;
    param = 0; % e.g. param = a if testing
    
    for run = 1:r

        % Generate random initial vector and ensure that it is feasible
        x = 10*rand(n,1);
        while ~ feasible(x)
            x = 10*rand(n,1);
        end

        tic
        x0 = x;
        f0 = KBF(x0);
        currentf = f0;
        c = 0; % Counter for new best solution

        maxf = f0; % Initial maximum f
        maxx = x0; % Store best x seen so far
        T = findT0(x,maxf); % Starting temperature
        etacount = 0; % Number of accepted trials at current temperature
        tempcount = 0; % Number of trials at current temperature
        D = 1*eye(length(x)); % Parks method of generating new solution
        R = (1/sqrt(n))*eye(length(x)); % Initial R such that dbar = 1

        for i = 1:iter
            trackx(:,i) = x;

            move = false;
            tempcount = tempcount + 1;
            trialx = newsolSA(x,g,D); % Generate new solution
            trialf = KBF(trialx);
            if trialf > maxf % New best solution
                maxf = trialf;
                maxx = trialx;
                move = true;
                c = 0;
            else
                c = c+1;
                if g == 1
                    p1 = exp(-(currentf-trialf)/T); % Acceptance probability
                else
                    dbar = sqrt(trace(R.^2));
                    p1 = exp(-(currentf-trialf)/(T*dbar));
                end
                p2 = rand(1);
                if p2 < p1
                    move = true;
                end
            end

            if move
                R = diag(abs(trialx - x));
                D = (1-alpha)*D + alpha*omega*R;
                x = trialx;
                currentf = trialf;
                etacount = etacount + 1;
            end

            if etacount == eta || tempcount == Lk
                T = a*T; % Decrease temperature
                etacount = 0;
                tempcount = 0;
            end

            if c >= cres % No new optmimum found in a while; restart
%                 disp('Restart!')
                x = maxx;
                currentf = maxf;
                D = 2*eye(length(x));
                R = (1/sqrt(n))*eye(length(x));
                c = 0;
            end
        end

        bx(:,run) = maxx;
        bf(run) = maxf;
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
    end
    scatter3(bx(1,:),bx(2,:),bf,20,'k','filled')
    if r == 1
        legend('KBF','Track','Start','Best')
    else
        legend('KBF','Best')
    end

end

