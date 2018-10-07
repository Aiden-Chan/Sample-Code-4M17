% PlotKBF.m
% Plot KBF for n=2
close all

x1 = (0:0.1:10)';
[X1,X2] = meshgrid(x1);
F = zeros(size(X1));

% Calculate F
F = abs((cos(X1).^4+cos(X2).^4-2*(cos(X1).^2).*(cos(X2).^2))./sqrt(X1.^2+2*(X2.^2)));

% Discard infeasible solutions
for u = 1:size(F,1)
    for v = 1:size(F,2)
        if ~ feasible([X1(u,v);X2(u,v)])
            F(u,v) = -1;
        end
    end
end

% Plot
set(0,'DefaultFigurePosition',[2 42 681 642])
surf(X1,X2,F,'EdgeColor','none')
xlabel('x1')
ylabel('x2')
zlabel('f')
axis([0 10 0 10 0 0.4])
view([160 30]); % Roughly isometric
% view([0 90]); % Bird's eye