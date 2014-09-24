lw = 1;
myfs = 12;
%addpath results/pso_vs_ga2/;

figure
st = 'a11';
mylegends = {'ABP', 'RESP', 'PLETH', 'II', 'V', 'AVR'};
load(['data/' st 'm.mat']);
sig = val(:,end-4999:end-3750);
numsigs = size(sig,1);
for i = 1:numsigs
    sig(i,:) = (sig(i,:) - mean(sig(i,:))) / std(sig(i,:));
end
for i = 1:numsigs
	subplot(numsigs,2,i*2-1)
    plot(linspace(0,10,1250), sig(i,:),'b','LineWidth', lw)
    set(gca,'xgrid', 'on', 'ygrid', 'on')
    ylabel(mylegends(i));

    if i == 1
        title('Prior')
    end

    if i == numsigs
        xlabel(gca, 'Time (s)');
    else
        set(gca,'xticklabel', {})
    end
    h = findobj(gca, 'type', 'text');
    set(h, 'FontSize', myfs);
end

subplot(6,2,[2 4 6])
%semilogy(gbests,'r', 'linewidth', lw)
plot(gbests,'r', 'linewidth', lw)
title('Learning curve')
set(gca,'xgrid', 'on', 'ygrid', 'on')
legend('PSO');
xlabel(gca, 'Generations' );
ylabel(gca, 'MSE');
h = findobj(gca, 'type', 'text');
set(h, 'FontSize', myfs);

tmppos = get(gca,'Position');
tmppos(4) = tmppos(4) - 0.05;
tmppos(2) = tmppos(2) + 0.05;
set(gca, 'Position', tmppos);


yy = load(['data/' st '.missing']);
subplot(6,2,[8 10 12])
xd = linspace(0, 30, 3750);
plot(xd, yy,'linewidth', lw)
hold on
plot(xd, ff, 'r--', 'linewidth', lw)
set(gca, 'xlim', [0 10]);

set(gca,'xgrid', 'on', 'ygrid', 'on')
legend('target', 'PSO', 'Orientation','horizontal');
set(gca, 'ylim', [1100 2300]);
xlabel(gca, 'Time (s)');
ylabel(gca, 'MSE');
title('Reconstruction')
h = findobj(gca, 'type', 'text');
set(h, 'FontSize', myfs);

tmppos = get(gca,'Position');
tmppos(4) = tmppos(4) - 0.05;
set(gca, 'Position', tmppos);

print('reconstruction', '-depsc');

