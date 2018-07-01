clear
clc
close all

resWin = [6 12 18 24 ...        % go up in 6 hour blocks to one day
    48:24:28*24 ...                 % go up in one day blocks to one month
    840:7*24:26*7*24];           % then go up in one week blocks to one year

szPerYear = [1:9 10:10:610];
nYears = 5;

N = 1e3;

upperBoundR = zeros(length(szPerYear), length(resWin));
maxRpoiss = zeros(length(szPerYear),N);
maxRnb = zeros(length(szPerYear),N);

for iRate = 1:length(szPerYear)
    
    lambda = szPerYear(iRate) / 8760;
    
    circVar = nan(2, length(resWin));
    for iSim = 1:N
        
        nSeizures = nYears * szPerYear(iRate);
        
        % generate a Poisson process
        SzTimes = -log(rand(1, nSeizures)) / lambda;
        SzTimes = cumsum(SzTimes)';
        
        % generate negbin data
        SzTimes2 = nbinrnd(1, lambda, 1, nSeizures);
        SzTimes2 = cumsum(SzTimes2)';
        
        % windows in hr
        count =1;
        for n = resWin
            
            % remove clusters
            %             if n > 48
            %                 leadTime = 24;
            %             else
            %                 leadTime = 0;
            %             end
            %             ISI = [leadTime ; diff(SzTimes)];
            %             SzTimesTemp = SzTimes(ISI > leadTime);
            SzTimesTemp = SzTimes;
            % actual variance
            SzPhase = 2 * pi * mod(SzTimesTemp, n) / n;
            circVar(1,count) = circ_r(SzPhase);
            
            SzTimesTemp = SzTimes2;
            % actual variance
            SzPhase = 2 * pi * mod(SzTimesTemp, n) / n;
            circVar(2,count) = circ_r(SzPhase);
            
            count = count + 1;
            
        end
        
        maxRpoiss(iRate,iSim) = max(circVar(1,:));
        maxRnb(iRate,iSim) = max(circVar(2,:));
    end
    
    %     upperBoundR(iRate, :) = prctile(circVar, 95);
    %     plot(resWin, prctile(circVar, 95), 'color', [1-iRate/60, 0.2, 0.2], 'linewidth', 1.5);
    %     hold on;
    %     drawnow;
    
    fprintf('tested %d of %d rates\n',iRate,length(szPerYear));
end
%%

plot(nYears*szPerYear,prctile(maxRpoiss',95),'k')
hold on;
plot(nYears*szPerYear,prctile(maxRnb',95),'b')

set(gca,'box','off')
set(gca,'ytick',[0.1,0.2,0.3])
set(gca,'xlim',[0 600],'xtick',[20, 50,100,500])
set(gca,'fontname','arial','fontsize',8)

xlabel('Total Seizures')
ylabel('r-value (95 percentile)')
hold on
line([0 700],[0.198 0.198],'color','r')

set(gcf,'paperunits','centimeters','paperposition',[0 0 6 5]);
print(gcf,'STsim','-dpng','-r300');